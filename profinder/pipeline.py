#!/usr/bin/env python3
"""
ProFinder — bacterial promoter identification pipeline.

Takes a single genome FASTA, annotates it with Prokka, extracts
intergenic regions, identifies operons, screens for HMM marker genes,
and outputs promoter sequences in 5'-to-3' orientation.

Usage
-----
    # Minimal (Prokka + IGR extraction + operon ID + promoter output):
    profinder -i genome.fasta -o results/

    # With HMM marker screening:
    profinder -i genome.fasta -o results/ \\
        --tigrfam /path/to/tigrfam.hmm --pfam /path/to/Pfam-A.hmm

    # Resume after interruption (skips completed steps):
    profinder -i genome.fasta -o results/

    # Force re-run from scratch:
    profinder -i genome.fasta -o results/ --force

    # Run specific steps:
    profinder -i genome.fasta -o results/ --start 3 --end 5

    # List steps:
    profinder --list

Steps
-----
    1.  Run Prokka (genome annotation)
    2.  Extract intergenic regions from GFF + FASTA
    3.  Identify operons from GFF
    4.  Run hmmsearch (TIGRfam + Pfam)       [bundled HMMs used by default]
    5.  Filter HMM output
    6.  Filter operons and add marker info
    7.  Match IGRs to marker operons
    8.  Extract promoter sequences (marker-filtered)
    9.  Extract all promoter sequences (orientation-based, no HMM filter)
    10. Annotate associated CDS (BLAST vs Swiss-Prot)
    11. Scan promoters for motifs (FIMO)
    12. Generate HTML report
"""

import argparse
import csv
import glob as globmod
import html as html_mod
import json
import os
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

from .config import Config
from .igr_extractor import extract_igrs


# =====================================================================
#  Operon identification (carried over from multi-genome pipeline)
# =====================================================================

def _parse_gff_for_operons(gff_path: str):
    """Parse a Prokka GFF into a DataFrame of CDS features."""
    rows = []
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9 or cols[2] != "CDS":
                continue
            rows.append({
                "seqid": cols[0],
                "start": int(cols[3]),
                "end": int(cols[4]),
                "strand": cols[6],
                "attributes": cols[8],
            })
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    df.sort_values(["seqid", "start"], inplace=True)
    df.reset_index(drop=True, inplace=True)
    return df


def _extract_gene_id(attributes: str) -> str:
    for attr in attributes.split(";"):
        if attr.startswith("ID="):
            return attr[3:]
    return ""


def _identify_operons(genes_df, max_internal_distance, min_flanking_distance):
    """Two-pass operon identification: proximity clustering then
    flanking-distance validation."""
    if genes_df.empty:
        return []

    # Pass 1: cluster consecutive genes within max_internal_distance
    raw_groups = []
    current_group = [0]
    for i in range(1, len(genes_df)):
        prev = genes_df.iloc[i - 1]
        curr = genes_df.iloc[i]
        same_contig = curr["seqid"] == prev["seqid"]
        close = (curr["start"] - prev["end"]) <= max_internal_distance
        if same_contig and close:
            current_group.append(i)
        else:
            raw_groups.append(current_group)
            current_group = [i]
    raw_groups.append(current_group)

    # Pass 2: validate flanking distances
    operons = []
    for group in raw_groups:
        first_i, last_i = group[0], group[-1]
        first_gene = genes_df.iloc[first_i]
        last_gene = genes_df.iloc[last_i]

        upstream_ok = True
        if first_i > 0:
            prev_gene = genes_df.iloc[first_i - 1]
            if prev_gene["seqid"] == first_gene["seqid"]:
                upstream_ok = (first_gene["start"] - prev_gene["end"]) >= min_flanking_distance

        downstream_ok = True
        if last_i < len(genes_df) - 1:
            next_gene = genes_df.iloc[last_i + 1]
            if next_gene["seqid"] == last_gene["seqid"]:
                downstream_ok = (next_gene["start"] - last_gene["end"]) >= min_flanking_distance

        if upstream_ok and downstream_ok:
            operons.append([genes_df.iloc[j] for j in group])

    return operons


# =====================================================================
#  HMM output parsing
# =====================================================================

_HMM_HEADER = [
    "target_name", "accession1", "query_name", "accession2",
    "full_sequence_evalue", "full_sequence_bitscore", "full_sequence_bias",
    "best_1_domain_evalue", "best_1_domain_score", "best_1_domain_bias",
    "exp", "reg", "glu", "ov", "env", "dom", "rep", "inc",
    "description_of_target",
]


def _parse_hmm_tblout(file_path: str) -> list:
    rows = []
    with open(file_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split(maxsplit=len(_HMM_HEADER) - 1)
            if len(parts) == len(_HMM_HEADER):
                rows.append(parts)
    return rows


# =====================================================================
#  Pipeline steps
# =====================================================================

def step01_run_prokka(cfg: Config, force: bool = False):
    """Annotate the input genome with Prokka."""
    if not force and cfg.gff_file.exists() and cfg.faa_file.exists():
        print("── Prokka output already exists, skipping ──")
        print("  Step 1 complete.\n")
        return

    print("── Running Prokka ──")
    cfg.prokka_dir.mkdir(parents=True, exist_ok=True)

    parts = []
    if cfg.conda_env_prokka:
        parts.append(
            f'eval "$(conda shell.bash hook 2>/dev/null)"; '
            f"conda activate {cfg.conda_env_prokka}; "
        )
    parts.append(
        f"{cfg.prokka_bin}"
        f" --outdir {cfg.prokka_dir}"
        f" --prefix {cfg.prokka_prefix}"
        f" --kingdom {cfg.prokka_kingdom}"
        f" --cpus {cfg.threads}"
        f" --force"
        f" {cfg.input_fasta}"
    )
    cmd = "set -euo pipefail; " + "".join(parts)
    result = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  STDERR: {result.stderr[-500:]}")
        sys.exit(f"Prokka failed with return code {result.returncode}")

    print(f"  GFF -> {cfg.gff_file}")
    print(f"  FAA -> {cfg.faa_file}")
    print("  Step 1 complete.\n")


def step02_extract_igrs(cfg: Config, force: bool = False):
    """Extract intergenic regions from the Prokka GFF and genome FASTA."""
    if not force and cfg.igr_summary.exists() and cfg.igr_fasta.exists():
        print("── IGR files already exist, skipping ──")
        print("  Step 2 complete.\n")
        return

    print("── Extracting intergenic regions ──")
    cfg.igr_dir.mkdir(parents=True, exist_ok=True)

    # Use the Prokka .fna (nucleotide FASTA with contig IDs matching the GFF)
    fasta_source = cfg.fna_file if cfg.fna_file.exists() else cfg.input_fasta
    igr_df = extract_igrs(cfg.gff_file, fasta_source,
                          size_min=cfg.igr_size_min, size_max=cfg.igr_size_max)

    igr_df.to_csv(cfg.igr_summary, sep="\t", index=False)
    print(f"  IGR summary ({len(igr_df)} regions) -> {cfg.igr_summary}")

    # Write IGR FASTA
    with open(cfg.igr_fasta, "w") as fh:
        for _, row in igr_df.iterrows():
            fh.write(f">{row['igr_id']}_{row['contig']}_{row['orientation']}"
                     f"_{row['left_gene']}_{row['right_gene']}\n")
            fh.write(f"{row['sequence']}\n")
    print(f"  IGR FASTA -> {cfg.igr_fasta}")
    print("  Step 2 complete.\n")


def step03_identify_operons(cfg: Config, force: bool = False):
    """Identify operons from the Prokka GFF."""
    if not force and cfg.operon_file.exists():
        print("── Operon file already exists, skipping ──")
        print("  Step 3 complete.\n")
        return

    print("── Identifying operons ──")
    genes_df = _parse_gff_for_operons(str(cfg.gff_file))
    if genes_df.empty:
        print("  WARNING: no CDS features found in GFF")
        pd.DataFrame().to_csv(cfg.operon_file, sep="\t", index=False)
        print("  Step 3 complete.\n")
        return

    operons = _identify_operons(genes_df, cfg.max_internal_distance,
                                cfg.min_flanking_distance)
    source = cfg.prokka_prefix

    rows = []
    for op_idx, operon in enumerate(operons, start=1):
        operon_id = f"operon{op_idx}"
        first_genes = set()
        if operon[0]["strand"] == "+":
            first_genes.add((0, "yes_1"))
        if operon[-1]["strand"] == "-":
            first_genes.add((len(operon) - 1, "yes_2"))

        for gi, gene in enumerate(operon):
            gene_id = _extract_gene_id(gene["attributes"])
            label = "No"
            for (idx, tag) in first_genes:
                if gi == idx:
                    label = tag
                    break
            rows.append({
                "Operon": operon_id,
                "SeqID": gene["seqid"],
                "Start": gene["start"],
                "End": gene["end"],
                "Attributes": gene["attributes"],
                "Gene": gene_id,
                "SourceFile": source,
                "FirstGene": label,
            })

    df = pd.DataFrame(rows)
    df.to_csv(cfg.operon_file, sep="\t", index=False)
    print(f"  Operons ({len(operons)} operons, {len(df)} gene rows) -> {cfg.operon_file}")
    print("  Step 3 complete.\n")


def step04_run_hmmsearch(cfg: Config, force: bool = False):
    """Run hmmsearch (TIGRfam + Pfam) on the Prokka .faa."""
    tigr_ok = cfg.tigrfam_hmm and str(cfg.tigrfam_hmm) != "." and cfg.tigrfam_hmm.is_file()
    pfam_ok = cfg.pfam_hmm and str(cfg.pfam_hmm) != "." and cfg.pfam_hmm.is_file()
    if not tigr_ok or not pfam_ok:
        print("── HMM database files not available ──")
        if not tigr_ok:
            print(f"  tigrfam: {cfg.tigrfam_hmm} (found={tigr_ok})")
        if not pfam_ok:
            print(f"  pfam:    {cfg.pfam_hmm} (found={pfam_ok})")
        print("  Provide --tigrfam and --pfam, or reinstall to restore bundled HMMs.")
        print("  Step 4 complete.\n")
        return

    print(f"  Using TIGRfam: {cfg.tigrfam_hmm}")
    print(f"  Using Pfam:    {cfg.pfam_hmm}")

    tigr_out = cfg.hmm_dir / "hmm_tigrfam.tblout"
    pfam_out = cfg.hmm_dir / "hmm_pfam.tblout"

    if not force and tigr_out.exists() and pfam_out.exists():
        print("── HMM output already exists, skipping ──")
        print("  Step 4 complete.\n")
        return

    print("── Running hmmsearch ──")
    cfg.hmm_dir.mkdir(parents=True, exist_ok=True)

    for db_label, db_path, out_file in [
        ("TIGRfam", cfg.tigrfam_hmm, tigr_out),
        ("Pfam", cfg.pfam_hmm, pfam_out),
    ]:
        log_file = out_file.with_suffix(".log")
        parts = []
        if cfg.conda_env_hmm:
            parts.append(
                f'eval "$(conda shell.bash hook 2>/dev/null)"; '
                f"conda activate {cfg.conda_env_hmm}; "
            )
        parts.append(
            f"{cfg.hmmsearch_bin}"
            f" --tblout {out_file}"
            f" -o {log_file}"
            f" --noali --cpu {cfg.threads}"
            f" {db_path} {cfg.faa_file}"
        )
        cmd = "set -euo pipefail; " + "".join(parts)
        print(f"  {db_label}...")
        result = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True)
        if result.returncode != 0:
            print(f"  WARNING: hmmsearch {db_label} failed: {result.stderr[-300:]}")
            return

    print("  Step 4 complete.\n")


def step05_filter_hmm(cfg: Config, force: bool = False):
    """Consolidate and filter HMM hits."""
    if not force and cfg.hmm_filtered.exists():
        print("── Filtered HMM file already exists, skipping ──")
        print("  Step 5 complete.\n")
        return

    tigr_out = cfg.hmm_dir / "hmm_tigrfam.tblout"
    pfam_out = cfg.hmm_dir / "hmm_pfam.tblout"

    if not tigr_out.exists() and not pfam_out.exists():
        print("── No HMM output to filter, skipping ──")
        print("  Step 5 complete.\n")
        return

    print("── Filtering HMM output ──")
    all_rows = []
    for f in [tigr_out, pfam_out]:
        if f.exists():
            all_rows.extend(_parse_hmm_tblout(str(f)))

    if not all_rows:
        print("  No HMM hits found.")
        pd.DataFrame(columns=_HMM_HEADER).to_csv(cfg.hmm_combined, sep="\t", index=False)
        pd.DataFrame(columns=_HMM_HEADER).to_csv(cfg.hmm_filtered, sep="\t", index=False)
        print("  Step 5 complete.\n")
        return

    df = pd.DataFrame(all_rows, columns=_HMM_HEADER)
    df.to_csv(cfg.hmm_combined, sep="\t", index=False)

    df["full_sequence_bitscore"] = pd.to_numeric(df["full_sequence_bitscore"], errors="coerce")
    df = df[df["full_sequence_bitscore"] >= cfg.hmm_bitscore_min]

    # Keep best hit per target+accession
    if not df.empty:
        best_idx = df.groupby(["target_name", "accession2"])["full_sequence_bitscore"].idxmax()
        df = df.loc[best_idx]

    df.to_csv(cfg.hmm_filtered, sep="\t", index=False)
    print(f"  Filtered HMM ({len(df)} rows) -> {cfg.hmm_filtered}")
    print("  Step 5 complete.\n")


def step06_filter_operons_add_markers(cfg: Config, force: bool = False):
    """Filter operons and merge with HMM marker gene info."""
    if not force and cfg.operon_filtered_markers.exists():
        print("── Filtered operons with markers already exists, skipping ──")
        print("  Step 6 complete.\n")
        return

    print("── Filtering operons ──")
    df = pd.read_csv(cfg.operon_file, sep="\t")
    if df.empty:
        print("  No operons to filter.")
        df.to_csv(cfg.operon_filtered_markers, sep="\t", index=False)
        print("  Step 6 complete.\n")
        return

    # Keep operons on a single contig
    single_contig = df.groupby("Operon").filter(
        lambda g: g["SeqID"].nunique() == 1
    )
    # Keep operons with at least one promoter-boundary gene
    has_first = single_contig.groupby("Operon").filter(
        lambda g: g["FirstGene"].isin(["yes_1", "yes_2"]).any()
    )
    # Remove operons with BOTH yes_1 and yes_2 (ambiguous)
    def _single_direction(g):
        labels = set(g["FirstGene"].unique())
        return not {"yes_1", "yes_2"}.issubset(labels)

    filtered = has_first.groupby("Operon").filter(_single_direction)
    filtered.to_csv(cfg.operon_filtered, sep="\t", index=False)
    print(f"  Filtered operons ({len(filtered)} rows) -> {cfg.operon_filtered}")

    # Merge with HMM markers if available
    if cfg.hmm_filtered.exists() and cfg.hmm_filtered.stat().st_size > 0:
        print("── Adding HMM marker info ──")
        try:
            hmm = pd.read_csv(cfg.hmm_filtered, sep="\t",
                               usecols=["target_name", "accession2"])
        except (pd.errors.EmptyDataError, ValueError):
            hmm = pd.DataFrame(columns=["target_name", "accession2"])

        merged = pd.merge(filtered, hmm, how="left",
                           left_on="Gene", right_on="target_name")
        valid_operons = merged.loc[
            merged["accession2"].notna() &
            (merged["accession2"] != "") &
            (merged["accession2"] != 0),
            "Operon"
        ].unique()
        merged = merged[merged["Operon"].isin(valid_operons)]
        merged.to_csv(cfg.operon_filtered_markers, sep="\t", index=False)
        print(f"  Operons with markers ({len(merged)} rows) "
              f"-> {cfg.operon_filtered_markers}")
    else:
        # No HMM data: use all filtered operons as-is
        print("  No HMM data available; using all filtered operons as markers")
        filtered["target_name"] = filtered["Gene"]
        filtered["accession2"] = "no_hmm"
        filtered.to_csv(cfg.operon_filtered_markers, sep="\t", index=False)

    print("  Step 6 complete.\n")


def step07_match_igrs_to_markers(cfg: Config, force: bool = False):
    """Match IGRs to marker operon genes."""
    if not force and cfg.promoter_markers.exists():
        print("── Promoter markers file already exists, skipping ──")
        print("  Step 7 complete.\n")
        return

    print("── Matching IGRs to marker operons ──")
    igr = pd.read_csv(cfg.igr_summary, sep="\t")
    try:
        markers = pd.read_csv(cfg.operon_filtered_markers, sep="\t")
    except (FileNotFoundError, pd.errors.EmptyDataError):
        print("  No marker data available.")
        pd.DataFrame().to_csv(cfg.promoter_markers, sep="\t", index=False)
        print("  Step 7 complete.\n")
        return

    marker_genes = set(markers["Gene"].dropna())

    # Match: an IGR is a promoter candidate if one of its flanking genes
    # is a marker gene, and the orientation is promoter-relevant.
    igr = igr.copy()
    igr["marker_match"] = "none"

    # CO_F: right_gene (downstream on + strand) is the marker
    co_f = (igr["orientation"] == "CO_F") & igr["right_gene"].isin(marker_genes)
    igr.loc[co_f, "marker_match"] = "CO_F"

    # CO_R: left_gene (downstream on - strand, since gene runs right to left) is the marker
    co_r = (igr["orientation"] == "CO_R") & igr["left_gene"].isin(marker_genes)
    igr.loc[co_r, "marker_match"] = "CO_R"

    # DP: either flanking gene could be a marker
    dp_left = (igr["orientation"] == "DP") & igr["left_gene"].isin(marker_genes)
    dp_right = (igr["orientation"] == "DP") & igr["right_gene"].isin(marker_genes)
    igr.loc[dp_left | dp_right, "marker_match"] = "DP"

    matched = igr[igr["marker_match"] != "none"].copy()
    matched.to_csv(cfg.promoter_markers, sep="\t", index=False)
    print(f"  Matched promoter markers ({len(matched)} rows) -> {cfg.promoter_markers}")
    print("  Step 7 complete.\n")


def _reverse_complement(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def _write_fasta(df, output_path, short_header=False):
    """Write promoter FASTA from a DataFrame with sequence_5p_to_3p."""
    with open(output_path, "w") as fh:
        for _, row in df.iterrows():
            seq = row["sequence_5p_to_3p"]
            igr_id = row["igr_id"]
            orient = row["orientation"]
            if short_header:
                header = f">{igr_id}_{orient}"
            else:
                header = (f">{igr_id}_{row['contig']}_{orient}"
                          f"_{row['left_gene']}_{row['right_gene']}")
            fh.write(f"{header}\n{seq}\n")
    print(f"  FASTA ({len(df)} seqs, short={short_header}) -> {output_path}")


def step08_extract_marker_promoters(cfg: Config, force: bool = False):
    """Extract marker-filtered promoter sequences in 5'→3' orientation."""
    if not force and cfg.promoter_fasta_short.exists():
        print("── Marker promoter FASTA already exists, skipping ──")
        print("  Step 8 complete.\n")
        return

    print("── Extracting marker-filtered promoter sequences ──")
    try:
        df = pd.read_csv(cfg.promoter_markers, sep="\t")
    except (FileNotFoundError, pd.errors.EmptyDataError):
        print("  No promoter markers to extract.")
        print("  Step 8 complete.\n")
        return

    if df.empty:
        print("  No promoter markers to extract.")
        print("  Step 8 complete.\n")
        return

    # Remove divergent promoters (ambiguous orientation)
    df = df[df["orientation"] != "DP"].copy()

    # Orient 5'→3': CO_R sequences need reverse complement
    df["sequence_5p_to_3p"] = df.apply(
        lambda r: _reverse_complement(r["sequence"]) if r["orientation"] == "CO_R"
        else r["sequence"],
        axis=1,
    )

    _write_fasta(df, cfg.promoter_fasta, short_header=False)
    _write_fasta(df, cfg.promoter_fasta_short, short_header=True)
    print("  Step 8 complete.\n")


def step09_extract_all_promoters(cfg: Config, force: bool = False):
    """Extract ALL promoter-orientation IGRs (no HMM filter)."""
    if not force and cfg.all_promoter_fasta_short.exists():
        print("── All-promoter FASTA already exists, skipping ──")
        print("  Step 9 complete.\n")
        return

    print("── Extracting all promoter-orientation IGRs ──")
    igr = pd.read_csv(cfg.igr_summary, sep="\t")

    # Keep only promoter-relevant orientations (exclude convergent/terminator)
    promoter_orientations = igr[igr["orientation"].isin(["CO_F", "CO_R", "DP"])].copy()
    n_dropped = len(igr) - len(promoter_orientations)
    if n_dropped:
        print(f"  Dropped {n_dropped} convergent (terminator) IGRs")

    # Remove divergent promoters (can't unambiguously orient)
    df = promoter_orientations[promoter_orientations["orientation"] != "DP"].copy()
    n_dp = len(promoter_orientations) - len(df)
    if n_dp:
        print(f"  Excluded {n_dp} divergent promoters (ambiguous orientation)")

    # Orient 5'→3'
    df["sequence_5p_to_3p"] = df.apply(
        lambda r: _reverse_complement(r["sequence"]) if r["orientation"] == "CO_R"
        else r["sequence"],
        axis=1,
    )

    _write_fasta(df, cfg.all_promoter_fasta, short_header=False)
    _write_fasta(df, cfg.all_promoter_fasta_short, short_header=True)
    print("  Step 9 complete.\n")


# =====================================================================
#  Step 10 — Annotate associated CDS via BLAST
# =====================================================================

_BLAST_COLS = [
    "qseqid", "sseqid", "stitle", "pident", "length", "mismatch",
    "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore",
]


def _run_blastp(cfg, query_fasta, output_tsv):
    """Run blastp against Swiss-Prot (remote or local)."""
    parts = []
    if cfg.conda_env_blast:
        parts.append(
            f'eval "$(conda shell.bash hook 2>/dev/null)"; '
            f"conda activate {cfg.conda_env_blast}; "
        )

    if cfg.blast_local_db is not None and cfg.blast_local_db.is_file():
        db_flag = f"-db {cfg.blast_local_db}"
    else:
        db_flag = f"-db {cfg.blast_db} -remote"

    parts.append(
        f"{cfg.blastp_bin}"
        f" -query {query_fasta}"
        f" {db_flag}"
        f" -evalue {cfg.blast_evalue}"
        f" -max_target_seqs {cfg.blast_max_targets}"
        f" -outfmt '6 {' '.join(_BLAST_COLS)}'"
        f" -out {output_tsv}"
    )
    if cfg.blast_local_db is not None and cfg.blast_local_db.is_file():
        parts[-1] += f" -num_threads {cfg.threads}"

    cmd = "set -euo pipefail; " + "".join(parts)
    return cmd


def _parse_swissprot_title(stitle):
    """Extract protein name and organism from Swiss-Prot stitle field.

    Swiss-Prot stitle format:
        sp|P12345|NAME_SPECIES Protein description OS=Organism Name OX=12345 ...
    """
    protein_name = stitle
    organism = ""

    os_idx = stitle.find(" OS=")
    if os_idx != -1:
        protein_name = stitle[:os_idx]
        # Strip the leading sp|...|... part
        sp_end = protein_name.find(" ")
        if sp_end != -1:
            protein_name = protein_name[sp_end + 1:]

        rest = stitle[os_idx + 4:]
        ox_idx = rest.find(" OX=")
        organism = rest[:ox_idx] if ox_idx != -1 else rest.strip()
    else:
        sp_end = stitle.find(" ")
        if sp_end != -1:
            protein_name = stitle[sp_end + 1:]

    return protein_name.strip(), organism.strip()


def step10_annotate_cds(cfg: Config, force: bool = False):
    """BLAST Prokka proteins against Swiss-Prot to annotate the CDS
    immediately associated with each promoter."""
    if not force and cfg.blast_annotated.exists():
        print("── CDS annotations already exist, skipping ──")
        print("  Step 10 complete.\n")
        return

    cfg.blast_dir.mkdir(parents=True, exist_ok=True)

    # Collect all unique CDS gene IDs associated with promoters.
    # For CO_F the associated CDS is right_gene; for CO_R it's left_gene.
    promoter_files = [cfg.promoter_markers, cfg.igr_summary]
    gene_ids = set()
    for pf in promoter_files:
        if pf.exists():
            try:
                df = pd.read_csv(pf, sep="\t")
                if "right_gene" in df.columns:
                    gene_ids.update(df.loc[df["orientation"] == "CO_F", "right_gene"].dropna())
                    gene_ids.update(df.loc[df["orientation"] == "CO_R", "left_gene"].dropna())
            except (pd.errors.EmptyDataError, KeyError):
                pass

    if not gene_ids:
        print("── No promoter-associated CDS to annotate ──")
        pd.DataFrame().to_csv(cfg.blast_annotated, sep="\t", index=False)
        print("  Step 10 complete.\n")
        return

    # Extract only those protein sequences from Prokka .faa
    subset_faa = cfg.blast_dir / "query_proteins.faa"
    kept = 0
    with open(subset_faa, "w") as out:
        for rec in SeqIO.parse(str(cfg.faa_file), "fasta"):
            if rec.id in gene_ids:
                SeqIO.write(rec, out, "fasta")
                kept += 1
    print(f"── BLASTing {kept} CDS proteins against Swiss-Prot ──")

    if kept == 0:
        pd.DataFrame().to_csv(cfg.blast_annotated, sep="\t", index=False)
        print("  Step 10 complete.\n")
        return

    # Run BLAST
    if not cfg.blast_raw.exists() or force:
        cmd = _run_blastp(cfg, subset_faa, cfg.blast_raw)
        is_remote = not (cfg.blast_local_db is not None
                         and cfg.blast_local_db.is_file())
        if is_remote:
            print("  Using NCBI remote BLAST (this may take several minutes)...")
        result = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True)
        if result.returncode != 0:
            print(f"  WARNING: blastp failed: {result.stderr[-500:]}")
            pd.DataFrame().to_csv(cfg.blast_annotated, sep="\t", index=False)
            print("  Step 10 complete.\n")
            return

    # Parse results
    try:
        blast_df = pd.read_csv(cfg.blast_raw, sep="\t", header=None,
                                names=_BLAST_COLS)
    except pd.errors.EmptyDataError:
        blast_df = pd.DataFrame(columns=_BLAST_COLS)

    if blast_df.empty:
        print("  No BLAST hits found.")
        blast_df.to_csv(cfg.blast_annotated, sep="\t", index=False)
        print("  Step 10 complete.\n")
        return

    # Keep best hit per query
    blast_df["evalue"] = pd.to_numeric(blast_df["evalue"], errors="coerce")
    blast_df["bitscore"] = pd.to_numeric(blast_df["bitscore"], errors="coerce")
    best_idx = blast_df.groupby("qseqid")["bitscore"].idxmax()
    best = blast_df.loc[best_idx].copy()

    # Parse Swiss-Prot stitle to extract protein name and organism
    parsed = best["stitle"].apply(_parse_swissprot_title)
    best["protein_name"] = [p[0] for p in parsed]
    best["organism"] = [p[1] for p in parsed]

    best.to_csv(cfg.blast_annotated, sep="\t", index=False)
    print(f"  Annotated {len(best)} CDS -> {cfg.blast_annotated}")
    print("  Step 10 complete.\n")


# =====================================================================
#  Step 11 — FIMO motif scanning
# =====================================================================

def step11_run_fimo(cfg: Config, force: bool = False):
    """Scan promoter sequences against motif databases with FIMO."""
    if not force and cfg.fimo_combined.exists():
        print("── FIMO results already exist, skipping ──")
        print("  Step 11 complete.\n")
        return

    # Determine which promoter FASTA to scan
    fasta_to_scan = None
    for candidate in [cfg.all_promoter_fasta, cfg.promoter_fasta]:
        if candidate.exists() and candidate.stat().st_size > 0:
            fasta_to_scan = candidate
            break
    if fasta_to_scan is None:
        print("── No promoter FASTA files found, skipping ──")
        print("  Step 11 complete.\n")
        return

    # Find .meme files
    motifs_dir = cfg.motifs_dir
    if motifs_dir is None or not motifs_dir.is_dir():
        # Fall back to bundled motifs alongside this script
        bundled = Path(__file__).parent / "motifs"
        if bundled.is_dir():
            motifs_dir = bundled
        else:
            print("── No motifs directory found, skipping ──")
            print("  Step 11 complete.\n")
            return

    meme_files = sorted(motifs_dir.glob("*.meme"))
    if not meme_files:
        print(f"── No .meme files in {motifs_dir}, skipping ──")
        print("  Step 11 complete.\n")
        return

    print(f"── Running FIMO against {len(meme_files)} motif databases ──")
    print(f"  Scanning: {fasta_to_scan.name}")
    cfg.fimo_dir.mkdir(parents=True, exist_ok=True)

    all_fimo_dfs = []
    for meme_file in meme_files:
        db_name = meme_file.stem
        out_subdir = cfg.fimo_dir / db_name
        tsv_out = out_subdir / "fimo.tsv"

        if not force and tsv_out.exists():
            print(f"  {db_name}: already done")
        else:
            out_subdir.mkdir(parents=True, exist_ok=True)
            parts = []
            if cfg.conda_env_meme:
                parts.append(
                    f'eval "$(conda shell.bash hook 2>/dev/null)"; '
                    f"conda activate {cfg.conda_env_meme}; "
                )
            parts.append(
                f"{cfg.fimo_bin}"
                f" --thresh {cfg.fimo_threshold}"
                f" --oc {out_subdir}"
                f" {meme_file}"
                f" {fasta_to_scan}"
            )
            cmd = "set -euo pipefail; " + "".join(parts)
            result = subprocess.run(["bash", "-c", cmd],
                                    capture_output=True, text=True)
            if result.returncode != 0:
                print(f"  WARNING: FIMO {db_name} failed: {result.stderr[-300:]}")
                continue
            print(f"  {db_name}: done")

        # Parse FIMO TSV output
        if tsv_out.exists():
            try:
                fimo_df = pd.read_csv(tsv_out, sep="\t", comment="#")
                if not fimo_df.empty:
                    fimo_df["motif_database"] = db_name
                    all_fimo_dfs.append(fimo_df)
            except pd.errors.EmptyDataError:
                pass

    if all_fimo_dfs:
        combined = pd.concat(all_fimo_dfs, ignore_index=True)
        combined.to_csv(cfg.fimo_combined, sep="\t", index=False)
        print(f"  Combined FIMO hits ({len(combined)} rows) -> {cfg.fimo_combined}")
    else:
        print("  No FIMO hits found across any database.")
        pd.DataFrame().to_csv(cfg.fimo_combined, sep="\t", index=False)

    print("  Step 11 complete.\n")


# =====================================================================
#  Step 12 — HTML report
# =====================================================================

def _get_associated_gene(row):
    """Return the gene ID of the CDS immediately downstream of this promoter."""
    if row["orientation"] == "CO_F":
        return row.get("right_gene", "")
    elif row["orientation"] == "CO_R":
        return row.get("left_gene", "")
    return ""


def _build_motif_diagram_svg(seq_len, motif_hits, width=700, height=50):
    """Build an inline SVG showing motif positions along a promoter sequence.

    Returns an SVG string.
    """
    if seq_len == 0:
        return ""

    margin = 10
    track_y = 25
    track_h = 10
    bar_w = width - 2 * margin
    scale = bar_w / seq_len

    # Colour palette for motif databases
    db_colours = {
        "collectf": "#2196F3",
        "prodoric_2021.9": "#FF9800",
        "prodoric_2021": "#FF9800",
        "regtransbase": "#4CAF50",
    }
    default_colour = "#9C27B0"

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'style="background:#fafafa;border:1px solid #ddd;border-radius:4px;">',
        # Sequence track (grey bar)
        f'<rect x="{margin}" y="{track_y}" width="{bar_w}" height="{track_h}" '
        f'fill="#e0e0e0" rx="3"/>',
        # 5' and 3' labels
        f'<text x="{margin}" y="16" font-size="11" fill="#666">5\'</text>',
        f'<text x="{width - margin - 10}" y="16" font-size="11" fill="#666">3\'</text>',
        # Length label
        f'<text x="{width / 2}" y="16" font-size="10" fill="#999" text-anchor="middle">'
        f'{seq_len} bp</text>',
    ]

    for hit in motif_hits:
        start = hit.get("start", 0)
        stop = hit.get("stop", start + 10)
        db = hit.get("motif_database", "").lower()
        motif_name = hit.get("motif_id", "")
        colour = db_colours.get(db, default_colour)

        x = margin + (start - 1) * scale
        w = max((stop - start + 1) * scale, 3)  # at least 3px wide

        parts.append(
            f'<rect x="{x:.1f}" y="{track_y - 2}" width="{w:.1f}" height="{track_h + 4}" '
            f'fill="{colour}" opacity="0.7" rx="2">'
            f'<title>{html_mod.escape(motif_name)} ({db}) pos {start}-{stop}</title>'
            f'</rect>'
        )

    parts.append('</svg>')
    return "\n".join(parts)


def step12_generate_report(cfg: Config, force: bool = False):
    """Generate an HTML report combining promoters, CDS annotations,
    and motif hits with visual sequence diagrams."""
    if not force and cfg.report_html.exists():
        print("── HTML report already exists, skipping ──")
        print("  Step 12 complete.\n")
        return

    print("── Generating HTML report ──")

    # Load promoter data (prefer all_promoters for completeness)
    promoter_df = None
    for candidate in [cfg.all_promoter_fasta]:
        igr_path = cfg.igr_summary
        if igr_path.exists():
            try:
                df = pd.read_csv(igr_path, sep="\t")
                # Filter to promoter orientations only
                df = df[df["orientation"].isin(["CO_F", "CO_R"])].copy()
                if not df.empty:
                    promoter_df = df
                    break
            except pd.errors.EmptyDataError:
                pass

    if promoter_df is None or promoter_df.empty:
        print("  No promoter data available for report.")
        print("  Step 12 complete.\n")
        return

    # Add associated gene column
    promoter_df["associated_gene"] = promoter_df.apply(_get_associated_gene, axis=1)

    # Orient sequences 5'→3'
    promoter_df["sequence_5p_to_3p"] = promoter_df.apply(
        lambda r: _reverse_complement(r["sequence"]) if r["orientation"] == "CO_R"
        else r["sequence"],
        axis=1,
    )

    # Load BLAST annotations
    blast_map = {}
    if cfg.blast_annotated.exists():
        try:
            blast_df = pd.read_csv(cfg.blast_annotated, sep="\t")
            for _, row in blast_df.iterrows():
                blast_map[row["qseqid"]] = {
                    "protein_name": row.get("protein_name", ""),
                    "organism": row.get("organism", ""),
                    "sseqid": row.get("sseqid", ""),
                    "pident": row.get("pident", ""),
                    "evalue": row.get("evalue", ""),
                    "bitscore": row.get("bitscore", ""),
                }
        except (pd.errors.EmptyDataError, KeyError):
            pass

    # Load FIMO hits grouped by promoter
    fimo_by_promoter = {}
    if cfg.fimo_combined.exists():
        try:
            fimo_df = pd.read_csv(cfg.fimo_combined, sep="\t")
            # FIMO output uses "sequence_name" for the FASTA header
            seq_col = "sequence_name"
            if seq_col not in fimo_df.columns:
                # Try alternative column names
                for alt in ["sequence name", "#pattern name"]:
                    if alt in fimo_df.columns:
                        seq_col = alt
                        break

            if seq_col in fimo_df.columns:
                for seq_name, grp in fimo_df.groupby(seq_col):
                    # Extract the igr_id from the FASTA header
                    # Headers look like: igr_000001_contigX_CO_F_geneA_geneB
                    igr_id = seq_name.split("_")[0] + "_" + seq_name.split("_")[1]
                    if igr_id not in fimo_by_promoter:
                        fimo_by_promoter[igr_id] = []
                    for _, hit in grp.iterrows():
                        fimo_by_promoter[igr_id].append({
                            "motif_id": str(hit.get("motif_id", hit.get("#pattern name", ""))),
                            "motif_alt_id": str(hit.get("motif_alt_id", "")),
                            "start": int(hit.get("start", 0)),
                            "stop": int(hit.get("stop", 0)),
                            "strand": str(hit.get("strand", "+")),
                            "p-value": float(hit.get("p-value", 1)),
                            "q-value": float(hit.get("q-value", 1)) if "q-value" in hit.index else None,
                            "matched_sequence": str(hit.get("matched_sequence", "")),
                            "motif_database": str(hit.get("motif_database", "")),
                        })
        except (pd.errors.EmptyDataError, KeyError):
            pass

    # Build HTML
    n_promoters = len(promoter_df)
    n_annotated = sum(1 for _, r in promoter_df.iterrows()
                      if r["associated_gene"] in blast_map)
    n_with_motifs = sum(1 for _, r in promoter_df.iterrows()
                        if r["igr_id"] in fimo_by_promoter)

    genome_name = cfg.input_fasta.stem

    html_parts = [_REPORT_HTML_HEAD.format(
        genome_name=html_mod.escape(genome_name),
        n_promoters=n_promoters,
        n_annotated=n_annotated,
        n_with_motifs=n_with_motifs,
        total_fimo=sum(len(v) for v in fimo_by_promoter.values()),
    )]

    # Legend
    html_parts.append("""
    <div class="legend">
        <strong>Motif database colours:</strong>
        <span class="legend-item"><span class="legend-swatch" style="background:#2196F3"></span>CollecTF</span>
        <span class="legend-item"><span class="legend-swatch" style="background:#FF9800"></span>PRODORIC</span>
        <span class="legend-item"><span class="legend-swatch" style="background:#4CAF50"></span>RegTransBase</span>
    </div>
    """)

    # Table
    html_parts.append("""
    <table>
    <thead>
    <tr>
        <th>Promoter ID</th>
        <th>Contig</th>
        <th>Position</th>
        <th>Length</th>
        <th>Orientation</th>
        <th>Associated CDS</th>
        <th>Protein (BLAST)</th>
        <th>Organism</th>
        <th>BLAST E-value</th>
        <th>Motif hits</th>
        <th>Sequence diagram</th>
    </tr>
    </thead>
    <tbody>
    """)

    for _, row in promoter_df.iterrows():
        igr_id = row["igr_id"]
        gene = row["associated_gene"]
        blast_info = blast_map.get(gene, {})
        motif_hits = fimo_by_promoter.get(igr_id, [])

        protein = html_mod.escape(str(blast_info.get("protein_name", "—")))
        organism = html_mod.escape(str(blast_info.get("organism", "—")))
        evalue = blast_info.get("evalue", "—")
        if isinstance(evalue, float):
            evalue = f"{evalue:.1e}"

        motif_count = len(motif_hits)
        motif_names = ", ".join(sorted(set(
            h.get("motif_alt_id") or h.get("motif_id", "") for h in motif_hits
        )))[:80]

        svg = _build_motif_diagram_svg(row["length"], motif_hits)

        orient_class = "co-f" if row["orientation"] == "CO_F" else "co-r"

        html_parts.append(f"""
        <tr>
            <td><code>{html_mod.escape(igr_id)}</code></td>
            <td>{html_mod.escape(str(row['contig']))}</td>
            <td>{row['start']:,}–{row['end']:,}</td>
            <td>{row['length']}</td>
            <td><span class="orient {orient_class}">{row['orientation']}</span></td>
            <td><code>{html_mod.escape(str(gene))}</code></td>
            <td>{protein if protein != '—' else '<span class="na">—</span>'}</td>
            <td><em>{organism if organism != '—' else '<span class="na">—</span>'}</em></td>
            <td>{evalue}</td>
            <td>{motif_count}{(' (' + html_mod.escape(motif_names) + ')') if motif_names else ''}</td>
            <td>{svg}</td>
        </tr>
        """)

    html_parts.append("</tbody></table>")

    # Motif detail section
    if fimo_by_promoter:
        html_parts.append("<h2>Motif hit details</h2>")
        html_parts.append("""
        <table class="motif-detail">
        <thead>
        <tr>
            <th>Promoter</th>
            <th>Motif ID</th>
            <th>Motif name</th>
            <th>Database</th>
            <th>Position</th>
            <th>Strand</th>
            <th>p-value</th>
            <th>Matched sequence</th>
        </tr>
        </thead>
        <tbody>
        """)
        for igr_id, hits in sorted(fimo_by_promoter.items()):
            for hit in sorted(hits, key=lambda h: h["start"]):
                pval = hit["p-value"]
                pval_str = f"{pval:.2e}" if isinstance(pval, float) else str(pval)
                html_parts.append(f"""
                <tr>
                    <td><code>{html_mod.escape(igr_id)}</code></td>
                    <td>{html_mod.escape(hit['motif_id'])}</td>
                    <td>{html_mod.escape(hit.get('motif_alt_id', ''))}</td>
                    <td>{html_mod.escape(hit['motif_database'])}</td>
                    <td>{hit['start']}–{hit['stop']}</td>
                    <td>{hit['strand']}</td>
                    <td>{pval_str}</td>
                    <td><code>{html_mod.escape(hit.get('matched_sequence', ''))}</code></td>
                </tr>
                """)
        html_parts.append("</tbody></table>")

    html_parts.append(_REPORT_HTML_FOOT)

    with open(cfg.report_html, "w") as fh:
        fh.write("\n".join(html_parts))
    print(f"  Report -> {cfg.report_html}")
    print("  Step 12 complete.\n")


# ── HTML template fragments ──────────────────────────────────────────

_REPORT_HTML_HEAD = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Promoter Report — {genome_name}</title>
<style>
    * {{ margin: 0; padding: 0; box-sizing: border-box; }}
    body {{
        font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
        color: #333; background: #fff; padding: 24px 32px;
        line-height: 1.5;
    }}
    h1 {{ font-size: 1.6rem; margin-bottom: 8px; }}
    h2 {{ font-size: 1.2rem; margin: 32px 0 12px; border-bottom: 2px solid #eee; padding-bottom: 6px; }}
    .summary {{
        display: flex; gap: 24px; margin: 16px 0 24px;
        flex-wrap: wrap;
    }}
    .summary .stat {{
        background: #f5f7fa; border-radius: 8px;
        padding: 12px 20px; min-width: 140px;
    }}
    .summary .stat .label {{ font-size: 0.8rem; color: #666; text-transform: uppercase; letter-spacing: 0.5px; }}
    .summary .stat .value {{ font-size: 1.5rem; font-weight: 600; color: #1a73e8; }}
    .legend {{
        margin: 12px 0 20px; font-size: 0.85rem; color: #555;
    }}
    .legend-item {{ margin-left: 16px; }}
    .legend-swatch {{
        display: inline-block; width: 14px; height: 14px;
        border-radius: 3px; vertical-align: middle; margin-right: 4px;
    }}
    table {{
        border-collapse: collapse; width: 100%; font-size: 0.85rem;
        margin-bottom: 24px;
    }}
    th, td {{
        padding: 8px 10px; text-align: left;
        border-bottom: 1px solid #e8e8e8;
    }}
    th {{ background: #f5f7fa; font-weight: 600; position: sticky; top: 0; }}
    tr:hover {{ background: #f9fbfd; }}
    code {{ font-size: 0.82rem; background: #f0f0f0; padding: 1px 5px; border-radius: 3px; }}
    .orient {{
        display: inline-block; padding: 2px 8px; border-radius: 4px;
        font-size: 0.78rem; font-weight: 600;
    }}
    .co-f {{ background: #e3f2fd; color: #1565c0; }}
    .co-r {{ background: #fce4ec; color: #c62828; }}
    .na {{ color: #bbb; }}
    .motif-detail td {{ font-size: 0.8rem; }}
    footer {{
        margin-top: 40px; padding-top: 16px;
        border-top: 1px solid #eee; font-size: 0.78rem; color: #999;
    }}
</style>
</head>
<body>
<h1>Promoter identification report</h1>
<p style="color:#666; margin-bottom:4px;">Genome: <strong>{genome_name}</strong></p>

<div class="summary">
    <div class="stat"><div class="label">Promoters</div><div class="value">{n_promoters}</div></div>
    <div class="stat"><div class="label">CDS annotated</div><div class="value">{n_annotated}</div></div>
    <div class="stat"><div class="label">With motif hits</div><div class="value">{n_with_motifs}</div></div>
    <div class="stat"><div class="label">Total motif hits</div><div class="value">{total_fimo}</div></div>
</div>

<h2>Promoter summary</h2>
"""

_REPORT_HTML_FOOT = """
<footer>
    Generated by ProFinder.
</footer>
</body>
</html>
"""

STEPS = [
    (1,  "Run Prokka",                      step01_run_prokka),
    (2,  "Extract intergenic regions",       step02_extract_igrs),
    (3,  "Identify operons",                 step03_identify_operons),
    (4,  "Run hmmsearch",                    step04_run_hmmsearch),
    (5,  "Filter HMM output",               step05_filter_hmm),
    (6,  "Filter operons + add markers",     step06_filter_operons_add_markers),
    (7,  "Match IGRs to marker operons",     step07_match_igrs_to_markers),
    (8,  "Extract marker promoters",         step08_extract_marker_promoters),
    (9,  "Extract all promoters",            step09_extract_all_promoters),
    (10, "Annotate CDS (BLAST)",             step10_annotate_cds),
    (11, "Scan motifs (FIMO)",               step11_run_fimo),
    (12, "Generate HTML report",             step12_generate_report),
]


def main():
    parser = argparse.ArgumentParser(
        description="ProFinder — bacterial promoter identification pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("-i", "--input", type=Path, required="--list" not in sys.argv,
                        help="Input genome FASTA file")
    parser.add_argument("-o", "--output", type=Path, default=Path("output"),
                        help="Output directory (default: output/)")
    parser.add_argument("--start", type=int, default=1,
                        help="First step to run (default: 1)")
    parser.add_argument("--end", type=int, default=max(s[0] for s in STEPS),
                        help="Last step to run (default: last)")
    parser.add_argument("--list", action="store_true",
                        help="List all steps and exit")
    parser.add_argument("--force", action="store_true",
                        help="Re-run steps even if output already exists")

    # Tool paths
    parser.add_argument("--prokka", default="prokka",
                        help="Path to prokka binary (default: prokka)")
    parser.add_argument("--hmmsearch", default="hmmsearch",
                        help="Path to hmmsearch binary (default: hmmsearch)")
    parser.add_argument("--tigrfam", type=Path, default=None,
                        help="Path to TIGRfam HMM database (default: bundled)")
    parser.add_argument("--pfam", type=Path, default=None,
                        help="Path to Pfam-A HMM database (default: bundled)")

    # BLAST
    parser.add_argument("--blastp", default="blastp",
                        help="Path to blastp binary (default: blastp)")
    parser.add_argument("--blast-db", default="swissprot",
                        help="NCBI remote BLAST database name (default: swissprot)")
    parser.add_argument("--blast-local-db", type=Path, default=None,
                        help="Path to local BLAST database (overrides --blast-db remote)")
    parser.add_argument("--blast-evalue", type=float, default=1e-5,
                        help="BLAST e-value threshold (default: 1e-5)")

    # FIMO / motifs
    parser.add_argument("--fimo", default="fimo",
                        help="Path to fimo binary (default: fimo)")
    parser.add_argument("--motifs-dir", type=Path, default=None,
                        help="Directory containing .meme motif files "
                             "(default: bundled motifs/)")
    parser.add_argument("--fimo-threshold", type=float, default=1e-4,
                        help="FIMO p-value threshold (default: 1e-4)")

    # Conda environments
    parser.add_argument("--conda-prokka", default="",
                        help="Conda env for Prokka (blank = use current env)")
    parser.add_argument("--conda-hmm", default="",
                        help="Conda env for hmmsearch (blank = use current env)")
    parser.add_argument("--conda-blast", default="",
                        help="Conda env for BLAST+ (blank = use current env)")
    parser.add_argument("--conda-meme", default="",
                        help="Conda env for MEME Suite (blank = use current env)")

    # Parameters
    parser.add_argument("--threads", type=int, default=4,
                        help="Threads for external tools (default: 4)")
    parser.add_argument("--kingdom", default="Bacteria",
                        help="Prokka --kingdom (default: Bacteria)")
    parser.add_argument("--prefix", default="genome",
                        help="Prokka --prefix (default: genome)")
    parser.add_argument("--igr-min", type=int, default=75,
                        help="Minimum IGR length (default: 75)")
    parser.add_argument("--igr-max", type=int, default=1000,
                        help="Maximum IGR length (default: 1000)")
    parser.add_argument("--max-internal-dist", type=int, default=25,
                        help="Max distance between genes in an operon (default: 25)")
    parser.add_argument("--min-flanking-dist", type=int, default=75,
                        help="Min flanking distance for operon boundaries (default: 75)")
    parser.add_argument("--hmm-bitscore", type=float, default=25.0,
                        help="Minimum HMM bitscore (default: 25.0)")

    args = parser.parse_args()

    if args.list:
        for num, name, _ in STEPS:
            print(f"  {num:2d}. {name}")
        sys.exit(0)

    if not args.input.exists():
        sys.exit(f"Input file not found: {args.input}")

    cfg = Config(
        input_fasta=args.input.resolve(),
        output_dir=args.output.resolve(),
        prokka_bin=args.prokka,
        hmmsearch_bin=args.hmmsearch,
        blastp_bin=args.blastp,
        fimo_bin=args.fimo,
        tigrfam_hmm=args.tigrfam,
        pfam_hmm=args.pfam,
        conda_env_prokka=args.conda_prokka,
        conda_env_hmm=args.conda_hmm,
        conda_env_blast=args.conda_blast,
        conda_env_meme=args.conda_meme,
        threads=args.threads,
        prokka_kingdom=args.kingdom,
        prokka_prefix=args.prefix,
        igr_size_min=args.igr_min,
        igr_size_max=args.igr_max,
        max_internal_distance=args.max_internal_dist,
        min_flanking_distance=args.min_flanking_dist,
        hmm_bitscore_min=args.hmm_bitscore,
        blast_db=args.blast_db,
        blast_evalue=args.blast_evalue,
        blast_local_db=args.blast_local_db,
        motifs_dir=args.motifs_dir,
        fimo_threshold=args.fimo_threshold,
    )
    cfg.ensure_dirs()

    if not args.force:
        print("Checkpoint mode: steps with existing output will be skipped.")
        print("Use --force to re-run all steps.\n")

    for num, name, func in STEPS:
        if args.start <= num <= args.end:
            print(f"\n{'=' * 60}")
            print(f"  STEP {num}: {name}")
            print(f"{'=' * 60}\n")
            func(cfg, force=args.force)

    print("\nPipeline complete.")


if __name__ == "__main__":
    main()
