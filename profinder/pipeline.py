#!/usr/bin/env python3
"""
ProFinder — bacterial and archaeal promoter identification pipeline.

Takes a single genome FASTA, annotates it with Prokka, extracts
intergenic regions, identifies operons, screens for HMM marker genes,
and outputs promoter sequences in 5'-to-3' orientation.

Use ``--domain bacteria`` (default) for σ70 promoter prediction via
PromoterLCNN, or ``--domain archaea`` for binary promoter prediction
via iProm-Archaea.

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
    10. Predict promoters (PromoterLCNN for bacteria, iProm-Archaea for archaea)
    11. Annotate associated CDS (Prokka product names, gene, locus tag)
    12. Scan promoters for motifs (FIMO)
    13. Build final output table (all promoters, all metadata)
    14. Generate HTML report
"""

import argparse
import html as html_mod
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd
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
    """Extract a gene identifier from GFF attribute string.

    Tries, in order: ID=, locus_tag=, gene=, Name= (stripping any
    trailing ' gene' or ' CDS' suffix from Geneious-style names).
    """
    for key in ("ID=", "locus_tag=", "gene=", "Name="):
        for attr in attributes.split(";"):
            if attr.startswith(key):
                val = attr[len(key):]
                # Geneious appends " gene" or " CDS" to Name values
                for suffix in (" gene", " CDS"):
                    if val.endswith(suffix):
                        val = val[:-len(suffix)]
                return val
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
    print(f"  GFF:   {cfg.gff_file}")
    print(f"  FASTA: {fasta_source}")
    igr_df = extract_igrs(cfg.gff_file, fasta_source,
                          size_min=cfg.igr_size_min, size_max=cfg.igr_size_max)

    # Warn if most sequences are empty (contig ID mismatch)
    if not igr_df.empty:
        n_empty = igr_df["sequence"].isna().sum() + (igr_df["sequence"] == "").sum()
        if n_empty > 0:
            print(f"  WARNING: {n_empty}/{len(igr_df)} IGRs have empty sequences.")
            print("    This usually means contig IDs in the GFF don't match "
                  "the FASTA headers.")
            print("    If you supplied an FNA file in the batch table, make "
                  "sure it contains contig sequences (not per-CDS nucleotides).")

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
    """Run hmmsearch for every individual .hmm profile in the profiles directory."""
    if cfg.hmm_profiles_dir is None or not cfg.hmm_profiles_dir.is_dir():
        print("── HMM profiles directory not available ──")
        print(f"  hmm_profiles_dir: {cfg.hmm_profiles_dir}")
        print("  Provide --hmm-dir or reinstall to restore bundled profiles.")
        print("  Step 4 complete.\n")
        return

    hmm_files = sorted(cfg.hmm_profiles_dir.glob("*.hmm"))
    if not hmm_files:
        print(f"── No .hmm files found in {cfg.hmm_profiles_dir} ──")
        print("  Step 4 complete.\n")
        return

    print(f"  Using {len(hmm_files)} HMM profiles from: {cfg.hmm_profiles_dir}")

    tblout_dir = cfg.hmm_dir / "tblout"
    tblout_dir.mkdir(parents=True, exist_ok=True)

    # Check if all tblout files already exist
    if not force and all(
        (tblout_dir / (f.stem + ".tblout")).exists() for f in hmm_files
    ):
        print("── HMM output already exists for all profiles, skipping ──")
        print("  Step 4 complete.\n")
        return

    print(f"── Running hmmsearch ({len(hmm_files)} profiles) ──")
    cfg.hmm_dir.mkdir(parents=True, exist_ok=True)

    failed = 0
    for i, hmm_path in enumerate(hmm_files, 1):
        out_file = tblout_dir / (hmm_path.stem + ".tblout")
        log_file = tblout_dir / (hmm_path.stem + ".log")

        if not force and out_file.exists():
            continue

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
            f" {hmm_path} {cfg.faa_file}"
        )
        cmd = "set -euo pipefail; " + "".join(parts)
        result = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True)
        if result.returncode != 0:
            print(f"  WARNING: hmmsearch failed for {hmm_path.name}: "
                  f"{result.stderr[-200:]}")
            failed += 1

        if i % 100 == 0 or i == len(hmm_files):
            print(f"  {i}/{len(hmm_files)} profiles searched...")

    if failed:
        print(f"  {failed} profile(s) failed during hmmsearch.")

    print("  Step 4 complete.\n")


def step05_filter_hmm(cfg: Config, force: bool = False):
    """Consolidate and filter HMM hits from individual profile searches."""
    if not force and cfg.hmm_filtered.exists():
        print("── Filtered HMM file already exists, skipping ──")
        print("  Step 5 complete.\n")
        return

    tblout_dir = cfg.hmm_dir / "tblout"
    tblout_files = sorted(tblout_dir.glob("*.tblout")) if tblout_dir.is_dir() else []

    if not tblout_files:
        print("── No HMM output to filter, skipping ──")
        print("  Step 5 complete.\n")
        return

    print(f"── Filtering HMM output ({len(tblout_files)} tblout files) ──")
    all_rows = []
    for f in tblout_files:
        all_rows.extend(_parse_hmm_tblout(str(f)))

    if not all_rows:
        print("  No HMM hits found.")
        pd.DataFrame(columns=_HMM_HEADER).to_csv(cfg.hmm_combined, sep="\t", index=False)
        pd.DataFrame(columns=_HMM_HEADER).to_csv(cfg.hmm_filtered, sep="\t", index=False)
        print("  Step 5 complete.\n")
        return

    df = pd.DataFrame(all_rows, columns=_HMM_HEADER)
    df.to_csv(cfg.hmm_combined, sep="\t", index=False)
    print(f"  Combined HMM hits: {len(df)}")

    df["full_sequence_bitscore"] = pd.to_numeric(df["full_sequence_bitscore"], errors="coerce")
    df = df[df["full_sequence_bitscore"] >= cfg.hmm_bitscore_min]

    # Keep ALL matches — a gene may match multiple profiles and all are retained.
    df.to_csv(cfg.hmm_filtered, sep="\t", index=False)
    n_genes = df["target_name"].nunique() if not df.empty else 0
    print(f"  Filtered HMM ({len(df)} hits across {n_genes} genes) -> {cfg.hmm_filtered}")
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

        if len(merged) == 0 and not hmm.empty and not filtered.empty:
            operon_genes = set(filtered["Gene"].dropna())
            hmm_genes = set(hmm["target_name"].dropna())
            overlap = operon_genes & hmm_genes
            print()
            print("  WARNING: 0 marker operons found despite HMM hits.")
            print(f"    Operon gene IDs (first 5): "
                  f"{sorted(operon_genes)[:5]}")
            print(f"    HMM target names (first 5): "
                  f"{sorted(hmm_genes)[:5]}")
            print(f"    Overlap: {len(overlap)} gene(s)")
            if not overlap:
                print("    The gene IDs in the GFF and the protein IDs in "
                      "the FAA appear to use different naming schemes.")
                print("    Check that both files come from the same "
                      "annotation run.")
    else:
        # No HMM data: use all filtered operons as-is
        print("  No HMM data available; using all filtered operons as markers")
        filtered["target_name"] = filtered["Gene"]
        filtered["accession2"] = "no_hmm"
        filtered.to_csv(cfg.operon_filtered_markers, sep="\t", index=False)

    print("  Step 6 complete.\n")


def step07_match_igrs_to_markers(cfg: Config, force: bool = False):
    """Match IGRs to marker operon genes."""
    if not force and cfg.promoter_markers.exists() and cfg.promoter_markers_hmm.exists():
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
        pd.DataFrame().to_csv(cfg.promoter_markers_hmm, sep="\t", index=False)
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

    # Build expanded variant with one row per IGR × HMM profile match.
    # Determine which gene column is the marker gene for each matched IGR.
    matched["marker_gene"] = matched.apply(
        lambda r: r["right_gene"] if r["marker_match"] == "CO_F" else r["left_gene"],
        axis=1,
    )

    # Load HMM filtered hits for the profile names
    hmm_profiles = pd.DataFrame()
    if cfg.hmm_filtered.exists():
        try:
            hmm_profiles = pd.read_csv(
                cfg.hmm_filtered, sep="\t",
                usecols=["target_name", "query_name"],
            )
        except (pd.errors.EmptyDataError, ValueError):
            pass

    if hmm_profiles.empty:
        matched["hmm_profile"] = "no_hmm"
        expanded = matched
    else:
        # Inner merge: only keep IGRs whose marker gene has a direct HMM
        # hit.  Step 6 keeps entire operons when any gene in the operon
        # matches, so some IGRs flank operon-member genes that have no HMM
        # hit themselves.  Those are excluded here.
        expanded = pd.merge(
            matched, hmm_profiles,
            how="inner",
            left_on="marker_gene", right_on="target_name",
        )
        expanded.rename(columns={"query_name": "hmm_profile"}, inplace=True)
        expanded.drop(columns=["target_name"], inplace=True, errors="ignore")

    expanded.drop(columns=["marker_gene"], inplace=True, errors="ignore")
    expanded.to_csv(cfg.promoter_markers_hmm, sep="\t", index=False)
    print(f"  Expanded markers with HMM profiles ({len(expanded)} rows) "
          f"-> {cfg.promoter_markers_hmm}")
    print("  Step 7 complete.\n")


def _reverse_complement(seq) -> str:
    if not isinstance(seq, str) or not seq:
        return ""
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


def _load_contigs(fasta_path):
    """Load contig sequences from a genome FASTA into a dict."""
    from Bio import SeqIO
    return {rec.id: str(rec.seq) for rec in SeqIO.parse(str(fasta_path), "fasta")}


def _extract_cds_start(row, contigs, n_bp):
    """Return the first *n_bp* nucleotides of the CDS downstream of a
    promoter, oriented 5'→3' on the coding strand.

    For CO_F the downstream CDS is the right gene, sitting on the +
    strand immediately after the IGR.  Its first *n_bp* nt run left-to-
    right starting at ``igr_end + 1``.

    For CO_R the downstream CDS is the left gene, sitting on the −
    strand immediately before the IGR.  Its coding sequence starts at
    ``igr_start - 1`` and runs right-to-left.  We extract the last
    *n_bp* nt of the forward strand up to that position and reverse-
    complement them.

    Returns an empty string if the contig is missing or there aren't
    enough nucleotides available.
    """
    contig_seq = contigs.get(row["contig"], "")
    if not contig_seq:
        return ""

    orient = row["orientation"]
    if orient == "CO_F":
        # CDS starts at igr_end + 1 (1-based), on the + strand.
        cds_start_0 = row["end"]            # 0-based start (= igr_end in 1-based)
        fragment = contig_seq[cds_start_0 : cds_start_0 + n_bp]
    elif orient == "CO_R":
        # CDS ends at igr_start - 1 (1-based) on the + strand, but the
        # gene is on the − strand, so its coding sequence begins at the
        # high coordinate and runs leftward.
        cds_end_0 = row["start"] - 1        # 0-based position (= igr_start - 1 in 1-based)
        start_0 = max(0, cds_end_0 - n_bp)
        fragment = _reverse_complement(contig_seq[start_0 : cds_end_0])
    else:
        return ""

    return fragment


def _add_cds_column(df, contigs, n_bp):
    """Add a ``cds_start_seq`` column with the first *n_bp* nt of CDS,
    and a ``sequence_5p_to_3p_cds`` column with promoter + CDS joined.
    """
    df["cds_start_seq"] = df.apply(
        lambda r: _extract_cds_start(r, contigs, n_bp), axis=1)
    df["sequence_5p_to_3p_cds"] = df["sequence_5p_to_3p"] + df["cds_start_seq"]
    return df


def _write_fasta_cds(df, output_path, short_header=False):
    """Write promoter+CDS FASTA from a DataFrame that has
    ``sequence_5p_to_3p_cds``.
    """
    with open(output_path, "w") as fh:
        for _, row in df.iterrows():
            seq = row["sequence_5p_to_3p_cds"]
            igr_id = row["igr_id"]
            orient = row["orientation"]
            if short_header:
                header = f">{igr_id}_{orient}"
            else:
                header = (f">{igr_id}_{row['contig']}_{orient}"
                          f"_{row['left_gene']}_{row['right_gene']}")
            fh.write(f"{header}\n{seq}\n")
    print(f"  FASTA+CDS ({len(df)} seqs, short={short_header}) -> {output_path}")


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

    # Drop rows with missing sequences
    df = df[df["sequence"].notna() & (df["sequence"] != "")].copy()

    # Orient 5'→3': CO_R sequences need reverse complement
    df["sequence_5p_to_3p"] = df.apply(
        lambda r: _reverse_complement(r["sequence"]) if r["orientation"] == "CO_R"
        else r["sequence"],
        axis=1,
    )

    _write_fasta(df, cfg.promoter_fasta, short_header=False)
    _write_fasta(df, cfg.promoter_fasta_short, short_header=True)

    if cfg.cds_bp > 0:
        contigs = _load_contigs(cfg.input_fasta)
        df = _add_cds_column(df, contigs, cfg.cds_bp)
        _write_fasta_cds(df, cfg.promoter_cds_fasta, short_header=False)
        _write_fasta_cds(df, cfg.promoter_cds_fasta_short, short_header=True)

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

    # Drop rows with missing sequences (can happen with degenerate IGRs)
    n_before = len(df)
    df = df[df["sequence"].notna() & (df["sequence"] != "")].copy()
    n_dropped_na = n_before - len(df)
    if n_dropped_na:
        print(f"  Dropped {n_dropped_na} IGRs with missing sequences")

    # Orient 5'→3'
    df["sequence_5p_to_3p"] = df.apply(
        lambda r: _reverse_complement(r["sequence"]) if r["orientation"] == "CO_R"
        else r["sequence"],
        axis=1,
    )

    _write_fasta(df, cfg.all_promoter_fasta, short_header=False)
    _write_fasta(df, cfg.all_promoter_fasta_short, short_header=True)

    if cfg.cds_bp > 0:
        contigs = _load_contigs(cfg.input_fasta)
        df = _add_cds_column(df, contigs, cfg.cds_bp)
        _write_fasta_cds(df, cfg.all_promoter_cds_fasta, short_header=False)
        _write_fasta_cds(df, cfg.all_promoter_cds_fasta_short, short_header=True)

    print("  Step 9 complete.\n")


# =====================================================================
#  Step 10 — Promoter prediction (final filter)
#
#  Bacteria  → PromoterLCNN  (two-stage: binary + σ-factor subtype)
#  Archaea   → iProm-Archaea (single-stage: binary only)
# =====================================================================

def _step10_bacteria(cfg: Config, all_igr, force: bool = False):
    """PromoterLCNN path — bacterial σ-factor classification."""

    # Check weights availability
    if cfg.lcnn_weights_dir is None or not cfg.lcnn_weights_dir.is_dir():
        print("── PromoterLCNN weights not available, skipping prediction ──")
        print("  Provide --lcnn-weights or reinstall to restore bundled weights.")
        if cfg.igr_summary.exists():
            pd.DataFrame(columns=["igr_id", "prediction", "sigma_type"]).to_csv(
                cfg.lcnn_predictions, sep="\t", index=False)
        if cfg.promoter_markers.exists():
            import shutil
            shutil.copy2(cfg.promoter_markers, cfg.promoter_markers_verified)
        return

    # Trim to 3'-terminal 81 nt (closest to TSS / CDS start)
    _LCNN_LEN = 81
    seqs_81 = []
    short_mask = []
    for seq in all_igr["sequence_5p_to_3p"]:
        if len(seq) >= _LCNN_LEN:
            seqs_81.append(seq[-_LCNN_LEN:])
            short_mask.append(False)
        else:
            seqs_81.append(None)
            short_mask.append(True)

    n_short = sum(short_mask)
    if n_short:
        print(f"  {n_short} sequence(s) shorter than {_LCNN_LEN} nt — "
              f"classified as NON_PROMOTER")

    from .promoter_lcnn import load_models, predict as lcnn_predict, PromoterType

    n_valid = len(seqs_81) - n_short
    print(f"── Running PromoterLCNN on {n_valid} sequences ──")
    models = load_models(
        cfg.lcnn_weights_dir / "IsPromoter_fold_5",
        cfg.lcnn_weights_dir / "PromotersOnly_fold_1",
    )

    valid_seqs = [s for s in seqs_81 if s is not None]
    preds = lcnn_predict(models, valid_seqs) if valid_seqs else []

    # Merge predictions back, filling short sequences as NON_PROMOTER
    full_preds = []
    pred_iter = iter(preds)
    for is_short in short_mask:
        if is_short:
            full_preds.append(PromoterType.NON_PROMOTER)
        else:
            full_preds.append(next(pred_iter))

    # Build predictions table
    pred_rows = []
    for i, (_, row) in enumerate(all_igr.iterrows()):
        pt = full_preds[i]
        pred_rows.append({
            "igr_id": row["igr_id"],
            "prediction": "PROMOTER" if pt != PromoterType.NON_PROMOTER else "NON_PROMOTER",
            "sigma_type": pt.name,
        })

    return pred_rows


def _step10_archaea(cfg: Config, all_igr, force: bool = False):
    """iProm-Archaea path — binary archaeal promoter classification."""

    # Check weights availability
    if cfg.ipromarchaea_weights is None or not cfg.ipromarchaea_weights.exists():
        print("── iProm-Archaea weights not available, skipping prediction ──")
        print("  Provide --ipromarchaea-weights or reinstall to restore bundled weights.")
        if cfg.igr_summary.exists():
            pd.DataFrame(columns=["igr_id", "prediction", "sigma_type"]).to_csv(
                cfg.lcnn_predictions, sep="\t", index=False)
        if cfg.promoter_markers.exists():
            import shutil
            shutil.copy2(cfg.promoter_markers, cfg.promoter_markers_verified)
        return

    # iProm-Archaea scans sequences in non-overlapping 100 bp windows,
    # classifying each window independently.  A sequence is called a
    # promoter if ANY window scores positive.  Windows shorter than
    # 80 bp are discarded, so sequences < 80 bp produce no windows and
    # are classified as non-promoters.  The windowing and any-positive
    # aggregation are handled inside arch_predict().
    from .promoter_ipromarchaea import load_model, predict as arch_predict

    # Sequences < 80 nt will yield no windows inside predict(), but we
    # count them here for the log message.
    _ARCHAEA_MIN = 80
    n_short = sum(1 for seq in all_igr["sequence_5p_to_3p"] if len(seq) < _ARCHAEA_MIN)
    if n_short:
        print(f"  {n_short} sequence(s) shorter than {_ARCHAEA_MIN} nt — "
              f"classified as NON_PROMOTER (no valid windows)")

    seqs = list(all_igr["sequence_5p_to_3p"])
    print(f"── Running iProm-Archaea on {len(seqs)} sequences "
          f"(sliding 100 bp windows) ──")
    model = load_model(cfg.ipromarchaea_weights)

    full_preds = arch_predict(model, seqs)

    # Build predictions table.  For archaea there are no sigma subtypes,
    # so sigma_type is "PROMOTER" or "NON_PROMOTER" to stay consistent
    # with the downstream column expectations.
    pred_rows = []
    for i, (_, row) in enumerate(all_igr.iterrows()):
        is_prom = full_preds[i]
        pred_rows.append({
            "igr_id": row["igr_id"],
            "prediction": "PROMOTER" if is_prom else "NON_PROMOTER",
            "sigma_type": "PROMOTER" if is_prom else "NON_PROMOTER",
        })

    return pred_rows


def step10_predict_promoters(cfg: Config, force: bool = False):
    """Run promoter prediction on ALL promoter-orientation IGRs.

    For bacteria, each sequence is trimmed to its 3'-terminal 81 nt and
    classified by PromoterLCNN (binary + σ-factor subtype).

    For archaea, each sequence is trimmed to its 3'-terminal 100 nt and
    classified by iProm-Archaea (binary promoter/non-promoter only).

    The step writes:

    * ``lcnn_predictions.tsv`` — per-sequence predictions for every
      promoter-orientation IGR.
    * ``promoter_markers_verified.tsv`` — subset of marker promoters
      confirmed by the CNN (used by the HTML report and downstream
      annotation steps).
    """
    if not force and cfg.lcnn_predictions.exists() and cfg.promoter_markers_verified.exists():
        classifier = "iProm-Archaea" if cfg.is_archaea else "PromoterLCNN"
        print(f"── {classifier} predictions already exist, skipping ──")
        print("  Step 10 complete.\n")
        return

    # Load ALL promoter-orientation IGRs
    try:
        all_igr = pd.read_csv(cfg.igr_summary, sep="\t")
        all_igr = all_igr[all_igr["orientation"].isin(["CO_F", "CO_R"])].copy()
    except (FileNotFoundError, pd.errors.EmptyDataError):
        all_igr = pd.DataFrame()

    if all_igr.empty:
        print("── No promoter-orientation IGRs to classify ──")
        pd.DataFrame(columns=["igr_id", "prediction", "sigma_type"]).to_csv(
            cfg.lcnn_predictions, sep="\t", index=False)
        pd.DataFrame().to_csv(cfg.promoter_markers_verified, sep="\t", index=False)
        print("  Step 10 complete.\n")
        return

    # Orient sequences 5'→3'
    all_igr["sequence_5p_to_3p"] = all_igr.apply(
        lambda r: _reverse_complement(r["sequence"]) if r["orientation"] == "CO_R"
        else r["sequence"],
        axis=1,
    )

    # Dispatch to domain-specific classifier
    if cfg.is_archaea:
        pred_rows = _step10_archaea(cfg, all_igr, force=force)
    else:
        pred_rows = _step10_bacteria(cfg, all_igr, force=force)

    if pred_rows is None:
        # Weights were missing; helper already wrote fallback files.
        print("  Step 10 complete.\n")
        return

    pred_df = pd.DataFrame(pred_rows)
    pred_df.to_csv(cfg.lcnn_predictions, sep="\t", index=False)

    n_promoter = sum(1 for r in pred_rows if r["prediction"] == "PROMOTER")
    n_non = len(pred_rows) - n_promoter
    print(f"  Predicted: {n_promoter} promoter(s), {n_non} non-promoter(s)")
    print(f"  Predictions -> {cfg.lcnn_predictions}")

    # Write verified markers: intersection of marker promoters and CNN-positive
    marker_igr_ids = set()
    try:
        markers_df = pd.read_csv(cfg.promoter_markers, sep="\t")
        markers_df = markers_df[markers_df["orientation"].isin(["CO_F", "CO_R"])]
        marker_igr_ids = set(markers_df["igr_id"])
    except (FileNotFoundError, pd.errors.EmptyDataError):
        pass

    cnn_positive = set(r["igr_id"] for r in pred_rows if r["prediction"] == "PROMOTER")
    verified_ids = marker_igr_ids & cnn_positive

    sigma_map = {r["igr_id"]: r["sigma_type"] for r in pred_rows}

    if marker_igr_ids:
        verified_df = markers_df[markers_df["igr_id"].isin(verified_ids)].copy()
        verified_df["sigma_type"] = verified_df["igr_id"].map(sigma_map)
    else:
        verified_df = pd.DataFrame()

    verified_df.to_csv(cfg.promoter_markers_verified, sep="\t", index=False)
    print(f"  Verified markers ({len(verified_df)} rows) -> {cfg.promoter_markers_verified}")
    print("  Step 10 complete.\n")


# =====================================================================
#  Step 11 — Annotate associated CDS from Prokka
# =====================================================================

def _parse_prokka_gff_annotations(gff_path, gene_ids=None):
    """Extract CDS annotations from a Prokka GFF.

    Parameters
    ----------
    gff_path : Path
        Prokka GFF3 file.
    gene_ids : set or None
        If provided, only return entries whose ID is in this set.
        If None, return all CDS entries.

    Returns
    -------
    dict
        Mapping gene_id -> {"product": str, "gene": str, "locus_tag": str}.
    """
    annotations = {}
    with open(str(gff_path)) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9 or cols[2] != "CDS":
                continue
            attrs = cols[8]
            product = "hypothetical protein"
            gene_name = ""
            locus_tag = ""
            raw_id = ""
            for attr in attrs.split(";"):
                if attr.startswith("ID=") and not raw_id:
                    raw_id = attr[3:]
                elif attr.startswith("product="):
                    product = attr[8:]
                elif attr.startswith("gene="):
                    gene_name = attr[5:]
                elif attr.startswith("locus_tag="):
                    locus_tag = attr[10:]
            gene_id = _extract_gene_id(attrs)
            if gene_ids is None or gene_id in gene_ids:
                annotations[gene_id] = {
                    "product": product,
                    "gene": gene_name,
                    "locus_tag": locus_tag,
                }
    return annotations


def step11_annotate_cds(cfg: Config, force: bool = False):
    """Extract CDS annotations (product, gene name, locus tag) from Prokka
    GFF for ALL promoter-orientation IGRs."""
    annotation_tsv = cfg.output_dir / "cds_annotations.tsv"
    if not force and annotation_tsv.exists():
        print("── CDS annotations already exist, skipping ──")
        print("  Step 11 complete.\n")
        return

    # Collect CDS gene IDs from ALL promoter-orientation IGRs.
    # For CO_F the associated CDS is right_gene; for CO_R it's left_gene.
    gene_ids = set()
    if cfg.igr_summary.exists():
        try:
            df = pd.read_csv(cfg.igr_summary, sep="\t")
            df = df[df["orientation"].isin(["CO_F", "CO_R"])]
            gene_ids.update(df.loc[df["orientation"] == "CO_F", "right_gene"].dropna())
            gene_ids.update(df.loc[df["orientation"] == "CO_R", "left_gene"].dropna())
        except (pd.errors.EmptyDataError, KeyError):
            pass

    if not gene_ids:
        print("── No promoter-associated CDS to annotate ──")
        pd.DataFrame(columns=["gene_id", "product", "gene", "locus_tag"]).to_csv(
            annotation_tsv, sep="\t", index=False)
        print("  Step 11 complete.\n")
        return

    print(f"── Annotating {len(gene_ids)} CDS from Prokka GFF ──")
    ann_map = _parse_prokka_gff_annotations(cfg.gff_file, gene_ids)

    rows = []
    for gid in sorted(gene_ids):
        info = ann_map.get(gid, {})
        rows.append({
            "gene_id": gid,
            "product": info.get("product", "hypothetical protein"),
            "gene": info.get("gene", ""),
            "locus_tag": info.get("locus_tag", ""),
        })
    result_df = pd.DataFrame(rows)
    result_df.to_csv(annotation_tsv, sep="\t", index=False)
    print(f"  Annotated {len(result_df)} CDS -> {annotation_tsv}")
    print("  Step 11 complete.\n")


# =====================================================================
#  Step 12 — FIMO motif scanning
# =====================================================================

def step12_run_fimo(cfg: Config, force: bool = False):
    """Scan promoter sequences against motif databases with FIMO."""
    if not force and cfg.fimo_combined.exists():
        print("── FIMO results already exist, skipping ──")
        print("  Step 12 complete.\n")
        return

    # Determine which promoter FASTA to scan
    fasta_to_scan = None
    for candidate in [cfg.all_promoter_fasta, cfg.promoter_fasta]:
        if candidate.exists() and candidate.stat().st_size > 0:
            fasta_to_scan = candidate
            break
    if fasta_to_scan is None:
        print("── No promoter FASTA files found, skipping ──")
        print("  Step 12 complete.\n")
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
            print("  Step 12 complete.\n")
            return

    meme_files = sorted(motifs_dir.glob("*.meme"))
    if not meme_files:
        print(f"── No .meme files in {motifs_dir}, skipping ──")
        print("  Step 12 complete.\n")
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

    print("  Step 12 complete.\n")


# =====================================================================
#  Step 14 — HTML report
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


def step14_generate_report(cfg: Config, force: bool = False):
    """Generate an HTML report combining promoters, CDS annotations,
    and motif hits with visual sequence diagrams."""
    if not force and cfg.report_html.exists():
        print("── HTML report already exists, skipping ──")
        print("  Step 14 complete.\n")
        return

    print("── Generating HTML report ──")

    # Load verified promoters (after CNN filter).
    # Falls back to promoter_markers if verified file doesn't exist.
    source_file = cfg.promoter_markers_verified if cfg.promoter_markers_verified.exists() else cfg.promoter_markers
    promoter_df = None
    if source_file.exists():
        try:
            df = pd.read_csv(source_file, sep="\t")
            df = df[df["orientation"].isin(["CO_F", "CO_R"])].copy()
            if not df.empty:
                promoter_df = df
        except pd.errors.EmptyDataError:
            pass

    if promoter_df is None or promoter_df.empty:
        print("  No verified promoter data available for report.")
        print("  Step 14 complete.\n")
        return

    # For bacteria, filter to σ70 promoters only.
    # For archaea, keep all CNN-confirmed promoters (no sigma subtypes).
    if not cfg.is_archaea and "sigma_type" in promoter_df.columns:
        promoter_df = promoter_df[promoter_df["sigma_type"] == "SIGMA_70"].copy()
    if promoter_df.empty:
        label = "promoters" if cfg.is_archaea else "σ70 promoters"
        print(f"  No {label} found for report.")
        print("  Step 14 complete.\n")
        return

    # Count all promoter-orientation IGRs for reference
    n_all_promoters = 0
    if cfg.igr_summary.exists():
        try:
            all_igr = pd.read_csv(cfg.igr_summary, sep="\t")
            n_all_promoters = len(all_igr[all_igr["orientation"].isin(["CO_F", "CO_R"])])
        except pd.errors.EmptyDataError:
            pass

    # Add associated gene column
    promoter_df["associated_gene"] = promoter_df.apply(_get_associated_gene, axis=1)

    # Orient sequences 5'→3'
    promoter_df["sequence_5p_to_3p"] = promoter_df.apply(
        lambda r: _reverse_complement(r["sequence"]) if r["orientation"] == "CO_R"
        else r["sequence"],
        axis=1,
    )

    # Load CDS annotations (Prokka product names)
    annotation_map = {}
    annotation_tsv = cfg.output_dir / "cds_annotations.tsv"
    if annotation_tsv.exists():
        try:
            ann_df = pd.read_csv(annotation_tsv, sep="\t")
            for _, row in ann_df.iterrows():
                annotation_map[row["gene_id"]] = str(row.get("product", ""))
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
                    m = re.match(r"(igr_\d+)", seq_name)
                    igr_id = m.group(1) if m else seq_name.split("_")[0] + "_" + seq_name.split("_")[1]
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

    # Filter FIMO hits to only marker-filtered promoter IDs
    marker_igr_ids = set(promoter_df["igr_id"])
    fimo_by_promoter = {k: v for k, v in fimo_by_promoter.items()
                        if k in marker_igr_ids}

    # Build HTML
    n_promoters = len(promoter_df)
    n_annotated = sum(1 for _, r in promoter_df.iterrows()
                      if r["associated_gene"] in annotation_map
                      and annotation_map[r["associated_gene"]] != "hypothetical protein")
    n_with_motifs = sum(1 for _, r in promoter_df.iterrows()
                        if r["igr_id"] in fimo_by_promoter)

    genome_name = cfg.input_fasta.stem

    # Paths for file links (relative to output dir)
    all_promo_name = cfg.all_promoter_fasta.name if cfg.all_promoter_fasta.exists() else ""
    verified_tsv_name = cfg.promoter_markers_verified.name if cfg.promoter_markers_verified.exists() else ""
    promoters_fasta_name = cfg.promoter_fasta.name if cfg.promoter_fasta.exists() else ""
    final_table_name = cfg.final_table.name if cfg.final_table.exists() else ""

    # Domain-aware labels
    promoter_type_label = "promoters" if cfg.is_archaea else "σ70 promoters"
    promoter_type_html = "promoters" if cfg.is_archaea else "&sigma;70 promoters"
    classifier_name = "iProm-Archaea" if cfg.is_archaea else "PromoterLCNN"
    domain_label = "Archaea" if cfg.is_archaea else "Bacteria"

    html_parts = [_REPORT_HTML_HEAD.format(
        genome_name=html_mod.escape(genome_name),
        n_marker_promoters=n_promoters,
        n_all_promoters=n_all_promoters,
        n_annotated=n_annotated,
        n_with_motifs=n_with_motifs,
        total_fimo=sum(len(v) for v in fimo_by_promoter.values()),
        all_promoters_fasta=html_mod.escape(all_promo_name),
        verified_tsv=html_mod.escape(verified_tsv_name),
        promoters_fasta=html_mod.escape(promoters_fasta_name),
        final_table=html_mod.escape(final_table_name),
        promoter_type_html=promoter_type_html,
        promoter_type_label=promoter_type_label,
        classifier_name=classifier_name,
        domain_label=domain_label,
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
    sigma_col_header = "Classification" if cfg.is_archaea else "Sigma type"
    html_parts.append(f"""
    <table>
    <thead>
    <tr>
        <th>Promoter ID</th>
        <th>Contig</th>
        <th>Position</th>
        <th>Length</th>
        <th>Orientation</th>
        <th>{sigma_col_header}</th>
        <th>Associated CDS</th>
        <th>Protein</th>
        <th>Sequence diagram</th>
        <th>Sequence (5'→3')</th>
    </tr>
    </thead>
    <tbody>
    """)

    for _, row in promoter_df.iterrows():
        igr_id = row["igr_id"]
        gene = row["associated_gene"]
        product = annotation_map.get(gene, "")
        motif_hits = fimo_by_promoter.get(igr_id, [])

        svg = _build_motif_diagram_svg(row["length"], motif_hits)

        orient_class = "co-f" if row["orientation"] == "CO_F" else "co-r"
        seq = row.get("sequence_5p_to_3p", "")
        sigma = str(row.get("sigma_type", "")) if pd.notna(row.get("sigma_type")) else ""
        # Format sigma type for display (e.g. SIGMA_70 -> σ70)
        sigma_display = sigma.replace("SIGMA_", "σ") if sigma.startswith("SIGMA_") else sigma

        product_cell = html_mod.escape(product) if product else '<span class="na">—</span>'

        html_parts.append(f"""
        <tr>
            <td><code>{html_mod.escape(igr_id)}</code></td>
            <td>{html_mod.escape(str(row['contig']))}</td>
            <td>{row['start']:,}–{row['end']:,}</td>
            <td>{row['length']}</td>
            <td><span class="orient {orient_class}">{row['orientation']}</span></td>
            <td>{html_mod.escape(sigma_display) if sigma_display else '<span class="na">—</span>'}</td>
            <td><code>{html_mod.escape(str(gene))}</code></td>
            <td>{product_cell}</td>
            <td>{svg}</td>
            <td><code class="seq">{html_mod.escape(str(seq))}</code></td>
        </tr>
        """)

    html_parts.append("</tbody></table>")

    html_parts.append(_REPORT_HTML_FOOT)

    with open(cfg.report_html, "w") as fh:
        fh.write("\n".join(html_parts))
    print(f"  Report -> {cfg.report_html}")
    print("  Step 14 complete.\n")


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
    code.seq {{
        display: inline-block; max-width: 280px; overflow-x: auto;
        white-space: nowrap; font-size: 0.72rem; letter-spacing: 0.5px;
    }}
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
    <div class="stat"><div class="label">{promoter_type_html}</div><div class="value">{n_marker_promoters}</div></div>
    <div class="stat"><div class="label">CDS annotated</div><div class="value">{n_annotated}</div></div>
    <div class="stat"><div class="label">With motif hits</div><div class="value">{n_with_motifs}</div></div>
    <div class="stat"><div class="label">Total motif hits</div><div class="value">{total_fimo}</div></div>
    <div class="stat"><div class="label">All promoter IGRs</div><div class="value">{n_all_promoters}</div></div>
</div>

<p style="font-size:0.85rem; color:#555; margin-bottom:20px;">
    Domain: <strong>{domain_label}</strong> · Classifier: <strong>{classifier_name}</strong><br/>
    This report shows the <strong>{n_marker_promoters}</strong> {promoter_type_label} (marker-filtered, CNN-confirmed).
    A total of {n_all_promoters} promoter-orientation IGRs were identified in the genome.
</p>
<p style="font-size:0.85rem; color:#555; margin-bottom:20px;">
    <strong>Downloads:</strong>
    <a href="{promoters_fasta}">promoters.fasta</a> ·
    <a href="{verified_tsv}">promoter_markers_verified.tsv</a> ·
    <a href="{all_promoters_fasta}">all_promoters.fasta</a> ·
    <a href="{final_table}">profinder_results.tsv</a>
</p>

<h2>{promoter_type_html}</h2>
"""

_REPORT_HTML_FOOT = """
<footer>
    Generated by ProFinder.
</footer>
</body>
</html>
"""

# =====================================================================
#  Step 13 — Final output table
# =====================================================================

def step13_final_table(cfg: Config, force: bool = False):
    """Build a comprehensive TSV table covering ALL promoter-orientation
    IGRs with every piece of metadata collected by the pipeline.

    Columns
    -------
    promoter_id           IGR identifier (igr_NNNNNN)
    contig                Source contig / scaffold
    start                 IGR start coordinate (1-based)
    end                   IGR end coordinate
    length                IGR length in bp
    orientation           CO_F or CO_R
    associated_cds        Prokka gene ID of downstream CDS
    gene_name             Gene name from Prokka GFF (gene= attribute)
    locus_tag             Locus tag from Prokka GFF (locus_tag= attribute)
    product               Protein product name from Prokka GFF
    is_marker             Whether this IGR flanks a marker-operon gene
    lcnn_prediction       PROMOTER or NON_PROMOTER
    sigma_type            Sigma-factor subtype from PromoterLCNN
    motif_hits            Semicolon-separated list of FIMO motif hits
    sequence_5p_to_3p     Full promoter sequence oriented 5'→3'
    """
    if not force and cfg.final_table.exists():
        print("── Final output table already exists, skipping ──")
        print("  Step 13 complete.\n")
        return

    print("── Building final output table ──")

    # 1. Start from ALL promoter-orientation IGRs
    try:
        igr = pd.read_csv(cfg.igr_summary, sep="\t")
        igr = igr[igr["orientation"].isin(["CO_F", "CO_R"])].copy()
    except (FileNotFoundError, pd.errors.EmptyDataError):
        igr = pd.DataFrame()

    if igr.empty:
        print("  No promoter-orientation IGRs found.")
        pd.DataFrame().to_csv(cfg.final_table, sep="\t", index=False)
        print("  Step 13 complete.\n")
        return

    # 2. Associated CDS
    igr["associated_cds"] = igr.apply(_get_associated_gene, axis=1)

    # 3. Orient sequences 5'→3'
    igr["sequence_5p_to_3p"] = igr.apply(
        lambda r: _reverse_complement(r["sequence"]) if r["orientation"] == "CO_R"
        else r["sequence"],
        axis=1,
    )

    # 4. CDS annotations (gene name, locus tag, product)
    annotation_tsv = cfg.output_dir / "cds_annotations.tsv"
    ann_map = {}   # gene_id -> {product, gene, locus_tag}
    if annotation_tsv.exists():
        try:
            ann_df = pd.read_csv(annotation_tsv, sep="\t")
            for _, row in ann_df.iterrows():
                ann_map[row["gene_id"]] = {
                    "product": str(row.get("product", "")),
                    "gene": str(row.get("gene", "")),
                    "locus_tag": str(row.get("locus_tag", "")),
                }
        except (pd.errors.EmptyDataError, KeyError):
            pass

    igr["gene_name"] = igr["associated_cds"].map(
        lambda g: ann_map.get(g, {}).get("gene", ""))
    igr["locus_tag"] = igr["associated_cds"].map(
        lambda g: ann_map.get(g, {}).get("locus_tag", ""))
    igr["product"] = igr["associated_cds"].map(
        lambda g: ann_map.get(g, {}).get("product", ""))

    # 5. Marker status
    marker_igr_ids = set()
    if cfg.promoter_markers.exists():
        try:
            m_df = pd.read_csv(cfg.promoter_markers, sep="\t")
            marker_igr_ids = set(m_df["igr_id"])
        except (pd.errors.EmptyDataError, KeyError):
            pass
    igr["is_marker"] = igr["igr_id"].isin(marker_igr_ids)

    # 6. PromoterLCNN predictions
    lcnn_map = {}   # igr_id -> {prediction, sigma_type}
    if cfg.lcnn_predictions.exists():
        try:
            lcnn_df = pd.read_csv(cfg.lcnn_predictions, sep="\t")
            for _, row in lcnn_df.iterrows():
                lcnn_map[row["igr_id"]] = {
                    "prediction": str(row.get("prediction", "")),
                    "sigma_type": str(row.get("sigma_type", "")),
                }
        except (pd.errors.EmptyDataError, KeyError):
            pass

    igr["lcnn_prediction"] = igr["igr_id"].map(
        lambda g: lcnn_map.get(g, {}).get("prediction", ""))
    igr["sigma_type"] = igr["igr_id"].map(
        lambda g: lcnn_map.get(g, {}).get("sigma_type", ""))

    # 7. FIMO motif hits (semicolon-separated summary per promoter)
    fimo_summary = {}   # igr_id -> semicolon-separated motif names
    if cfg.fimo_combined.exists():
        try:
            fimo_df = pd.read_csv(cfg.fimo_combined, sep="\t")
            seq_col = "sequence_name"
            if seq_col not in fimo_df.columns:
                for alt in ["sequence name", "#pattern name"]:
                    if alt in fimo_df.columns:
                        seq_col = alt
                        break
            if seq_col in fimo_df.columns:
                for seq_name, grp in fimo_df.groupby(seq_col):
                    m = re.match(r"(igr_\d+)", seq_name)
                    igr_id = m.group(1) if m else seq_name.split("_")[0] + "_" + seq_name.split("_")[1]
                    names = []
                    for _, hit in grp.iterrows():
                        alt = str(hit.get("motif_alt_id", ""))
                        mid = str(hit.get("motif_id", hit.get("#pattern name", "")))
                        db = str(hit.get("motif_database", ""))
                        pval = hit.get("p-value", "")
                        label = alt if alt and alt != "nan" else mid
                        names.append(f"{label}({db},p={pval})")
                    fimo_summary[igr_id] = "; ".join(names)
        except (pd.errors.EmptyDataError, KeyError):
            pass

    igr["motif_hits"] = igr["igr_id"].map(lambda g: fimo_summary.get(g, ""))

    # 8. CDS-extended sequences (optional)
    columns = [
        "igr_id", "contig", "start", "end", "length", "orientation",
        "associated_cds", "gene_name", "locus_tag", "product",
        "is_marker", "lcnn_prediction", "sigma_type",
        "motif_hits", "sequence_5p_to_3p",
    ]

    if cfg.cds_bp > 0:
        contigs = _load_contigs(cfg.input_fasta)
        igr = _add_cds_column(igr, contigs, cfg.cds_bp)
        columns.append("sequence_5p_to_3p_cds")

    # 9. Select and order final columns
    out = igr[columns].copy()
    out.rename(columns={"igr_id": "promoter_id"}, inplace=True)

    out.to_csv(cfg.final_table, sep="\t", index=False)
    print(f"  Final table ({len(out)} rows) -> {cfg.final_table}")
    print("  Step 13 complete.\n")


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
    (10, "Predict promoters (CNN classifier)", step10_predict_promoters),
    (11, "Annotate CDS (Prokka)",            step11_annotate_cds),
    (12, "Scan motifs (FIMO)",               step12_run_fimo),
    (13, "Build final output table",         step13_final_table),
    (14, "Generate HTML report",             step14_generate_report),
]


def main():
    # Suppress TensorFlow C++ runtime warnings (CUDA probing, oneDNN,
    # cudart_stub, TensorRT, etc.) BEFORE any TF import can happen.
    # Level 3 = FATAL only; these must be set before the C++ runtime
    # initialises, which happens at the first ``import tensorflow``.
    import os
    os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
    os.environ["TF_ENABLE_ONEDNN_OPTS"] = "0"
    os.environ["CUDA_VISIBLE_DEVICES"] = os.environ.get("CUDA_VISIBLE_DEVICES", "")

    parser = argparse.ArgumentParser(
        description="ProFinder — bacterial and archaeal promoter identification pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    input_group = parser.add_mutually_exclusive_group(
        required="--list" not in sys.argv,
    )
    input_group.add_argument("-i", "--input", type=Path,
                             help="Input genome FASTA file (single-sample mode)")
    input_group.add_argument("--batch", type=Path,
                             help="TSV table for batch mode. Required columns: "
                                  "sample_id, fasta. Optional columns: gff, faa, "
                                  "fna. When gff and faa are provided, Prokka "
                                  "(step 1) is skipped for that sample.")
    parser.add_argument("-o", "--output", type=Path, default=Path("output"),
                        help="Output directory (default: output/). In batch mode "
                             "each sample gets a subdirectory.")
    parser.add_argument("--domain", choices=["bacteria", "archaea"],
                        default="bacteria",
                        help="Target domain: 'bacteria' uses PromoterLCNN "
                             "(σ-factor classification); 'archaea' uses "
                             "iProm-Archaea (binary promoter/non-promoter). "
                             "Default: bacteria")
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
    parser.add_argument("--hmm-dir", type=Path, default=None,
                        help="Directory containing individual .hmm profile files "
                             "(default: bundled profiles)")

    # PromoterLCNN (bacteria)
    parser.add_argument("--lcnn-weights", type=Path, default=None,
                        help="Directory containing PromoterLCNN SavedModel "
                             "subdirectories (default: bundled weights/)")

    # iProm-Archaea (archaea)
    parser.add_argument("--ipromarchaea-weights", type=Path, default=None,
                        help="Path to iProm-Archaea .h5 weights file "
                             "(default: bundled weights/iPromArchaea/)")

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
    parser.add_argument("--conda-meme", default="",
                        help="Conda env for MEME Suite (blank = use current env)")

    # Parameters
    parser.add_argument("--threads", type=int, default=4,
                        help="Threads for external tools (default: 4)")
    parser.add_argument("--kingdom", default="Bacteria",
                        help="Prokka --kingdom (default: Bacteria)")
    parser.add_argument("--prefix", default="genome",
                        help="Prokka --prefix (default: genome)")
    parser.add_argument("--igr-min", type=int, default=81,
                        help="Minimum IGR length (default: 81)")
    parser.add_argument("--igr-max", type=int, default=1000,
                        help="Maximum IGR length (default: 1000)")
    parser.add_argument("--max-internal-dist", type=int, default=25,
                        help="Max distance between genes in an operon (default: 25)")
    parser.add_argument("--min-flanking-dist", type=int, default=75,
                        help="Min flanking distance for operon boundaries (default: 75)")
    parser.add_argument("--hmm-bitscore", type=float, default=25.0,
                        help="Minimum HMM bitscore (default: 25.0)")
    parser.add_argument("--cds-bp", type=int, default=0,
                        help="Number of CDS-start nucleotides to append to "
                             "promoter sequences in additional FASTA outputs "
                             "and the final table. 0 = disabled (default: 0)")

    args = parser.parse_args()

    if args.list:
        for num, name, _ in STEPS:
            print(f"  {num:2d}. {name}")
        sys.exit(0)

    # When running in archaea mode, default Prokka kingdom to Archaea
    # unless the user explicitly provided --kingdom.
    kingdom = args.kingdom
    if args.domain == "archaea" and kingdom == "Bacteria":
        kingdom = "Archaea"

    # ── Batch mode ──────────────────────────────────────────────────
    if args.batch is not None:
        if not args.batch.exists():
            sys.exit(f"Batch file not found: {args.batch}")
        _run_batch(args, kingdom)
        return

    # ── Single-sample mode ──────────────────────────────────────────
    if not args.input.exists():
        sys.exit(f"Input file not found: {args.input}")

    cfg = _build_config(args, kingdom,
                        input_fasta=args.input.resolve(),
                        output_dir=args.output.resolve())
    _run_pipeline(cfg, args)

    print("\nPipeline complete.")


def _build_config(args, kingdom, *, input_fasta, output_dir, prokka_prefix=None):
    """Create a Config from CLI args, overriding paths as needed."""
    return Config(
        input_fasta=input_fasta,
        output_dir=output_dir,
        domain=args.domain,
        prokka_bin=args.prokka,
        hmmsearch_bin=args.hmmsearch,
        fimo_bin=args.fimo,
        hmm_profiles_dir=args.hmm_dir,
        lcnn_weights_dir=args.lcnn_weights,
        ipromarchaea_weights=args.ipromarchaea_weights,
        conda_env_prokka=args.conda_prokka,
        conda_env_hmm=args.conda_hmm,
        conda_env_meme=args.conda_meme,
        threads=args.threads,
        prokka_kingdom=kingdom,
        prokka_prefix=prokka_prefix or args.prefix,
        igr_size_min=args.igr_min,
        igr_size_max=args.igr_max,
        max_internal_distance=args.max_internal_dist,
        min_flanking_distance=args.min_flanking_dist,
        hmm_bitscore_min=args.hmm_bitscore,
        motifs_dir=args.motifs_dir,
        fimo_threshold=args.fimo_threshold,
        cds_bp=args.cds_bp,
    )


def _run_pipeline(cfg, args):
    """Run pipeline steps on a single Config."""
    cfg.ensure_dirs()

    if cfg.cds_bp > 0 and cfg.cds_bp % 3 != 0:
        print()
        print("!" * 60)
        print("  WARNING: --cds-bp %d is not divisible by 3." % cfg.cds_bp)
        print()
        print("  The CDS extension will be %d nt, which is not a whole" % cfg.cds_bp)
        print("  number of codons. If you concatenate these promoter+CDS")
        print("  sequences upstream of a coding sequence, the downstream")
        print("  reading frame will be shifted. Use a multiple of 3")
        print("  (e.g. 90, 150, 300) to keep the CDS fragment in-frame.")
        print("!" * 60)
        print()

    if not args.force:
        print("Checkpoint mode: steps with existing output will be skipped.")
        print("Use --force to re-run all steps.\n")

    for num, name, func in STEPS:
        if args.start <= num <= args.end:
            print(f"\n{'=' * 60}")
            print(f"  STEP {num}: {name}")
            print(f"{'=' * 60}\n")
            func(cfg, force=args.force)


def _link_prokka_files(cfg, row):
    """Symlink user-supplied Prokka files into the expected output layout.

    Returns the step number to start from: 2 if Prokka files were linked
    (skip step 1), or 1 if they were not provided.
    """
    gff_path = row.get("gff", "")
    faa_path = row.get("faa", "")

    if not gff_path or not faa_path:
        return 1  # no Prokka files supplied — run step 1 normally

    gff_src = Path(str(gff_path)).resolve()
    faa_src = Path(str(faa_path)).resolve()

    for label, src in [("gff", gff_src), ("faa", faa_src)]:
        if not src.exists():
            sys.exit(f"Batch table references missing {label} file: {src}")

    cfg.prokka_dir.mkdir(parents=True, exist_ok=True)

    # Derive the prefix from the GFF filename so downstream paths resolve.
    prefix = gff_src.stem  # e.g. "genome" from "genome.gff"
    cfg.prokka_prefix = prefix

    # Symlink GFF and FAA
    for src, suffix in [(gff_src, ".gff"), (faa_src, ".faa")]:
        dest = cfg.prokka_dir / f"{prefix}{suffix}"
        if dest.exists() or dest.is_symlink():
            dest.unlink()
        dest.symlink_to(src)

    # Optional: FNA (nucleotide FASTA from Prokka)
    fna_path = row.get("fna", "")
    if fna_path:
        fna_src = Path(str(fna_path)).resolve()
        if fna_src.exists():
            dest = cfg.prokka_dir / f"{prefix}.fna"
            if dest.exists() or dest.is_symlink():
                dest.unlink()
            dest.symlink_to(fna_src)

    return 2  # skip step 1


def _run_batch(args, kingdom):
    """Run the pipeline once per row in a batch TSV table."""
    batch_df = pd.read_csv(args.batch, sep="\t")

    required_cols = {"sample_id", "fasta"}
    missing = required_cols - set(batch_df.columns)
    if missing:
        sys.exit(f"Batch table is missing required columns: {', '.join(sorted(missing))}")

    n_samples = len(batch_df)
    print(f"Batch mode: {n_samples} sample(s) from {args.batch}\n")

    for idx, row in batch_df.iterrows():
        sample_id = str(row["sample_id"])
        fasta_path = Path(str(row["fasta"])).resolve()

        if not fasta_path.exists():
            sys.exit(f"FASTA not found for sample '{sample_id}': {fasta_path}")

        sample_output = args.output.resolve() / sample_id

        print("\n" + "#" * 60)
        print(f"  SAMPLE {idx + 1}/{n_samples}: {sample_id}")
        print("#" * 60 + "\n")

        cfg = _build_config(
            args, kingdom,
            input_fasta=fasta_path,
            output_dir=sample_output,
        )

        # If Prokka files are pre-supplied, symlink them and skip step 1.
        min_step = _link_prokka_files(cfg, row)
        effective_start = max(args.start, min_step)

        # Temporarily override start so _run_pipeline respects it.
        saved_start = args.start
        args.start = effective_start
        _run_pipeline(cfg, args)
        args.start = saved_start

        print(f"\n  Sample '{sample_id}' complete.\n")

    print(f"\nBatch complete. {n_samples} sample(s) processed.")


if __name__ == "__main__":
    main()
