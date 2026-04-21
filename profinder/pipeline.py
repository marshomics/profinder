#!/usr/bin/env python3
"""
ProFinder — bacterial and archaeal promoter identification pipeline.

Takes a single genome FASTA, annotates it with Prokka, extracts
intergenic regions, identifies operons, screens for HMM marker genes,
and outputs promoter sequences in 5'-to-3' orientation.

Use ``--domain bacteria`` (default) or ``--domain archaea``.
Promoter verification is performed by scanning for -10/-35 hexamer
motifs using position weight matrices from the bundled
``all_unique_subgroups.meme`` file. Each -10 hit is classified as
Path A (linked -10/-35 same subgroup, 15–19 bp spacer), Path B
(extended -10 with variable or absent -35), or Path C (unlinked
-10/-35, 15–19 bp spacer). Anything else is no hit.

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
    8.  Scan promoter motifs (-10/-35 hexamer verification, A/B/C)
    9.  Annotate associated CDS (Prokka product names, gene, locus tag)
    10. Build final output table (all promoters, all metadata)
    11. Generate HTML report
"""

import argparse
import html as html_mod
import math
import re
import subprocess
import sys
from itertools import product as iprod
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

from .config import Config
from .igr_extractor import extract_igrs


# Shared reverse-complement translation table.  Used by the motif-scan
# inner loop which builds two strands per sequence — recreating the
# translation table per row is measurable on large IGR sets.
_REVCOMP_TRANS = str.maketrans("ACGTacgt", "TGCAtgca")


def _revcomp_simple(s: str) -> str:
    """Fast reverse-complement using the module-level translation table."""
    return s.translate(_REVCOMP_TRANS)[::-1]


# =====================================================================
#  Motif scanning helpers (adapted from meme_scan.py)
# =====================================================================

_MOTIF_WIDTH = 6
_SPACER_MIN = 15
_SPACER_MAX = 19
_BASE_IDX = {"A": 0, "C": 1, "G": 2, "T": 3}


def _parse_meme_file(filepath):
    """Parse a MEME-format file. Returns {motif_name: [[pA, pC, pG, pT], ...]}.

    The motif name is the first token following ``MOTIF``. For
    ``all_unique_subgroups.meme`` these look like ``M001_m35``,
    ``M001_m10``, ``M002_m35`` etc — the prefix before the underscore
    identifies the subgroup, and the ``_m10`` / ``_m35`` suffix the
    element type.
    """
    motifs = {}
    current = None
    matrix = []
    in_matrix = False

    with open(filepath) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("MOTIF "):
                if current and matrix:
                    motifs[current] = matrix
                current = line.split()[1]
                matrix = []
                in_matrix = False
            elif line.startswith("letter-probability matrix"):
                in_matrix = True
            elif in_matrix:
                if not line:
                    in_matrix = False
                    continue
                vals = line.split()
                if len(vals) == 4:
                    try:
                        matrix.append([float(v) for v in vals])
                    except ValueError:
                        in_matrix = False
                else:
                    in_matrix = False

    if current and matrix:
        motifs[current] = matrix
    return motifs


def _freq_to_log_odds(freq_matrix, bg=0.25, pseudocount=1e-6):
    """Convert frequency matrix to log2-odds scoring matrix."""
    return [[math.log2(max(f, pseudocount) / bg) for f in row]
            for row in freq_matrix]


def _compute_score_threshold(log_odds_matrix, p_value):
    """Compute score threshold for a given p-value by enumerating all 4^w k-mers."""
    w = len(log_odds_matrix)
    scores = sorted(
        [sum(log_odds_matrix[i][b] for i, b in enumerate(kmer))
         for kmer in iprod(range(4), repeat=w)],
        reverse=True,
    )
    rank = max(0, min(int(math.ceil(p_value * len(scores))) - 1, len(scores) - 1))
    return scores[rank]


def _score_kmer(seq, pos, log_odds_matrix):
    """Score a 6-mer at the given position. Returns None if non-ACGT bases present."""
    s = 0.0
    for j in range(_MOTIF_WIDTH):
        idx = _BASE_IDX.get(seq[pos + j])
        if idx is None:
            return None
        s += log_odds_matrix[j][idx]
    return s


class _MotifSet:
    """A collection of subgroup motifs for one element type (-10 or -35),
    each with its own log-odds matrix and score threshold. Subgroup IDs
    (e.g. ``M001``) are used to test for "linked" hits where a -10 and
    a -35 come from the same subgroup."""

    def __init__(self):
        self.entries = []  # list of (subgroup, log_odds_matrix, threshold)

    def add(self, subgroup, freq_matrix, p_value):
        lom = _freq_to_log_odds(freq_matrix)
        thresh = _compute_score_threshold(lom, p_value)
        self.entries.append((subgroup, lom, thresh))

    def best_hit(self, seq, pos):
        """Return (score, subgroup) of the best-scoring motif that exceeds
        its threshold at this position, or None."""
        best = None
        for subgroup, lom, thresh in self.entries:
            s = _score_kmer(seq, pos, lom)
            if s is not None and s >= thresh:
                if best is None or s > best[0]:
                    best = (s, subgroup)
        return best


def _load_motifs_from_file(meme_path, p10, p35):
    """Load paired subgroup motifs from a single MEME file.

    Expects motifs named ``<subgroup>_m10`` and ``<subgroup>_m35`` —
    e.g. ``M001_m10`` / ``M001_m35``. Returns (m10, m35) MotifSets
    where each entry's source is the subgroup ID.
    """
    m10 = _MotifSet()
    m35 = _MotifSet()

    meme_path = Path(meme_path)
    if not meme_path.is_file():
        print(f"  WARNING: MEME file not found: {meme_path}", file=sys.stderr)
        return m10, m35

    motifs = _parse_meme_file(meme_path)
    if not motifs:
        print(f"  WARNING: no motifs parsed from {meme_path}", file=sys.stderr)
        return m10, m35

    for name, matrix in motifs.items():
        # Expected naming: "<subgroup>_m10" or "<subgroup>_m35"
        if name.endswith("_m10"):
            subgroup = name[:-4]
            m10.add(subgroup, matrix, p10)
        elif name.endswith("_m35"):
            subgroup = name[:-4]
            m35.add(subgroup, matrix, p35)

    print(f"  Loaded {len(m10.entries)} -10 and {len(m35.entries)} -35 "
          f"subgroup motifs from {meme_path.name}")
    return m10, m35


def _find_promoters_on_strand(seq, m10, m35):
    """Scan one strand for promoter candidates. Returns a list of result dicts.

    Each -10 hit is classified into one of three paths, in priority
    order:

        A  Linked -10 and -35 (same subgroup) with 15–19 bp spacer.
        B  Extended -10 (TG dinucleotide immediately upstream of the
           hexamer), regardless of any -35 hit.
        C  Unlinked -10 and -35 (different subgroups) with 15–19 bp
           spacer.

    Anything else is no hit and is dropped.
    """
    w = _MOTIF_WIDTH
    max_pos = len(seq) - w

    # Pre-scan all -10 hits
    hits_10 = []
    for i in range(max_pos + 1):
        hit = m10.best_hit(seq, i)
        if hit:
            hits_10.append((i, hit[0], hit[1]))

    if not hits_10:
        return []

    # Pre-scan all -35 hits
    hits_35_by_pos = {}
    for i in range(max_pos + 1):
        hit = m35.best_hit(seq, i)
        if hit:
            hits_35_by_pos[i] = hit

    results = []
    for pos_10, score_10, subgroup_10 in hits_10:
        has_ext10 = pos_10 >= 2 and seq[pos_10 - 2:pos_10] == "TG"

        # Collect -35 hits within the 15–19 bp spacer window. Split into
        # linked (same subgroup) vs unlinked (different subgroup).
        best_linked = None    # (score, subgroup, pos, spacer)
        best_unlinked = None
        for spacer in range(_SPACER_MIN, _SPACER_MAX + 1):
            p35 = pos_10 - w - spacer
            if p35 < 0:
                continue
            hit = hits_35_by_pos.get(p35)
            if hit is None:
                continue
            s35, sub35 = hit
            if sub35 == subgroup_10:
                if best_linked is None or s35 > best_linked[0]:
                    best_linked = (s35, sub35, p35, spacer)
            else:
                if best_unlinked is None or s35 > best_unlinked[0]:
                    best_unlinked = (s35, sub35, p35, spacer)

        # Classify this -10 hit: A > B > C, else drop.
        if best_linked is not None:
            path = "A"
            chosen_35 = best_linked
        elif has_ext10:
            path = "B"
            # For Path B the -35 is "variable or absent"; keep the best
            # -35 in the 15–19 bp window if there is one, preferring a
            # linked hit (there isn't one here by construction) then the
            # best unlinked hit.
            chosen_35 = best_unlinked
        elif best_unlinked is not None:
            path = "C"
            chosen_35 = best_unlinked
        else:
            continue

        r = {
            "pos_10": pos_10,
            "seq_10": seq[pos_10:pos_10 + w],
            "score_10": round(score_10, 3),
            "source_10": subgroup_10,
            "has_ext10": has_ext10,
            "path": path,
        }

        if chosen_35 is not None:
            s35, sub35, p35, spacer = chosen_35
            r["pos_35"] = p35
            r["seq_35"] = seq[p35:p35 + w]
            r["score_35"] = round(s35, 3)
            r["source_35"] = sub35
            r["spacer_len"] = spacer
            r["spacer_seq"] = seq[p35 + w:pos_10]
        else:
            r["pos_35"] = "."
            r["seq_35"] = "."
            r["score_35"] = "."
            r["source_35"] = "."
            r["spacer_len"] = "."
            r["spacer_seq"] = "."

        results.append(r)

    return results


_MOTIF_COLUMNS = [
    "strand", "pos_10", "seq_10", "score_10", "source_10",
    "has_ext10", "pos_35", "seq_35", "score_35", "source_35",
    "spacer_len", "spacer_seq", "path",
]

_PATH_RANK = {"A": 0, "B": 1, "C": 2}


def _scan_sequences_for_motifs(df, m10, m35,
                               seq_col="sequence_5p_to_3p",
                               id_col="igr_id"):
    """Scan a DataFrame of sequences for -10/-35 promoter motifs.

    Returns two DataFrames:
        all_hits  — every motif hit found (multiple rows per sequence possible)
        best_hits — one row per sequence, keeping the best hit
                    (Path A > B > C, then highest -10 score)
    """
    all_rows = []
    best_per_seq = {}  # igr_id -> (path_rank, neg_score_10, row_dict)

    for _, row in df.iterrows():
        raw = row[seq_col]
        if not isinstance(raw, str) or not raw:
            continue
        raw = raw.upper()
        igr_id = row[id_col]

        for strand, seq in [("+", raw), ("-", _revcomp_simple(raw))]:
            for r in _find_promoters_on_strand(seq, m10, m35):
                out = {id_col: igr_id, "strand": strand}
                for k in _MOTIF_COLUMNS:
                    if k != "strand":
                        out[k] = r.get(k, ".")
                all_rows.append(out)

                path_rank = _PATH_RANK.get(r["path"], 99)
                s10 = r["score_10"] if isinstance(r["score_10"], float) else 0.0
                prev = best_per_seq.get(igr_id)
                if prev is None or (path_rank, -s10) < (prev[0], prev[1]):
                    best_per_seq[igr_id] = (path_rank, -s10, out)

    all_hits = pd.DataFrame(all_rows) if all_rows else pd.DataFrame(columns=[id_col] + _MOTIF_COLUMNS)
    best_rows = [v[2] for v in best_per_seq.values()]
    best_hits = pd.DataFrame(best_rows) if best_rows else pd.DataFrame(columns=[id_col] + _MOTIF_COLUMNS)

    return all_hits, best_hits


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
    # but only if it actually contains contig-level sequences. Per-CDS
    # nucleotide files (.ffn/.fnn) use gene locus tags as headers rather
    # than contig accessions, so they won't match the GFF seqid column.
    fasta_source = cfg.input_fasta
    if cfg.fna_file.exists():
        # Read the contig IDs from the GFF (first column of data lines)
        gff_contigs = set()
        with open(str(cfg.gff_file)) as fh:
            for line in fh:
                if line.startswith("##FASTA"):
                    break
                if line.startswith("#"):
                    continue
                gff_contigs.add(line.split("\t", 1)[0])
        # Check if any GFF contig appears in the FNA headers
        fna_ids = {rec.id for rec in SeqIO.parse(str(cfg.fna_file), "fasta")}
        if gff_contigs & fna_ids:
            fasta_source = cfg.fna_file
        else:
            print(f"  NOTE: FNA file ({cfg.fna_file.name}) does not contain "
                  f"contig sequences matching the GFF — using genome FASTA instead.")

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

    # Add 5'→3' oriented sequence: reverse-complement CO_R, keep others as-is
    if not igr_df.empty:
        igr_df["sequence_5p_to_3p"] = igr_df.apply(
            lambda r: _reverse_complement(r["sequence"]) if r["orientation"] == "CO_R"
            else r["sequence"],
            axis=1,
        )

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


# =====================================================================
#  Step 8 — Motif-based promoter filtering
#
#  Scans promoter sequences for -10 and -35 hexamer motifs using PWMs
#  from the bundled all_unique_subgroups.meme file, which contains
#  paired M###_m10 / M###_m35 subgroups.
#
#  Each -10 hit is classified into one of three paths:
#    A — linked -10 and -35 (same subgroup) with 15–19 bp spacer
#    B — extended -10 (TG dinucleotide upstream) with variable or
#        absent -35
#    C — unlinked -10 and -35 (different subgroups) with 15–19 bp
#        spacer
#  Anything else is regarded as no hit.
# =====================================================================

def step08_scan_motifs(cfg: Config, force: bool = False):
    """Scan promoter sequences for -10/-35 motifs and classify as A/B/C.

    Writes into ``motifs/``:
    * ``motif_hits_all.tsv``       — all hits for all promoter-orientation IGRs
    * ``motif_best_all.tsv``       — best hit per IGR (all promoters)
    * ``motif_hits_markers.tsv``   — all hits for marker promoters only
    * ``motif_best_markers.tsv``   — best hit per IGR (markers only)
    * ``promoter_markers_verified.tsv`` — marker promoters confirmed by motif scan

    Writes at the top level of ``output_dir``:
    * ``all_promoters_verified.fasta``    — all motif-confirmed promoters,
      CO_R sequences reverse-complemented so every record is 5'→3'
    * ``marker_promoters_verified.fasta`` — marker + motif-confirmed promoters,
      CO_R sequences reverse-complemented so every record is 5'→3'
    """
    if not force and cfg.motif_best_all.exists() and cfg.promoter_markers_verified.exists():
        print("── Motif scan results already exist, skipping ──")
        print("  Step 8 complete.\n")
        return

    # Locate MEME file
    meme_file = cfg.meme_file
    if meme_file is None or not Path(meme_file).is_file():
        bundled = Path(__file__).parent / "all_unique_subgroups.meme"
        if bundled.is_file():
            meme_file = bundled
        else:
            print("── No MEME file found, skipping motif scan ──")
            pd.DataFrame().to_csv(cfg.promoter_markers_verified, sep="\t", index=False)
            print("  Step 8 complete.\n")
            return

    # Load paired subgroup motifs
    p10 = cfg.motif_p10
    p35 = cfg.motif_p35
    print(f"── Loading motifs from {Path(meme_file).name} "
          f"(p10={p10}, p35={p35}) ──")
    m10, m35 = _load_motifs_from_file(meme_file, p10, p35)

    if not m10.entries:
        print("  No -10 motifs loaded — cannot scan.")
        pd.DataFrame().to_csv(cfg.promoter_markers_verified, sep="\t", index=False)
        print("  Step 8 complete.\n")
        return

    # Report threshold ranges (one per subgroup is noisy to print, so
    # summarise across all subgroups).
    t10 = [thresh for _, _, thresh in m10.entries]
    t35 = [thresh for _, _, thresh in m35.entries]
    if t10:
        print(f"  -10 score thresholds: {min(t10):.3f} – {max(t10):.3f} "
              f"(p < {p10}, n={len(t10)} subgroups)")
    if t35:
        print(f"  -35 score thresholds: {min(t35):.3f} – {max(t35):.3f} "
              f"(p < {p35}, n={len(t35)} subgroups)")

    # --- Scan ALL promoter-orientation IGRs ---
    try:
        all_igr = pd.read_csv(cfg.igr_summary, sep="\t")
        all_igr = all_igr[all_igr["orientation"].isin(["CO_F", "CO_R"])].copy()
    except (FileNotFoundError, pd.errors.EmptyDataError):
        all_igr = pd.DataFrame()

    if all_igr.empty:
        print("── No promoter-orientation IGRs to scan ──")
        for p in [cfg.motif_hits_all, cfg.motif_best_all,
                  cfg.motif_hits_markers, cfg.motif_best_markers]:
            pd.DataFrame().to_csv(p, sep="\t", index=False)
        pd.DataFrame().to_csv(cfg.promoter_markers_verified, sep="\t", index=False)
        print("  Step 8 complete.\n")
        return

    # Orient sequences 5'→3'
    all_igr["sequence_5p_to_3p"] = all_igr.apply(
        lambda r: _reverse_complement(r["sequence"]) if r["orientation"] == "CO_R"
        else r["sequence"],
        axis=1,
    )

    print(f"\n  Scanning {len(all_igr)} promoter-orientation IGRs...")
    all_hits, best_all = _scan_sequences_for_motifs(all_igr, m10, m35)

    all_hits.to_csv(cfg.motif_hits_all, sep="\t", index=False)
    best_all.to_csv(cfg.motif_best_all, sep="\t", index=False)
    n_confirmed_all = len(best_all)
    print(f"  All promoters: {n_confirmed_all}/{len(all_igr)} confirmed by motif scan")
    print(f"    All hits -> {cfg.motif_hits_all}")
    print(f"    Best hits -> {cfg.motif_best_all}")

    # Write verified all-promoter FASTA
    confirmed_ids = set(best_all["igr_id"])
    confirmed_all_df = all_igr[all_igr["igr_id"].isin(confirmed_ids)].copy()
    _write_fasta(confirmed_all_df, cfg.all_promoters_verified_fasta, short_header=False)

    # --- Scan MARKER promoters ---
    marker_igr_ids = set()
    markers_df = pd.DataFrame()
    try:
        markers_df = pd.read_csv(cfg.promoter_markers, sep="\t")
        markers_df = markers_df[markers_df["orientation"].isin(["CO_F", "CO_R"])].copy()
        marker_igr_ids = set(markers_df["igr_id"])
    except (FileNotFoundError, pd.errors.EmptyDataError):
        pass

    if markers_df.empty:
        print("  No marker promoters to scan.")
        for p in [cfg.motif_hits_markers, cfg.motif_best_markers]:
            pd.DataFrame().to_csv(p, sep="\t", index=False)
        pd.DataFrame().to_csv(cfg.promoter_markers_verified, sep="\t", index=False)
        print("  Step 8 complete.\n")
        return

    # Orient marker sequences 5'→3'
    if "sequence_5p_to_3p" not in markers_df.columns:
        markers_df["sequence_5p_to_3p"] = markers_df.apply(
            lambda r: _reverse_complement(r["sequence"]) if r["orientation"] == "CO_R"
            else r["sequence"],
            axis=1,
        )

    print(f"\n  Scanning {len(markers_df)} marker promoters...")
    marker_hits, best_markers = _scan_sequences_for_motifs(markers_df, m10, m35)

    marker_hits.to_csv(cfg.motif_hits_markers, sep="\t", index=False)
    best_markers.to_csv(cfg.motif_best_markers, sep="\t", index=False)

    # Write verified markers: marker promoters that passed motif scan
    verified_ids = set(best_markers["igr_id"])
    verified_df = markers_df[markers_df["igr_id"].isin(verified_ids)].copy()

    # Merge best motif hit info onto verified markers
    motif_info = best_markers.set_index("igr_id")[
        ["strand", "pos_10", "seq_10", "score_10", "source_10",
         "has_ext10", "pos_35", "seq_35", "score_35", "source_35",
         "spacer_len", "spacer_seq", "path"]
    ].rename(columns=lambda c: f"motif_{c}")
    verified_df = verified_df.merge(motif_info, left_on="igr_id",
                                     right_index=True, how="left")

    verified_df.to_csv(cfg.promoter_markers_verified, sep="\t", index=False)
    n_confirmed_markers = len(verified_df)
    print(f"  Marker promoters: {n_confirmed_markers}/{len(markers_df)} confirmed")
    print(f"    Marker hits -> {cfg.motif_hits_markers}")
    print(f"    Best marker hits -> {cfg.motif_best_markers}")
    print(f"    Verified markers -> {cfg.promoter_markers_verified}")

    # Write verified marker promoter FASTA
    _write_fasta(verified_df, cfg.marker_promoters_verified_fasta, short_header=False)

    print("  Step 8 complete.\n")


# =====================================================================
#  Step 9 — Annotate associated CDS from Prokka
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


def step09_annotate_cds(cfg: Config, force: bool = False):
    """Extract CDS annotations (product, gene name, locus tag) from Prokka
    GFF for ALL promoter-orientation IGRs."""
    annotation_tsv = cfg.cds_annotations
    if not force and annotation_tsv.exists():
        print("── CDS annotations already exist, skipping ──")
        print("  Step 9 complete.\n")
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
        print("  Step 9 complete.\n")
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
    print("  Step 9 complete.\n")


# =====================================================================
#  Step 11 — HTML report
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

    # Colour palette for motif elements
    db_colours = {
        "minus35": "#1565c0",
        "minus10": "#c62828",
        "ext10": "#ff8f00",
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


def step11_generate_report(cfg: Config, force: bool = False):
    """Generate an HTML report combining promoters, CDS annotations,
    and motif hit positions with visual sequence diagrams."""
    if not force and cfg.report_html.exists():
        print("── HTML report already exists, skipping ──")
        print("  Step 11 complete.\n")
        return

    print("── Generating HTML report ──")

    # Load verified promoters (after motif scan).  Step 8 may have
    # written an empty placeholder file on a bailout path, so fall back
    # to promoter_markers.tsv if the verified file is empty or lacks the
    # expected columns.
    candidates = []
    if cfg.promoter_markers_verified.exists():
        candidates.append(cfg.promoter_markers_verified)
    if cfg.promoter_markers.exists():
        candidates.append(cfg.promoter_markers)

    promoter_df = None
    for source_file in candidates:
        try:
            df = pd.read_csv(source_file, sep="\t")
        except pd.errors.EmptyDataError:
            continue
        if "orientation" not in df.columns:
            continue
        df = df[df["orientation"].isin(["CO_F", "CO_R"])].copy()
        if not df.empty:
            promoter_df = df
            break

    if promoter_df is None or promoter_df.empty:
        print("  No verified promoter data available for report.")
        print("  Step 11 complete.\n")
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
    annotation_tsv = cfg.cds_annotations
    if annotation_tsv.exists():
        try:
            ann_df = pd.read_csv(annotation_tsv, sep="\t")
            for _, row in ann_df.iterrows():
                annotation_map[row["gene_id"]] = str(row.get("product", ""))
        except (pd.errors.EmptyDataError, KeyError):
            pass

    # Build motif hit info from verified markers (motif_pos_10, motif_pos_35, etc.)
    # These columns were added by step 8 when writing promoter_markers_verified.tsv
    has_motif_cols = "motif_pos_10" in promoter_df.columns

    # Build HTML
    n_promoters = len(promoter_df)
    n_annotated = sum(1 for _, r in promoter_df.iterrows()
                      if r["associated_gene"] in annotation_map
                      and annotation_map[r["associated_gene"]] != "hypothetical protein")
    n_with_motifs = sum(1 for _, r in promoter_df.iterrows()
                        if has_motif_cols and pd.notna(r.get("motif_pos_10"))
                        and str(r.get("motif_pos_10")) != ".")

    genome_name = cfg.input_fasta.stem

    # Paths for file links (relative to output_dir).
    all_verified_name = (cfg.all_promoters_verified_fasta.name
                         if cfg.all_promoters_verified_fasta.exists() else "")
    marker_verified_name = (cfg.marker_promoters_verified_fasta.name
                            if cfg.marker_promoters_verified_fasta.exists() else "")
    cds_ann_name = cfg.cds_annotations.name if cfg.cds_annotations.exists() else ""
    final_table_name = cfg.final_table.name if cfg.final_table.exists() else ""

    domain_label = "Archaea" if cfg.is_archaea else "Bacteria"

    html_parts = [_REPORT_HTML_HEAD.format(
        genome_name=html_mod.escape(genome_name),
        n_marker_promoters=n_promoters,
        n_all_promoters=n_all_promoters,
        n_annotated=n_annotated,
        n_with_motifs=n_with_motifs,
        all_verified_fasta=html_mod.escape(all_verified_name),
        marker_verified_fasta=html_mod.escape(marker_verified_name),
        cds_annotations_tsv=html_mod.escape(cds_ann_name),
        final_table=html_mod.escape(final_table_name),
        domain_label=domain_label,
    )]

    # Legend for motif diagram
    html_parts.append("""
    <div class="legend">
        <strong>Motif diagram:</strong>
        <span class="legend-item"><span class="legend-swatch" style="background:#1565c0"></span>&minus;35 element</span>
        <span class="legend-item"><span class="legend-swatch" style="background:#c62828"></span>&minus;10 element</span>
        <span class="legend-item"><span class="legend-swatch" style="background:#ff8f00"></span>Extended &minus;10 (TG)</span>
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
        <th>Motif path</th>
        <th>&minus;10</th>
        <th>&minus;35</th>
        <th>Spacer</th>
        <th>Associated CDS</th>
        <th>Protein</th>
        <th>Sequence diagram</th>
        <th>Sequence (5'&rarr;3')</th>
    </tr>
    </thead>
    <tbody>
    """)

    for _, row in promoter_df.iterrows():
        igr_id = row["igr_id"]
        gene = row["associated_gene"]
        product = annotation_map.get(gene, "")

        # Build motif diagram from -10 and -35 positions
        motif_hits = []
        if has_motif_cols:
            pos_10 = row.get("motif_pos_10", ".")
            pos_35 = row.get("motif_pos_35", ".")
            if pos_10 != "." and pd.notna(pos_10):
                try:
                    p10 = int(float(pos_10))
                    motif_hits.append({
                        "motif_id": "-10",
                        "start": p10 + 1,  # 0-based to 1-based
                        "stop": p10 + _MOTIF_WIDTH,
                        "motif_database": "minus10",
                    })
                    # Extended -10 TG.  The column may be a Python bool
                    # (unread DataFrame), the string "True"/"False" (TSV
                    # round-trip), or missing altogether.
                    ext10_raw = row.get("motif_has_ext10", False)
                    if isinstance(ext10_raw, bool):
                        has_ext10_flag = ext10_raw
                    else:
                        has_ext10_flag = str(ext10_raw).strip().lower() == "true"
                    if has_ext10_flag:
                        motif_hits.append({
                            "motif_id": "ext-10 (TG)",
                            "start": p10 - 1,  # TG is 2 nt before -10
                            "stop": p10,
                            "motif_database": "ext10",
                        })
                except (ValueError, TypeError):
                    pass

            if pos_35 != "." and pd.notna(pos_35):
                try:
                    p35 = int(float(pos_35))
                    motif_hits.append({
                        "motif_id": "-35",
                        "start": p35 + 1,
                        "stop": p35 + _MOTIF_WIDTH,
                        "motif_database": "minus35",
                    })
                except (ValueError, TypeError):
                    pass

        svg = _build_motif_diagram_svg(row["length"], motif_hits)

        orient_class = "co-f" if row["orientation"] == "CO_F" else "co-r"
        seq = row.get("sequence_5p_to_3p", "")

        # Motif details
        motif_path = str(row.get("motif_path", "")) if has_motif_cols else ""
        seq_10 = str(row.get("motif_seq_10", "")) if has_motif_cols else ""
        score_10 = row.get("motif_score_10", "") if has_motif_cols else ""
        source_10 = str(row.get("motif_source_10", "")) if has_motif_cols else ""
        seq_35 = str(row.get("motif_seq_35", "")) if has_motif_cols else ""
        score_35 = row.get("motif_score_35", "") if has_motif_cols else ""
        source_35 = str(row.get("motif_source_35", "")) if has_motif_cols else ""
        spacer_len = row.get("motif_spacer_len", "") if has_motif_cols else ""

        minus10_cell = f"{html_mod.escape(seq_10)}" if seq_10 and seq_10 != "." else '<span class="na">—</span>'
        if score_10 and str(score_10) != ".":
            minus10_cell += f' <small>({score_10}, {html_mod.escape(source_10)})</small>'

        minus35_cell = f"{html_mod.escape(seq_35)}" if seq_35 and seq_35 != "." else '<span class="na">—</span>'
        if score_35 and str(score_35) != ".":
            minus35_cell += f' <small>({score_35}, {html_mod.escape(source_35)})</small>'

        spacer_cell = str(spacer_len) if spacer_len and str(spacer_len) != "." else '<span class="na">—</span>'

        product_cell = html_mod.escape(product) if product else '<span class="na">—</span>'

        html_parts.append(f"""
        <tr>
            <td><code>{html_mod.escape(igr_id)}</code></td>
            <td>{html_mod.escape(str(row['contig']))}</td>
            <td>{row['start']:,}–{row['end']:,}</td>
            <td>{row['length']}</td>
            <td><span class="orient {orient_class}">{row['orientation']}</span></td>
            <td>{html_mod.escape(motif_path)}</td>
            <td>{minus10_cell}</td>
            <td>{minus35_cell}</td>
            <td>{spacer_cell}</td>
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
    print("  Step 11 complete.\n")


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
    small {{ color: #888; }}
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
    <div class="stat"><div class="label">Verified marker promoters</div><div class="value">{n_marker_promoters}</div></div>
    <div class="stat"><div class="label">CDS annotated</div><div class="value">{n_annotated}</div></div>
    <div class="stat"><div class="label">With &minus;10/&minus;35 motifs</div><div class="value">{n_with_motifs}</div></div>
    <div class="stat"><div class="label">All promoter IGRs</div><div class="value">{n_all_promoters}</div></div>
</div>

<p style="font-size:0.85rem; color:#555; margin-bottom:20px;">
    Domain: <strong>{domain_label}</strong> · Verification: <strong>&minus;10/&minus;35 motif scan</strong><br/>
    This report shows the <strong>{n_marker_promoters}</strong> marker-filtered promoters confirmed by &minus;10/&minus;35 motif scanning.
    A total of {n_all_promoters} promoter-orientation IGRs were identified in the genome.
</p>
<p style="font-size:0.85rem; color:#555; margin-bottom:20px;">
    <strong>Downloads:</strong>
    <a href="{marker_verified_fasta}">marker_promoters_verified.fasta</a> &middot;
    <a href="{all_verified_fasta}">all_promoters_verified.fasta</a> &middot;
    <a href="{cds_annotations_tsv}">cds_annotations.tsv</a> &middot;
    <a href="{final_table}">profinder_results.tsv</a>
</p>

<h2>Verified marker promoters</h2>
"""

_REPORT_HTML_FOOT = """
<footer>
    Generated by ProFinder.
</footer>
</body>
</html>
"""

# =====================================================================
#  Step 10 — Final output table
# =====================================================================

def step10_final_table(cfg: Config, force: bool = False):
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
    motif_confirmed       Whether a -10/-35 motif combination was found
    motif_path            Motif path classification (A, B, or C)
    motif_strand          Strand on which best motif hit was found
    motif_pos_10          Position of -10 element in scanned sequence
    motif_seq_10          Sequence of -10 element
    motif_score_10        Log-odds score of -10 hit
    motif_source_10       Subgroup ID of the best -10 hit (e.g. M001)
    motif_has_ext10       Whether extended -10 TG dinucleotide is present
    motif_pos_35          Position of -35 element
    motif_seq_35          Sequence of -35 element
    motif_score_35        Log-odds score of -35 hit
    motif_source_35       Subgroup ID of the best -35 hit (e.g. M001)
    motif_spacer_len      Spacer length between -35 and -10
    motif_spacer_seq      Spacer sequence
    sequence_5p_to_3p     Full promoter sequence oriented 5'→3'
    """
    if not force and cfg.final_table.exists():
        print("── Final output table already exists, skipping ──")
        print("  Step 10 complete.\n")
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
        print("  Step 10 complete.\n")
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
    annotation_tsv = cfg.cds_annotations
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

    # 6. Motif scan results (best hit per IGR from step 8)
    motif_cols = ["strand", "pos_10", "seq_10", "score_10", "source_10",
                  "has_ext10", "pos_35", "seq_35", "score_35", "source_35",
                  "spacer_len", "spacer_seq", "path"]
    motif_map = {}   # igr_id -> dict of motif columns
    if cfg.motif_best_all.exists():
        try:
            motif_df = pd.read_csv(cfg.motif_best_all, sep="\t")
            for _, row in motif_df.iterrows():
                motif_map[row["igr_id"]] = {c: row.get(c, "") for c in motif_cols}
        except (pd.errors.EmptyDataError, KeyError):
            pass

    igr["motif_confirmed"] = igr["igr_id"].isin(set(motif_map.keys()))
    for c in motif_cols:
        col_name = f"motif_{c}"
        igr[col_name] = igr["igr_id"].map(
            lambda g, _c=c: motif_map.get(g, {}).get(_c, ""))

    # 7. CDS-extended sequences (optional)
    columns = [
        "igr_id", "contig", "start", "end", "length", "orientation",
        "associated_cds", "gene_name", "locus_tag", "product",
        "is_marker", "motif_confirmed", "motif_path", "motif_strand",
        "motif_pos_10", "motif_seq_10", "motif_score_10", "motif_source_10",
        "motif_has_ext10",
        "motif_pos_35", "motif_seq_35", "motif_score_35", "motif_source_35",
        "motif_spacer_len", "motif_spacer_seq",
        "sequence_5p_to_3p",
    ]

    if cfg.cds_bp > 0:
        contigs = _load_contigs(cfg.input_fasta)
        igr = _add_cds_column(igr, contigs, cfg.cds_bp)
        columns.append("sequence_5p_to_3p_cds")

    # 8. Select and order final columns
    out = igr[columns].copy()
    out.rename(columns={"igr_id": "promoter_id"}, inplace=True)

    out.to_csv(cfg.final_table, sep="\t", index=False)
    print(f"  Final table ({len(out)} rows) -> {cfg.final_table}")
    print("  Step 10 complete.\n")


STEPS = [
    (1,  "Run Prokka",                       step01_run_prokka),
    (2,  "Extract intergenic regions",       step02_extract_igrs),
    (3,  "Identify operons",                 step03_identify_operons),
    (4,  "Run hmmsearch",                    step04_run_hmmsearch),
    (5,  "Filter HMM output",                step05_filter_hmm),
    (6,  "Filter operons + add markers",     step06_filter_operons_add_markers),
    (7,  "Match IGRs to marker operons",     step07_match_igrs_to_markers),
    (8,  "Scan promoter motifs (-10/-35)",   step08_scan_motifs),
    (9,  "Annotate CDS (Prokka)",            step09_annotate_cds),
    (10, "Build final output table",         step10_final_table),
    (11, "Generate HTML report",             step11_generate_report),
]


def main():
    parser = argparse.ArgumentParser(
        description="ProFinder — bacterial and archaeal promoter identification pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    # Single-sample and batch inputs are mutually exclusive.  Requirement
    # is enforced manually after parse_args so --list can run without
    # either (an argparse `required=` flag here would need to know about
    # --list before the parser has run).
    input_group = parser.add_mutually_exclusive_group()
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
                        help="Target domain (affects Prokka --kingdom). "
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

    # Motif scanning
    parser.add_argument("--meme-file", type=Path, default=None,
                        help="Single MEME file with paired M###_m10 / "
                             "M###_m35 subgroup motifs "
                             "(default: bundled all_unique_subgroups.meme)")
    parser.add_argument("--p10", type=float, default=2.5e-3,
                        help="p-value threshold for -10 motif hits "
                             "(default: 0.0025)")
    parser.add_argument("--p35", type=float, default=2.5e-3,
                        help="p-value threshold for -35 motif hits "
                             "(default: 0.0025)")

    # Conda environments
    parser.add_argument("--conda-prokka", default="",
                        help="Conda env for Prokka (blank = use current env)")
    parser.add_argument("--conda-hmm", default="",
                        help="Conda env for hmmsearch (blank = use current env)")

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

    # One of --input / --batch is required for any real run. We enforce
    # this here (rather than via argparse `required=`) so --list can run
    # without either.
    if args.input is None and args.batch is None:
        parser.error("one of the arguments -i/--input --batch is required")

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
        hmm_profiles_dir=args.hmm_dir,
        conda_env_prokka=args.conda_prokka,
        conda_env_hmm=args.conda_hmm,
        threads=args.threads,
        prokka_kingdom=kingdom,
        prokka_prefix=prokka_prefix or args.prefix,
        igr_size_min=args.igr_min,
        igr_size_max=args.igr_max,
        max_internal_distance=args.max_internal_dist,
        min_flanking_distance=args.min_flanking_dist,
        hmm_bitscore_min=args.hmm_bitscore,
        meme_file=args.meme_file,
        motif_p10=args.p10,
        motif_p35=args.p35,
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
