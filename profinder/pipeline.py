#!/usr/bin/env python3
"""
ProFinder — bacterial and archaeal promoter identification pipeline.

Takes a single genome FASTA, annotates it with Prokka, extracts
intergenic regions, identifies operons, screens for HMM marker genes,
and outputs promoter sequences in 5'-to-3' orientation.

Usage
-----
    # Minimal (Prokka + IGR extraction + operon ID + promoter output):
    profinder -i genome.fasta -o results/

    # With a custom HMM profile directory (default: bundled 98 profiles):
    profinder -i genome.fasta -o results/ --hmm-dir /path/to/hmms/

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
    4.  Run hmmsearch                        [bundled HMM profiles by default]
    5.  Filter HMM output
    6.  Filter operons and add marker info
    7.  Match IGRs to marker operons
    8.  Scan promoter motifs (-10/-35 hexamer verification, A/B/C)
    9.  Annotate associated CDS (Prokka product names, gene, locus tag)
    10. Build final output table (all promoters, all metadata)
    11. Generate HTML report

Diagnostic mode
---------------
    --diagnose-p-sweep replaces steps 8, 10, 11 with a p-value sweep.
    The pipeline runs steps 1-7 and step 9 (Prokka, IGR extraction,
    operons, HMM, marker matching, CDS annotation — all p-value
    independent), then scans every CO_F/CO_R IGR at a range of
    (--p10, --p35) settings in three modes — matched, p10_only,
    p35_only — and writes:

      diagnostics/p_sweep.tsv         per-sweep summary counts
      diagnostics/p_sweep.png         3-panel figure (if matplotlib)
      diagnostics/sweep_tables/       one profinder_results.tsv per
                                      (mode, p10, p35), pre-filtered
                                      to motif_confirmed = True and
                                      motif_path in {A, B}

    Each per-sweep-point table has the same column layout as the
    canonical profinder_results.tsv produced by step 10, so it can be
    consumed directly by downstream benchmarking (e.g. Section 5).
    Useful for choosing defensible p-value thresholds for a given
    genome before committing to a full run.
"""

import argparse
import hashlib
import html as html_mod
import json
import math
import os
import re
import shlex
import subprocess
import sys
from collections import Counter
from pathlib import Path
import gzip
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

from .config import Config
from .igr_extractor import extract_igrs,_open_maybe_gzip


_REVCOMP_TRANS = str.maketrans("ACGTacgt", "TGCAtgca")


def _revcomp_simple(s: str) -> str:
    """Fast reverse-complement using the module-level translation table."""
    return s.translate(_REVCOMP_TRANS)[::-1]


# =====================================================================
#  Input, tool and prerequisite guards
# =====================================================================

_PROKKA_MIN_VERSION = (1, 14, 6)
_PROKKA_MAX_CONTIG_ID = 37          # Prokka refuses longer contig IDs
_IGR_COLUMNS = ["igr_id", "contig", "start", "end", "length", "orientation",
                "left_gene", "right_gene", "left_type", "right_type", "sequence",
                "sequence_5p_to_3p"]
_OPERON_COLUMNS = ["Operon", "SeqID", "Start", "End", "Attributes", "Gene",
                   "SourceFile", "FirstGene", "Strand"]


def _require(step, *needed):
    """Stop with a clear message when a step's inputs are missing.
    ``needed`` holds (path, producing_step) pairs."""
    missing = [(Path(p), s) for p, s in needed if not Path(p).exists()]
    if missing:
        lines = "\n".join(f"    {p}   (written by step {s})" for p, s in missing)
        first = min(s for _, s in missing)
        sys.exit(f"Step {step} cannot run: required input(s) missing:\n{lines}\n"
                 f"Run the earlier steps first (e.g. --start {first}), or check that "
                 f"-o points at the output directory of a previous run.")


def _run_tool(conda_env, command):
    """Run a shell command under bash, optionally inside a conda env."""
    prefix = (f'eval "$(conda shell.bash hook 2>/dev/null)"; conda activate {conda_env}; '
              if conda_env else "")
    return subprocess.run(["bash", "-c", "set -euo pipefail; " + prefix + command],
                          capture_output=True, text=True)


def _check_prokka(cfg):
    """Fail fast when Prokka is missing or too old to annotate."""
    r = _run_tool(cfg.conda_env_prokka, f"{cfg.prokka_bin} --version")
    text = (r.stdout + r.stderr).strip()
    m = re.search(r"prokka\s+(\d+)\.(\d+)(?:\.(\d+))?", text, re.I)
    where = f" in conda env {cfg.conda_env_prokka}" if cfg.conda_env_prokka else ""
    if r.returncode != 0 or not m:
        sys.exit(f"Prokka is not runnable as '{cfg.prokka_bin}'{where}: "
                 f"{text[-300:] or 'no output'}\n"
                 f"Install it with: conda install -c conda-forge -c bioconda 'prokka>=1.14.6' "
                 f"(or pass --prokka / --conda-prokka).")
    ver = tuple(int(x or 0) for x in m.groups())
    if ver < _PROKKA_MIN_VERSION:
        sys.exit(f"Prokka {'.'.join(map(str, ver))} found{where}, but ProFinder needs Prokka >= 1.14.6: "
                 f"older releases misread BLAST >= 2.10 as 'too old' and cannot annotate. "
                 f"Upgrade with: conda install -c conda-forge -c bioconda 'prokka>=1.14.6'.")
    return ver


def _validate_genome_fasta(path):
    """Check the input is a nucleotide FASTA with unique IDs before spending a
    Prokka run on it. Returns the contig IDs."""
    ids, n_bases, n_acgtn = [], 0, 0
    try:
        with _open_maybe_gzip(path) as fh:
            for rec in SeqIO.parse(fh, "fasta"):
                seq = str(rec.seq).upper()
                ids.append(rec.id)
                n_bases += len(seq)
                n_acgtn += sum(seq.count(b) for b in "ACGTN")
    except (OSError, ValueError, UnicodeDecodeError, EOFError) as e:
        sys.exit(f"Cannot read input FASTA {path}: {e}")
    if not ids or n_bases == 0:
        sys.exit(f"Input FASTA {path} contains no sequences (empty file, or not FASTA format).")
    frac = n_acgtn / n_bases
    if frac < 0.9:
        sys.exit(f"Input FASTA {path} does not look like a nucleotide genome: only "
                 f"{100 * frac:.0f}% of its characters are A/C/G/T/N. ProFinder needs a "
                 f"genome assembly (nucleotide FASTA), not protein sequences.")
    dup = [i for i, k in Counter(ids).items() if k > 1]
    if dup:
        sys.exit(f"Input FASTA {path} has duplicate contig IDs: {', '.join(dup[:5])}"
                 f"{' ...' if len(dup) > 5 else ''}. Make every FASTA header unique.")
    return ids


def _shorten_contig_ids(cfg, src):
    """Prokka refuses contig IDs longer than 37 characters. Write a copy of
    ``src`` with short IDs and a map back to the originals."""
    mapping = {}
    out = cfg.prokka_dir / "input.short_ids.fna"
    with _open_maybe_gzip(src) as fh, open(out, "w") as fo:
        for k, rec in enumerate(SeqIO.parse(fh, "fasta"), 1):
            short = f"contig_{k:06d}"
            mapping[short] = rec.id
            fo.write(f">{short}\n{str(rec.seq)}\n")
    with open(cfg.prokka_dir / "contig_id_map.tsv", "w") as fo:
        fo.write("prokka_id\toriginal_id\n")
        for short, orig in mapping.items():
            fo.write(f"{short}\t{orig}\n")
    return out, mapping


def _restore_contig_ids(cfg, mapping):
    """Rewrite Prokka's GFF (column 1, ##sequence-region, embedded FASTA) and
    FNA headers from the short IDs back to the original contig IDs, so every
    ProFinder output carries the user's contig names."""
    out, in_fasta = [], False
    for line in cfg.gff_file.read_text().splitlines(keepends=True):
        if line.startswith("##FASTA"):
            in_fasta = True
        elif in_fasta:
            if line.startswith(">"):
                sid = line[1:].split()[0]
                line = f">{mapping.get(sid, sid)}\n"
        elif line.startswith("##sequence-region"):
            parts = line.split()
            if len(parts) >= 2 and parts[1] in mapping:
                parts[1] = mapping[parts[1]]
                line = " ".join(parts) + "\n"
        elif not line.startswith("#"):
            cols = line.split("\t", 1)
            if len(cols) == 2 and cols[0] in mapping:
                line = mapping[cols[0]] + "\t" + cols[1]
        out.append(line)
    cfg.gff_file.write_text("".join(out))
    if cfg.fna_file.exists():
        fna = []
        for line in cfg.fna_file.read_text().splitlines(keepends=True):
            if line.startswith(">"):
                sid = line[1:].split()[0]
                line = f">{mapping.get(sid, sid)}\n"
            fna.append(line)
        cfg.fna_file.write_text("".join(fna))


# =====================================================================
#  Motif scanning helpers (adapted from meme_scan.py)
# =====================================================================

_MOTIF_WIDTH = 6
_SPACER_MIN = 16
_SPACER_MAX = 18
_BASE_IDX = {"A": 0, "C": 1, "G": 2, "T": 3}


def _parse_meme_file(filepath):
    """Parse a MEME-format file. Returns {motif_name: [[pA, pC, pG, pT], ...]}.
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


_UNIFORM_BG = (0.25, 0.25, 0.25, 0.25)


def _compute_background_from_fasta(fasta_path):
    """Compute genome-wide A/C/G/T frequencies from a FASTA.
    """
    counts = {"A": 0, "C": 0, "G": 0, "T": 0}
    try:
        with _open_maybe_gzip(fasta_path) as fh:
            for rec in SeqIO.parse(fh, "fasta"):
                s = str(rec.seq).upper()
                for b in "ACGT":
                    counts[b] += s.count(b)
    except Exception:
        return _UNIFORM_BG
    total = sum(counts.values())
    if total == 0:
        return _UNIFORM_BG
    return (counts["A"] / total, counts["C"] / total,
            counts["G"] / total, counts["T"] / total)


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


OPERON_RULE = "same-strand runs; promoter-side flank"   # recorded with operons.tsv (step state)


def _identify_operons(genes_df, max_internal_distance, min_flanking_distance):
    """Strand-aware operon identification.

    Pass 1 groups consecutive CDS on the same contig and the same strand
    whose gap is at most ``max_internal_distance`` bp. The gap is measured
    from the furthest end the group has reached, so a gene nested inside
    another cannot split it. A change of strand always ends a group, however
    close the genes are: genes on opposite strands cannot share a transcript.

    Pass 2 tests only the promoter side of each group, the gap upstream of
    its first gene in the direction of transcription (left of a + strand
    group, right of a - strand group). If the neighbouring group on that side
    lies on the same strand and less than ``min_flanking_distance`` bp away,
    the group may continue an operon that starts further upstream, and it is
    dropped. An opposite-strand neighbour, a contig end or a gap of at least
    ``min_flanking_distance`` bp marks an unambiguous operon start. The gap
    after the last gene is not tested: it bounds the operon's 3' end, not its
    promoter.

    A CDS whose strand is neither '+' nor '-' ends any group and is never
    placed in an operon. ``genes_df`` must be sorted by contig and start
    (as returned by _parse_gff_for_operons). Returns a list of operons, each
    a list of gene rows in genome order.
    """
    if genes_df.empty:
        return []
    seq = genes_df["seqid"].tolist()
    beg = genes_df["start"].tolist()
    end = genes_df["end"].tolist()
    strand = genes_df["strand"].tolist()

    # Pass 1: same-strand proximity runs
    groups, cur, reach = [], [0], end[0]
    for i in range(1, len(seq)):
        j = cur[-1]
        if (seq[i] == seq[j] and strand[i] == strand[j] and strand[i] in ("+", "-")
                and beg[i] - reach <= max_internal_distance):
            cur.append(i)
            reach = max(reach, end[i])
        else:
            groups.append(cur)
            cur, reach = [i], end[i]
    groups.append(cur)

    # Pass 2: the promoter-side flank
    lo = [beg[g[0]] for g in groups]
    hi = [max(end[j] for j in g) for g in groups]
    operons = []
    for k, g in enumerate(groups):
        s = strand[g[0]]
        if s not in ("+", "-"):
            continue
        nb = k - 1 if s == "+" else k + 1
        if (0 <= nb < len(groups) and seq[groups[nb][0]] == seq[g[0]]
                and strand[groups[nb][0]] == s):
            gap = lo[k] - hi[nb] if s == "+" else lo[nb] - hi[k]
            if gap < min_flanking_distance:
                continue
        operons.append([genes_df.iloc[j] for j in g])
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

# ── Rebuild on change: the operon and marker steps ───────────────────────
# Steps 3, 6, 7, the marker shortlist of step 8, steps 10 and 11 and the
# p-value sweep are quick and all depend on the operon rule. Each records a
# fingerprint of its input files and settings in <output>/.profinder_state.json
# and is skipped on a resumed run only when its outputs exist and that
# fingerprint still matches. A changed operon rule, operon setting or upstream
# file therefore rebuilds every dependent output, and a directory without the
# state file (an older ProFinder, or a copy made without dot-files) rebuilds
# them once. The genome-wide motif scan of step 8 is never repeated for this.
# Bump _STATE_FORMAT whenever the logic of one of these steps changes.
_STATE_FORMAT = 1
_STATE_NAME = ".profinder_state.json"
_LEGACY_WARNED = set()


def _fingerprint(paths, params=None):
    """sha256 over the contents of ``paths`` (None entries skipped; a missing
    file is recorded as missing) and the JSON form of ``params``."""
    h = hashlib.sha256(f"profinder-state-{_STATE_FORMAT}".encode())
    for p in paths:
        if p is None:
            continue
        p = Path(p)
        h.update(b"\x00" + p.name.encode() + b"\x00")
        if p.is_file():
            h.update(str(p.stat().st_size).encode() + b"\x00")
            with open(p, "rb") as fh:
                for chunk in iter(lambda: fh.read(1 << 20), b""):
                    h.update(chunk)
        else:
            h.update(b"missing")
    h.update(json.dumps(params or {}, sort_keys=True, default=str).encode())
    return h.hexdigest()


def _state_path(cfg):
    return Path(cfg.output_dir) / _STATE_NAME


def _load_state(cfg):
    try:
        with open(_state_path(cfg)) as fh:
            state = json.load(fh)
    except (OSError, ValueError):
        return {}
    return state if isinstance(state, dict) else {}


def _save_state(cfg, state):
    path = _state_path(cfg)
    tmp = path.with_name(path.name + ".tmp")
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(tmp, "w") as fh:
            json.dump(state, fh, indent=1, sort_keys=True)
            fh.write("\n")
        os.replace(tmp, path)
    except OSError as e:
        print(f"  WARNING: could not write {path} ({e}); a resumed run will rebuild this step")


def _is_current(cfg, key, outputs, inputs, params=None):
    """True when every output exists and step ``key`` was last built from
    exactly these input files and settings."""
    if not all(Path(o).exists() for o in outputs):
        return False
    return _load_state(cfg).get(key) == _fingerprint(inputs, params)


def _mark_current(cfg, key, inputs, params=None):
    """Record that step ``key`` has just been built from these inputs."""
    state = _load_state(cfg)
    state[key] = _fingerprint(inputs, params)
    _save_state(cfg, state)


def _forget(cfg, key):
    """Drop ``key`` before rebuilding it, so an interrupted rebuild is redone."""
    state = _load_state(cfg)
    if key in state:
        del state[key]
        _save_state(cfg, state)


def _operon_params(cfg):
    return {"rule": OPERON_RULE, "max_internal_distance": cfg.max_internal_distance,
            "min_flanking_distance": cfg.min_flanking_distance}


def _hmm_profile_names(cfg):
    d = cfg.hmm_profiles_dir
    return sorted(p.name for p in Path(d).glob("*.hmm")) if d is not None and Path(d).is_dir() else []


def _step10_inputs(cfg):
    return [cfg.igr_summary, cfg.promoter_markers, cfg.motif_best_all, cfg.cds_annotations,
            cfg.marker_best_hits, cfg.operon_file, cfg.operon_filtered_markers, cfg.gff_file,
            cfg.input_fasta]


def _warn_legacy_operons(cfg):
    """Warn once per output directory when operons.tsv predates the strand-aware rule."""
    key = str(Path(cfg.output_dir).resolve())
    if key in _LEGACY_WARNED or not cfg.operon_file.exists() or _has_column(cfg.operon_file, "Strand"):
        return
    _LEGACY_WARNED.add(key)
    print("  WARNING: operons/operons.tsv was made by an older ProFinder whose operons ignore "
          "strand; the marker results use it as it is. Re-run from step 3 to apply the "
          "current operon rule.")


def step01_run_prokka(cfg: Config, force: bool = False):
    """Annotate the input genome with Prokka."""
    if not force and cfg.gff_file.exists() and cfg.faa_file.exists():
        print("── Prokka output already exists, skipping ──")
        print("  Step 1 complete.\n")
        return

    print("── Running Prokka ──")
    cfg.prokka_dir.mkdir(parents=True, exist_ok=True)
    ids = _validate_genome_fasta(cfg.input_fasta)
    _check_prokka(cfg)

    input_for_prokka = cfg.input_fasta
    mapping = None
    long_ids = [i for i in ids if len(i) > _PROKKA_MAX_CONTIG_ID]
    if long_ids:
        input_for_prokka, mapping = _shorten_contig_ids(cfg, cfg.input_fasta)
        print(f"  NOTE: {len(long_ids)} contig ID(s) exceed Prokka's {_PROKKA_MAX_CONTIG_ID}-character "
              f"limit (e.g. {long_ids[0][:60]}). Prokka runs on short IDs; the original IDs are "
              f"restored in prokka/{cfg.gff_file.name} and prokka/{cfg.fna_file.name} "
              f"(map: prokka/contig_id_map.tsv).")
    elif str(cfg.input_fasta).endswith(".gz"):
        input_for_prokka = cfg.output_dir / f"{cfg.input_fasta.stem}.decompressed.fna"
        with gzip.open(cfg.input_fasta, "rt") as fin, open(input_for_prokka, "w") as fout:
            fout.write(fin.read())

    q = shlex.quote
    result = _run_tool(cfg.conda_env_prokka,
                       f"{cfg.prokka_bin} --outdir {q(str(cfg.prokka_dir))} --prefix {q(cfg.prokka_prefix)}"
                       f" --kingdom {q(cfg.prokka_kingdom)} --cpus {cfg.threads} --rfam --force"
                       f" {q(str(input_for_prokka))}")
    if result.returncode != 0:
        print(f"  STDERR: {result.stderr[-500:]}")
        sys.exit(f"Prokka failed with return code {result.returncode}")
    if mapping:
        _restore_contig_ids(cfg, mapping)

    print(f"  GFF -> {cfg.gff_file}")
    print(f"  FAA -> {cfg.faa_file}")
    print("  Step 1 complete.\n")


def step02_extract_igrs(cfg: Config, force: bool = False):
    """Extract intergenic regions from the Prokka GFF and genome FASTA."""
    if not force and cfg.igr_summary.exists() and cfg.igr_fasta.exists():
        print("── IGR files already exist, skipping ──")
        print("  Step 2 complete.\n")
        return

    _require(2, (cfg.gff_file, 1))
    print("── Extracting intergenic regions ──")
    cfg.igr_dir.mkdir(parents=True, exist_ok=True)
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

    if igr_df.empty:
        igr_df = pd.DataFrame(columns=_IGR_COLUMNS)     # header-only: later steps can read it
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
    """Identify operons from the Prokka GFF (rule: see _identify_operons)."""
    s3_in, s3_par = [cfg.gff_file], _operon_params(cfg)
    if not force and _is_current(cfg, "step3", [cfg.operon_file], s3_in, s3_par):
        print("── Operon file already exists, skipping ──")
        print("  Step 3 complete.\n")
        return
    if not force and cfg.operon_file.exists():
        print("── Operon file was made from another annotation, other operon settings or an "
              "older operon rule; rebuilding ──")

    _require(3, (cfg.gff_file, 1))
    _forget(cfg, "step3")
    print("── Identifying operons ──")
    genes_df = _parse_gff_for_operons(str(cfg.gff_file))
    if genes_df.empty:
        print("  WARNING: no CDS features found in GFF")
        pd.DataFrame(columns=_OPERON_COLUMNS).to_csv(cfg.operon_file, sep="\t", index=False)
        _mark_current(cfg, "step3", s3_in, s3_par)
        print("  Step 3 complete.\n")
        return

    operons = _identify_operons(genes_df, cfg.max_internal_distance,
                                cfg.min_flanking_distance)
    source = cfg.prokka_prefix

    rows = []
    for op_idx, operon in enumerate(operons, start=1):
        operon_id = f"operon{op_idx}"
        # An operon lies on one strand. Its first gene in the direction of
        # transcription is the leftmost on the + strand (yes_1) and the one
        # that reaches furthest right on the - strand (yes_2).
        if operon[0]["strand"] == "+":
            first, tag = 0, "yes_1"
        else:
            first = max(range(len(operon)), key=lambda j: (operon[j]["end"], j))
            tag = "yes_2"
        for gi, gene in enumerate(operon):
            rows.append({
                "Operon": operon_id,
                "SeqID": gene["seqid"],
                "Start": gene["start"],
                "End": gene["end"],
                "Attributes": gene["attributes"],
                "Gene": _extract_gene_id(gene["attributes"]),
                "SourceFile": source,
                "FirstGene": tag if gi == first else "No",
                "Strand": gene["strand"],
            })

    df = pd.DataFrame(rows, columns=_OPERON_COLUMNS)
    df.to_csv(cfg.operon_file, sep="\t", index=False)
    _mark_current(cfg, "step3", s3_in, s3_par)
    print(f"  Operons ({len(operons)} operons, {len(df)} gene rows) -> {cfg.operon_file}")
    print("  Step 3 complete.\n")


def step04_run_hmmsearch(cfg: Config, force: bool = False):
    """Run hmmsearch for every individual .hmm profile in the profiles directory."""
    # Without the marker search there is no marker set: stop rather than
    # letting later steps fall back to labelling every operon a "marker".
    if cfg.hmm_profiles_dir is None or not cfg.hmm_profiles_dir.is_dir():
        sys.exit(f"HMM profile directory not available ({cfg.hmm_profiles_dir}). Reinstall "
                 f"ProFinder to restore the bundled marker profiles, or pass --hmm-dir.")

    hmm_files = sorted(cfg.hmm_profiles_dir.glob("*.hmm"))
    if not hmm_files:
        sys.exit(f"No .hmm profiles found in {cfg.hmm_profiles_dir}; pass --hmm-dir with the "
                 f"marker profiles.")

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

    _require(4, (cfg.faa_file, 1))
    print(f"── Running hmmsearch ({len(hmm_files)} profiles) ──")
    cfg.hmm_dir.mkdir(parents=True, exist_ok=True)

    # Fail once, clearly, when hmmsearch itself is unavailable.
    probe = _run_tool(cfg.conda_env_hmm, f"{cfg.hmmsearch_bin} -h")
    if probe.returncode != 0:
        where = f" in conda env {cfg.conda_env_hmm}" if cfg.conda_env_hmm else ""
        sys.exit(f"hmmsearch is not runnable as '{cfg.hmmsearch_bin}'{where}: "
                 f"{(probe.stdout + probe.stderr).strip()[-300:] or 'no output'}\n"
                 f"Install HMMER (conda install -c conda-forge -c bioconda hmmer) or pass "
                 f"--hmmsearch / --conda-hmm.")
    with open(cfg.faa_file) as fh:
        n_proteins = sum(1 for line in fh if line.startswith(">"))

    q = shlex.quote
    failed = []
    for i, hmm_path in enumerate(hmm_files, 1):
        out_file = tblout_dir / (hmm_path.stem + ".tblout")
        log_file = tblout_dir / (hmm_path.stem + ".log")
        part = out_file.with_name(out_file.name + ".part")

        if not force and out_file.exists():
            continue

        if n_proteins == 0:
            # No proteins predicted: record an empty, completed search.
            part.write_text("# no proteins predicted in this genome\n# [ok]\n")
            os.replace(part, out_file)
            continue
        # Write to a .part file and rename only on success, so an interrupted
        # or failed search never looks complete to the checkpoint test.
        result = _run_tool(cfg.conda_env_hmm,
                           f"{cfg.hmmsearch_bin} --tblout {q(str(part))} -o {q(str(log_file))}"
                           f" --noali --cpu {cfg.threads} {q(str(hmm_path))} {q(str(cfg.faa_file))}")
        if result.returncode != 0:
            part.unlink(missing_ok=True)
            failed.append((hmm_path.name, (result.stderr or result.stdout).strip()[-200:]))
            print(f"  WARNING: hmmsearch failed for {hmm_path.name}: {failed[-1][1]}")
        else:
            os.replace(part, out_file)

        if i % 100 == 0 or i == len(hmm_files):
            print(f"  {i}/{len(hmm_files)} profiles searched...")

    if failed:
        sys.exit(f"hmmsearch failed for {len(failed)} of {len(hmm_files)} profiles (first: "
                 f"{failed[0][0]}: {failed[0][1] or 'no message'}). Per-profile logs are in "
                 f"{tblout_dir}/. Re-running the same command retries only the failed profiles.")

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
        sys.exit(f"Step 5 cannot run: no hmmsearch output in {tblout_dir} (written by step 4). "
                 f"Run step 4 first (e.g. --start 4).")

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
    s6_out, s6_in = [cfg.operon_filtered, cfg.operon_filtered_markers], [cfg.operon_file, cfg.hmm_filtered]
    if not force and _is_current(cfg, "step6", s6_out, s6_in):
        print("── Filtered operons with markers already exists, skipping ──")
        print("  Step 6 complete.\n")
        return

    _require(6, (cfg.operon_file, 3), (cfg.hmm_filtered, 5))
    _forget(cfg, "step6")
    print("── Filtering operons ──")
    try:
        df = pd.read_csv(cfg.operon_file, sep="\t")
    except pd.errors.EmptyDataError:                # header-less file from an older run
        df = pd.DataFrame(columns=_OPERON_COLUMNS)
    _warn_legacy_operons(cfg)
    if df.empty:
        print("  No operons to filter.")
        df.to_csv(cfg.operon_filtered, sep="\t", index=False)
        df.to_csv(cfg.operon_filtered_markers, sep="\t", index=False)
        _mark_current(cfg, "step6", s6_in)
        print("  Step 6 complete.\n")
        return

    # Keep operons on a single contig
    single_contig = df.groupby("Operon").filter(
        lambda g: g["SeqID"].nunique() == 1
    )
    # Step 3 puts every operon on one strand; drop a mixed-strand operon from a
    # hand-edited or older operons.tsv rather than misorder its genes.
    if "Strand" in single_contig.columns:
        single_contig = single_contig.groupby("Operon").filter(
            lambda g: g["Strand"].nunique() == 1
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
        # An empty hmm_filtered.tsv means step 5 was interrupted. Never fall
        # back to treating every operon as a marker operon.
        sys.exit(f"{cfg.hmm_filtered} is empty - re-run step 5 (e.g. --start 5 --force).")

    _mark_current(cfg, "step6", s6_in)
    print("  Step 6 complete.\n")


# ── Marker evidence: best match first ────────────────────────────────────
# Every gene with an HMM hit >= --hmm-bitscore is a marker candidate, but each
# profile's best-scoring gene is the likeliest true single-copy marker; lower-
# scoring hits are mostly paralogs. A promoter-orientation IGR counts as a
# marker region only if it can drive a hit gene: the gene immediately
# downstream of it (in the direction of transcription) or a later gene of the
# same operon. Labels, in marker-shortlist order:
#   best_match               the driven gene is some profile's best hit
#   best_match_via_operon    a later gene of the operon is some profile's best hit
#   weaker_match             the driven gene has only lower-scoring hits
#   weaker_match_via_operon  only lower-scoring hits lie later in the operon
MARKER_EVIDENCE_RANK = {"best_match": 1, "best_match_via_operon": 2,
                        "weaker_match": 3, "weaker_match_via_operon": 4}
MARKER_EVIDENCE_LABEL = {"best_match": "best match",
                         "best_match_via_operon": "best match (via operon)",
                         "weaker_match": "weaker match",
                         "weaker_match_via_operon": "weaker match (via operon)"}
_MARKER_EVIDENCE_COLS = ["marker_evidence", "marker_rank", "marker_gene",
                         "marker_profiles", "marker_bitscore", "marker_genes_downstream"]
_MARKER_BEST_COLUMNS = ["profile", "best_gene", "best_bitscore", "n_genes_hit",
                        "n_tied_best", "runner_up_gene", "runner_up_bitscore", "margin_bits"]
_MARKER_SUMMARY_COLUMNS = ["profile", "best_gene", "gene_name", "product", "best_bitscore",
                           "n_genes_hit", "n_tied_best", "runner_up_gene", "runner_up_bitscore",
                           "margin_bits", "promoter_id", "link", "genes_between",
                           "marker_evidence", "motif_confirmed", "motif_path",
                           "motif_shortlist_score", "note"]


def _hmm_hit_table(hmm):
    """(profile, gene, bit score) per hit, with ``is_best`` marking each
    profile's highest-scoring gene(s); tied genes are all kept."""
    cols = ["query_name", "target_name", "full_sequence_bitscore"]
    if hmm is None or len(hmm) == 0 or not set(cols) <= set(hmm.columns):
        return pd.DataFrame(columns=cols + ["is_best"])
    h = hmm[cols].copy()
    h["full_sequence_bitscore"] = pd.to_numeric(h["full_sequence_bitscore"], errors="coerce")
    h = h.dropna(subset=["full_sequence_bitscore"])
    if h.empty:
        return pd.DataFrame(columns=cols + ["is_best"])
    h = h.groupby(["query_name", "target_name"], as_index=False)["full_sequence_bitscore"].max()
    top = h.groupby("query_name")["full_sequence_bitscore"].transform("max")
    h["is_best"] = h["full_sequence_bitscore"] == top
    return h


def _marker_best_hits_table(h):
    """One row per profile with at least one hit: best gene(s), runner-up, margin."""
    rows = []
    for q, g in h.groupby("query_name", sort=True):
        g = g.sort_values(["full_sequence_bitscore", "target_name"], ascending=[False, True])
        best, rest = g[g["is_best"]], g[~g["is_best"]]
        b = float(best["full_sequence_bitscore"].iloc[0])
        r = float(rest["full_sequence_bitscore"].iloc[0]) if len(rest) else np.nan
        rows.append({"profile": q, "best_gene": ",".join(best["target_name"].astype(str)),
                     "best_bitscore": b, "n_genes_hit": int(len(g)), "n_tied_best": int(len(best)),
                     "runner_up_gene": str(rest["target_name"].iloc[0]) if len(rest) else "",
                     "runner_up_bitscore": r, "margin_bits": (b - r) if len(rest) else np.nan})
    return pd.DataFrame(rows, columns=_MARKER_BEST_COLUMNS)


def _operon_order(ops):
    """gene -> (operon, position, genes in transcription order). With a
    Strand column (operons.tsv from this version) a + strand operon runs left
    to right by Start and a - strand operon right to left by End. Without it,
    ``yes_1`` marks an operon transcribed from its leftmost gene and ``yes_2``
    one transcribed from its rightmost (step 6 drops operons carrying both)."""
    order = {}
    if ops is None or len(ops) == 0 or not {"Operon", "Gene", "Start", "FirstGene"} <= set(ops.columns):
        return order
    o = ops.drop_duplicates(["Operon", "Gene"])
    stranded = {"Strand", "End"} <= set(o.columns)
    for op, g in o.groupby("Operon", sort=False):
        strands = set(g["Strand"]) if stranded else set()
        if strands == {"+"}:
            genes = list(g.sort_values("Start", kind="mergesort")["Gene"])
        elif strands == {"-"}:
            genes = list(g.sort_values("End", ascending=False, kind="mergesort")["Gene"])
        else:
            fwd = bool((g["FirstGene"] == "yes_1").any())
            genes = list(g.sort_values("Start", ascending=fwd, kind="mergesort")["Gene"])
        for k, gene in enumerate(genes):
            order[gene] = (op, k, genes)
    return order


def _driven_label(d, order, best_prof, best_bits, any_prof, any_bits):
    """Marker evidence for a region whose downstream gene is ``d``:
    (label, marker gene, profiles, bit score, genes downstream) or None."""
    if not isinstance(d, str) or d not in order:
        return None
    _op, k, genes = order[d]
    later = genes[k + 1:]
    if d in best_prof:
        return ("best_match", d, best_prof[d], best_bits[d], 0)
    for j, g in enumerate(later, 1):
        if g in best_prof:
            return ("best_match_via_operon", g, best_prof[g], best_bits[g], j)
    if d in any_prof:
        return ("weaker_match", d, any_prof[d], any_bits[d], 0)
    for j, g in enumerate(later, 1):
        if g in any_prof:
            return ("weaker_match_via_operon", g, any_prof[g], any_bits[g], j)
    return None


def _marker_evidence_table(igr, ops, h):
    """One row per IGR that can drive a gene with an HMM hit (CO_F, CO_R, and
    DP regions for information), with its marker evidence."""
    order = _operon_order(ops)
    join = lambda s: ",".join(sorted(set(s)))
    hb = h[h["is_best"].astype(bool)] if len(h) else h
    best_prof = hb.groupby("target_name")["query_name"].agg(join).to_dict() if len(hb) else {}
    best_bits = hb.groupby("target_name")["full_sequence_bitscore"].max().to_dict() if len(hb) else {}
    any_prof = h.groupby("target_name")["query_name"].agg(join).to_dict() if len(h) else {}
    any_bits = h.groupby("target_name")["full_sequence_bitscore"].max().to_dict() if len(h) else {}
    recs = []
    for _, r in igr.iterrows():
        o = r.get("orientation")
        if o == "CO_F":
            cands = [r.get("right_gene")]
        elif o == "CO_R":
            cands = [r.get("left_gene")]
        elif o == "DP":
            cands = [r.get("left_gene"), r.get("right_gene")]
        else:
            continue
        labs = [x for x in (_driven_label(d, order, best_prof, best_bits, any_prof, any_bits)
                            for d in cands) if x]
        if not labs:
            continue
        lab = min(labs, key=lambda x: (MARKER_EVIDENCE_RANK[x[0]], x[4]))
        recs.append({"igr_id": r["igr_id"], "marker_evidence": lab[0],
                     "marker_rank": MARKER_EVIDENCE_RANK[lab[0]], "marker_gene": lab[1],
                     "marker_profiles": lab[2], "marker_bitscore": float(lab[3]),
                     "marker_genes_downstream": int(lab[4])})
    return pd.DataFrame(recs, columns=["igr_id"] + _MARKER_EVIDENCE_COLS)


def _has_column(path, col):
    """True if the TSV at ``path`` has a header containing ``col``."""
    try:
        with open(path) as fh:
            return col in fh.readline().rstrip("\n").split("\t")
    except OSError:
        return False


def _load_marker_evidence(cfg):
    """promoter_markers.tsv -> evidence columns indexed by igr_id (CO_F/CO_R)."""
    empty = pd.DataFrame(columns=_MARKER_EVIDENCE_COLS)
    try:
        m = pd.read_csv(cfg.promoter_markers, sep="\t")
    except (FileNotFoundError, pd.errors.EmptyDataError):
        return empty
    if m.empty or not {"igr_id", "orientation"} <= set(m.columns):
        return empty
    m = m[m["orientation"].isin(["CO_F", "CO_R"])].copy()
    if "marker_evidence" not in m.columns:          # written by an older version
        for c, v in zip(_MARKER_EVIDENCE_COLS, ["", 9, "", "", np.nan, np.nan]):
            m[c] = v
    return m.drop_duplicates("igr_id").set_index("igr_id")[_MARKER_EVIDENCE_COLS]


def _attach_marker_evidence(igr, ev):
    """Add marker_evidence / marker_gene / marker_profiles (blank when none)."""
    for c in ("marker_evidence", "marker_gene", "marker_profiles"):
        igr[c] = (igr["igr_id"].map(ev[c]).fillna("") if ev is not None and len(ev)
                  else "")
    return igr


def _marker_order(df, prefix=""):
    """Marker shortlist order: marker evidence (best match first), then the
    shortlist order within each evidence class."""
    if df is None or len(df) == 0:
        return df
    out = _shortlist_order(df, prefix=prefix)
    rank_col = out[prefix + "marker_rank"] if (prefix + "marker_rank") in out.columns else None
    if rank_col is None:
        return out
    rank = pd.to_numeric(rank_col, errors="coerce").fillna(9).to_numpy()
    return out.iloc[np.lexsort((np.arange(len(out)), rank))].reset_index(drop=True)


def _write_promoter_outputs(cfg, calls, n_prom=None):
    """Write the promoter FASTAs and the marker shortlist from per-IGR calls.
    all_promoters_verified.fasta is in shortlist order; the marker shortlist
    (motif_best_markers.tsv, promoter_markers_verified.tsv,
    marker_promoters_verified.fasta) is best-match first, then shortlist order."""
    have = calls is not None and len(calls) > 0
    ranked = _shortlist_order(calls) if have else calls
    if have:
        _write_fasta(ranked, cfg.all_promoters_verified_fasta, short_header=False)
    else:
        open(cfg.all_promoters_verified_fasta, "w").close()
    if n_prom is not None:
        print(f"  All promoters: {len(ranked) if have else 0}/{n_prom} intervals called (path A/B/C)")
    ev = _load_marker_evidence(cfg)
    if have and len(ev):
        mcalls = ranked[ranked["igr_id"].isin(ev.index)].copy()
        for c in _MARKER_EVIDENCE_COLS:
            mcalls[c] = mcalls["igr_id"].map(ev[c])
        mcalls = _marker_order(mcalls)
    else:
        base_cols = list(calls.columns) if calls is not None else []
        mcalls = pd.DataFrame(columns=base_cols + [c for c in _MARKER_EVIDENCE_COLS if c not in base_cols])
    mcalls.to_csv(cfg.motif_best_markers, sep="\t", index=False)
    mcalls.to_csv(cfg.promoter_markers_verified, sep="\t", index=False)
    if len(mcalls):
        _write_fasta(mcalls, cfg.marker_promoters_verified_fasta, short_header=False)
    else:
        open(cfg.marker_promoters_verified_fasta, "w").close()
    by = mcalls["marker_evidence"].value_counts().to_dict() if len(mcalls) else {}
    print(f"  Marker promoters: {len(mcalls)} called ("
          + ", ".join(f"{MARKER_EVIDENCE_LABEL[k]} {by.get(k, 0)}" for k in MARKER_EVIDENCE_RANK)
          + f") -> {cfg.promoter_markers_verified}")


def _gff_gene_names(gff_path):
    """locus tag / ID -> (gene name, product) from a Prokka GFF."""
    names = {}
    try:
        with open(gff_path) as fh:
            for line in fh:
                if line.startswith("##FASTA"):
                    break
                f = line.rstrip("\n").split("\t")
                if len(f) != 9 or f[2] in ("gene", "region"):
                    continue
                a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
                for key in (a.get("locus_tag"), a.get("ID")):   # operon gene IDs prefer ID=
                    if key and key not in names:
                        names[key] = (a.get("gene", ""), a.get("product", ""))
    except OSError:
        pass
    return names


def _lost_marker_note(cfg, gene, names, gene_op, op_labels, op_contigs,
                      op_strands=None, gene_strand=None, legacy=False):
    """Why a profile's best gene lies outside every operon kept by step 6.
    Step 3 places a + or - strand gene in a run of same-strand genes and
    drops the run only when a gene on the same strand lies within
    (--max-internal-dist, --min-flanking-dist) bp upstream of its first gene.
    Step 6 removes operons that span contigs, mix strands or have an
    ambiguous first gene; only a hand-edited or older operons.tsv has them."""
    if gene not in names:
        return "best gene not found in the annotation (GFF)"
    op = gene_op.get(gene)
    if op is None:
        if legacy:
            return ("best gene is in no operon (step 3; operons/operons.tsv was made by an "
                    "older ProFinder whose operons ignore strand)")
        if gene_strand is not None and gene_strand.get(gene) not in ("+", "-"):
            return "best gene has no strand in the annotation, so it is in no operon (step 3)"
        return (f"best gene is in no operon (step 3): a gene on the same strand lies "
                f"{cfg.max_internal_distance + 1}-{cfg.min_flanking_distance - 1} bp upstream "
                f"of its gene cluster, so its operon may start further upstream")
    labels = op_labels.get(op, set())
    if len(op_contigs.get(op, set())) > 1:
        why = "it spans contigs"
    elif op_strands is not None and len(op_strands.get(op, set())) > 1:
        why = "it has genes on both strands"
    elif {"yes_1", "yes_2"} <= labels:
        why = f"it joins convergent genes {cfg.max_internal_distance} bp or less apart"
    elif not labels & {"yes_1", "yes_2"}:
        why = f"it joins divergent genes {cfg.max_internal_distance} bp or less apart"
    else:
        return "best gene is not in a retained operon (step 6)"
    return f"best gene's operon was removed (step 6): {why}"


def _write_marker_summary(cfg, final):
    """marker_summary.tsv: one row per marker profile with its best-scoring
    gene, the nearest promoter-orientation region upstream of that gene in its
    operon, and that region's call; ``note`` says why no region was scanned."""
    profiles = []
    hdir = cfg.hmm_profiles_dir
    if hdir is not None and Path(hdir).is_dir():
        for f in sorted(Path(hdir).glob("*.hmm")):
            name = f.stem
            with open(f) as fh:
                for line in fh:
                    if line.startswith("NAME"):
                        name = line.split(None, 1)[1].strip()
                        break
            profiles.append(name)
    try:
        best = pd.read_csv(cfg.marker_best_hits, sep="\t")
    except (FileNotFoundError, pd.errors.EmptyDataError):
        best = pd.DataFrame(columns=_MARKER_BEST_COLUMNS)
    try:
        ops = pd.read_csv(cfg.operon_filtered_markers, sep="\t")
    except (FileNotFoundError, pd.errors.EmptyDataError):
        ops = pd.DataFrame()
    try:
        allops = pd.read_csv(cfg.operon_file, sep="\t")
    except (FileNotFoundError, pd.errors.EmptyDataError):
        allops = pd.DataFrame(columns=_OPERON_COLUMNS)
    legacy = "Strand" not in allops.columns
    gene_op, op_labels, op_contigs, op_strands = {}, {}, {}, {}
    for o, gene, lab, sid in zip(allops["Operon"], allops["Gene"], allops["FirstGene"], allops["SeqID"]):
        gene_op[gene] = o
        op_labels.setdefault(o, set()).add(lab)
        op_contigs.setdefault(o, set()).add(sid)
    if not legacy:
        for o, s in zip(allops["Operon"], allops["Strand"]):
            op_strands.setdefault(o, set()).add(s)
    gene_strand = {}
    if cfg.gff_file.exists():
        gdf = _parse_gff_for_operons(str(cfg.gff_file))
        for a, s in zip(gdf.get("attributes", []), gdf.get("strand", [])):
            gene_strand.setdefault(_extract_gene_id(a), s)
    try:
        igr = pd.read_csv(cfg.igr_summary, sep="\t")
    except (FileNotFoundError, pd.errors.EmptyDataError):
        igr = pd.DataFrame(columns=_IGR_COLUMNS)
    order = _operon_order(ops)
    drive = {}
    for _, r in igr.iterrows():
        if r.get("orientation") == "CO_F":
            drive[r.get("right_gene")] = r["igr_id"]
        elif r.get("orientation") == "CO_R":
            drive[r.get("left_gene")] = r["igr_id"]
    dp = igr[igr["orientation"] == "DP"] if "orientation" in igr.columns else igr.iloc[0:0]
    dp_flank = set(dp.get("left_gene", pd.Series(dtype=str))) | set(dp.get("right_gene", pd.Series(dtype=str)))
    names = _gff_gene_names(cfg.gff_file)
    fin = final.set_index("promoter_id") if final is not None and len(final) else pd.DataFrame()
    bt = best.set_index("profile") if len(best) else best
    rows = []
    for p in sorted(set(profiles) | set(best["profile"].astype(str))):
        rec = {c: "" for c in _MARKER_SUMMARY_COLUMNS}
        rec["profile"] = p
        if p not in getattr(bt, "index", []):
            rec["note"] = f"no hit >= {cfg.hmm_bitscore_min:g} bits"
            rows.append(rec)
            continue
        b = bt.loc[p]
        g = str(b["best_gene"]).split(",")[0]
        gn, pr = names.get(g, ("", ""))
        rec.update(best_gene=b["best_gene"], gene_name=gn, product=pr,
                   best_bitscore=b["best_bitscore"], n_genes_hit=int(b["n_genes_hit"]),
                   n_tied_best=int(b["n_tied_best"]), runner_up_gene=b["runner_up_gene"],
                   runner_up_bitscore=b["runner_up_bitscore"], margin_bits=b["margin_bits"])
        if g not in order:
            rec["note"] = _lost_marker_note(cfg, g, names, gene_op, op_labels, op_contigs,
                                            op_strands, gene_strand, legacy)
        else:
            _op, k, og = order[g]
            link = next(((drive[og[j]], k - j) for j in range(k, -1, -1) if og[j] in drive), None)
            if link is None:
                rec["note"] = ("upstream region is divergent (not scanned)" if og[0] in dp_flank
                               else "no intergenic region upstream within the length limits")
            else:
                pid, between = link
                rec.update(promoter_id=pid, link="direct" if between == 0 else "via_operon",
                           genes_between=int(between))
                if len(fin) and pid in fin.index:
                    row = fin.loc[pid]
                    rec.update(marker_evidence=row.get("marker_evidence", ""),
                               motif_confirmed=bool(row.get("motif_confirmed", False)),
                               motif_path=row.get("motif_path", ""),
                               motif_shortlist_score=row.get("motif_shortlist_score", ""))
                    if not bool(row.get("motif_confirmed", False)):
                        rec["note"] = "region scanned; no -10/-35 call at the current thresholds"
        rows.append(rec)
    pd.DataFrame(rows, columns=_MARKER_SUMMARY_COLUMNS).to_csv(cfg.marker_summary, sep="\t", index=False)
    n_called = sum(1 for r in rows if r.get("motif_confirmed") is True)
    print(f"  Marker summary ({len(rows)} profiles; {n_called} with a called promoter) -> {cfg.marker_summary}")


def step07_match_igrs_to_markers(cfg: Config, force: bool = False):
    """Match IGRs to marker operon genes."""
    s7_out = [cfg.promoter_markers, cfg.promoter_markers_hmm, cfg.marker_best_hits]
    s7_in = [cfg.igr_summary, cfg.operon_filtered_markers, cfg.hmm_filtered]
    if not force and _is_current(cfg, "step7", s7_out, s7_in):
        print("── Promoter markers file already exists, skipping ──")
        print("  Step 7 complete.\n")
        return

    _require(7, (cfg.igr_summary, 2), (cfg.operon_filtered_markers, 6), (cfg.hmm_filtered, 5))
    _forget(cfg, "step7")
    _warn_legacy_operons(cfg)
    print("── Matching IGRs to marker operons ──")
    try:
        igr = pd.read_csv(cfg.igr_summary, sep="\t")
    except pd.errors.EmptyDataError:                # header-less file from an older run
        igr = pd.DataFrame(columns=_IGR_COLUMNS)
    try:
        markers = pd.read_csv(cfg.operon_filtered_markers, sep="\t")
    except (FileNotFoundError, pd.errors.EmptyDataError):
        print("  No marker data available.")
        pd.DataFrame().to_csv(cfg.promoter_markers, sep="\t", index=False)
        pd.DataFrame().to_csv(cfg.promoter_markers_hmm, sep="\t", index=False)
        pd.DataFrame(columns=_MARKER_BEST_COLUMNS).to_csv(cfg.marker_best_hits, sep="\t", index=False)
        _mark_current(cfg, "step7", s7_in)
        print("  Step 7 complete.\n")
        return

    try:
        hmm = pd.read_csv(cfg.hmm_filtered, sep="\t")
    except pd.errors.EmptyDataError:
        hmm = pd.DataFrame()
    h = _hmm_hit_table(hmm)
    best_tab = _marker_best_hits_table(h)
    best_tab.to_csv(cfg.marker_best_hits, sep="\t", index=False)
    print(f"  Best-scoring gene per profile ({len(best_tab)} profiles with a hit) -> {cfg.marker_best_hits}")

    # A region is a marker region only if it can drive a gene with an HMM hit:
    # the gene immediately downstream of it or a later gene of the same operon.
    # Regions that lie downstream of every hit gene in an operon cannot drive
    # any of them and are not marker regions.
    ev = _marker_evidence_table(igr, markers, h)
    igr = igr.copy()
    igr["marker_match"] = "none"
    hit_ids = set(ev["igr_id"])
    for o in ("CO_F", "CO_R", "DP"):
        igr.loc[(igr["orientation"] == o) & igr["igr_id"].isin(hit_ids), "marker_match"] = o
    matched = igr[igr["marker_match"] != "none"].merge(ev, on="igr_id", how="left")
    matched.to_csv(cfg.promoter_markers, sep="\t", index=False)
    by = matched["marker_evidence"].value_counts().to_dict() if len(matched) else {}
    print(f"  Matched promoter markers ({len(matched)} rows: "
          + ", ".join(f"{MARKER_EVIDENCE_LABEL[k]} {by.get(k, 0)}" for k in MARKER_EVIDENCE_RANK)
          + f") -> {cfg.promoter_markers}")

    # Expanded variant: one row per marker region x HMM profile hit by the
    # DRIVEN gene (the gene immediately downstream), flagged is_best_hit when
    # that gene is the profile's best-scoring hit. Regions labelled via_operon
    # have no direct hit and are absent here.
    matched["driven_gene"] = np.where(matched["marker_match"] == "CO_F",
                                      matched["right_gene"], matched["left_gene"])
    if h.empty:
        matched["hmm_profile"] = "no_hmm"
        matched["is_best_hit"] = False
        expanded = matched
    else:
        hp = h.rename(columns={"target_name": "driven_gene", "query_name": "hmm_profile",
                               "is_best": "is_best_hit"})[["driven_gene", "hmm_profile", "is_best_hit"]]
        expanded = pd.merge(matched, hp, how="inner", on="driven_gene")
    expanded = expanded.drop(columns=["driven_gene"], errors="ignore")
    expanded.to_csv(cfg.promoter_markers_hmm, sep="\t", index=False)
    print(f"  Expanded markers with HMM profiles ({len(expanded)} rows) "
          f"-> {cfg.promoter_markers_hmm}")
    _mark_current(cfg, "step7", s7_in)
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
    # return {rec.id: str(rec.seq) for rec in SeqIO.parse(str(fasta_path), "fasta")}
    with _open_maybe_gzip(fasta_path) as fh:
        return {rec.id: str(rec.seq) for rec in SeqIO.parse(fh, "fasta")}
    


def _extract_cds_start(row, contigs, n_bp):
    """Return the first *n_bp* nucleotides of the CDS downstream of a
    promoter, oriented 5'→3' on the coding strand.
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
#  Each -10 hit is classified into one of three paths. Paths A–C
#  require a -35 hit with a 16-18 bp spacer to the -10 hit (canonical
#  σ⁷⁰ spacer is 17 ± 1). -10
#  motifs whose right edge sits more than cfg.max_dist_to_cds_start
#  bp from the downstream CDS start (default: see Config) are dropped
#  before the spacer search.
#    A — linked -10 and strict -35 (same subgroup)
#    B — extended -10 ("TGn" immediately 5' of the hexamer) and a
#        linked relaxed -35 (same subgroup)
#    C — unlinked -10 and strict -35 (different subgroups)
#  Anything else is regarded as no hit.
#
#  Each element carries an empirical quantile (emp_p_10 / emp_p_35)
#  read from the per-genome dinucleotide-shuffle null. Within an IGR
#  the reported call minimises (path rank, -[-log10(p10) - log10(p35)]).
#  Output tables and the HTML report are ordered by shortlist_score
#  (see _shortlist_order).
# =====================================================================

def step08_scan_motifs(cfg: Config, force: bool = False):
    """Score promoter-orientation IGRs with the interval-level empirical
    method (replaces the analytic per-position p-value).
    """
    from . import empirical as emp

    s8_out = [cfg.motif_best_markers, cfg.promoter_markers_verified,
              cfg.all_promoters_verified_fasta, cfg.marker_promoters_verified_fasta]
    s8_in = [cfg.promoter_markers, cfg.motif_best_all]
    if not force and cfg.motif_best_all.exists() and cfg.promoter_markers_verified.exists():
        if _is_current(cfg, "step8_markers", s8_out, s8_in):
            print("── Empirical motif scan already done, skipping ──")
        elif not _has_column(cfg.motif_best_all, "igr_id"):
            # The scan found no promoter-orientation IGRs: nothing to rebuild.
            print("── Empirical motif scan already done, skipping ──")
            _mark_current(cfg, "step8_markers", s8_in)
        elif not _has_column(cfg.promoter_markers, "marker_evidence"):
            print("── Empirical motif scan already done, skipping ──")
            print("  WARNING: operons/promoter_markers.tsv has no marker evidence (made by an "
                  "older ProFinder); re-run from step 7 to rebuild the marker shortlist.")
        else:
            # The marker regions changed since the shortlist was written (or an
            # older version wrote it): rebuild the shortlist outputs from the
            # existing calls; no rescan needed.
            _warn_legacy_operons(cfg)
            print("── Empirical motif scan already done; rebuilding the marker shortlist "
                  "from the current marker regions ──")
            _forget(cfg, "step8_markers")
            # round_trip parsing: the rebuilt shortlist must carry the exact
            # floats a full scan writes, not the default parser's rounding
            _write_promoter_outputs(cfg, pd.read_csv(cfg.motif_best_all, sep="\t",
                                                     float_precision="round_trip"))
            _mark_current(cfg, "step8_markers", s8_in)
        print("  Step 8 complete.\n")
        return

    _require(8, (cfg.igr_summary, 2), (cfg.promoter_markers, 7))
    _forget(cfg, "step8_markers")
    meme_file = cfg.meme_file
    if meme_file is None or not Path(meme_file).is_file():
        bundled = Path(__file__).parent / "all_unique_subgroups.meme"
        meme_file = bundled if bundled.is_file() else None
    if meme_file is None:
        sys.exit("Motif file not found: reinstall ProFinder to restore the bundled "
                 "all_unique_subgroups.meme, or pass --meme-file.")

    bg_source = cfg.fna_file if cfg.fna_file.exists() else cfg.input_fasta
    bg = _compute_background_from_fasta(bg_source)
    genome_gc = float(bg[1] + bg[2])
    subs, m10f, m35f = emp.paired_subgroups(emp.parse_meme(str(meme_file)))
    pwm10 = emp.freq_to_logodds(m10f, bg)
    pwm35 = emp.freq_to_logodds(m35f, bg)
    print(f"── Empirical scan: {len(subs)} subgroups; background GC={genome_gc:.3f}; "
          f"thresholds p10={cfg.motif_p10} p35={cfg.motif_p35} "
          f"p35_relaxed={cfg.motif_p35_relaxed}; {cfg.n_shuffles} shuffles/interval ──")
    print("  note: thresholds are applied to the per-architecture p-value, but the REPORTED "
          "emp_p_35 for Path C is")
    print(f"    Šidák-adjusted over the {len(subs)-1} other architectures it searches, so a "
          f"Path C call passing p35={cfg.motif_p35:g}")
    print(f"    is reported at up to "
          f"{emp.sidak(cfg.motif_p35, len(subs)-1):.3g}. That is not an inconsistency: it is "
          f"the same evidence on a")
    print("    single-test scale, comparable with Paths A and B. "
          "emp_p_35_uncorrected keeps the raw value.")

    if cfg.motif_p35_relaxed <= cfg.motif_p35:
        print(f"  WARNING: p35_relaxed ({cfg.motif_p35_relaxed}) <= p35 "
              f"({cfg.motif_p35}). Path B (extended -10 with a relaxed same-architecture "
              f"-35) loses the weak -35 matches it was added to capture: with the relaxed "
              f"threshold no looser than the strict one, a qualifying -35 is normally claimed "
              f"by Path A first. Path B is not silenced outright — it still fires where a "
              f"different architecture outscores the -10's own at the -35 site, which routes "
              f"the candidate past Path A — but most calls will be split across Paths A and "
              f"C. Set p35_relaxed > p35 for the intended behaviour.")

    cal = emp.Calibration()
    cal_path = cfg.ppv_calibration
    if cal_path is None:
        bundled_cal = Path(__file__).parent / "ppv_calibration.json"
        cal_path = str(bundled_cal) if bundled_cal.is_file() else None
    if cal_path is not None and Path(cal_path).is_file():
        cal = emp.Calibration.from_json(str(cal_path))
        print(f"  calibration: {Path(cal_path).name} "
              f"(paths {sorted(cal.per_path)}; PPV against prior={cfg.ppv_prior})")
    else:
        print("  calibration: none bundled — reporting score + path structure; "
              "lr_pos / ppv_interval left blank (not fabricated).")
        print("    Blank is deliberate: a fit on 49 held-out genomes was rejected because "
              "within a genome the score's")
        print("    AUC against TSS truth is 0.5007 (Wilcoxon P = 0.74) — it does not rank "
              "true against false among")
        print("    the calls. Use the path (A/B/C) as the confidence statement, not the "
              "score within a path.")

    try:
        igr = pd.read_csv(cfg.igr_summary, sep="\t")
    except (FileNotFoundError, pd.errors.EmptyDataError):
        igr = pd.DataFrame()
    n_prom = int(igr["orientation"].isin(["CO_F", "CO_R"]).sum()) if not igr.empty else 0
    if n_prom == 0:
        print("── No promoter-orientation IGRs to scan ──")
        for p in [cfg.motif_best_all, cfg.motif_best_markers, cfg.promoter_markers_verified]:
            pd.DataFrame().to_csv(p, sep="\t", index=False)
        open(cfg.all_promoters_verified_fasta, "w").close()
        open(cfg.marker_promoters_verified_fasta, "w").close()
        _mark_current(cfg, "step8_markers", s8_in)
        print("  Step 8 complete.\n")
        return

    calls = emp.scan_igrs(igr, pwm10, pwm35, subs,
                          p10=cfg.motif_p10, p35=cfg.motif_p35,
                          p35_relaxed=cfg.motif_p35_relaxed,
                          max_dist=cfg.max_dist_to_cds_start,
                          n_shuffles=cfg.n_shuffles, seed=0,
                          calibration=(cal if cal.calibrated else None),
                          prior=cfg.ppv_prior, genome_gc=genome_gc)
    if not calls.empty and not cal.calibrated:
        for c in ("lr_pos", "ppv_interval", "in_domain"):
            calls[c] = ""

    calls.to_csv(cfg.motif_best_all, sep="\t", index=False)
    _write_promoter_outputs(cfg, calls, n_prom)
    _mark_current(cfg, "step8_markers", s8_in)
    print("  Step 8 complete.\n")


def _logspace_p_values(p_min, p_max, n):
    """Return `n` p-values log-spaced inclusive of both endpoints. For
    `n == 1`, returns just `p_max`."""
    if n <= 1:
        return [p_max]
    log_min = math.log10(p_min)
    log_max = math.log10(p_max)
    step = (log_max - log_min) / (n - 1)
    return [10.0 ** (log_min + i * step) for i in range(n)]

_SCAN_MOTIF_COLS = [
    "path", "source_10", "source_35", "has_ext10",
    "emp_p_10", "emp_p_35", "emp_p_35_uncorrected", "p35_multiplicity_k",
    "seq_10", "seq_35",
    "pos_10", "pos_35", "spacer", "score_std", "n_pairs",
    "shortlist_score",
    "lr_pos", "ppv_interval", "in_domain",
]
_FINAL_TABLE_COLUMNS = [
    "igr_id", "contig", "start", "end", "length", "orientation",
    "associated_cds", "gene_name", "locus_tag", "product",
    "is_marker", "marker_evidence", "marker_gene", "marker_profiles",
    "motif_confirmed", "motif_path",
    "motif_source_10", "motif_source_35", "motif_has_ext10",
    "motif_emp_p_10", "motif_emp_p_35",
    "motif_emp_p_35_uncorrected", "motif_p35_multiplicity_k",
    "motif_seq_10", "motif_seq_35",
    "motif_pos_10", "motif_pos_35", "motif_spacer",
    "motif_score_std", "motif_n_pairs", "motif_shortlist_score",
    "motif_lr_pos", "motif_ppv_interval", "motif_in_domain",
    "sequence_5p_to_3p",
]


def _shortlist_order(df, prefix=""):
    """Order a calls/results table for presentation.
    """
    if df is None or len(df) == 0:
        return df
    n = len(df)
    def col(name, default):
        return df[prefix + name] if (prefix + name) in df.columns else pd.Series(
            [default] * n, index=df.index)
    path = col("path", "").astype(str).str.strip()
    called = path.isin(["A", "B", "C"]).to_numpy()
    prank = path.map({"A": 0, "B": 1, "C": 2}).fillna(3).to_numpy(dtype=float)
    score = pd.to_numeric(col("shortlist_score", np.nan), errors="coerce").to_numpy(dtype=float)
    p10 = pd.to_numeric(col("emp_p_10", np.nan), errors="coerce").to_numpy(dtype=float)
    p35 = pd.to_numeric(col("emp_p_35", np.nan), errors="coerce").to_numpy(dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        comb = -np.log10(np.clip(p10, 1e-300, None)) - np.log10(np.clip(p35, 1e-300, None))
    comb = np.where(np.isfinite(comb), comb, -np.inf)
    primary = (~called).astype(int)                      # called rows first
    if np.isfinite(score).any():
        k2 = -np.where(np.isfinite(score), score, -np.inf)   # score desc
        k3 = np.zeros(n)
    else:
        k2 = prank
        k3 = -comb
    order = np.lexsort((np.arange(n), k3, k2, primary))
    return df.iloc[order].reset_index(drop=True)


def _assemble_final_table(igr_df, calls_df, marker_ids, ann_map,
                          cds_bp=0, contigs=None, marker_ev=None):
    """Build a ``profinder_results.tsv``-shaped DataFrame (the canonical
    ``_FINAL_TABLE_COLUMNS`` layout) from pre-loaded IGRs and a per-IGR calls
    DataFrame (the empirical scan output — ``scan_igrs`` / ``motif_best_all``
    columns). Shared by ``step10_final_table`` and the p-sweep diagnostic so
    both emit one identical schema.

    ``igr_df`` needs ``igr_id, contig, start, end, length, orientation,
    left_gene, right_gene, sequence`` (``associated_cds`` and
    ``sequence_5p_to_3p`` are derived if absent). ``ann_map`` maps gene IDs to
    ``{product, gene, locus_tag}``. ``marker_ids`` flags ``is_marker``.
    ``calls_df`` is one row per motif-confirmed IGR (empty / None -> nothing
    confirmed). ``motif_*`` cells are blank for unconfirmed IGRs, and
    ``motif_lr_pos`` / ``motif_ppv_interval`` are blank unless the scan ran
    with a calibration (never fabricated)."""
    igr = igr_df.copy()
    if "associated_cds" not in igr.columns:
        igr["associated_cds"] = igr.apply(_get_associated_gene, axis=1)
    if "sequence_5p_to_3p" not in igr.columns:
        igr["sequence_5p_to_3p"] = igr.apply(
            lambda r: _reverse_complement(r["sequence"])
            if r["orientation"] == "CO_R" else r["sequence"], axis=1)
    igr["gene_name"] = igr["associated_cds"].map(
        lambda g: (ann_map or {}).get(g, {}).get("gene", ""))
    igr["locus_tag"] = igr["associated_cds"].map(
        lambda g: (ann_map or {}).get(g, {}).get("locus_tag", ""))
    igr["product"] = igr["associated_cds"].map(
        lambda g: (ann_map or {}).get(g, {}).get("product", ""))
    igr["is_marker"] = (igr["igr_id"].isin(marker_ids) if marker_ids else False)
    igr = _attach_marker_evidence(igr, marker_ev)

    motif_map = {}
    if calls_df is not None and len(calls_df) > 0:
        for _, row in calls_df.iterrows():
            motif_map[row["igr_id"]] = {c: row.get(c, "") for c in _SCAN_MOTIF_COLS}
    igr["motif_confirmed"] = igr["igr_id"].isin(set(motif_map.keys()))
    for c in _SCAN_MOTIF_COLS:
        igr[f"motif_{c}"] = igr["igr_id"].map(
            lambda g, _c=c: motif_map.get(g, {}).get(_c, ""))

    columns = list(_FINAL_TABLE_COLUMNS)
    if cds_bp > 0 and contigs is not None:
        igr = _add_cds_column(igr, contigs, cds_bp)
        columns.append("sequence_5p_to_3p_cds")
    out = igr[columns].copy()
    out.rename(columns={"igr_id": "promoter_id"}, inplace=True)
    return _shortlist_order(out, prefix="motif_")


def _load_annotation_map(annotations_tsv):
    """Read ``cds_annotations.tsv`` into a {gene_id -> {gene, locus_tag,
    product}} dict. Returns empty dict if the file is missing or empty."""
    if not annotations_tsv.exists():
        return {}
    try:
        df = pd.read_csv(annotations_tsv, sep="\t")
    except (pd.errors.EmptyDataError, KeyError):
        return {}
    out = {}
    for _, row in df.iterrows():
        out[row["gene_id"]] = {
            "product": str(row.get("product", "")),
            "gene": str(row.get("gene", "")),
            "locus_tag": str(row.get("locus_tag", "")),
        }
    return out


def _sweep_table_filename(mode, p10, p35, p35_relaxed):
    """Return a deterministic filename for a sweep-point TSV that
    encodes mode + (p10, p35, p35_relaxed). Filename-safe and round-
    trips through ``str.split``: e.g.
    ``sweep_matched__p10_1.00e-03__p35_1.00e-03__p35r_2.00e-02.tsv``."""
    return (
        f"sweep_{mode}"
        f"__p10_{p10:.2e}"
        f"__p35_{p35:.2e}"
        f"__p35r_{p35_relaxed:.2e}"
        ".tsv"
    )


def _scan_counts(best_hits_df, marker_ids):
    """Summarize a best-hits DataFrame into per-path counts, both
    overall and restricted to marker-anchored IGRs."""
    out = {"n_total": 0, "n_path_A": 0, "n_path_B": 0,
           "n_path_C": 0,
           "n_marker": 0, "n_marker_path_A": 0, "n_marker_path_B": 0,
           "n_marker_path_C": 0}
    if best_hits_df is None or len(best_hits_df) == 0:
        return out
    out["n_total"] = len(best_hits_df)
    path_counts = best_hits_df["path"].value_counts().to_dict()
    for k in ("A", "B", "C"):
        out[f"n_path_{k}"] = int(path_counts.get(k, 0))
    if marker_ids:
        marker_mask = best_hits_df["igr_id"].isin(marker_ids)
        marker_df = best_hits_df[marker_mask]
        out["n_marker"] = int(len(marker_df))
        mpath = marker_df["path"].value_counts().to_dict()
        for k in ("A", "B", "C"):
            out[f"n_marker_path_{k}"] = int(mpath.get(k, 0))
    return out


def _render_p_sweep_figure(sweep_df, out_path):
    """Optional 3-panel figure: one subplot per mode, x = swept p,
    y = stacked counts of marker-anchored hits by path (A / B / C).
    Returns False (and skips) if matplotlib is not available."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"  matplotlib unavailable ({e}); skipping figure")
        return False

    modes = list(sweep_df["mode"].drop_duplicates())
    n_modes = len(modes)
    fig, axes = plt.subplots(1, max(1, n_modes),
                              figsize=(4.2 * max(1, n_modes), 3.5),
                              sharey=True, squeeze=False)
    axes = axes[0]
    path_colors = {"A": "#2e7d32", "B": "#e65100",
                    "C": "#c62828"}
    for ax, mode in zip(axes, modes):
        sub = sweep_df[sweep_df["mode"] == mode].sort_values("sweep_p")
        if len(sub) == 0:
            ax.set_title(f"{mode} (no data)")
            continue
        x = sub["sweep_p"].values
        bottom = [0] * len(sub)
        for path in ("A", "B", "C"):
            y = sub[f"n_marker_path_{path}"].values
            ax.bar(range(len(sub)), y, bottom=bottom,
                   color=path_colors[path], label=f"path {path}",
                   width=0.85, edgecolor="white", linewidth=0.5)
            bottom = [b + v for b, v in zip(bottom, y)]
        ax.set_xticks(range(len(sub)))
        ax.set_xticklabels([f"{p:.1e}" for p in x],
                            rotation=45, ha="right", fontsize=7)
        ax.set_title(mode, fontsize=10)
        ax.set_xlabel("swept p-value")
        ax.set_xscale("linear")
    axes[0].set_ylabel("marker-anchored hits (count)")
    axes[-1].legend(loc="upper right", fontsize=8, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  Figure -> {out_path}")
    return True


def diagnose_p_sweep(cfg: Config, sweep_min=1e-5, sweep_max=1e-1,
                      sweep_n=13, modes=("matched", "p10_only", "p35_only"),
                      p35_relaxed_ratio=20.0, force: bool = False):
    """Sweep the empirical element thresholds (as null quantiles) and record
    per-path promoter counts at each sweep point, using the SAME per-genome
    dinucleotide-shuffle null as step 8. The expensive null is built ONCE and
    re-thresholded at every sweep point (no re-shuffling). Reuses cached
    outputs from steps 2 (igr_summary.tsv) and 7 (promoter_markers.tsv).
    Writes ``diagnostics/p_sweep.tsv``, an optional 3-panel figure
    ``diagnostics/p_sweep.png``, and one profinder_results.tsv-shaped table
    per sweep point (canonical schema, filtered to motif_confirmed).
    Outputs are deterministic for the same inputs and CLI flags."""
    import time
    from . import empirical as emp

    out_tsv = cfg.diagnostics_dir / "p_sweep.tsv"
    out_png = cfg.diagnostics_dir / "p_sweep.png"
    bundled_meme = Path(__file__).parent / "all_unique_subgroups.meme"
    ps_in = [cfg.igr_summary, cfg.promoter_markers, cfg.cds_annotations,
             cfg.fna_file if cfg.fna_file.exists() else cfg.input_fasta,
             cfg.meme_file if cfg.meme_file is not None and Path(cfg.meme_file).is_file()
             else bundled_meme]
    ps_par = {"sweep": [sweep_min, sweep_max, sweep_n, list(modes), p35_relaxed_ratio],
              "cfg": [cfg.motif_p10, cfg.motif_p35, cfg.n_shuffles,
                      cfg.max_dist_to_cds_start, cfg.cds_bp]}
    if not force and _is_current(cfg, "p_sweep", [out_tsv], ps_in, ps_par):
        print("── p-value sweep diagnostic already exists, skipping ──")
        print(f"  -> {out_tsv}")
        print("  Diagnostic complete.\n")
        return
    _forget(cfg, "p_sweep")

    cfg.diagnostics_dir.mkdir(parents=True, exist_ok=True)
    tables_dir = cfg.diagnostics_dir / "sweep_tables"
    tables_dir.mkdir(parents=True, exist_ok=True)

    meme_file = cfg.meme_file
    if meme_file is None or not Path(meme_file).is_file():
        bundled = Path(__file__).parent / "all_unique_subgroups.meme"
        if bundled.is_file():
            meme_file = bundled
        else:
            sys.exit("Motif file not found: reinstall ProFinder to restore the bundled "
                     "all_unique_subgroups.meme, or pass --meme-file.")
    _require("p-sweep", (cfg.igr_summary, 2), (cfg.promoter_markers, 7))

    try:
        all_igr = pd.read_csv(cfg.igr_summary, sep="\t")
        all_igr = all_igr[all_igr["orientation"].isin(["CO_F", "CO_R"])].copy()
    except (FileNotFoundError, pd.errors.EmptyDataError):
        all_igr = pd.DataFrame()
    if all_igr.empty:
        print("── No promoter-orientation IGRs to scan; run steps 1-2 first ──")
        return

    if "sequence_5p_to_3p" not in all_igr.columns:
        all_igr["sequence_5p_to_3p"] = all_igr.apply(
            lambda r: _reverse_complement(r["sequence"])
            if r["orientation"] == "CO_R" else r["sequence"], axis=1)
    if "associated_cds" not in all_igr.columns:
        all_igr["associated_cds"] = all_igr.apply(_get_associated_gene, axis=1)

    marker_ids = set()
    try:
        mdf = pd.read_csv(cfg.promoter_markers, sep="\t")
        mdf = mdf[mdf["orientation"].isin(["CO_F", "CO_R"])]
        marker_ids = set(mdf["igr_id"])
    except (FileNotFoundError, pd.errors.EmptyDataError):
        print("  WARNING: no marker IGR file (run step 7 first); "
              "marker counts will be 0.")

    ann_map = _load_annotation_map(cfg.cds_annotations)
    contigs = (_load_contigs(cfg.input_fasta) if cfg.cds_bp > 0 else None)

    # Per-genome background + PWMs (identical to step 8), computed once.
    bg_source = cfg.fna_file if cfg.fna_file.exists() else cfg.input_fasta
    bg = _compute_background_from_fasta(bg_source)
    genome_gc = float(bg[1] + bg[2])
    subs, m10f, m35f = emp.paired_subgroups(emp.parse_meme(str(meme_file)))
    pwm10 = emp.freq_to_logodds(m10f, bg)
    pwm35 = emp.freq_to_logodds(m35f, bg)

    sweep_ps = _logspace_p_values(sweep_min, sweep_max, sweep_n)
    print("── p-value sweep diagnostic (empirical null) ──")
    print(f"  IGRs to scan: {len(all_igr)}  (markers: {len(marker_ids)})")
    marker_ev = _load_marker_evidence(cfg)
    print(f"  Background ({Path(bg_source).name}): "
          f"A={bg[0]:.3f} C={bg[1]:.3f} G={bg[2]:.3f} T={bg[3]:.3f}")
    print(f"  Sweep: {sweep_n} log-spaced points in "
          f"[{sweep_min:.1e}, {sweep_max:.1e}]")
    print(f"  Modes: {', '.join(modes)}")
    print(f"  p35_relaxed = min(1.0, {p35_relaxed_ratio} x p35)")
    print(f"  max_dist_to_cds_start = {cfg.max_dist_to_cds_start} bp; "
          f"{cfg.n_shuffles} shuffles/interval")
    print(f"  Per-sweep-point tables (motif_confirmed) -> {tables_dir}/")

    t0 = time.perf_counter()
    co, seqs, strata, engines = emp.build_scan_engine(
        all_igr, pwm10, pwm35, subs, max_dist=cfg.max_dist_to_cds_start,
        n_shuffles=cfg.n_shuffles, seed=0)
    if co.empty:
        print("── No promoter-orientation IGRs after engine build ──")
        return
    print(f"  Built shuffle null once ({len(engines)} strata) in "
          f"{time.perf_counter() - t0:.0f}s; re-thresholding per sweep point.")

    records = []
    for mode in modes:
        print(f"\n  -- mode: {mode} --")
        for p in sweep_ps:
            if mode == "matched":
                p10, p35 = p, p
            elif mode == "p10_only":
                p10, p35 = p, cfg.motif_p35
            elif mode == "p35_only":
                p10, p35 = cfg.motif_p10, p
            else:
                print(f"    unknown mode {mode!r}, skipping")
                continue
            p35_relaxed = min(1.0, p35 * p35_relaxed_ratio)

            t0 = time.perf_counter()
            calls = emp.reclassify(co, seqs, engines, p10=p10, p35=p35,
                                   p35_relaxed=p35_relaxed, genome_gc=genome_gc)
            elapsed = time.perf_counter() - t0
            counts = _scan_counts(calls, marker_ids)

            final = _assemble_final_table(all_igr, calls, marker_ids, ann_map,
                                          cds_bp=cfg.cds_bp, contigs=contigs,
                                          marker_ev=marker_ev)
            final_conf = final[final["motif_confirmed"]].copy()
            sweep_table_path = (tables_dir /
                                _sweep_table_filename(mode, p10, p35, p35_relaxed))
            final_conf.to_csv(sweep_table_path, sep="\t", index=False)

            row = {"mode": mode, "sweep_p": p,
                   "p10": p10, "p35": p35, "p35_relaxed": p35_relaxed,
                   "bg_a": bg[0], "bg_c": bg[1], "bg_g": bg[2], "bg_t": bg[3],
                   "n_igrs_scanned": len(all_igr),
                   "elapsed_s": elapsed,
                   "sweep_table_path": str(sweep_table_path.relative_to(cfg.output_dir)),
                   "n_in_sweep_table": int(len(final_conf))}
            row.update(counts)
            records.append(row)
            print(f"    p={p:.2e}  p10={p10:.2e} p35={p35:.2e} p35r={p35_relaxed:.2e}  "
                  f"total={counts['n_total']} "
                  f"(A={counts['n_path_A']} B={counts['n_path_B']} "
                  f"C={counts['n_path_C']})  marker={counts['n_marker']}  "
                  f"{elapsed:.1f}s")

    sweep_df = pd.DataFrame(records)
    cols = ["mode", "sweep_p", "p10", "p35", "p35_relaxed",
            "bg_a", "bg_c", "bg_g", "bg_t", "n_igrs_scanned",
            "n_total", "n_path_A", "n_path_B", "n_path_C",
            "n_marker", "n_marker_path_A", "n_marker_path_B", "n_marker_path_C",
            "n_in_sweep_table", "sweep_table_path", "elapsed_s"]
    sweep_df = sweep_df[cols]
    sweep_df.to_csv(out_tsv, sep="\t", index=False)
    print(f"\n  Sweep summary ({len(sweep_df)} rows) -> {out_tsv}")
    print(f"  Per-sweep-point tables ({len(sweep_df)} files) in {tables_dir}/")

    try:
        _render_p_sweep_figure(sweep_df, out_png)
    except Exception as e:
        print(f"  Figure rendering failed: {e}")

    _mark_current(cfg, "p_sweep", ps_in, ps_par)
    print("  Diagnostic complete.\n")


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

    _require(9, (cfg.igr_summary, 2), (cfg.gff_file, 1))
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


def _build_motif_diagram_svg(seq_len, motif_hits, width=700, height=50,
                             max_len=None):
    """Build an inline SVG showing motif positions along a promoter sequence.

    Returns an SVG string.
    """
    if seq_len == 0:
        return ""

    if not isinstance(max_len, (int, float)) or max_len <= 0:
        max_len = seq_len

    margin = 10
    track_y = 25
    track_h = 10
    full_w = width - 2 * margin             # drawing area between margins
    scale = full_w / max_len                # shared px per bp across rows
    bar_w = seq_len * scale                  # this promoter's bar width
    bar_left = margin + (full_w - bar_w)     # right-align: 3' end at right margin
    bar_right = width - margin

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
        # Sequence track (grey bar): proportional to length, right-aligned
        f'<rect x="{bar_left:.1f}" y="{track_y}" width="{bar_w:.1f}" height="{track_h}" '
        f'fill="#e0e0e0" rx="3"/>',
        # 5' label at this bar's (variable) left edge; 3' at the right margin
        f'<text x="{bar_left:.1f}" y="16" font-size="11" fill="#666">5\'</text>',
        f'<text x="{bar_right - 10}" y="16" font-size="11" fill="#666">3\'</text>',
        # Length label centered over this bar
        f'<text x="{bar_left + bar_w / 2:.1f}" y="16" font-size="10" fill="#999" '
        f'text-anchor="middle">{seq_len} bp</text>',
    ]

    for hit in motif_hits:
        start = hit.get("start", 0)
        stop = hit.get("stop", start + 10)
        db = hit.get("motif_database", "").lower()
        motif_name = hit.get("motif_id", "")
        colour = db_colours.get(db, default_colour)

        # Positions are 1-based along the 5'->3' sequence; place them
        # relative to this bar's left edge using the shared scale.
        x = bar_left + (start - 1) * scale
        w = max((stop - start + 1) * scale, 2)  # at least 2px wide

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
    s11_in = [cfg.final_table, cfg.promoter_markers_verified, cfg.all_promoters_verified_fasta,
              cfg.marker_promoters_verified_fasta, cfg.cds_annotations, cfg.igr_summary,
              cfg.input_fasta]
    s11_par = {"archaea": bool(cfg.is_archaea)}
    if not force and _is_current(cfg, "step11", [cfg.report_html], s11_in, s11_par):
        print("── HTML report already exists, skipping ──")
        print("  Step 11 complete.\n")
        return

    _require(11, (cfg.igr_summary, 2), (cfg.promoter_markers_verified, 8), (cfg.cds_annotations, 9))
    _forget(cfg, "step11")
    print("── Generating HTML report ──")

    # Only motif-verified marker promoters (step 8) are reported. The
    # unverified marker list (step 7) is deliberately NOT used as a
    # fallback: the report is headed "Verified marker promoters".
    promoter_df = None
    if cfg.promoter_markers_verified.exists():
        try:
            df = pd.read_csv(cfg.promoter_markers_verified, sep="\t")
            if "orientation" in df.columns:
                promoter_df = df[df["orientation"].isin(["CO_F", "CO_R"])].copy()
        except pd.errors.EmptyDataError:
            pass

    if promoter_df is None or promoter_df.empty:
        print("  No verified marker promoters available for report "
              "(run step 8 first, or the scan called none).")
        if cfg.report_html.exists():
            cfg.report_html.unlink()
            print("  Removed the report of an earlier run, which no longer matches these results.")
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

    need_orient = ("sequence_5p_to_3p" not in promoter_df.columns
                   or promoter_df["sequence_5p_to_3p"].isna().any())
    if need_orient and "sequence" in promoter_df.columns:
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

    has_motif_cols = "pos_10" in promoter_df.columns

    if has_motif_cols:
        promoter_df = _marker_order(promoter_df)
    _ranked_by_score = bool(
        has_motif_cols and "shortlist_score" in promoter_df.columns
        and pd.to_numeric(promoter_df["shortlist_score"], errors="coerce").notna().any())

    n_promoters = len(promoter_df)
    n_annotated = sum(1 for _, r in promoter_df.iterrows()
                      if r["associated_gene"] in annotation_map
                      and annotation_map[r["associated_gene"]] != "hypothetical protein")
    n_with_motifs = sum(1 for _, r in promoter_df.iterrows()
                        if has_motif_cols and pd.notna(r.get("pos_10"))
                        and str(r.get("pos_10")) != ".")

    top_df = promoter_df.head(10).copy()
    n_shown = len(top_df)
    try:
        max_len_display = int(top_df["length"].max())
    except (ValueError, TypeError):
        max_len_display = 0

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
        n_shown=n_shown,
        n_all_promoters=n_all_promoters,
        n_annotated=n_annotated,
        n_with_motifs=n_with_motifs,
        all_verified_fasta=html_mod.escape(all_verified_name),
        marker_verified_fasta=html_mod.escape(marker_verified_name),
        cds_annotations_tsv=html_mod.escape(cds_ann_name),
        final_table=html_mod.escape(final_table_name),
        domain_label=domain_label,
    )]

    # Legend for motif diagram + path colour code
    html_parts.append("""
    <div class="legend">
        <strong>Motif path:</strong>
        <span class="legend-item"><span class="path path-a">A</span> linked &minus;10/&minus;35 (same subgroup), 16&ndash;18 bp</span>
        <span class="legend-item"><span class="path path-b">B</span> extended &minus;10 + linked relaxed &minus;35 (same subgroup), 16&ndash;18 bp</span>
        <span class="legend-item"><span class="path path-c">C</span> unlinked &minus;10/&minus;35 (different subgroups), 16&ndash;18 bp</span>
    </div>
    <div class="legend">
        <strong>Motif diagram:</strong>
        <span class="legend-item"><span class="legend-swatch" style="background:#1565c0"></span>&minus;35 element</span>
        <span class="legend-item"><span class="legend-swatch" style="background:#c62828"></span>&minus;10 element</span>
        <span class="legend-item"><span class="legend-swatch" style="background:#ff8f00"></span>Extended &minus;10 (TG)</span>
    </div>
    <p style="font-size:0.8rem; color:#777; margin-bottom:16px;">
        """ + (
        "Showing the top 10 marker promoters, best marker matches first (a profile's "
        "highest-scoring gene, then one reached through its operon, then weaker matches) and "
        "within each class by <code>shortlist_score</code>, a ranker fitted on 49 "
        "genomes and validated by leave-one-genome-out to put real promoters at the top of the "
        "list: about 11 of the top 20 calls carry an annotated intergenic transcription start, "
        "against about 9 under the ordering this replaces. The score orders calls &mdash; it is "
        "not a calibrated probability that a given promoter will work in a construct."
        if _ranked_by_score else
        "Showing the top 10 marker promoters, best marker matches first, then ordered by motif "
        "path (A&rarr;B&rarr;C) and combined "
        "&minus;log<sub>10</sub>(p); no ranker spec was available for <code>shortlist_score</code>.") + """
        Sequence diagrams are drawn
        proportional to promoter length and aligned on the 3' (gene-proximal) end.
    </p>
    """)

    def _col_has_values(col):
        # A calibration value is present iff the column has any real
        # (non-NaN) number; blank/uncalibrated reads back as NaN.
        if not has_motif_cols or col not in promoter_df.columns:
            return False
        return bool(pd.to_numeric(promoter_df[col], errors="coerce").notna().any())
    _has_ppv = _col_has_values("ppv_interval")
    _has_lr = _col_has_values("lr_pos")
    _sig_label = "PPV (interval)" if _has_ppv else ("LR+" if _has_lr else "PPV (interval)")
    _sig_title = ("Calibrated probability the interval carries a promoter, at the stated prior (--prior)."
                  if _has_ppv else
                  "Positive likelihood ratio (base-rate-free), from the held-out calibration."
                  if _has_lr else
                  "Requires a calibration JSON (--ppv-calibration); blank when uncalibrated.")

    # Table
    html_parts.append(f"""
    <table>
    <thead>
    <tr>
        <th>Promoter ID</th>
        <th>Contig</th>
        <th>Position</th>
        <th>Length</th>
        <th>Orientation</th>
        <th>Motif path</th>
        <th title="How strongly the region is tied to a marker profile; best matches are listed first">Marker evidence</th>
        <th>&minus;10</th>
        <th>&minus;35</th>
        <th title="Length-standardised empirical score (score_std): best combined &minus;log10 empirical quantile of a &minus;10/&minus;35 pair over all spacer-linked positions, minus the shuffle-null location fitted on log(n_pairs). A ranker input, not the sort key.">Score</th>
        <th title="{html_mod.escape(_sig_title)}">{_sig_label}</th>
        <th>Spacer</th>
        <th>Associated CDS</th>
        <th>Protein</th>
        <th>Sequence diagram</th>
        <th>Sequence (5'&rarr;3')</th>
    </tr>
    </thead>
    <tbody>
    """)

    for _, row in top_df.iterrows():
        igr_id = row["igr_id"]
        gene = row["associated_gene"]
        product = annotation_map.get(gene, "")

        # Build motif diagram from -10 and -35 positions
        motif_hits = []
        if has_motif_cols:
            pos_10 = row.get("pos_10", ".")
            pos_35 = row.get("pos_35", ".")
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
                    ext10_raw = row.get("has_ext10", False)
                    if isinstance(ext10_raw, bool):
                        has_ext10_flag = ext10_raw
                    else:
                        has_ext10_flag = str(ext10_raw).strip().lower() == "true"
                    if has_ext10_flag:
                        motif_hits.append({
                            "motif_id": "ext-10 (TG)",
                            # TGn: TG occupies 0-based p10-3, p10-2
                            # (1-based p10-2, p10-1), one base 5' of
                            # the hexamer.
                            "start": p10 - 2,
                            "stop": p10 - 1,
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

        svg = _build_motif_diagram_svg(row["length"], motif_hits,
                                       max_len=max_len_display)

        orient_class = "co-f" if row["orientation"] == "CO_F" else "co-r"
        seq = row.get("sequence_5p_to_3p", "")

        # Motif details (raw empirical-scan columns from step 8)
        motif_path = str(row.get("path", "")).strip() if has_motif_cols else ""
        path_key = motif_path if motif_path in ("A", "B", "C") else ""
        if path_key:
            path_cell = f'<span class="path path-{path_key.lower()}">{path_key}</span>'
            row_class = f"row-{path_key.lower()}"
        else:
            path_cell = '<span class="path path-none">&mdash;</span>'
            row_class = ""
        ev_key = row.get("marker_evidence", "")
        ev_key = ev_key if isinstance(ev_key, str) else ""
        if ev_key in MARKER_EVIDENCE_LABEL:
            ev_title = html_mod.escape(f"{row.get('marker_profiles', '')}: {row.get('marker_gene', '')}")
            evidence_cell = f'<span title="{ev_title}">{MARKER_EVIDENCE_LABEL[ev_key]}</span>'
        else:
            evidence_cell = "&mdash;"

        def _to_float(v):
            try:
                f = float(v)
            except (TypeError, ValueError):
                return None
            return None if f != f else f  # NaN -> None

        def _fmt_p(v):
            f = _to_float(v)
            return f"{f:.1e}" if (f is not None and f > 0) else ""

        seq_10 = str(row.get("seq_10", "")) if has_motif_cols else ""
        source_10 = str(row.get("source_10", "")) if has_motif_cols else ""
        emp_p_10 = _fmt_p(row.get("emp_p_10")) if has_motif_cols else ""
        seq_35 = str(row.get("seq_35", "")) if has_motif_cols else ""
        source_35 = str(row.get("source_35", "")) if has_motif_cols else ""
        emp_p_35 = _fmt_p(row.get("emp_p_35")) if has_motif_cols else ""
        spacer_len = row.get("spacer", "") if has_motif_cols else ""

        def _element_cell(seq_x, source_x, emp_p_x):
            if not seq_x or seq_x in (".", "nan"):
                return '<span class="na">—</span>'
            cell = html_mod.escape(seq_x)
            annot = ", ".join(x for x in (
                html_mod.escape(source_x) if source_x and source_x != "nan" else "",
                f"p={emp_p_x}" if emp_p_x else "") if x)
            if annot:
                cell += f' <small>({annot})</small>'
            return cell
        minus10_cell = _element_cell(seq_10, source_10, emp_p_10)
        minus35_cell = _element_cell(seq_35, source_35, emp_p_35)

        spacer_cell = (str(int(_to_float(spacer_len)))
                       if str(spacer_len) not in ("", ".", "nan")
                       and _to_float(spacer_len) is not None
                       else '<span class="na">—</span>')

        # Length-standardised empirical score (always shown) + the
        # calibration-dependent significance value: PPV (interval) at the
        # stated prior, else base-rate-free LR+, else blank (uncalibrated
        # values are never fabricated).
        score_std_f = _to_float(row.get("score_std")) if has_motif_cols else None
        score_cell = f"{score_std_f:.2f}" if score_std_f is not None else '<span class="na">—</span>'
        if _has_ppv:
            _sv = _to_float(row.get("ppv_interval"))
            sig_cell = f"{_sv * 100:.0f}%" if _sv is not None else '<span class="na">—</span>'
        elif _has_lr:
            _sv = _to_float(row.get("lr_pos"))
            sig_cell = f"{_sv:.2f}" if _sv is not None else '<span class="na">—</span>'
        else:
            sig_cell = '<span class="na">—</span>'

        product_cell = html_mod.escape(product) if product else '<span class="na">—</span>'

        html_parts.append(f"""
        <tr class="{row_class}">
            <td><code>{html_mod.escape(igr_id)}</code></td>
            <td>{html_mod.escape(str(row['contig']))}</td>
            <td>{row['start']:,}–{row['end']:,}</td>
            <td>{row['length']}</td>
            <td><span class="orient {orient_class}">{row['orientation']}</span></td>
            <td>{path_cell}</td>
            <td>{evidence_cell}</td>
            <td>{minus10_cell}</td>
            <td>{minus35_cell}</td>
            <td>{score_cell}</td>
            <td>{sig_cell}</td>
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
    _mark_current(cfg, "step11", s11_in, s11_par)
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
    .path {{
        display: inline-block; padding: 2px 10px; border-radius: 10px;
        font-size: 0.78rem; font-weight: 700; letter-spacing: 0.3px;
        min-width: 22px; text-align: center;
    }}
    .path-a {{ background: #e8f5e9; color: #2e7d32; }}
    .path-b {{ background: #fff8e1; color: #e65100; }}
    .path-c {{ background: #ffebee; color: #c62828; }}
    .path-none {{ background: #f5f5f5; color: #999; }}
    tr.row-a {{ background: #f4faf4; }}
    tr.row-b {{ background: #fffcf2; }}
    tr.row-c {{ background: #fdf4f4; }}
    tr.row-a:hover, tr.row-b:hover, tr.row-c:hover {{ background: #eef4ef; }}
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
    This table lists the <strong>top {n_shown}</strong> of {n_marker_promoters} marker-filtered promoters confirmed by &minus;10/&minus;35 motif scanning, ordered by <code>shortlist_score</code> (see note below the legend).
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
    is_marker             Whether this IGR can drive a gene with a marker
                          HMM hit (the gene immediately downstream, or a
                          later gene of the same operon)
    marker_evidence       best_match | best_match_via_operon |
                          weaker_match | weaker_match_via_operon (blank
                          when not a marker region); see README
    marker_gene           The hit gene the label refers to
    marker_profiles       HMM profile(s) behind the label
    motif_confirmed       Whether a -10/-35 motif combination was found
    motif_path            Motif path classification (A, B, or C)
    motif_source_10       Subgroup ID of the best -10 hit (e.g. M001)
    motif_source_35       Subgroup ID of the best -35 hit (e.g. M001)
    motif_has_ext10       Whether extended -10 TG dinucleotide is present
    motif_emp_p_10        Empirical quantile of the -10 hit in the
                          per-genome dinucleotide-shuffle null
    motif_emp_p_35        Empirical quantile of the -35 hit in the
                          per-genome dinucleotide-shuffle null. For
                          Path C this is Šidák-adjusted for the number
                          of architectures searched, because Path C's
                          -35 is the best match among all architectures
                          OTHER than the one supplying the -10, and the
                          minimum of k comparisons overstates the
                          evidence if read as a single test. Paths A
                          and B fix the architecture from the -10 and
                          are reported unadjusted. Comparable across
                          paths; use with motif_path in view
    motif_emp_p_35_uncorrected
                          The same quantity before adjustment, so the
                          raw per-architecture match is not lost. Equal
                          to motif_emp_p_35 for Paths A and B
    motif_p35_multiplicity_k
                          Number of comparisons the -35 p-value is a
                          minimum over: 1 for Paths A and B, and the
                          count of other architectures for Path C
    motif_seq_10          Sequence of -10 element
    motif_seq_35          Sequence of -35 element
    motif_pos_10          Position of -10 element in scanned sequence
    motif_pos_35          Position of -35 element
    motif_spacer          Spacer length between -35 and -10
    motif_score_std       Length-standardised empirical score: best
                          combined -10/-35 log-odds, standardised
                          against the shuffle null and normalised for
                          the number of tested positions.
                          NOT a ranking quantity: sorting on it puts
                          WORSE calls first (leave-one-genome-out
                          precision@20 41.3% against a 47.96% base
                          rate), and the ranker's fitted weight on it
                          is negative. Retained because it is the
                          calibration input and the threshold statistic.
    motif_n_pairs         Number of (architecture, spacer) pairs that
                          passed in this interval. Higher means the
                          call rests on several independent matches
                          rather than one; a ranker input
    motif_shortlist_score Shortlist ranking score in (0,1), descending
                          = most trustworthy first. This is what the
                          final table and the HTML report are sorted
                          by. Fitted on 49 genomes against "interval
                          contains an annotated INTERGENIC TSS" and
                          validated leave-one-genome-out: precision@20
                          55.9% vs 46.9% for the ordering it replaced
                          (+9.0 pp, better in 32/49 genomes,
                          Wilcoxon P = 1.1e-4), precision@10 58.8%.
                          Its ORDER is validated; the value is not a
                          calibrated probability that a promoter will
                          work in a construct, and the target is
                          promoter PRESENCE, not strength. See
                          profinder/ranker.py
    motif_lr_pos          Positive likelihood ratio (base-rate-free)
                          from the held-out calibration; blank when
                          uncalibrated (never fabricated)
    motif_ppv_interval    Calibrated probability the INTERVAL carries a
                          promoter, at the stated prior (--prior); blank
                          when uncalibrated or no prior given.
                          Blank in the shipped configuration BY DESIGN:
                          a fit on 49 held-out genomes was rejected
                          because within a genome the score's AUC against
                          TSS truth is 0.5007 (Wilcoxon P = 0.74), so a
                          per-call PPV from it would rank calls no better
                          than chance. Use motif_path (A/B/C), whose
                          precision differences are measured, as the
                          confidence statement instead.
    motif_in_domain       Whether the genome's GC sits in the
                          calibration's applicability domain
    sequence_5p_to_3p     Full promoter sequence oriented 5'→3'
    """
    s10_out, s10_in = [cfg.final_table, cfg.marker_summary], _step10_inputs(cfg)
    s10_par = {"cds_bp": cfg.cds_bp, "hmm_bitscore_min": cfg.hmm_bitscore_min,
               "profiles": _hmm_profile_names(cfg), "operons": _operon_params(cfg)}
    if not force and _is_current(cfg, "step10", s10_out, s10_in, s10_par):
        print("── Final output table already exists, skipping ──")
        print("  Step 10 complete.\n")
        return

    _require(10, (cfg.igr_summary, 2), (cfg.promoter_markers, 7),
             (cfg.motif_best_all, 8), (cfg.cds_annotations, 9))
    _forget(cfg, "step10")
    _warn_legacy_operons(cfg)
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
        _write_marker_summary(cfg, pd.DataFrame())
        _mark_current(cfg, "step10", s10_in, s10_par)
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
    igr = _attach_marker_evidence(igr, _load_marker_evidence(cfg))

    # 6. Motif scan results (best hit per IGR from step 8). Column layout is
    # shared with the p-sweep diagnostic via the module-level constants.
    motif_cols = _SCAN_MOTIF_COLS
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
    columns = list(_FINAL_TABLE_COLUMNS)

    if cfg.cds_bp > 0:
        contigs = _load_contigs(cfg.input_fasta)
        igr = _add_cds_column(igr, contigs, cfg.cds_bp)
        columns.append("sequence_5p_to_3p_cds")

    # 8. Select and order final columns; rows ordered by shortlist_score
    #    (called IGRs first), see _shortlist_order.
    out = igr[columns].copy()
    out.rename(columns={"igr_id": "promoter_id"}, inplace=True)
    out = _shortlist_order(out, prefix="motif_")

    out.to_csv(cfg.final_table, sep="\t", index=False)
    print(f"  Final table ({len(out)} rows) -> {cfg.final_table}")
    _write_marker_summary(cfg, out)
    _mark_current(cfg, "step10", s10_in, s10_par)
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
    # These three read their defaults from Config rather than restating them. When they
    # carried literals they silently diverged from it: the literals stayed at an older
    # operating point while Config moved, so a command-line run used thresholds the config
    # no longer specified -- and the help text quoted a third set of numbers again. %(default)s
    # keeps the advertised value and the applied value the same object.
    parser.add_argument("--p10", type=float,
                        default=Config.__dataclass_fields__["motif_p10"].default,
                        help="p-value threshold for -10 motif hits. The default is the "
                             "swept value that gets the most genomes to a full top-10 "
                             "candidate list; see profinder/config.py for the sweep. "
                             "(default: %(default)g)")
    parser.add_argument("--p35", type=float,
                        default=Config.__dataclass_fields__["motif_p35"].default,
                        help="Strict p-value threshold for -35 motif hits, "
                             "used for Path A and Path C (default: %(default)g)")
    parser.add_argument("--p35-relaxed", type=float,
                        default=Config.__dataclass_fields__["motif_p35_relaxed"].default,
                        help="Relaxed p-value threshold for -35 motif hits, "
                             "used only for Path B (extended -10). Should exceed --p35; at or "
                             "below it Path B loses the weak -35 matches it exists to "
                             "capture. (default: %(default)g)")
    parser.add_argument("--max-dist-to-cds-start", type=int,
                        default=Config.__dataclass_fields__[
                            "max_dist_to_cds_start"].default,
                        help="Maximum distance (bp) between the right edge "
                             "of the -10 motif and the downstream CDS start. "
                             "-10 hits further out are dropped before the "
                             "spacer search. The default is the swept elbow "
                             "on the 49-genome panel: it recovers 94.8%% of "
                             "the recall attainable with no limit, gains in "
                             "49/49 genomes over 200 bp, and costs no "
                             "precision (see config.py for the sweep). "
                             "(default: %(default)s)")
    parser.add_argument("--n-shuffles", type=int, default=10,
                        help="Dinucleotide shuffles per interval for the "
                             "per-genome empirical null (default: 10)")
    parser.add_argument("--ppv-calibration", type=Path, default=None,
                        help="Per-path PPV calibration JSON (default: "
                             "bundled; if none, LR+/ppv_interval are blank)")
    parser.add_argument("--prior", type=float, default=None,
                        help="Prevalence prior for ppv_interval; without "
                             "it only the base-rate-free LR+ and the "
                             "calibration IQR are reported")

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
                        help="Max gap between consecutive same-strand genes of an "
                             "operon in bp (default: 25)")
    parser.add_argument("--min-flanking-dist", type=int, default=75,
                        help="Min gap to a same-strand gene upstream of an operon's "
                             "first gene in bp; a closer one means the operon may start "
                             "further upstream, so it is dropped. The gap after the last "
                             "gene is not tested (default: 75)")
    parser.add_argument("--hmm-bitscore", type=float, default=25.0,
                        help="Minimum HMM bitscore (default: 25.0)")
    parser.add_argument("--cds-bp", type=int, default=0,
                        help="Number of CDS-start nucleotides to append to "
                             "promoter sequences in additional FASTA outputs "
                             "and the final table. 0 = disabled (default: 0)")

    # Diagnostic: -10 / -35 p-value sweep
    parser.add_argument("--diagnose-p-sweep", action="store_true",
                        help="Run the p-value sweep diagnostic after step 7 "
                             "instead of running step 8 onward. Sweeps -10 "
                             "and -35 PWM thresholds in three modes "
                             "(matched, p10_only, p35_only) and writes "
                             "diagnostics/p_sweep.tsv plus a figure.")
    parser.add_argument("--sweep-min", type=float, default=1e-5,
                        help="Smallest p-value in the sweep "
                             "(default: 1e-5)")
    parser.add_argument("--sweep-max", type=float, default=1e-1,
                        help="Largest p-value in the sweep "
                             "(default: 1e-1)")
    parser.add_argument("--sweep-n", type=int, default=13,
                        help="Number of log-spaced sweep points "
                             "(default: 13)")
    parser.add_argument("--sweep-modes", default="matched,p10_only,p35_only",
                        help="Comma-separated list of sweep modes "
                             "(default: matched,p10_only,p35_only). "
                             "Valid modes: matched | p10_only | p35_only.")
    parser.add_argument("--sweep-p35-relaxed-ratio", type=float,
                        default=(Config.__dataclass_fields__[
                                     "motif_p35_relaxed"].default
                                 / Config.__dataclass_fields__[
                                     "motif_p35"].default),
                        help="In every sweep point, p35_relaxed is set to "
                             "this ratio x p35 (capped at 1.0). The default "
                             "is derived from the shipped pair so that the "
                             "sweep passes through the operating point the "
                             "tool actually uses; hardcoding it let the two "
                             "drift apart when the defaults changed. "
                             "(default: %(default)g)")

    args = parser.parse_args()

    if args.list:
        for num, name, _ in STEPS:
            print(f"  {num:2d}. {name}")
        sys.exit(0)

    last = max(num for num, _, _ in STEPS)
    if not 1 <= args.start <= args.end <= last:
        parser.error(f"--start/--end must satisfy 1 <= start <= end <= {last}")
    if args.threads < 1:
        parser.error("--threads must be >= 1")
    if args.n_shuffles < 1:
        parser.error("--n-shuffles must be >= 1 (the empirical null needs at least one shuffle)")
    for flag, val in (("--p10", args.p10), ("--p35", args.p35), ("--p35-relaxed", args.p35_relaxed)):
        if not 0.0 < val < 1.0:
            parser.error(f"{flag} must lie strictly between 0 and 1 (got {val})")
    if args.prior is not None and not 0.0 < args.prior < 1.0:
        parser.error(f"--prior must lie strictly between 0 and 1 (got {args.prior})")
    if args.igr_min < 1 or args.igr_min > args.igr_max:
        parser.error("--igr-min must be >= 1 and <= --igr-max")
    if args.cds_bp < 0 or args.max_dist_to_cds_start < 0:
        parser.error("--cds-bp and --max-dist-to-cds-start must be >= 0")

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
        motif_p35_relaxed=args.p35_relaxed,
        max_dist_to_cds_start=args.max_dist_to_cds_start,
        cds_bp=args.cds_bp,
        n_shuffles=args.n_shuffles,
        ppv_calibration=args.ppv_calibration,
        ppv_prior=args.prior,
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

    # Diagnostic mode: run the p-value-independent prerequisites
    # (Prokka, IGR extraction, operons, HMM, marker matching, CDS
    # annotation — steps 1-7 plus 9), skip the default-threshold motif
    # scan (step 8), and replace steps 10-11 with the p-value sweep.
    # Step 9 is included because the per-sweep-point profinder_results.tsv
    # files carry gene_name / locus_tag / product columns, which are
    # p-value-independent and worth computing once. Cached outputs from
    # any prior run (full or diagnostic) are reused.
    if getattr(args, "diagnose_p_sweep", False):
        DIAGNOSTIC_PREREQ_STEPS = {1, 2, 3, 4, 5, 6, 7, 9}
        for num, name, func in STEPS:
            if num not in DIAGNOSTIC_PREREQ_STEPS:
                continue
            if args.start <= num <= args.end:
                print(f"\n{'=' * 60}")
                print(f"  STEP {num}: {name}")
                print(f"{'=' * 60}\n")
                func(cfg, force=args.force)
        modes = tuple(m.strip() for m in args.sweep_modes.split(",") if m.strip())
        print(f"\n{'=' * 60}")
        print("  DIAGNOSTIC: -10 / -35 p-value sweep")
        print(f"{'=' * 60}\n")
        diagnose_p_sweep(
            cfg,
            sweep_min=args.sweep_min,
            sweep_max=args.sweep_max,
            sweep_n=args.sweep_n,
            modes=modes,
            p35_relaxed_ratio=args.sweep_p35_relaxed_ratio,
            force=args.force,
        )
        return

    for num, name, func in STEPS:
        if args.start <= num <= args.end:
            print(f"\n{'=' * 60}")
            print(f"  STEP {num}: {name}")
            print(f"{'=' * 60}\n")
            func(cfg, force=args.force)


def _batch_cell(row, key):
    """Return a batch-table cell as a stripped string; empty cells (pandas
    NaN / None / blank) become "" so an unfilled optional column is treated
    as absent rather than as the literal path "nan"."""
    v = row.get(key, "")
    if v is None or (isinstance(v, float) and v != v):
        return ""
    v = str(v).strip()
    return "" if v.lower() in ("", "nan", "none") else v


def _link_prokka_files(cfg, row):
    """Symlink user-supplied Prokka files into the expected output layout.

    Returns the step number to start from: 2 if Prokka files were linked
    (skip step 1), or 1 if they were not provided.
    """
    gff_path = _batch_cell(row, "gff")
    faa_path = _batch_cell(row, "faa")

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
    fna_path = _batch_cell(row, "fna")
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
        sample_id = _batch_cell(row, "sample_id")
        fasta_cell = _batch_cell(row, "fasta")
        if not sample_id or not fasta_cell:
            sys.exit(f"Batch table row {idx + 1}: sample_id and fasta are required")
        fasta_path = Path(fasta_cell).resolve()

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
