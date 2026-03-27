"""
Configuration for ProFinder — bacterial promoter identification pipeline.

All paths are resolved at runtime so they work whether passed via CLI
or set in the dataclass defaults.
"""

from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class Config:
    # ── Required inputs (typically set via CLI) ──────────────────────
    input_fasta: Path = Path("genome.fasta")
    output_dir: Path = Path("output")

    # ── Tool paths ───────────────────────────────────────────────────
    prokka_bin: str = "prokka"
    hmmsearch_bin: str = "hmmsearch"
    fimo_bin: str = "fimo"
    conda_env_prokka: str = ""   # leave blank to skip conda activate
    conda_env_hmm: str = ""
    conda_env_meme: str = ""

    # ── HMM databases (bundled by default) ──────────────────────────
    tigrfam_hmm: Path = None
    pfam_hmm: Path = None

    def __post_init__(self):
        """Resolve bundled HMM and weight paths if none were provided."""
        bundled_hmms = Path(__file__).parent / "hmms"
        if self.tigrfam_hmm is None:
            candidate = bundled_hmms / "tigrfam.hmm"
            self.tigrfam_hmm = candidate if candidate.exists() else None
        if self.pfam_hmm is None:
            candidate = bundled_hmms / "Pfam-A.hmm"
            self.pfam_hmm = candidate if candidate.exists() else None
        if self.lcnn_weights_dir is None:
            candidate = Path(__file__).parent / "weights" / "PromoterLCNN"
            self.lcnn_weights_dir = candidate if candidate.is_dir() else None

    # ── Parallelism ──────────────────────────────────────────────────
    threads: int = 4
    num_workers: int = 4

    # ── Prokka parameters ────────────────────────────────────────────
    prokka_kingdom: str = "Bacteria"
    prokka_prefix: str = "genome"

    # ── IGR extraction parameters ────────────────────────────────────
    igr_size_min: int = 75
    igr_size_max: int = 1000

    # ── Operon identification parameters ─────────────────────────────
    max_internal_distance: int = 25
    min_flanking_distance: int = 75

    # ── HMM filtering ────────────────────────────────────────────────
    hmm_bitscore_min: float = 25.0

    # ── PromoterLCNN parameters ────────────────────────────────────────
    lcnn_weights_dir: Path = None     # parent dir containing IsPromoter_fold_5/ and PromotersOnly_fold_1/

    # ── FIMO parameters ──────────────────────────────────────────────
    motifs_dir: Path = None           # directory containing .meme files
    fimo_threshold: float = 1e-4     # p-value threshold for motif hits

    # ── Derived paths ────────────────────────────────────────────────
    @property
    def prokka_dir(self) -> Path:
        return self.output_dir / "prokka"

    @property
    def gff_file(self) -> Path:
        return self.prokka_dir / f"{self.prokka_prefix}.gff"

    @property
    def faa_file(self) -> Path:
        return self.prokka_dir / f"{self.prokka_prefix}.faa"

    @property
    def fna_file(self) -> Path:
        return self.prokka_dir / f"{self.prokka_prefix}.fna"

    @property
    def igr_dir(self) -> Path:
        return self.output_dir / "igr"

    @property
    def igr_fasta(self) -> Path:
        return self.igr_dir / "intergenic_regions.fasta"

    @property
    def igr_summary(self) -> Path:
        return self.igr_dir / "igr_summary.tsv"

    @property
    def operon_file(self) -> Path:
        return self.output_dir / "operons.tsv"

    @property
    def hmm_dir(self) -> Path:
        return self.output_dir / "hmm"

    @property
    def hmm_combined(self) -> Path:
        return self.hmm_dir / "hmm_combined.tsv"

    @property
    def hmm_filtered(self) -> Path:
        return self.hmm_dir / "hmm_filtered.tsv"

    @property
    def operon_filtered(self) -> Path:
        return self.output_dir / "operons_filtered.tsv"

    @property
    def operon_filtered_markers(self) -> Path:
        return self.output_dir / "operons_filtered_markers.tsv"

    @property
    def promoter_markers(self) -> Path:
        return self.output_dir / "promoter_markers.tsv"

    @property
    def promoter_fasta(self) -> Path:
        return self.output_dir / "promoters.fasta"

    @property
    def promoter_fasta_short(self) -> Path:
        return self.output_dir / "promoters_short.fasta"

    @property
    def all_promoter_fasta(self) -> Path:
        return self.output_dir / "all_promoters.fasta"

    @property
    def all_promoter_fasta_short(self) -> Path:
        return self.output_dir / "all_promoters_short.fasta"

    # ── PromoterLCNN output ────────────────────────────────────────────
    @property
    def lcnn_predictions(self) -> Path:
        return self.output_dir / "lcnn_predictions.tsv"

    @property
    def promoter_markers_verified(self) -> Path:
        return self.output_dir / "promoter_markers_verified.tsv"

    # ── FIMO output ───────────────────────────────────────────────────
    @property
    def fimo_dir(self) -> Path:
        return self.output_dir / "fimo"

    @property
    def fimo_combined(self) -> Path:
        return self.fimo_dir / "fimo_combined.tsv"

    # ── Report ────────────────────────────────────────────────────────
    @property
    def report_html(self) -> Path:
        return self.output_dir / "promoter_report.html"

    def ensure_dirs(self):
        """Create all output directories."""
        for d in [self.output_dir, self.prokka_dir, self.igr_dir,
                  self.hmm_dir, self.fimo_dir]:
            d.mkdir(parents=True, exist_ok=True)
