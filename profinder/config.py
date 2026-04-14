"""
Configuration for ProFinder — bacterial and archaeal promoter identification
pipeline.

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

    # ── Domain selection ─────────────────────────────────────────────
    domain: str = "bacteria"  # "bacteria" or "archaea"

    # ── Tool paths ───────────────────────────────────────────────────
    prokka_bin: str = "prokka"
    hmmsearch_bin: str = "hmmsearch"
    conda_env_prokka: str = ""   # leave blank to skip conda activate
    conda_env_hmm: str = ""

    # ── HMM profiles directory (bundled by default) ────────────────
    hmm_profiles_dir: Path = None

    def __post_init__(self):
        """Resolve bundled HMM path if none was provided."""
        if self.hmm_profiles_dir is None:
            candidate = Path(__file__).parent / "hmms"
            self.hmm_profiles_dir = candidate if candidate.is_dir() else None

        # Normalise domain to lowercase
        self.domain = self.domain.lower()


    # ── Parallelism ──────────────────────────────────────────────────
    threads: int = 4
    num_workers: int = 4

    # ── Prokka parameters ────────────────────────────────────────────
    prokka_kingdom: str = "Bacteria"
    prokka_prefix: str = "genome"

    # ── IGR extraction parameters ────────────────────────────────────
    igr_size_min: int = 81
    igr_size_max: int = 1000

    # ── Operon identification parameters ─────────────────────────────
    max_internal_distance: int = 25
    min_flanking_distance: int = 75

    # ── HMM filtering ────────────────────────────────────────────────
    hmm_bitscore_min: float = 25.0

    # ── CDS extension ────────────────────────────────────────────────
    cds_bp: int = 0                   # number of CDS-start nt to append (0 = disabled)

    # ── Motif scanning parameters ────────────────────────────────────
    motifs_dir: Path = None           # directory containing .meme files
    motif_p10: float = 2.5e-3        # p-value threshold for -10 hits
    motif_p35: float = 2.5e-3        # p-value threshold for -35 hits (strict)
    motif_p35_relaxed: float = 0.05  # relaxed -35 threshold when ext -10 present

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
    def promoter_markers_hmm(self) -> Path:
        return self.output_dir / "promoter_markers_hmm.tsv"

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

    # ── CDS-extended FASTA outputs ──────────────────────────────────
    @property
    def promoter_cds_fasta(self) -> Path:
        return self.output_dir / "promoters_cds_bp.fasta"

    @property
    def promoter_cds_fasta_short(self) -> Path:
        return self.output_dir / "promoters_short_cds_bp.fasta"

    @property
    def all_promoter_cds_fasta(self) -> Path:
        return self.output_dir / "all_promoters_cds_bp.fasta"

    @property
    def all_promoter_cds_fasta_short(self) -> Path:
        return self.output_dir / "all_promoters_short_cds_bp.fasta"

    @property
    def is_archaea(self) -> bool:
        return self.domain == "archaea"

    # ── Motif scan output ────────────────────────────────────────────
    @property
    def motif_hits_all(self) -> Path:
        return self.output_dir / "motif_hits_all.tsv"

    @property
    def motif_best_all(self) -> Path:
        return self.output_dir / "motif_best_all.tsv"

    @property
    def motif_hits_markers(self) -> Path:
        return self.output_dir / "motif_hits_markers.tsv"

    @property
    def motif_best_markers(self) -> Path:
        return self.output_dir / "motif_best_markers.tsv"

    @property
    def promoter_markers_verified(self) -> Path:
        return self.output_dir / "promoter_markers_verified.tsv"

    @property
    def all_promoters_verified_fasta(self) -> Path:
        return self.output_dir / "all_promoters_verified.fasta"

    @property
    def marker_promoters_verified_fasta(self) -> Path:
        return self.output_dir / "marker_promoters_verified.fasta"

    # ── Final output table ─────────────────────────────────────────────
    @property
    def final_table(self) -> Path:
        return self.output_dir / "profinder_results.tsv"

    # ── Report ────────────────────────────────────────────────────────
    @property
    def report_html(self) -> Path:
        return self.output_dir / "promoter_report.html"

    def ensure_dirs(self):
        """Create all output directories."""
        for d in [self.output_dir, self.prokka_dir, self.igr_dir,
                  self.hmm_dir]:
            d.mkdir(parents=True, exist_ok=True)
