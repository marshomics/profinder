# ProFinder

Identify and extract bacterial promoter sequences from a single genome FASTA file, annotate their associated coding sequences, scan for known transcription factor binding motifs, and produce an HTML report.

The pipeline annotates a genome with [Prokka](https://github.com/tseemann/prokka), extracts intergenic regions (IGRs) between annotated genes, identifies operons using a two-pass proximity/flanking-distance algorithm, screens for marker genes via hmmsearch against TIGRfam and Pfam (HMM profiles are bundled), and outputs promoter sequences oriented 5'→3'. Candidate promoters are then verified with [PromoterLCNN](https://github.com/occasumlux/Promoters), a two-stage CNN that classifies each sequence as promoter or non-promoter and assigns a sigma-factor subtype (σ70, σ24, σ28, σ38, σ32, σ54). CDS functional annotations come from Prokka's gene product predictions. Promoters are scanned for known motifs with FIMO, and a visual HTML report is generated.

All databases ship with the pipeline. No additional downloads are needed beyond installing the external tools.

## Requirements

**Python ≥ 3.7** with:

- pandas ≥ 1.3
- biopython ≥ 1.79
- tensorflow ≥ 2.6

**External tools** (must be on `$PATH` or specified via CLI flags):

- [Prokka](https://github.com/tseemann/prokka)
- [MEME Suite](https://meme-suite.org/) (specifically `fimo`, for motif scanning)
- [HMMER](http://hmmer.org/) (for marker gene screening; HMM profiles are bundled)

## Installation

```bash
pip install .
```

Or for development:

```bash
pip install -e .
```

## Quick start

Minimal run (bundled HMMs and motifs):

```bash
profinder -i my_genome.fasta -o results/
```

Override bundled HMMs with your own:

```bash
profinder -i my_genome.fasta -o results/ \
    --tigrfam /path/to/custom_tigrfam.hmm \
    --pfam /path/to/custom_Pfam-A.hmm
```

Custom motifs directory:

```bash
profinder -i my_genome.fasta -o results/ \
    --motifs-dir /path/to/my_motifs/
```

The pipeline produces checkpoint files at each step. Re-running the same command skips completed steps automatically. Use `--force` to re-run everything from scratch.

## Pipeline steps

| Step | Description | Key output |
|------|-------------|------------|
| 1 | Run Prokka | `prokka/genome.gff`, `prokka/genome.faa` |
| 2 | Extract intergenic regions | `igr/igr_summary.tsv`, `igr/intergenic_regions.fasta` |
| 3 | Identify operons | `operons.tsv` |
| 4 | Run hmmsearch | `hmm/hmm_tigrfam.tblout`, `hmm/hmm_pfam.tblout` |
| 5 | Filter HMM output | `hmm/hmm_filtered.tsv` |
| 6 | Filter operons + add markers | `operons_filtered_markers.tsv` |
| 7 | Match IGRs to marker operons | `promoter_markers.tsv` |
| 8 | Extract marker-filtered promoters | `promoters.fasta`, `promoters_short.fasta` |
| 9 | Extract all promoter-orientation IGRs | `all_promoters.fasta`, `all_promoters_short.fasta` |
| 10 | Predict promoters (PromoterLCNN) | `lcnn_predictions.tsv`, `promoter_markers_verified.tsv` |
| 11 | Annotate CDS (Prokka) | `cds_annotations.tsv` |
| 12 | Scan motifs with FIMO | `fimo/fimo_combined.tsv` |
| 13 | Generate HTML report | `promoter_report.html` |

Steps 4–5 use bundled TIGRfam and Pfam HMM profiles by default. You can override them with `--tigrfam` and `--pfam` if you have custom profiles.

Step 10 runs PromoterLCNN on each marker-filtered promoter. The 3'-terminal 81 nt of each sequence (closest to the transcription start site) is passed through a two-stage CNN: a binary classifier filters non-promoters, and a sigma-factor classifier assigns one of six subtypes to each confirmed promoter. Sequences shorter than 81 nt are classified as non-promoters. Pre-trained weights are bundled in `weights/PromoterLCNN/`.

Step 11 extracts protein product names directly from Prokka's GFF annotations for each verified promoter's associated CDS.

Step 12 uses FIMO from the MEME Suite. Three motif databases are bundled with the pipeline (CollecTF, PRODORIC, RegTransBase). You can point to a different directory with `--motifs-dir`.

## CDS annotation

Each promoter is associated with its immediately downstream CDS based on orientation. For CO_F promoters (both flanking genes on the + strand), the right gene is the associated CDS. For CO_R (both on the − strand), the left gene is the associated CDS, because that gene's transcription start site faces the IGR.

Functional annotation uses the `product=` field from Prokka's GFF output for each associated CDS.

## Motif scanning

FIMO scans each promoter sequence (oriented 5'→3') against position weight matrices from three bacterial transcription factor databases:

- **CollecTF** — experimentally validated TF binding sites
- **PRODORIC** — prokaryotic gene regulation database
- **RegTransBase** — regulatory interactions in prokaryotes

The HTML report includes an inline SVG diagram for each promoter showing where motif hits fall along the sequence, colour-coded by database.

## IGR orientation logic

Each intergenic region is classified by the orientation of its flanking genes:

- **CO_F** (co-oriented forward): both flanking genes on the + strand. The IGR is a promoter candidate for the downstream gene. Sequence is already 5'→3'.
- **CO_R** (co-oriented reverse): both on the − strand. The IGR is a promoter candidate for the downstream gene (which runs right-to-left). Sequence is reverse-complemented to achieve 5'→3'.
- **DP** (divergent promoter): genes point away from each other (← IGR →). Contains promoters for both genes, so orientation is ambiguous. These are excluded from the final FASTA output.
- **CONV** (convergent): genes point toward each other (→ IGR ←). This is a terminator region, not a promoter. Excluded.

## Operon identification

Operons are identified in two passes:

1. **Proximity clustering**: consecutive CDS on the same contig within `--max-internal-dist` bp (default 25) are grouped.
2. **Flanking validation**: groups are kept only if both upstream and downstream gaps exceed `--min-flanking-dist` bp (default 75).

## CLI reference

```
profinder -i FASTA -o DIR [options]

Required:
  -i, --input            Input genome FASTA
  -o, --output           Output directory (default: output/)

Step control:
  --start N              First step to run (default: 1)
  --end N                Last step to run (default: 13)
  --list                 List steps and exit
  --force                Re-run all steps, ignoring checkpoints

Tool paths:
  --prokka PATH          Prokka binary (default: prokka)
  --hmmsearch PATH       hmmsearch binary (default: hmmsearch)
  --tigrfam PATH         TIGRfam HMM database
  --pfam PATH            Pfam-A HMM database
  --lcnn-weights DIR     PromoterLCNN weights directory (default: bundled)
  --fimo PATH            fimo binary (default: fimo)
  --motifs-dir DIR       Directory of .meme files (default: bundled)
  --fimo-threshold F     FIMO p-value threshold (default: 1e-4)

Conda environments:
  --conda-prokka ENV     Conda env for Prokka (blank = current env)
  --conda-hmm ENV        Conda env for hmmsearch (blank = current env)
  --conda-meme ENV       Conda env for MEME Suite (blank = current env)

Parameters:
  --threads N            Threads for external tools (default: 4)
  --kingdom STR          Prokka kingdom (default: Bacteria)
  --prefix STR           Prokka output prefix (default: genome)
  --igr-min N            Minimum IGR length in bp (default: 81)
  --igr-max N            Maximum IGR length in bp (default: 1000)
  --max-internal-dist N  Max gap within an operon in bp (default: 25)
  --min-flanking-dist N  Min gap at operon boundaries in bp (default: 75)
  --hmm-bitscore F       Minimum HMM bitscore (default: 25.0)
```

## Output structure

```
output/
├── prokka/                      # Prokka annotation files
├── igr/                         # Intergenic region extraction
│   ├── igr_summary.tsv
│   └── intergenic_regions.fasta
├── operons.tsv                  # Operon identification
├── hmm/                         # HMM marker screening (optional)
├── promoters.fasta              # Marker-filtered promoters
├── all_promoters.fasta          # All promoter-orientation IGRs
├── lcnn_predictions.tsv         # PromoterLCNN per-sequence predictions
├── promoter_markers_verified.tsv # Promoters confirmed by CNN
├── cds_annotations.tsv          # CDS product names from Prokka
├── fimo/                        # Motif scanning
│   ├── collectf/
│   ├── prodoric_2021.9/
│   ├── regtransbase/
│   └── fimo_combined.tsv
└── promoter_report.html         # Visual HTML report
```

## License

MIT
