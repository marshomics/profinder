# ProFinder

A bacterial promoter identification pipeline that extracts high-confidence constitutive (σ70) promoter candidates from any genome assembly. Given a FASTA file, ProFinder returns a curated shortlist of promoter sequences upstream of conserved housekeeping genes, ready for integration into genetic system development workflows.

## What ProFinder does

ProFinder is built for a specific use case: identifying constitutive promoters that are likely to drive reliable expression in a microbe of interest. It is not a genome-wide promoter annotation tool. Instead of returning thousands of candidate sites with varying confidence, ProFinder applies a series of biologically motivated filters to produce a short, high-confidence list of σ70 (constitutive) promoter sequences upstream of single-copy phylogenetic marker genes (ribosomal proteins, tRNA synthetases, DNA replication components, and other essential housekeeping functions). These marker genes are known to be expressed across growth conditions because the cell cannot afford to silence them, and their promoters are strong candidates for use in expression constructs, reporter systems, and metabolic engineering.

The pipeline annotates a genome with [Prokka](https://github.com/tseemann/prokka), extracts intergenic regions (IGRs) between annotated genes, identifies operons using a two-pass proximity/flanking-distance algorithm, and screens for marker genes via hmmsearch against TIGRfam and Pfam (HMM profiles are bundled). Candidate promoters are classified with [PromoterLCNN](https://github.com/occasumlux/Promoters), a two-stage CNN that distinguishes promoter from non-promoter sequences and assigns a sigma-factor subtype (σ70, σ24, σ28, σ38, σ32, σ54). Promoters are scanned for known transcription factor binding motifs with FIMO, and a visual HTML report is generated. All databases ship with the pipeline.

## Why an *E. coli*-trained classifier works across species

PromoterLCNN was trained on experimentally verified *E. coli* K-12 promoter sequences from RegulonDB. The *E. coli* σ70 promoter is one of the best-characterised regulatory elements in bacterial biology: the −10 (TATAAT) and −35 (TTGACA) hexamers, the 17 ± 1 bp spacer, and the AT-rich upstream element have been validated by decades of mutagenesis, footprinting, and structural studies. A classifier trained on these features is anchored to known biochemistry rather than statistical patterns of uncertain biological meaning.

Because the core σ70 recognition logic is conserved across bacteria, these learned features transfer well to other species. In benchmarking across eleven species spanning seven phyla, ProFinder returned useful promoter shortlists (13–33 σ70 marker promoters) for every species with intergenic GC content below 50%, including organisms as distant from *E. coli* as *Bacillus subtilis* (Firmicutes), *Synechocystis* sp. PCC 6803 (Cyanobacteria), and *Vibrio cholerae* (Vibrionaceae). The critical variable is base composition, not phylogenetic distance: the AT-rich motifs that define σ70 promoters are present in low-to-moderate-GC genomes regardless of taxonomic position.

## GC content limitation

For organisms with intergenic GC content above approximately 60%, the pipeline will return fewer σ70 marker promoters. This is because the AT-rich features the classifier relies on are sparse in high-GC intergenic regions. The classifier rejects these sequences rather than misclassifying them, so the output remains reliable (no false positives).

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
| 13 | Build final output table | `profinder_results.tsv` |
| 14 | Generate HTML report | `promoter_report.html` |

Steps 4–5 use bundled TIGRfam and Pfam HMM profiles by default. You can override them with `--tigrfam` and `--pfam` if you have custom profiles.

Step 10 runs PromoterLCNN on all promoter-orientation IGRs (not just marker-filtered ones). The 3'-terminal 81 nt of each sequence (closest to the transcription start site) is passed through a two-stage CNN: a binary classifier filters non-promoters, and a sigma-factor classifier assigns one of six subtypes to each confirmed promoter. Sequences shorter than 81 nt are classified as non-promoters. Pre-trained weights are bundled in `weights/PromoterLCNN/`.

Step 11 extracts gene names, locus tags, and protein product names from Prokka's GFF annotations for each promoter's associated CDS.

Step 12 uses FIMO from the MEME Suite. Three motif databases are bundled with the pipeline (CollecTF, PRODORIC, RegTransBase). You can point to a different directory with `--motifs-dir`.

Step 13 joins data from all previous steps into a single TSV (`profinder_results.tsv`) covering every promoter-orientation IGR. Columns include promoter ID, contig coordinates, associated CDS annotation (gene name, locus tag, product), marker status, LCNN prediction and sigma-factor subtype, FIMO motif hits, and the full-length 5'→3' sequence.

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
  --end N                Last step to run (default: 14)
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
├── profinder_results.tsv        # Comprehensive results table
└── promoter_report.html         # Visual HTML report
```

## License

MIT
