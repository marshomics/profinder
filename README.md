# ProFinder

ProFinder extracts high-confidence constitutive promoter candidates from bacterial or archaeal genome assemblies. Given a FASTA file, it returns a curated shortlist of promoter sequences upstream of conserved housekeeping genes, ready for use in expression constructs, reporter systems, or metabolic engineering.

<p align="center">
  <img src="pipeline.png" width="100%">
</p>

## How the pipeline works

ProFinder is not a genome-wide promoter annotation tool. It applies a series of biologically motivated filters to produce a short, high-confidence list of promoter sequences upstream of single-copy phylogenetic marker genes: ribosomal proteins, tRNA synthetases, DNA replication components, and similar housekeeping functions. These genes are constitutively expressed because the cell cannot afford to silence them, so their promoters are strong candidates for driving reliable expression.

The pipeline annotates a genome with [Prokka](https://github.com/tseemann/prokka), extracts intergenic regions between annotated CDS features, identifies operons using a two-pass proximity/flanking-distance algorithm, and screens for marker genes via hmmsearch against TIGRfam and Pfam (profiles are bundled). Candidate promoters are then classified by the domain-appropriate CNN. Promoters are scanned for transcription factor binding motifs with FIMO, and a visual HTML report is generated.

## Bacteria vs. archaea

The `--domain` flag controls which promoter classifier is used.

`--domain bacteria` (the default) runs [PromoterLCNN](https://github.com/occasumlux/Promoters), a two-stage CNN trained on experimentally verified *E. coli* K-12 promoters from RegulonDB. Stage 1 separates promoters from non-promoters. Stage 2 assigns a sigma-factor subtype (σ70, σ24, σ28, σ32, σ38, σ54). Only the 3'-terminal 81 nt of each intergenic region (closest to the transcription start site) is classified. The HTML report filters to σ70 promoters.

`--domain archaea` runs [iProm-Archaea](https://github.com/PromoterTools/iPromArchaea), a CNN that uses 6-mer frequency encoding to perform binary promoter/non-promoter classification. Archaea lack bacterial sigma factors, so there is no subtype assignment. Each intergenic region is scanned in non-overlapping 100 bp windows from 5' to 3', matching the original tool's design. An IGR is called a promoter if any window scores positive. This means longer intergenic regions are fully scanned rather than truncated to a fixed-length suffix.

When `--domain archaea` is specified, Prokka's `--kingdom` is automatically set to `Archaea` (override with `--kingdom` if needed).

All other pipeline steps (Prokka annotation, IGR extraction, operon identification, HMM marker screening, FIMO motif scanning, report generation) run identically for both domains.

## Why an *E. coli*-trained classifier works across bacterial species

PromoterLCNN was trained on *E. coli* K-12 σ70 promoters, but the core recognition logic (the −10 TATAAT and −35 TTGACA hexamers, the 17 ± 1 bp spacer, the AT-rich upstream element) is conserved across bacteria. In benchmarking across eleven species spanning seven phyla, ProFinder returned useful promoter shortlists (13 to 33 σ70 marker promoters) for every species with intergenic GC content below 50%, including organisms as distant from *E. coli* as *Bacillus subtilis* (Firmicutes), *Synechocystis* sp. PCC 6803 (Cyanobacteria), and *Vibrio cholerae* (Vibrionaceae).

For organisms with intergenic GC content above approximately 60% (*Pseudomonas aeruginosa*, *Caulobacter crescentus*, *Mycobacterium tuberculosis*), the pipeline returns very few or no σ70 marker promoters. The AT-rich motifs that define σ70 promoters are genuinely sparse in high-GC intergenic regions. The classifier rejects these sequences rather than misclassifying them, so output remains reliable but incomplete.

## Requirements

**Python ≥ 3.7** with pandas ≥ 1.3, biopython ≥ 1.79, and tensorflow ≥ 2.6.

**External tools** (must be on `$PATH` or specified via CLI flags):

- [Prokka](https://github.com/tseemann/prokka)
- [MEME Suite](https://meme-suite.org/) (specifically `fimo`)
- [HMMER](http://hmmer.org/)

HMM profiles, motif databases, and CNN weights for both classifiers are bundled with the pipeline.

## Installation

```bash
pip install .
```

Or for development:

```bash
pip install -e .
```

## Quick start

Bacterial genome (default):

```bash
profinder -i my_genome.fasta -o results/
```

Archaeal genome:

```bash
profinder -i my_genome.fasta -o results/ --domain archaea
```

Append the first 90 nt of each downstream CDS to promoter sequences:

```bash
profinder -i my_genome.fasta -o results/ --cds-bp 90
```

Override bundled HMMs:

```bash
profinder -i my_genome.fasta -o results/ \
    --tigrfam /path/to/custom_tigrfam.hmm \
    --pfam /path/to/custom_Pfam-A.hmm
```

The pipeline checkpoints at each step. Re-running the same command skips completed steps. Use `--force` to re-run everything.

## Pipeline steps

| Step | Description | Key output |
|------|-------------|------------|
| 1 | Run Prokka | `prokka/genome.gff`, `prokka/genome.faa` |
| 2 | Extract intergenic regions | `igr/igr_summary.tsv`, `igr/intergenic_regions.fasta` |
| 3 | Identify operons | `operons.tsv` |
| 4 | Run hmmsearch (TIGRfam + Pfam) | `hmm/hmm_tigrfam.tblout`, `hmm/hmm_pfam.tblout` |
| 5 | Filter HMM output | `hmm/hmm_filtered.tsv` |
| 6 | Filter operons + add markers | `operons_filtered_markers.tsv` |
| 7 | Match IGRs to marker operons | `promoter_markers.tsv` |
| 8 | Extract marker-filtered promoters | `promoters.fasta`, `promoters_short.fasta` |
| 9 | Extract all promoter-orientation IGRs | `all_promoters.fasta`, `all_promoters_short.fasta` |
| 10 | Predict promoters (CNN classifier) | `lcnn_predictions.tsv`, `promoter_markers_verified.tsv` |
| 11 | Annotate CDS (Prokka) | `cds_annotations.tsv` |
| 12 | Scan motifs (FIMO) | `fimo/fimo_combined.tsv` |
| 13 | Build final output table | `profinder_results.tsv` |
| 14 | Generate HTML report | `promoter_report.html` |

Step 10 runs the domain-appropriate CNN on all promoter-orientation IGRs (not just marker-filtered ones). For bacteria, PromoterLCNN classifies the 3'-terminal 81 nt through its two-stage cascade. For archaea, iProm-Archaea scans the full IGR in 100 bp non-overlapping windows and calls it a promoter if any window scores positive. Pre-trained weights for both classifiers are bundled in `weights/`.

Steps 8 and 9 also produce CDS-extended FASTA files (`promoters_cds_bp.fasta`, `promoters_short_cds_bp.fasta`, `all_promoters_cds_bp.fasta`, `all_promoters_short_cds_bp.fasta`) when `--cds-bp` is set to a value greater than 0.

## CDS extension (`--cds-bp`)

By default, FASTA outputs contain promoter sequences only. The `--cds-bp N` option appends the first N nucleotides of the downstream CDS to each promoter sequence, producing four additional FASTA files alongside the standard ones. The final results table (`profinder_results.tsv`) also gains a `sequence_5p_to_3p_cds` column with the extended sequences.

Strand handling: for CO_F promoters the CDS sits on the + strand immediately after the IGR, so the first N bp are read directly from the forward strand starting at the CDS start coordinate. For CO_R promoters the CDS is on the − strand immediately before the IGR, so the first N bp of the coding sequence are extracted by taking N bp upstream of the IGR start on the forward strand and reverse-complementing them. In both cases the result is the first N coding nucleotides in the 5'→3' reading direction, appended to the 5'→3' promoter sequence.

Set `--cds-bp 0` (the default) to disable this feature and produce no additional files.

## IGR orientation logic

Each intergenic region is classified by the orientation of its flanking genes:

- **CO_F** (co-oriented forward): both flanking genes on the + strand. The IGR is a promoter candidate for the right (downstream) gene. Sequence is already 5'→3'.
- **CO_R** (co-oriented reverse): both flanking genes on the − strand. The IGR is a promoter candidate for the left (downstream) gene, which runs right to left. Sequence is reverse-complemented to 5'→3'.
- **DP** (divergent promoter): genes point away from each other (← IGR →). Contains promoters for both flanking genes, so orientation is ambiguous. Excluded from FASTA output.
- **CONV** (convergent): genes point toward each other (→ IGR ←). This is a terminator region. Excluded.

## CDS annotation

Each promoter is associated with its immediately downstream CDS based on orientation: the right gene for CO_F, the left gene for CO_R. Functional annotation uses the `product=`, `gene=`, and `locus_tag=` fields from Prokka's GFF output.

## Motif scanning

FIMO scans each promoter sequence (oriented 5'→3') against position weight matrices from three transcription factor databases: CollecTF (experimentally validated binding sites), PRODORIC (prokaryotic gene regulation), and RegTransBase (regulatory interactions). The HTML report includes inline SVG diagrams showing motif positions along each promoter, colour-coded by database.

## Operon identification

Operons are identified in two passes. First, consecutive CDS on the same contig within `--max-internal-dist` bp (default 25) are clustered. Second, clusters are kept only if both upstream and downstream gaps exceed `--min-flanking-dist` bp (default 75).

## CLI reference

```
profinder -i FASTA -o DIR [options]

Required:
  -i, --input              Input genome FASTA
  -o, --output             Output directory (default: output/)

Domain:
  --domain {bacteria,archaea}
                           Classifier to use (default: bacteria)

Step control:
  --start N                First step to run (default: 1)
  --end N                  Last step to run (default: 14)
  --list                   List steps and exit
  --force                  Re-run all steps, ignoring checkpoints

Tool paths:
  --prokka PATH            Prokka binary (default: prokka)
  --hmmsearch PATH         hmmsearch binary (default: hmmsearch)
  --tigrfam PATH           TIGRfam HMM database (default: bundled)
  --pfam PATH              Pfam-A HMM database (default: bundled)
  --lcnn-weights DIR       PromoterLCNN weights directory (default: bundled)
  --ipromarchaea-weights PATH
                           iProm-Archaea .h5 weights file (default: bundled)
  --fimo PATH              fimo binary (default: fimo)
  --motifs-dir DIR         Directory of .meme files (default: bundled)
  --fimo-threshold F       FIMO p-value threshold (default: 1e-4)

Conda environments:
  --conda-prokka ENV       Conda env for Prokka
  --conda-hmm ENV          Conda env for hmmsearch
  --conda-meme ENV         Conda env for MEME Suite

Parameters:
  --threads N              Threads for external tools (default: 4)
  --kingdom STR            Prokka kingdom (default: Bacteria; auto-set to
                           Archaea when --domain archaea)
  --prefix STR             Prokka output prefix (default: genome)
  --igr-min N              Minimum IGR length in bp (default: 81)
  --igr-max N              Maximum IGR length in bp (default: 1000)
  --max-internal-dist N    Max gap within an operon in bp (default: 25)
  --min-flanking-dist N    Min gap at operon boundaries in bp (default: 75)
  --hmm-bitscore F         Minimum HMM bitscore (default: 25.0)
  --cds-bp N               Append first N nt of downstream CDS to promoter
                           sequences in additional FASTA files and the final
                           table (default: 0, disabled)
```

## Output structure

```
output/
├── prokka/                          # Prokka annotation files
├── igr/                             # Intergenic region extraction
│   ├── igr_summary.tsv
│   └── intergenic_regions.fasta
├── operons.tsv                      # Operon identification
├── hmm/                             # HMM marker screening
├── promoters.fasta                  # Marker-filtered promoters (5'→3')
├── promoters_short.fasta            # Same, short headers
├── all_promoters.fasta              # All promoter-orientation IGRs (5'→3')
├── all_promoters_short.fasta        # Same, short headers
├── promoters_cds_bp.fasta           # Promoter + CDS extension (if --cds-bp > 0)
├── promoters_short_cds_bp.fasta     # Same, short headers
├── all_promoters_cds_bp.fasta       # All IGRs + CDS extension (if --cds-bp > 0)
├── all_promoters_short_cds_bp.fasta # Same, short headers
├── lcnn_predictions.tsv             # Per-sequence CNN predictions
├── promoter_markers_verified.tsv    # Marker promoters confirmed by CNN
├── cds_annotations.tsv              # CDS product names from Prokka
├── fimo/                            # Motif scanning results
│   ├── collectf/
│   ├── prodoric_2021.9/
│   ├── regtransbase/
│   └── fimo_combined.tsv
├── profinder_results.tsv            # Comprehensive results table
└── promoter_report.html             # Visual HTML report
```

The `profinder_results.tsv` table contains one row per promoter-orientation IGR with columns for contig coordinates, associated CDS annotation (gene name, locus tag, product), marker status, CNN prediction, sigma-factor subtype (bacteria) or binary classification (archaea), FIMO motif hits, the full-length 5'→3' sequence, and (when `--cds-bp > 0`) the promoter+CDS extended sequence.

## Bundled data

```
profinder/
├── hmms/
│   ├── tigrfam.hmm          # TIGRfam marker gene profiles
│   └── Pfam-A.hmm           # Pfam-A marker gene profiles
├── motifs/
│   ├── collectf.meme         # CollecTF TF binding sites
│   ├── prodoric_2021.9.meme  # PRODORIC motifs
│   └── regtransbase.meme     # RegTransBase motifs
└── weights/
    ├── PromoterLCNN/         # Bacterial CNN (TensorFlow SavedModel)
    │   ├── IsPromoter_fold_5/
    │   └── PromotersOnly_fold_1/
    └── iPromArchaea/         # Archaeal CNN (Keras weights)
        └── model_cnn.weights.h5
```

## License

MIT (pipeline code). The iProm-Archaea model weights are from [PromoterTools/iPromArchaea](https://github.com/PromoterTools/iPromArchaea) under MIT (code) / CC-BY-NC-ND (original training data).
