# ProFinder

Identify and extract bacterial promoter sequences from a single genome FASTA file, annotate their associated coding sequences via BLAST, scan for known transcription factor binding motifs, and produce an HTML report.

The pipeline annotates a genome with [Prokka](https://github.com/tseemann/prokka), extracts intergenic regions (IGRs) between annotated genes, identifies operons using a two-pass proximity/flanking-distance algorithm, screens for marker genes via hmmsearch against TIGRfam and Pfam (HMM profiles are bundled), and outputs promoter sequences oriented 5'→3'. It then BLASTs the associated CDS proteins against Swiss-Prot for functional annotation and organism identification, scans promoters for known motifs with FIMO, and generates a visual HTML report.

All databases ship with the pipeline. No additional downloads are needed beyond installing the external tools.

## Requirements

**Python ≥ 3.7** with:

- pandas ≥ 1.3
- biopython ≥ 1.79

**External tools** (must be on `$PATH` or specified via CLI flags):

- [Prokka](https://github.com/tseemann/prokka)
- [BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) (for CDS annotation; uses NCBI remote BLAST by default, no local database required)
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

Minimal run (bundled HMMs and motifs, remote BLAST for CDS annotation):

```bash
profinder -i my_genome.fasta -o results/
```

With a local BLAST database (faster than remote):

```bash
profinder -i my_genome.fasta -o results/ \
    --blast-local-db /path/to/swissprot
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
| 10 | Annotate CDS via BLAST | `blast/cds_annotations.tsv` |
| 11 | Scan motifs with FIMO | `fimo/fimo_combined.tsv` |
| 12 | Generate HTML report | `promoter_report.html` |

Steps 4–5 use bundled TIGRfam and Pfam HMM profiles by default. You can override them with `--tigrfam` and `--pfam` if you have custom profiles.

Step 10 uses NCBI remote BLAST against Swiss-Prot by default. This requires an internet connection and can take several minutes depending on the number of CDS queries. Provide `--blast-local-db` to use a pre-formatted local database instead.

Step 11 uses FIMO from the MEME Suite. Three motif databases are bundled with the pipeline (CollecTF, PRODORIC, RegTransBase). You can point to a different directory with `--motifs-dir`.

## CDS annotation

Each promoter is associated with its immediately downstream CDS based on orientation. For CO_F promoters (both flanking genes on the + strand), the right gene is the associated CDS. For CO_R (both on the − strand), the left gene is the associated CDS, because that gene's transcription start site faces the IGR.

Prokka annotations are used only for gene coordinates. Functional annotation comes from BLASTing the translated CDS protein against Swiss-Prot, which provides a curated protein name and source organism for each hit. The pipeline extracts protein name and organism from the Swiss-Prot title format (`OS=` and `OX=` fields).

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
  --end N                Last step to run (default: 12)
  --list                 List steps and exit
  --force                Re-run all steps, ignoring checkpoints

Tool paths:
  --prokka PATH          Prokka binary (default: prokka)
  --hmmsearch PATH       hmmsearch binary (default: hmmsearch)
  --tigrfam PATH         TIGRfam HMM database
  --pfam PATH            Pfam-A HMM database
  --blastp PATH          blastp binary (default: blastp)
  --blast-db NAME        NCBI remote database (default: swissprot)
  --blast-local-db PATH  Local BLAST database (overrides remote)
  --blast-evalue F       BLAST e-value threshold (default: 1e-5)
  --fimo PATH            fimo binary (default: fimo)
  --motifs-dir DIR       Directory of .meme files (default: bundled)
  --fimo-threshold F     FIMO p-value threshold (default: 1e-4)

Conda environments:
  --conda-prokka ENV     Conda env for Prokka (blank = current env)
  --conda-hmm ENV        Conda env for hmmsearch (blank = current env)
  --conda-blast ENV      Conda env for BLAST+ (blank = current env)
  --conda-meme ENV       Conda env for MEME Suite (blank = current env)

Parameters:
  --threads N            Threads for external tools (default: 4)
  --kingdom STR          Prokka kingdom (default: Bacteria)
  --prefix STR           Prokka output prefix (default: genome)
  --igr-min N            Minimum IGR length in bp (default: 75)
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
├── blast/                       # CDS annotation
│   ├── query_proteins.faa
│   ├── blastp_results.tsv
│   └── cds_annotations.tsv
├── fimo/                        # Motif scanning
│   ├── collectf/
│   ├── prodoric_2021.9/
│   ├── regtransbase/
│   └── fimo_combined.tsv
└── promoter_report.html         # Visual HTML report
```

## License

MIT
