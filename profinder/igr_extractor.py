"""
Extract intergenic regions (IGRs) directly from a Prokka GFF and the
corresponding genome FASTA.

This replaces PIGGY for the single-genome case. PIGGY is designed for
pangenome analysis across multiple genomes and requires Roary output;
here we parse gene coordinates from the GFF, compute intergenic gaps,
classify each by the orientation of the flanking genes, and pull the
sequence from the FASTA.

Orientation labels mirror PIGGY's convention:

    CO_F   → ← IGR → →   co-oriented forward (downstream gene on + strand)
    CO_R   ← ← IGR ← ←   co-oriented reverse (downstream gene on − strand)
    DP     ← IGR →        divergent promoter  (genes point away from IGR)
    CONV   → IGR ←        convergent / terminator (genes point into IGR)
"""

import pandas as pd
from Bio import SeqIO
import gzip



def _open_maybe_gzip(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def _classify_orientation(left_strand: str, right_strand: str) -> str:
    """Return PIGGY-compatible orientation label for a gene pair."""
    if left_strand == "-" and right_strand == "+":
        return "DP"
    if left_strand == "+" and right_strand == "-":
        return "CONV"
    if left_strand == "+" and right_strand == "+":
        return "CO_F"
    # left == "-" and right == "-"
    return "CO_R"


def _parse_gene_id(attributes: str) -> str:
    """Extract a gene identifier from GFF attributes.

    Tries, in order: ID=, locus_tag=, gene=, Name= (stripping any
    trailing ' gene' or ' CDS' suffix from Geneious-style names).
    """
    for key in ("ID=", "locus_tag=", "gene=", "Name="):
        for attr in attributes.split(";"):
            if attr.startswith(key):
                val = attr[len(key):]
                for suffix in (" gene", " CDS"):
                    if val.endswith(suffix):
                        val = val[:-len(suffix)]
                return val
    return ""


def extract_igrs(gff_path, fasta_path, size_min=75, size_max=1000):
    """Parse a Prokka GFF and extract intergenic regions.

    Parameters
    ----------
    gff_path : str or Path
        Prokka GFF3 file.
    fasta_path : str or Path
        Genome FASTA (the .fna from Prokka, or the original input).
    size_min, size_max : int
        Keep IGRs whose length falls in [size_min, size_max].

    Returns
    -------
    pd.DataFrame
        Columns: igr_id, contig, start, end, length, orientation,
                 left_gene, right_gene, sequence
    """
    # Load contig sequences
    # contigs = {rec.id: str(rec.seq) for rec in SeqIO.parse(str(fasta_path), "fasta")}
    with _open_maybe_gzip(fasta_path) as fh:
        contigs = {rec.id: str(rec.seq) for rec in SeqIO.parse(fh, "fasta")}
        
    # Parse gene-like entries from the GFF. Boundary features are CDS
    # plus all non-coding RNA classes that Prokka emits when invoked with
    # --rfam (Infernal/Rfam covariance models): tRNA, rRNA, tmRNA, ncRNA,
    # misc_RNA. Treating these as boundaries prevents IGRs that overlap
    # or sit inside a non-coding gene from reaching the motif scan.
    BOUNDARY_FEATURES = {
        "CDS", "tRNA", "rRNA", "tmRNA", "ncRNA", "misc_RNA",
    }
    genes = []
    feature_spans = []   # for overlap exclusion of IGR windows
    with open(str(gff_path)) as fh:
        for line in fh:
            if line.startswith("#"):
                if line.startswith("##FASTA"):
                    break          # Prokka GFF embeds FASTA after this marker
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9:
                continue
            ftype = cols[2]
            if ftype not in BOUNDARY_FEATURES:
                continue
            entry = {
                "contig": cols[0],
                "start": int(cols[3]),
                "end": int(cols[4]),
                "strand": cols[6],
                "gene_id": _parse_gene_id(cols[8]),
                "feature_type": ftype,
            }
            genes.append(entry)
            feature_spans.append(entry)

    if not genes:
        return pd.DataFrame()

    df = pd.DataFrame(genes).sort_values(["contig", "start"]).reset_index(drop=True)

    # Build per-contig sorted span lists for quick overlap rejection.
    spans_by_contig: dict[str, list[tuple[int, int]]] = {}
    for f in feature_spans:
        spans_by_contig.setdefault(f["contig"], []).append((f["start"], f["end"]))
    for k in spans_by_contig:
        spans_by_contig[k].sort()

    # Walk consecutive gene pairs on each contig.  We track the farthest
    # end seen so far (``max_end``) rather than the previous gene's end,
    # so nested/contained genes don't produce a bogus IGR inside an
    # enclosing gene's body.
    rows = []
    igr_counter = 0
    for contig_id, grp in df.groupby("contig", sort=False):
        grp = grp.sort_values("start").reset_index(drop=True)
        contig_seq = contigs.get(contig_id, "")

        max_end = -1                         # 1-based; -1 means no gene seen yet
        left_gene_for_gap = None             # the gene whose end == max_end
        for i in range(len(grp)):
            curr = grp.iloc[i]

            # If the current gene is entirely inside the span already
            # covered by an earlier gene, skip it — no IGR to emit and
            # we don't want to update max_end downward.
            if curr["end"] <= max_end:
                continue

            if left_gene_for_gap is not None:
                igr_start = max_end + 1        # 1-based, inclusive
                igr_end = curr["start"] - 1
                igr_len = igr_end - igr_start + 1

                # Defensive: reject IGR windows that overlap any annotated
                # feature on the same contig. With --rfam on, this catches
                # Rfam-only features (sRNAs, riboswitches, etc.) whose
                # coordinates fall between two CDS without bounding them.
                contig_spans = spans_by_contig.get(contig_id, [])
                overlaps_feature = any(
                    not (f_end < igr_start or f_start > igr_end)
                    and not (f_start == left_gene_for_gap["start"]
                             and f_end == left_gene_for_gap["end"])
                    and not (f_start == curr["start"]
                             and f_end == curr["end"])
                    for f_start, f_end in contig_spans
                )

                if size_min <= igr_len <= size_max and not overlaps_feature:
                    left = left_gene_for_gap
                    right = curr
                    orientation = _classify_orientation(left["strand"], right["strand"])

                    # Extract sequence (GFF is 1-based; Python slicing is 0-based)
                    seq = contig_seq[igr_start - 1 : igr_end] if contig_seq else ""

                    igr_counter += 1
                    rows.append({
                        "igr_id": f"igr_{igr_counter:06d}",
                        "contig": contig_id,
                        "start": igr_start,
                        "end": igr_end,
                        "length": igr_len,
                        "orientation": orientation,
                        "left_gene": left["gene_id"],
                        "right_gene": right["gene_id"],
                        "sequence": seq,
                    })

            # Advance the right-boundary tracker
            max_end = curr["end"]
            left_gene_for_gap = curr

    return pd.DataFrame(rows)
