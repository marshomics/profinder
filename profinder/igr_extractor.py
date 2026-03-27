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
    """Extract ID= value from GFF attributes column."""
    for attr in attributes.split(";"):
        if attr.startswith("ID="):
            return attr[3:]
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
    contigs = {rec.id: str(rec.seq) for rec in SeqIO.parse(str(fasta_path), "fasta")}

    # Parse CDS entries from the GFF
    genes = []
    with open(str(gff_path)) as fh:
        for line in fh:
            if line.startswith("#"):
                if line.startswith("##FASTA"):
                    break          # Prokka GFF embeds FASTA after this marker
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9 or cols[2] != "CDS":
                continue
            genes.append({
                "contig": cols[0],
                "start": int(cols[3]),
                "end": int(cols[4]),
                "strand": cols[6],
                "gene_id": _parse_gene_id(cols[8]),
            })

    if not genes:
        return pd.DataFrame()

    df = pd.DataFrame(genes).sort_values(["contig", "start"]).reset_index(drop=True)

    # Walk consecutive gene pairs on each contig
    rows = []
    igr_counter = 0
    for contig_id, grp in df.groupby("contig", sort=False):
        grp = grp.sort_values("start").reset_index(drop=True)
        contig_seq = contigs.get(contig_id, "")

        for i in range(len(grp) - 1):
            left = grp.iloc[i]
            right = grp.iloc[i + 1]

            igr_start = left["end"] + 1        # 1-based, inclusive
            igr_end = right["start"] - 1
            igr_len = igr_end - igr_start + 1

            if igr_len < size_min or igr_len > size_max:
                continue

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

    return pd.DataFrame(rows)
