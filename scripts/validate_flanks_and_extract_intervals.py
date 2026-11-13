import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from pathlib import Path

base = Path("results/synteny")
data = Path("data/raw/escherichia_albertii/unpacked")

candidates = pd.read_csv(base / "nearby_flank_candidates.tsv", sep="\t")
annotations = pd.read_csv(base / "albertii_cds_annotations.tsv", sep="\t")

genome_file = next(data.rglob("*genomic.fna"))
genomes = {r.id: r.seq for r in SeqIO.parse(genome_file, "fasta")}

validated = []
intervals = []

for _, r in candidates.iterrows():
    selected = {}

    for side in ["left", "right"]:
        seqid = r[f"{side}_target_seqid"]
        start = int(r[f"{side}_target_start"])
        end = int(r[f"{side}_target_end"])

        hits = annotations[
            (annotations.seqid == seqid) &
            (annotations.start <= end) &
            (annotations.end >= start)
        ].copy()

        if hits.empty:
            selected[side] = None
            continue

        hits["overlap"] = (
            hits[["end"]].min(axis=1).clip(upper=end) -
            hits[["start"]].max(axis=1).clip(lower=start) + 1
        )

        best = hits.sort_values("overlap", ascending=False).iloc[0]
        selected[side] = best

    if selected["left"] is None or selected["right"] is None:
        continue

    left = selected["left"]
    right = selected["right"]

    same_replicon = left.seqid == right.seqid
    same_strand = left.strand == right.strand

    if not same_replicon:
        status = "different_replicons"
        distance = None
    else:
        distance = max(0, max(left.start, right.start) -
                       min(left.end, right.end) - 1)

        if left.end < right.start:
            order = "left_before_right"
        elif right.end < left.start:
            order = "right_before_left"
        else:
            order = "overlapping"

        status = "same_order_same_strand" if (
            order == "left_before_right" and same_strand
        ) else "rearranged_or_opposite_orientation"

        genome = genomes[left.seqid]

        if order == "left_before_right":
            start = int(left.end)
            end = int(right.start)
            strand = left.strand
        elif order == "right_before_left":
            start = int(right.end)
            end = int(left.start)
            strand = right.strand
        else:
            start = end = None
            strand = None

        if start is not None and end > start:
            seq = genome[start:end]

            if strand == "-":
                seq = seq.reverse_complement()

            intervals.append({
                "candidate_id": r.candidate_id,
                "sample": r["sample"],
                "seqid": left.seqid,
                "start": start + 1,
                "end": end,
                "length": len(seq),
                "strand": strand,
                "sequence": str(seq)
            })

    validated.append({
        "candidate_id": r.candidate_id,
        "sample": r["sample"],
        "left_protein": left.protein_id,
        "left_gene": left.gene,
        "left_product": left["product"],
        "left_start": left.start,
        "left_end": left.end,
        "left_strand": left.strand,
        "right_protein": right.protein_id,
        "right_gene": right.gene,
        "right_product": right["product"],
        "right_start": right.start,
        "right_end": right.end,
        "right_strand": right.strand,
        "same_replicon": same_replicon,
        "same_strand": same_strand,
        "distance": distance,
        "status": status
    })

pd.DataFrame(validated).to_csv(
    base / "validated_flanking_synteny.tsv", sep="\t", index=False
)

pd.DataFrame(intervals).to_csv(
    base / "albertii_intervals.fasta.tsv", sep="\t", index=False
)

with open("results/orf/albertii_intervals.fna", "w") as f:
    for i, x in enumerate(intervals, 1):
        f.write(f">{x['candidate_id']}|{x['seqid']}|{x['start']}-{x['end']}|{x['strand']}\n")
        f.write(x["sequence"] + "\n")

print("Validated candidates:", len(validated))
print("Extracted intervals:", len(intervals))
print("\nSynteny statuses:")
print(pd.Series([x["status"] for x in validated]).value_counts().to_string())
print("\nInterval lengths:")
print(pd.Series([x["length"] for x in intervals]).describe().to_string())
