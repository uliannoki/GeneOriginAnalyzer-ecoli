
import pandas as pd

orfs = pd.read_csv(
    "results/orf/albertii_orfs.tsv", sep="\t"
)

intervals = pd.read_csv(
    "results/synteny/albertii_intervals.fasta.tsv", sep="\t"
)

annotations = pd.read_csv(
    "results/synteny/albertii_cds_annotations.tsv", sep="\t"
)

selected_ids = [
    "ORF405|ecoli_REL606__ACT39985.1|+|1522237-1522428",
    "ORF467|ecoli_Sakai__NP_308259.1|-|4007735-4007911"
]

seqids = dict(zip(
    intervals["candidate_id"],
    intervals["seqid"]
))

orfs["orf_id"] = [
    f"ORF{i+1}|{row.candidate_id}|{row.strand}|"
    f"{row.genomic_start}-{row.genomic_end}"
    for i, row in orfs.iterrows()
]

results = []

for orf_id in selected_ids:
    matches = orfs[orfs["orf_id"] == orf_id]

    if matches.empty:
        print("ORF not found:", orf_id)
        continue

    o = matches.iloc[0]

    seqid = seqids.get(o["candidate_id"])
    start = int(o["genomic_start"])
    end = int(o["genomic_end"])
    strand = o["strand"]

    nearby = annotations[
        (annotations["seqid"] == seqid)
        & (annotations["start"] <= end + 5000)
        & (annotations["end"] >= start - 5000)
    ].copy()

    nearby["overlap_bp"] = (
        nearby["end"].clip(upper=end)
        - nearby["start"].clip(lower=start)
        + 1
    ).clip(lower=0)

    nearby["same_strand"] = (
        nearby["strand"] == strand
    )

    def relation(row):
        if row["overlap_bp"] > 0:
            return "overlapping"
        if row["end"] < start:
            return "upstream"
        return "downstream"

    nearby["relation"] = nearby.apply(relation, axis=1)

    def distance(row):
        if row["overlap_bp"] > 0:
            return 0
        if row["end"] < start:
            return start - int(row["end"]) - 1
        return int(row["start"]) - end - 1

    nearby["distance_bp"] = nearby.apply(distance, axis=1)

    nearby = nearby.sort_values(
        ["distance_bp", "start"]
    )

    print("\n" + "=" * 80)
    print("ORF:", orf_id)
    print("Replicon:", seqid)
    print("Coordinates:", start, "-", end)
    print("Strand:", strand)
    print("Length:", end - start + 1, "bp")

    print("\nNearby annotated CDS:")

    if nearby.empty:
        print("No annotated CDS within 5 kb")
        continue

    columns = [
        "protein_id", "gene", "product",
        "start", "end", "strand",
        "relation", "distance_bp",
        "overlap_bp", "same_strand"
    ]

    print(nearby[columns].to_string(index=False))

    for _, a in nearby.iterrows():
        results.append({
            "orf_id": orf_id,
            "candidate_id": o["candidate_id"],
            "seqid": seqid,
            "orf_start": start,
            "orf_end": end,
            "orf_strand": strand,
            "protein_id": a["protein_id"],
            "gene": a["gene"],
            "product": a["product"],
            "cds_start": a["start"],
            "cds_end": a["end"],
            "cds_strand": a["strand"],
            "relation": a["relation"],
            "distance_bp": a["distance_bp"],
            "overlap_bp": a["overlap_bp"],
            "same_strand": a["same_strand"]
        })

pd.DataFrame(results).to_csv(
    "results/orf/interesting_orf_neighborhoods.tsv",
    sep="\t",
    index=False
)

print("\nSaved:")
print("results/orf/interesting_orf_neighborhoods.tsv")

