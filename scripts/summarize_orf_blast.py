import pandas as pd

blast_file = "results/orf/albertii_orfs_vs_ecoli_candidates.tsv"
interval_file = "results/synteny/albertii_intervals.fasta.tsv"
synteny_file = "results/synteny/validated_flanking_synteny.tsv"

cols = [
    "orf_id", "hit_candidate", "identity", "alignment_length",
    "orf_length", "candidate_length", "qstart", "qend",
    "evalue", "bitscore"
]

df = pd.read_csv(blast_file, sep="\t", names=cols)

df["origin_candidate"] = df["orf_id"].str.split("|").str[1]
df["query_coverage"] = df["alignment_length"] / df["orf_length"]
df["target_coverage"] = df["alignment_length"] / df["candidate_length"]

df["match_type"] = (
    df["origin_candidate"] == df["hit_candidate"]
).map({
    True: "same_candidate",
    False: "other_candidate"
})

best = (
    df.sort_values("bitscore", ascending=False)
    .drop_duplicates("orf_id")
)

intervals = pd.read_csv(interval_file, sep="\t")
synteny = pd.read_csv(synteny_file, sep="\t")

same = df[df["match_type"] == "same_candidate"].copy()

summary = (
    same.sort_values("bitscore", ascending=False)
    .groupby("origin_candidate")
    .agg(
        number_of_ORFs=("orf_id", "nunique"),
        best_identity=("identity", "max"),
        best_alignment=("alignment_length", "max"),
        best_evalue=("evalue", "min"),
        best_bitscore=("bitscore", "max")
    )
    .reset_index()
)

best_same = (
    same.sort_values("bitscore", ascending=False)
    .drop_duplicates("origin_candidate")
)

summary = summary.merge(
    best_same[[
        "origin_candidate", "orf_id", "orf_length",
        "candidate_length", "query_coverage", "target_coverage",
        "qstart", "qend"
    ]],
    on="origin_candidate"
)

summary = summary.merge(
    intervals[["candidate_id", "length"]],
    left_on="origin_candidate",
    right_on="candidate_id",
    how="left"
)

summary = summary.merge(
    synteny[["candidate_id", "status"]],
    on="candidate_id",
    how="left"
)

summary.to_csv(
    "results/orf/same_candidate_hits_summary.tsv",
    sep="\t",
    index=False
)

print("Distinct candidates with same-candidate hits:", len(summary))
print("\nSummary:")
print(summary.to_string(index=False))

print("\nBest hits per ORF — other candidates:")
other = best[best["match_type"] == "other_candidate"]
print(
    other.sort_values("bitscore", ascending=False)[
        ["origin_candidate", "orf_id", "hit_candidate",
         "identity", "alignment_length", "orf_length",
         "candidate_length", "evalue", "bitscore"]
    ].head(15).to_string(index=False)
)
