from pathlib import Path
import pandas as pd

INPUT = Path("results/homology/final_candidate_proteome_search.tsv")
OUTPUT = Path("results/homology/final_candidate_proteome_search_annotated.tsv")

columns = [
    "qseqid", "sseqid", "pident", "length",
    "qlen", "slen", "qstart", "qend",
    "sstart", "send", "evalue", "bitscore"
]

df = pd.read_csv(INPUT, sep="\t", names=columns)

if df.empty:
    print("BLAST returned no hits.")
    raise SystemExit(0)

df["query_coverage_pct"] = (
    (df["qend"] - df["qstart"] + 1) / df["qlen"] * 100
)

df["subject_coverage_pct"] = (
    (df["send"] - df["sstart"] + 1) / df["slen"] * 100
)

df["query_sample"] = df["qseqid"].str.split("__").str[0]
df["hit_sample"] = df["sseqid"].str.split("__").str[0]

df["relationship"] = "other_sample"

df.loc[df["qseqid"] == df["sseqid"], "relationship"] = "self_hit"

df.loc[
    (df["query_sample"] == df["hit_sample"])
    & (df["qseqid"] != df["sseqid"]),
    "relationship"
] = "same_sample"

df["coverage_class"] = "limited_coverage"

df.loc[
    (df["query_coverage_pct"] >= 50)
    & (df["subject_coverage_pct"] >= 50),
    "coverage_class"
] = "both_coverage_ge_50"

df = df.sort_values(
    ["qseqid", "bitscore"],
    ascending=[True, False]
)

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(OUTPUT, sep="\t", index=False)

for query in sorted(df["qseqid"].unique()):
    subset = df[df["qseqid"] == query].copy()

    print("\n" + "=" * 90)
    print("QUERY:", query)
    print("Query length:", int(subset.iloc[0]["qlen"]), "aa")
    print("Total BLAST hits:", len(subset))

    print("\nHIT COUNTS:")
    print(subset["relationship"].value_counts().to_string())

    nonself = subset[subset["relationship"] != "self_hit"]

    print("\nTOP NON-SELF HITS:")
    display_columns = [
        "sseqid", "hit_sample", "pident", "length",
        "qlen", "slen",
        "query_coverage_pct", "subject_coverage_pct",
        "evalue", "bitscore", "relationship"
    ]

    if nonself.empty:
        print("No non-self hits.")
    else:
        print(
            nonself[display_columns]
            .head(15)
            .to_string(index=False)
        )

    substantial = nonself[
        (nonself["query_coverage_pct"] >= 50)
        & (nonself["subject_coverage_pct"] >= 50)
        & (nonself["evalue"] <= 1e-5)
    ]

    print("\nHITS WITH BOTH COVERAGES >=50% AND E <= 1e-5:")
    if substantial.empty:
        print("None.")
    else:
        print(
            substantial[display_columns]
            .to_string(index=False)
        )

print("\nSaved:", OUTPUT)
