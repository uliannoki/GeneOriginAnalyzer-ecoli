from pathlib import Path
import pandas as pd

ANNOTATION_FILE = Path("results/reports/all_ecoli_cds_annotations.tsv")

annotations = pd.read_csv(ANNOTATION_FILE, sep="\t")

TARGETS = {
    "NP_310248.1": [
        "ecoli_Sakai__NP_310248.1",
        "ecoli_Sakai__NP_309366.1",
        "ecoli_Sakai__NP_311896.1",
        "ecoli_Sakai__NP_311893.2",
        "ecoli_Sakai__NP_310818.1",
        "ecoli_Sakai__NP_310251.1",
    ],
    "NP_303550.1": [
        "ecoli_Sakai__NP_303550.1",
        "ecoli_Sakai__NP_309554.1",
        "ecoli_Sakai__NP_313020.1",
    ],
}

WINDOW = 5000

OUTPUT = Path("results/reports/candidate_homolog_neighborhoods.tsv")

records = []

for group, protein_ids in TARGETS.items():
    for protein_id in protein_ids:
        matches = annotations[
            annotations["protein_id"] == protein_id
        ]

        if matches.empty:
            print(f"NOT FOUND: {protein_id}")
            continue

        gene = matches.iloc[0]

        nearby = annotations[
            (annotations["sample"] == gene["sample"])
            & (annotations["seqid"] == gene["seqid"])
            & (annotations["start"] <= gene["end"] + WINDOW)
            & (annotations["end"] >= gene["start"] - WINDOW)
        ].copy()

        nearby["distance_from_target"] = nearby.apply(
            lambda row: (
                row["start"] - gene["end"]
                if row["start"] > gene["end"]
                else gene["start"] - row["end"]
                if row["end"] < gene["start"]
                else 0
            ),
            axis=1
        )

        nearby = nearby.sort_values("start")

        print(f"\n=== {group}: {protein_id} ===")
        print(
            nearby[
                [
                    "protein_id", "start", "end", "strand",
                    "gene", "locus_tag", "product",
                    "distance_from_target"
                ]
            ].to_string(index=False)
        )

        for _, row in nearby.iterrows():
            record = row.to_dict()
            record["candidate_group"] = group
            record["target_protein"] = protein_id
            records.append(record)

result = pd.DataFrame(records)

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
result.to_csv(OUTPUT, sep="\t", index=False)

print("\nSaved:", OUTPUT)
print("Neighborhood records:", len(result))
