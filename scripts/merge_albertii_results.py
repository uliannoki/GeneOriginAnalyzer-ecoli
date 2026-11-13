import csv
from pathlib import Path

protein_file = Path(
    "results/homology/albertii_hit_summary.tsv"
)

genome_file = Path(
    "results/homology/albertii_genome_hit_summary.tsv"
)

output_file = Path(
    "results/homology/albertii_combined_summary.tsv"
)

with open(protein_file) as f:
    protein_rows = {
        row["protein_id"]: row
        for row in csv.DictReader(f, delimiter="\t")
    }

with open(genome_file) as f:
    genome_rows = {
        row["protein_id"]: row
        for row in csv.DictReader(f, delimiter="\t")
    }

columns = [
    "protein_id",
    "original_id",
    "sample",
    "length",
    "protein_hit",
    "genome_hit",
    "combined_category",
]

with open(output_file, "w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=columns,
        delimiter="\t",
    )

    writer.writeheader()

    for protein_id, protein in protein_rows.items():
        genome = genome_rows[protein_id]

        protein_hit = protein["albertii_hit"] == "True"
        genome_hit = genome["genome_hit"] == "True"

        if protein_hit:
            category = "protein_level_hit"
        elif genome_hit:
            category = "genome_only_hit"
        else:
            category = "no_detected_hit"

        writer.writerow({
            "protein_id": protein_id,
            "original_id": protein["original_id"],
            "sample": protein["sample"],
            "length": protein["length"],
            "protein_hit": protein_hit,
            "genome_hit": genome_hit,
            "combined_category": category,
        })

print("Total:", len(protein_rows))
print("Saved:", output_file)
