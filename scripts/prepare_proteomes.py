from Bio import SeqIO
from pathlib import Path
import csv

samples = [
    "ecoli_MG1655",
    "ecoli_W3110",
    "ecoli_REL606",
    "ecoli_Sakai",
]

base = Path("data/raw")
output_dir = Path("data/processed")

output_dir.mkdir(parents=True, exist_ok=True)

combined_records = []
mapping = []

for sample in samples:
    sample_dir = base / sample / "unpacked"
    protein_files = list(sample_dir.rglob("*.faa"))

    if not protein_files:
        print(f"WARNING: No protein file for {sample}")
        continue

    if len(protein_files) > 1:
        print(f"WARNING: Multiple protein files for {sample}")

    protein_file = protein_files[0]
    records = list(SeqIO.parse(protein_file, "fasta"))

    print(f"{sample}: {len(records)} proteins")

    for record in records:
        original_id = record.id
        new_id = f"{sample}__{original_id}"

        record.id = new_id
        record.name = new_id
        record.description = ""

        combined_records.append(record)

        mapping.append({
            "protein_id": new_id,
            "sample": sample,
            "original_id": original_id,
            "length": len(record.seq),
        })

fasta_path = output_dir / "all_proteins.faa"

SeqIO.write(combined_records, fasta_path, "fasta")

mapping_path = output_dir / "protein_mapping.tsv"

with open(mapping_path, "w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=[
            "protein_id",
            "sample",
            "original_id",
            "length",
        ],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(mapping)

print("\nTOTAL")
print("Proteins:", len(combined_records))
print("FASTA:", fasta_path)
print("Mapping:", mapping_path)
