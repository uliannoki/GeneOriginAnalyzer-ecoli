import csv
from pathlib import Path
from Bio import SeqIO

PROJECT = Path(__file__).resolve().parent.parent

NEIGHBORHOODS = PROJECT / "results/synteny/candidate_neighborhoods.tsv"
PROTEOME = PROJECT / "data/processed/all_proteins.faa"
OUTPUT = PROJECT / "data/processed/flanking_proteins.faa"

wanted = set()

with NEIGHBORHOODS.open(newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")

    for row in reader:
        if row["is_candidate"] == "True":
            continue

        original_id = row["neighbor_protein_id"]
        sample = row["sample"]
        protein_id = f"{sample}__{original_id}"
        wanted.add(protein_id)

found = set()

with OUTPUT.open("w") as out:
    for record in SeqIO.parse(PROTEOME, "fasta"):
        if record.id in wanted:
            SeqIO.write(record, out, "fasta")
            found.add(record.id)

missing = wanted - found

print("Unique flanking proteins requested:", len(wanted))
print("Proteins extracted:", len(found))
print("Missing proteins:", len(missing))

if missing:
    print("First missing IDs:")
    for protein_id in sorted(missing)[:20]:
        print(protein_id)

print("Output:", OUTPUT)
