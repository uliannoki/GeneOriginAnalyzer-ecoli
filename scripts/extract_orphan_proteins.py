import csv
from pathlib import Path
from Bio import SeqIO

candidate_file = Path(
    "results/homology/orphan_candidates_50aa.tsv"
)

protein_file = Path(
    "data/processed/all_proteins.faa"
)

output_file = Path(
    "data/processed/orphan_candidates.faa"
)

with open(candidate_file) as f:
    reader = csv.DictReader(f, delimiter="\t")
    candidate_ids = {
        row["protein_id"]
        for row in reader
    }

found = set()
written = 0

with open(output_file, "w") as out:
    for record in SeqIO.parse(protein_file, "fasta"):
        if record.id in candidate_ids:
            SeqIO.write(record, out, "fasta")
            found.add(record.id)
            written += 1

missing = candidate_ids - found

print("Candidates requested:", len(candidate_ids))
print("Sequences extracted:", written)
print("Missing sequences:", len(missing))
print("Output:", output_file)

if missing:
    print("First missing IDs:", list(missing)[:10])
