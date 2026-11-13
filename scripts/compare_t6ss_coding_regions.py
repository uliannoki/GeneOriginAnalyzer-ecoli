from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
import csv

INPUT = Path("results/synteny/t6ss_intervals_aligned.fasta")
OUTPUT = Path("results/synteny/t6ss_coding_comparison.tsv")

records = {}

for record in SeqIO.parse(INPUT, "fasta"):
    sequence = str(record.seq).upper()
    records[record.id] = sequence

print("=== SEQUENCE VALIDATION ===")

rows = []

for sample, aligned in records.items():
    ungapped = aligned.replace("-", "")

    valid_chars = set("ATGC")
    invalid = sorted(set(ungapped) - valid_chars)

    divisible_by_3 = len(ungapped) % 3 == 0

    if invalid:
        status = "invalid_nucleotides"
        protein = ""
    elif not divisible_by_3:
        status = "length_not_multiple_of_3"
        protein = ""
    else:
        protein = str(Seq(ungapped).translate())
        internal_stops = protein[:-1].count("*")

        if internal_stops:
            status = "internal_stop_codons"
        else:
            status = "translation_possible"

    print(
        f"\n{sample}\n"
        f"Aligned length: {len(aligned)}\n"
        f"Ungapped length: {len(ungapped)}\n"
        f"Length mod 3: {len(ungapped) % 3}\n"
        f"Status: {status}\n"
        f"Protein length: {len(protein)}"
    )

    rows.append({
        "sample": sample,
        "aligned_length": len(aligned),
        "ungapped_length": len(ungapped),
        "length_mod_3": len(ungapped) % 3,
        "status": status,
        "protein_length": len(protein),
        "internal_stops": protein[:-1].count("*") if protein else "",
    })

OUTPUT.parent.mkdir(parents=True, exist_ok=True)

with OUTPUT.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "sample",
            "aligned_length",
            "ungapped_length",
            "length_mod_3",
            "status",
            "protein_length",
            "internal_stops",
        ],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(rows)

print(f"\nSaved: {OUTPUT}")
