import csv
from pathlib import Path
from urllib.parse import unquote

PROJECT = Path(__file__).resolve().parent.parent

GFF = PROJECT / "data/raw/escherichia_albertii/unpacked/ncbi_dataset/data/GCF_028622335.1/genomic.gff"
OUTPUT = PROJECT / "results/synteny/albertii_cds_annotations.tsv"


def parse_attributes(text):
    result = {}

    for item in text.strip().split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = unquote(value)

    return result


rows = []

with GFF.open() as handle:
    for line in handle:
        if line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")

        if len(fields) != 9 or fields[2] != "CDS":
            continue

        attrs = parse_attributes(fields[8])

        rows.append({
            "seqid": fields[0],
            "start": int(fields[3]),
            "end": int(fields[4]),
            "strand": fields[6],
            "protein_id": attrs.get("protein_id", ""),
            "gene": attrs.get("gene", ""),
            "locus_tag": attrs.get("locus_tag", ""),
            "product": attrs.get("product", ""),
        })

with OUTPUT.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "seqid", "start", "end", "strand",
            "protein_id", "gene", "locus_tag", "product"
        ],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(rows)

print("CDS records:", len(rows))
print("Unique protein IDs:", len(set(r["protein_id"] for r in rows)))
print("Output:", OUTPUT)
