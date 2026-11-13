from pathlib import Path
import pandas as pd

RAW = Path("data/raw")
OUTPUT = Path("results/reports/all_ecoli_cds_annotations.tsv")

rows = []

for sample_dir in sorted(RAW.glob("ecoli_*")):
    sample = sample_dir.name

    for gff in sample_dir.rglob("*.gff"):
        with open(gff) as handle:
            for line in handle:
                if line.startswith("#"):
                    continue

                fields = line.rstrip("\n").split("\t")

                if len(fields) != 9 or fields[2] != "CDS":
                    continue

                seqid, source, feature, start, end, score, strand, phase, attr = fields

                attributes = {}

                for item in attr.split(";"):
                    if "=" in item:
                        key, value = item.split("=", 1)
                        attributes[key] = value

                protein_id = attributes.get("protein_id")

                if not protein_id:
                    continue

                rows.append({
                    "sample": sample,
                    "protein_id": f"{sample}__{protein_id}",
                    "original_id": protein_id,
                    "seqid": seqid,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                    "gene": attributes.get("gene", ""),
                    "locus_tag": attributes.get("locus_tag", ""),
                    "product": attributes.get("product", ""),
                    "gff_file": str(gff),
                })

result = pd.DataFrame(rows)

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
result.to_csv(OUTPUT, sep="\t", index=False)

print("CDS records:", len(result))
print("Unique protein IDs:", result["protein_id"].nunique())
print("Samples:")
print(result.groupby("sample")["protein_id"].nunique())
print("Saved:", OUTPUT)
