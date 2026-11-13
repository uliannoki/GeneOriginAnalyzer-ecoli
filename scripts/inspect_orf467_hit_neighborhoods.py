from pathlib import Path
import csv

GFF_FILES = {
    "ecoli_Sakai": Path(
        "data/raw/ecoli_Sakai/unpacked/ncbi_dataset/data/GCF_000008865.2/genomic.gff"
    ),
    "ecoli_W3110": Path(
        "data/raw/ecoli_W3110/unpacked/ncbi_dataset/data/GCA_000010245.1/genomic.gff"
    ),
    "ecoli_MG1655": Path(
        "data/raw/ecoli_MG1655/unpacked/ncbi_dataset/data/GCF_000005845.2/genomic.gff"
    ),
    "ecoli_REL606": Path(
        "data/raw/ecoli_REL606/unpacked/ncbi_dataset/data/GCA_000017985.1/genomic.gff"
    ),
}

TARGET_IDS = {
    "NP_308259.1",
    "NP_309765.1",
    "ACT38889.1",
    "BAA36104.1",
    "NP_415752.1",
}

OUTPUT = Path("results/orf/orf467_hit_neighborhoods.tsv")


def parse_attributes(text):
    result = {}
    for item in text.strip().split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
    return result


def read_gff(path):
    records = []

    with path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) != 9 or fields[2] != "CDS":
                continue

            attrs = parse_attributes(fields[8])

            record = {
                "seqid": fields[0],
                "start": int(fields[3]),
                "end": int(fields[4]),
                "strand": fields[6],
                "protein_id": attrs.get("protein_id", ""),
                "gene": attrs.get("gene", ""),
                "locus_tag": attrs.get("locus_tag", ""),
                "product": attrs.get("product", ""),
            }

            records.append(record)

    records.sort(key=lambda x: (x["seqid"], x["start"]))
    return records


def main():
    output_rows = []

    for sample, path in GFF_FILES.items():
        print(f"\nReading {sample}...")

        if not path.exists():
            print(f"WARNING: file not found: {path}")
            continue

        records = read_gff(path)

        print(f"CDS records: {len(records)}")

        for target_id in sorted(TARGET_IDS):
            indices = [
                i for i, record in enumerate(records)
                if record["protein_id"] == target_id
            ]

            for index in indices:
                target = records[index]

                print(
                    f"\n{sample}: {target_id} "
                    f"{target['gene']} {target['product']}"
                )

                same_replicon = [
                    (i, record)
                    for i, record in enumerate(records)
                    if record["seqid"] == target["seqid"]
                ]

                local_position = next(
                    j for j, (i, _) in enumerate(same_replicon)
                    if i == index
                )

                left = max(0, local_position - 5)
                right = min(len(same_replicon), local_position + 6)

                for j in range(left, right):
                    genome_index, record = same_replicon[j]

                    if genome_index == index:
                        relation = "TARGET"
                    elif genome_index < index:
                        relation = "LEFT"
                    else:
                        relation = "RIGHT"

                    print(
                        f"{relation:6} "
                        f"{record['start']}-{record['end']} "
                        f"{record['strand']:2} "
                        f"{record['protein_id']:20} "
                        f"{record['gene']:12} "
                        f"{record['product']}"
                    )

                    output_rows.append({
                        "sample": sample,
                        "target_id": target_id,
                        "relation": relation,
                        **record,
                    })

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)

    columns = [
        "sample", "target_id", "relation", "seqid",
        "start", "end", "strand", "protein_id",
        "gene", "locus_tag", "product"
    ]

    with OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(output_rows)

    print("\nDONE")
    print(f"Neighborhood records: {len(output_rows)}")
    print(f"Output: {OUTPUT}")


if __name__ == "__main__":
    main()
