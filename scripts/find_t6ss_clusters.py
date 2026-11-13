from pathlib import Path
from urllib.parse import unquote
import csv

GFF_FILES = {
    "ecoli_Sakai": "data/raw/ecoli_Sakai/unpacked/ncbi_dataset/data/GCF_000008865.2/genomic.gff",
    "ecoli_W3110": "data/raw/ecoli_W3110/unpacked/ncbi_dataset/data/GCA_000010245.1/genomic.gff",
    "ecoli_MG1655": "data/raw/ecoli_MG1655/unpacked/ncbi_dataset/data/GCF_000005845.2/genomic.gff",
    "ecoli_REL606": "data/raw/ecoli_REL606/unpacked/ncbi_dataset/data/GCA_000017985.1/genomic.gff",
    "escherichia_albertii": "data/raw/escherichia_albertii/unpacked/ncbi_dataset/data/GCF_028622335.1/genomic.gff",
}

KEYWORDS = [
    "type vi secretion",
    "type 6 secretion",
    "tss",
    "imp",
    "hcp",
    "vgr",
    "paar",
    "clpv",
    "tagh",
    "vas",
]

MAX_GAP = 30000
OUTPUT = Path("results/synteny/t6ss_related_regions.tsv")


def parse_attributes(text):
    result = {}
    for item in text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = unquote(value)
    return result


def read_cds(path):
    records = []

    with open(path) as handle:
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

            text = (
                record["gene"] + " " +
                record["product"] + " " +
                record["locus_tag"]
            ).lower()

            if any(keyword in text for keyword in KEYWORDS):
                records.append(record)

    return records


def make_clusters(records):
    clusters = []

    by_seqid = {}
    for record in records:
        by_seqid.setdefault(record["seqid"], []).append(record)

    for seqid, items in by_seqid.items():
        items.sort(key=lambda x: x["start"])

        current = []

        for record in items:
            if not current:
                current = [record]
                continue

            gap = record["start"] - current[-1]["end"]

            if gap <= MAX_GAP:
                current.append(record)
            else:
                clusters.append((seqid, current))
                current = [record]

        if current:
            clusters.append((seqid, current))

    return clusters


def main():
    all_rows = []

    for sample, filename in GFF_FILES.items():
        path = Path(filename)
        print(f"\n{'=' * 70}")
        print(sample)

        if not path.exists():
            print("GFF NOT FOUND")
            continue

        records = read_cds(path)
        clusters = make_clusters(records)

        print(f"Related CDS: {len(records)}")
        print(f"Regions: {len(clusters)}")

        for number, (seqid, cluster) in enumerate(clusters, 1):
            if len(cluster) < 2:
                continue

            start = min(x["start"] for x in cluster)
            end = max(x["end"] for x in cluster)

            print(
                f"\nRegion {number}: {seqid}:{start}-{end} "
                f"({len(cluster)} CDS)"
            )

            for record in cluster:
                print(
                    f"{record['start']}-{record['end']} "
                    f"{record['strand']} "
                    f"{record['protein_id']} "
                    f"{record['gene']} "
                    f"{record['product']}"
                )

                all_rows.append({
                    "sample": sample,
                    "region": number,
                    **record,
                })

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)

    columns = [
        "sample", "region", "seqid", "start", "end",
        "strand", "protein_id", "gene", "locus_tag", "product"
    ]

    with OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(all_rows)

    print("\nDONE")
    print(f"Output: {OUTPUT}")


if __name__ == "__main__":
    main()
