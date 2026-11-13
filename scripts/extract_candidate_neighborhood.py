import csv
from pathlib import Path
from urllib.parse import unquote

PROJECT = Path(__file__).resolve().parent.parent

GFF_FILES = {
    "ecoli_MG1655": PROJECT / "data/raw/ecoli_MG1655/unpacked/ncbi_dataset/data/GCF_000005845.2/genomic.gff",
    "ecoli_W3110": PROJECT / "data/raw/ecoli_W3110/unpacked/ncbi_dataset/data/GCA_000010245.1/genomic.gff",
    "ecoli_REL606": PROJECT / "data/raw/ecoli_REL606/unpacked/ncbi_dataset/data/GCA_000017985.1/genomic.gff",
    "ecoli_Sakai": PROJECT / "data/raw/ecoli_Sakai/unpacked/ncbi_dataset/data/GCF_000008865.2/genomic.gff",
}

CANDIDATES_FILE = PROJECT / "results/homology/no_albertii_hits.tsv"
OUT_DIR = PROJECT / "results/synteny"

ATTRIBUTES = [
    "protein_id",
    "gene",
    "locus_tag",
    "product",
    "Name",
    "ID",
    "Parent",
]


def parse_attributes(text):
    result = {}

    for item in text.strip().split(";"):
        if "=" not in item:
            continue

        key, value = item.split("=", 1)
        result[key] = unquote(value)

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

            seqid, source, feature, start, end, score, strand, phase, attr_text = fields
            attrs = parse_attributes(attr_text)

            record = {
                "seqid": seqid,
                "start": int(start),
                "end": int(end),
                "strand": strand,
                "phase": phase,
            }

            for key in ATTRIBUTES:
                record[key] = attrs.get(key, "")

            if record["protein_id"]:
                records.append(record)

    records.sort(key=lambda x: (x["seqid"], x["start"], x["end"]))
    return records


def read_candidates(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, columns):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    candidates = read_candidates(CANDIDATES_FILE)

    coordinates = []
    neighborhoods = []
    missing = []

    for sample, gff_path in GFF_FILES.items():
        print(f"Reading {sample}...")

        if not gff_path.exists():
            print(f"WARNING: GFF not found: {gff_path}")
            continue

        records = read_gff(gff_path)

        by_protein = {}
        for index, record in enumerate(records):
            by_protein.setdefault(record["protein_id"], []).append(index)

        sample_candidates = [
            c for c in candidates if c["sample"] == sample
        ]

        for candidate in sample_candidates:
            original_id = candidate["original_id"]
            matches = by_protein.get(original_id, [])

            if not matches:
                missing.append({
                    "sample": sample,
                    "protein_id": candidate["protein_id"],
                    "original_id": original_id,
                    "length": candidate["length"],
                    "reason": "protein_id_not_found_in_GFF",
                })
                continue

            for match_number, index in enumerate(matches, start=1):
                record = records[index]

                candidate_row = {
                    "sample": sample,
                    "protein_id": candidate["protein_id"],
                    "original_id": original_id,
                    "length": candidate["length"],
                    "seqid": record["seqid"],
                    "start": record["start"],
                    "end": record["end"],
                    "strand": record["strand"],
                    "gene": record["gene"],
                    "locus_tag": record["locus_tag"],
                    "product": record["product"],
                    "match_number": match_number,
                    "number_of_gff_matches": len(matches),
                }
                coordinates.append(candidate_row)

                same_sequence = [
                    (i, r) for i, r in enumerate(records)
                    if r["seqid"] == record["seqid"]
                ]

                local_position = next(
                    j for j, (global_index, _) in enumerate(same_sequence)
                    if global_index == index
                )

                left = max(0, local_position - 5)
                right = min(len(same_sequence), local_position + 6)

                for neighbor_position in range(left, right):
                    global_index, neighbor = same_sequence[neighbor_position]
                    relative_order = neighbor_position - local_position

                    if relative_order < 0:
                        gap = record["start"] - neighbor["end"] - 1
                    elif relative_order > 0:
                        gap = neighbor["start"] - record["end"] - 1
                    else:
                        gap = 0

                    neighborhoods.append({
                        "sample": sample,
                        "candidate_protein_id": candidate["protein_id"],
                        "candidate_original_id": original_id,
                        "candidate_start": record["start"],
                        "candidate_end": record["end"],
                        "candidate_strand": record["strand"],
                        "relative_order": relative_order,
                        "is_candidate": relative_order == 0,
                        "neighbor_protein_id": neighbor["protein_id"],
                        "neighbor_start": neighbor["start"],
                        "neighbor_end": neighbor["end"],
                        "neighbor_strand": neighbor["strand"],
                        "neighbor_gene": neighbor["gene"],
                        "neighbor_locus_tag": neighbor["locus_tag"],
                        "neighbor_product": neighbor["product"],
                        "gap_bp": gap,
                        "seqid": neighbor["seqid"],
                    })

    coordinate_columns = [
        "sample", "protein_id", "original_id", "length",
        "seqid", "start", "end", "strand", "gene",
        "locus_tag", "product", "match_number",
        "number_of_gff_matches",
    ]

    neighborhood_columns = [
        "sample", "candidate_protein_id", "candidate_original_id",
        "candidate_start", "candidate_end", "candidate_strand",
        "relative_order", "is_candidate",
        "neighbor_protein_id", "neighbor_start", "neighbor_end",
        "neighbor_strand", "neighbor_gene", "neighbor_locus_tag",
        "neighbor_product", "gap_bp", "seqid",
    ]

    missing_columns = [
        "sample", "protein_id", "original_id", "length", "reason"
    ]

    write_tsv(
        OUT_DIR / "candidate_coordinates.tsv",
        coordinates,
        coordinate_columns,
    )

    write_tsv(
        OUT_DIR / "candidate_neighborhoods.tsv",
        neighborhoods,
        neighborhood_columns,
    )

    write_tsv(
        OUT_DIR / "candidates_missing_from_gff.tsv",
        missing,
        missing_columns,
    )

    print("\nDONE")
    print("Candidates in input:", len(candidates))
    print("Coordinate records:", len(coordinates))
    print("Neighborhood records:", len(neighborhoods))
    print("Candidates not found in GFF:", len(missing))
    print("Output directory:", OUT_DIR)


if __name__ == "__main__":
    main()
