import csv
from pathlib import Path

PROJECT = Path(__file__).resolve().parent.parent
INPUT = PROJECT / "results/synteny/candidate_neighborhoods.tsv"
OUT_DIR = PROJECT / "results/synteny"

MOBILE_KEYWORDS = [
    "integrase",
    "transposase",
    "insertion sequence",
    "IS element",
    "prophage",
    "phage",
    "bacteriophage",
    "conjugative",
    "plasmid",
    "mobile element",
    "recombinase",
    "transposon",
    "integron",
    "restriction-modification",
]

def read_tsv(path):
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
    rows = read_tsv(INPUT)
    groups = {}

    for row in rows:
        key = (row["sample"], row["candidate_protein_id"])
        groups.setdefault(key, []).append(row)

    summaries = []

    for (sample, candidate_id), group in groups.items():
        group.sort(key=lambda x: int(x["relative_order"]))

        candidate = next(
            row for row in group if row["is_candidate"] == "True"
        )

        neighbors = [
            row for row in group if row["is_candidate"] != "True"
        ]

        annotated_neighbors = [
            row for row in neighbors if row["neighbor_product"].strip()
        ]

        mobile_neighbors = []

        for row in neighbors:
            annotation = " ".join([
                row["neighbor_product"],
                row["neighbor_gene"],
                row["neighbor_locus_tag"],
            ]).lower()

            matched = [
                keyword for keyword in MOBILE_KEYWORDS
                if keyword.lower() in annotation
            ]

            if matched:
                mobile_neighbors.append({
                    "relative_order": row["relative_order"],
                    "protein_id": row["neighbor_protein_id"],
                    "gene": row["neighbor_gene"],
                    "product": row["neighbor_product"],
                    "keywords": ";".join(matched),
                })

        left = [
            row for row in neighbors
            if int(row["relative_order"]) < 0
        ]
        right = [
            row for row in neighbors
            if int(row["relative_order"]) > 0
        ]

        nearest_left = max(
            left, key=lambda x: int(x["relative_order"])
        ) if left else {}

        nearest_right = min(
            right, key=lambda x: int(x["relative_order"])
        ) if right else {}

        mobile_description = " | ".join(
            f'{x["relative_order"]}:{x["gene"] or x["protein_id"]} '
            f'({x["product"]})'
            for x in mobile_neighbors
        )

        summaries.append({
            "sample": sample,
            "candidate_protein_id": candidate_id,
            "original_id": candidate["candidate_original_id"],
            "seqid": candidate["seqid"],
            "start": candidate["candidate_start"],
            "end": candidate["candidate_end"],
            "strand": candidate["candidate_strand"],
            "candidate_product": candidate["neighbor_product"],
            "annotated_neighbor_count": len(annotated_neighbors),
            "mobile_element_neighbor_count": len(mobile_neighbors),
            "mobile_element_flag": bool(mobile_neighbors),
            "mobile_element_annotations": mobile_description,
            "nearest_left_gene": nearest_left.get("neighbor_gene", ""),
            "nearest_left_product": nearest_left.get("neighbor_product", ""),
            "nearest_right_gene": nearest_right.get("neighbor_gene", ""),
            "nearest_right_product": nearest_right.get("neighbor_product", ""),
        })

    columns = [
        "sample", "candidate_protein_id", "original_id",
        "seqid", "start", "end", "strand", "candidate_product",
        "annotated_neighbor_count", "mobile_element_neighbor_count",
        "mobile_element_flag", "mobile_element_annotations",
        "nearest_left_gene", "nearest_left_product",
        "nearest_right_gene", "nearest_right_product",
    ]

    summaries.sort(
        key=lambda x: (
            x["mobile_element_flag"] != True,
            x["sample"],
            x["start"],
        )
    )

    write_tsv(
        OUT_DIR / "candidate_context_summary.tsv",
        summaries,
        columns,
    )

    flagged = [
        row for row in summaries
        if row["mobile_element_flag"]
    ]

    write_tsv(
        OUT_DIR / "mobile_context_candidates.tsv",
        flagged,
        columns,
    )

    print("Candidates summarized:", len(summaries))
    print("Candidates near possible mobile elements:", len(flagged))
    print("Output: candidate_context_summary.tsv")
    print("Output: mobile_context_candidates.tsv")

    print("\nCounts by sample:")
    for sample in sorted(set(row["sample"] for row in summaries)):
        sample_rows = [
            row for row in summaries if row["sample"] == sample
        ]
        sample_flagged = sum(
            row["mobile_element_flag"] for row in sample_rows
        )
        print(
            sample,
            "total:", len(sample_rows),
            "mobile context:", sample_flagged,
        )

if __name__ == "__main__":
    main()
