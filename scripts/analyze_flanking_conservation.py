import csv
from pathlib import Path

PROJECT = Path(__file__).resolve().parent.parent
BASE = PROJECT / "results/synteny"

NEIGHBORHOODS = BASE / "candidate_neighborhoods.tsv"
ANNOTATIONS = BASE / "albertii_cds_annotations.tsv"

PROTEIN_HITS = BASE / "flanking_proteins_vs_albertii_coords.tsv"
GENOME_HITS = BASE / "flanking_proteins_vs_albertii_genome_coords.tsv"

OUTPUT = BASE / "flanking_conservation_summary.tsv"
HITS_OUTPUT = BASE / "best_flanking_hits.tsv"

MIN_IDENTITY = 25.0
MIN_QUERY_COVERAGE = 50.0

# --- ИСПРАВЛЕНИЕ 1: Добавлены колонки ---
PROTEIN_COLUMNS = [
    "qseqid", "sseqid", "pident", "length",
    "qlen", "slen", "qstart", "qend",
    "sstart", "send", "evalue", "bitscore"
]

GENOME_COLUMNS = [
    "qseqid", "sseqid", "pident", "length",
    "qlen", "qstart", "qend", "sstart",
    "send", "sframe", "evalue", "bitscore"
]


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

# --- ИСПРАВЛЕНИЕ 2: Добавлена функция чтения BLAST ---
def read_blast_tsv(path, columns):
    with path.open(newline="") as handle:
        return list(
            csv.DictReader(
                handle,
                fieldnames=columns,
                delimiter="\t"
            )
        )


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


def clean_accession(value):
    parts = value.split("|")
    if len(parts) >= 3 and parts[1]:
        return parts[1]
    return value


def load_annotations():
    annotations = read_tsv(ANNOTATIONS)
    by_protein = {}

    for row in annotations:
        by_protein[row["protein_id"]] = row

    return annotations, by_protein


def best_protein_hits():
    best = {}

    # --- ИСПРАВЛЕНИЕ 3.1: Заменен вызов ---
    for row in read_blast_tsv(PROTEIN_HITS, PROTEIN_COLUMNS):
        qlen = int(row["qlen"])
        aligned = int(row["length"])
        coverage = 100 * aligned / qlen if qlen else 0
        identity = float(row["pident"])

        if identity < MIN_IDENTITY or coverage < MIN_QUERY_COVERAGE:
            continue

        query = row["qseqid"]
        subject = clean_accession(row["sseqid"])

        hit = {
            "query": query,
            "method": "BLASTP",
            "subject": subject,
            "identity": identity,
            "coverage": coverage,
            "evalue": float(row["evalue"]),
            "bitscore": float(row["bitscore"]),
        }

        if query not in best or hit["bitscore"] > best[query]["bitscore"]:
            best[query] = hit

    return best


def best_genome_hits(annotations):
    best = {}

    # --- ИСПРАВЛЕНИЕ 3.2: Заменен вызов ---
    for row in read_blast_tsv(GENOME_HITS, GENOME_COLUMNS):
        qlen = int(row["qlen"])
        aligned = int(row["length"])
        coverage = 100 * aligned / qlen if qlen else 0
        identity = float(row["pident"])

        if identity < MIN_IDENTITY or coverage < MIN_QUERY_COVERAGE:
            continue

        query = row["qseqid"]
        seqid = clean_accession(row["sseqid"])
        hit_start = min(int(row["sstart"]), int(row["send"]))
        hit_end = max(int(row["sstart"]), int(row["send"]))

        overlapping = []

        for feature in annotations:
            if feature["seqid"] != seqid:
                continue

            overlap = max(
                0,
                min(hit_end, int(feature["end"]))
                - max(hit_start, int(feature["start"])) + 1,
            )

            if overlap:
                overlapping.append((overlap, feature))

        matched_feature = (
            max(overlapping, key=lambda x: x[0])[1]
            if overlapping else None
        )

        hit = {
            "query": query,
            "method": "TBLASTN",
            "subject": seqid,
            "start": hit_start,
            "end": hit_end,
            "identity": identity,
            "coverage": coverage,
            "evalue": float(row["evalue"]),
            "bitscore": float(row["bitscore"]),
            "gene": matched_feature["gene"] if matched_feature else "",
            "protein_id": matched_feature["protein_id"] if matched_feature else "",
            "feature_start": matched_feature["start"] if matched_feature else "",
            "feature_end": matched_feature["end"] if matched_feature else "",
        }

        if query not in best or hit["bitscore"] > best[query]["bitscore"]:
            best[query] = hit

    return best


def main():
    neighborhoods = read_tsv(NEIGHBORHOODS)
    annotations, annotation_by_protein = load_annotations()

    protein_best = best_protein_hits()
    genome_best = best_genome_hits(annotations)

    candidate_groups = {}

    for row in neighborhoods:
        key = (row["sample"], row["candidate_protein_id"])
        candidate_groups.setdefault(key, []).append(row)

    summaries = []
    hit_rows = []

    for (sample, candidate_id), group in candidate_groups.items():
        candidate = next(
            row for row in group if row["is_candidate"] == "True"
        )

        left_rows = [
            row for row in group
            if int(row["relative_order"]) == -1
        ]
        right_rows = [
            row for row in group
            if int(row["relative_order"]) == 1
        ]

        flank_results = {}

        for side, flank_rows in [
            ("left", left_rows),
            ("right", right_rows),
        ]:
            if not flank_rows:
                flank_results[side] = None
                continue

            flank = flank_rows[0]
            query = f'{sample}__{flank["neighbor_protein_id"]}'

            protein_hit = protein_best.get(query)
            genome_hit = genome_best.get(query)

            annotation = None
            method = ""
            target_seqid = ""
            target_start = ""
            target_end = ""
            target_gene = ""

            if protein_hit:
                annotation = annotation_by_protein.get(
                    protein_hit["subject"]
                )

                if annotation:
                    method = "BLASTP"
                    target_seqid = annotation["seqid"]
                    target_start = int(annotation["start"])
                    target_end = int(annotation["end"])
                    target_gene = annotation["gene"]

            if annotation is None and genome_hit:
                method = "TBLASTN"
                target_seqid = genome_hit["subject"]
                target_start = genome_hit["start"]
                target_end = genome_hit["end"]
                target_gene = genome_hit["gene"]

            result = {
                "query": query,
                "method": method,
                "seqid": target_seqid,
                "start": target_start,
                "end": target_end,
                "gene": target_gene,
                "protein_hit": bool(protein_hit),
                "genome_hit": bool(genome_hit),
                "identity": (
                    protein_hit["identity"] if protein_hit
                    else genome_hit["identity"] if genome_hit
                    else ""
                ),
                "coverage": (
                    protein_hit["coverage"] if protein_hit
                    else genome_hit["coverage"] if genome_hit
                    else ""
                ),
            }

            flank_results[side] = result

            hit_rows.append({
                "sample": sample,
                "candidate_id": candidate_id,
                "side": side,
                **result,
            })

        left = flank_results["left"]
        right = flank_results["right"]

        same_replicon = False
        distance = ""
        flank_order = ""

        if left and right and left["seqid"] and right["seqid"]:
            same_replicon = left["seqid"] == right["seqid"]

            if same_replicon:
                left_mid = (int(left["start"]) + int(left["end"])) / 2
                right_mid = (int(right["start"]) + int(right["end"])) / 2

                flank_order = (
                    "left_before_right"
                    if left_mid < right_mid
                    else "right_before_left"
                )

                distance = max(
                    0,
                    max(int(left["start"]), int(right["start"]))
                    - min(int(left["end"]), int(right["end"])) - 1,
                )

        if left and right:
            if not left["seqid"] or not right["seqid"]:
                status = "both_flanks_not_mapped"
            elif not same_replicon:
                status = "different_replicons"
            elif distance <= 10000:
                status = "both_flanks_nearby"
            else:
                status = "both_flanks_far_apart"
        elif left and left["seqid"]:
            status = "left_flank_only"
        elif right and right["seqid"]:
            status = "right_flank_only"
        else:
            status = "no_flanks_mapped"

        summaries.append({
            "sample": sample,
            "candidate_id": candidate_id,
            "original_id": candidate["candidate_original_id"],
            "candidate_start": candidate["candidate_start"],
            "candidate_end": candidate["candidate_end"],
            "left_query": left["query"] if left else "",
            "left_method": left["method"] if left else "",
            "left_target_seqid": left["seqid"] if left else "",
            "left_target_start": left["start"] if left else "",
            "left_target_end": left["end"] if left else "",
            "left_target_gene": left["gene"] if left else "",
            "right_query": right["query"] if right else "",
            "right_method": right["method"] if right else "",
            "right_target_seqid": right["seqid"] if right else "",
            "right_target_start": right["start"] if right else "",
            "right_target_end": right["end"] if right else "",
            "right_target_gene": right["gene"] if right else "",
            "same_replicon": same_replicon,
            "flank_order": flank_order,
            "distance_between_flanks_bp": distance,
            "conservation_status": status,
        })

    summary_columns = list(summaries[0].keys()) if summaries else []
    hit_columns = list(hit_rows[0].keys()) if hit_rows else []

    write_tsv(OUTPUT, summaries, summary_columns)
    write_tsv(HITS_OUTPUT, hit_rows, hit_columns)

    print("Candidates analyzed:", len(summaries))
    print("Best flank hits saved:", len(hit_rows))
    print("\nConservation status:")

    counts = {}
    for row in summaries:
        status = row["conservation_status"]
        counts[status] = counts.get(status, 0) + 1

    for status, count in sorted(counts.items()):
        print(status, count)

    print("\nOutput:", OUTPUT)
    print("Output:", HITS_OUTPUT)


if __name__ == "__main__":
    main()
