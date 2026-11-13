import csv
from collections import defaultdict
from pathlib import Path

blast_file = Path(
    "results/homology/orphans_vs_albertii_genome.tsv"
)

candidate_file = Path(
    "results/homology/orphan_candidates_50aa.tsv"
)

output_file = Path(
    "results/homology/albertii_genome_hit_summary.tsv"
)

hits = defaultdict(list)

with open(blast_file) as f:
    for line in f:
        fields = line.strip().split("\t")

        if len(fields) != 11:
            continue

        (
            query,
            subject,
            identity,
            alignment_length,
            qlen,
            qstart,
            qend,
            sstart,
            send,
            evalue,
            bitscore,
        ) = fields

        identity = float(identity)
        alignment_length = int(alignment_length)
        qlen = int(qlen)
        evalue = float(evalue)
        bitscore = float(bitscore)

        query_coverage = alignment_length / qlen

        if evalue <= 1e-5 and query_coverage >= 0.50:
            hits[query].append({
                "subject": subject,
                "identity": identity,
                "alignment_length": alignment_length,
                "query_coverage": query_coverage,
                "qstart": qstart,
                "qend": qend,
                "sstart": sstart,
                "send": send,
                "evalue": evalue,
                "bitscore": bitscore,
            })

with open(candidate_file) as f:
    candidates = list(csv.DictReader(f, delimiter="\t"))

columns = [
    "protein_id",
    "original_id",
    "sample",
    "length",
    "genome_hit",
    "number_of_hits",
    "best_identity",
    "best_coverage",
    "best_evalue",
    "best_bitscore",
    "best_subject",
    "subject_start",
    "subject_end",
]

with open(output_file, "w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=columns,
        delimiter="\t",
    )

    writer.writeheader()

    for candidate in candidates:
        protein_id = candidate["protein_id"]
        candidate_hits = hits[protein_id]

        best = (
            max(candidate_hits, key=lambda x: x["bitscore"])
            if candidate_hits
            else None
        )

        writer.writerow({
            "protein_id": protein_id,
            "original_id": candidate["original_id"],
            "sample": candidate["sample"],
            "length": candidate["length"],
            "genome_hit": bool(candidate_hits),
            "number_of_hits": len(candidate_hits),
            "best_identity": best["identity"] if best else "",
            "best_coverage": round(best["query_coverage"], 3) if best else "",
            "best_evalue": best["evalue"] if best else "",
            "best_bitscore": best["bitscore"] if best else "",
            "best_subject": best["subject"] if best else "",
            "subject_start": best["sstart"] if best else "",
            "subject_end": best["send"] if best else "",
        })

with_hits = sum(bool(hits[x["protein_id"]]) for x in candidates)

print("Total candidates:", len(candidates))
print("With genomic hits:", with_hits)
print("Without genomic hits:", len(candidates) - with_hits)
print("Saved:", output_file)
