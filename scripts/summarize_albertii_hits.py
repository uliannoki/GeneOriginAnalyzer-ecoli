import csv
from collections import defaultdict
from pathlib import Path

blast_file = Path(
    "results/homology/orphans_vs_albertii.tsv"
)

candidate_file = Path(
    "results/homology/orphan_candidates_50aa.tsv"
)

output_file = Path(
    "results/homology/albertii_hit_summary.tsv"
)

hits = defaultdict(list)

with open(blast_file) as f:
    for line in f:
        fields = line.strip().split("\t")

        if len(fields) != 8:
            continue

        query, subject, identity, aln_len, qlen, slen, evalue, bitscore = fields

        identity = float(identity)
        aln_len = int(aln_len)
        qlen = int(qlen)
        slen = int(slen)
        evalue = float(evalue)

        qcov = aln_len / qlen
        scov = aln_len / slen

        if (
            evalue <= 1e-5
            and qcov >= 0.50
            and scov >= 0.50
        ):
            hits[query].append({
                "subject": subject,
                "identity": identity,
                "evalue": evalue,
                "qcov": qcov,
                "scov": scov,
            })

with open(candidate_file) as f:
    candidates = list(csv.DictReader(f, delimiter="\t"))

with open(output_file, "w", newline="") as f:
    columns = [
        "protein_id",
        "original_id",
        "sample",
        "length",
        "albertii_hit",
        "number_of_hits",
        "best_identity",
        "best_evalue",
        "best_subject",
    ]

    writer = csv.DictWriter(
        f,
        fieldnames=columns,
        delimiter="\t",
    )

    writer.writeheader()

    for candidate in candidates:
        protein_id = candidate["protein_id"]
        candidate_hits = hits[protein_id]

        best = max(
            candidate_hits,
            key=lambda x: -x["evalue"],
        ) if candidate_hits else None

        writer.writerow({
            "protein_id": protein_id,
            "original_id": candidate["original_id"],
            "sample": candidate["sample"],
            "length": candidate["length"],
            "albertii_hit": bool(candidate_hits),
            "number_of_hits": len(candidate_hits),
            "best_identity": best["identity"] if best else "",
            "best_evalue": best["evalue"] if best else "",
            "best_subject": best["subject"] if best else "",
        })

print("Candidates:", len(candidates))
print("With Albertii hits:", sum(bool(hits[x["protein_id"]]) for x in candidates))
print("Without Albertii hits:", sum(not hits[x["protein_id"]] for x in candidates))
print("Saved:", output_file)
