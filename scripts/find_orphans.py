import csv
from pathlib import Path
from collections import defaultdict

blast_file = Path("results/homology/all_vs_all.tsv")
mapping_file = Path("data/processed/protein_mapping.tsv")

output_dir = Path("results/homology")
output_dir.mkdir(parents=True, exist_ok=True)

# Загружаем информацию о белках
protein_info = {}

with open(mapping_file) as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        protein_info[row["protein_id"]] = {
            "sample": row["sample"],
            "length": int(row["length"]),
            "original_id": row["original_id"],
        }

# Список всех штаммов
all_samples = {
    info["sample"]
    for info in protein_info.values()
}

# Для каждого белка сохраняем штаммы,
# в которых найдены убедительные совпадения
cross_hits = defaultdict(set)

with open(blast_file) as f:
    for line in f:
        fields = line.strip().split("\t")

        if len(fields) != 8:
            continue

        (
            query,
            subject,
            identity,
            alignment_length,
            qlen,
            slen,
            evalue,
            bitscore,
        ) = fields

        if query not in protein_info or subject not in protein_info:
            continue

        query_sample = protein_info[query]["sample"]
        subject_sample = protein_info[subject]["sample"]

        # Не учитываем совпадения внутри одного штамма
        if query_sample == subject_sample:
            continue

        qlen = int(qlen)
        slen = int(slen)
        alignment_length = int(alignment_length)
        evalue = float(evalue)

        if qlen == 0 or slen == 0:
            continue

        query_coverage = alignment_length / qlen
        subject_coverage = alignment_length / slen

        if (
            evalue <= 1e-5
            and query_coverage >= 0.50
            and subject_coverage >= 0.50
        ):
            cross_hits[query].add(subject_sample)

# Отбираем предварительных кандидатов
candidates = []

for protein_id, info in protein_info.items():

    other_samples = all_samples - {info["sample"]}
    samples_with_hits = cross_hits[protein_id]

    missing_samples = other_samples - samples_with_hits

    candidates.append({
        "protein_id": protein_id,
        "original_id": info["original_id"],
        "sample": info["sample"],
        "length": info["length"],
        "number_of_other_samples_with_hits": len(samples_with_hits),
        "matched_samples": ",".join(sorted(samples_with_hits)),
        "number_of_other_samples_without_hits": len(missing_samples),
        "preliminary_orphan": len(samples_with_hits) == 0,
        "short_protein_warning": info["length"] < 50,
    })

output_file = output_dir / "orphan_candidates.tsv"

with open(output_file, "w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=candidates[0].keys(),
        delimiter="\t",
    )

    writer.writeheader()
    writer.writerows(candidates)

orphans = [
    protein for protein in candidates
    if protein["preliminary_orphan"]
]

short_orphans = [
    protein for protein in orphans
    if protein["short_protein_warning"]
]

print("Total proteins:", len(candidates))
print("Preliminary orphan candidates:", len(orphans))
print("Candidates shorter than 50 aa:", len(short_orphans))
print(
    "Candidates at least 50 aa:",
    len(orphans) - len(short_orphans),
)
print("Output:", output_file)
