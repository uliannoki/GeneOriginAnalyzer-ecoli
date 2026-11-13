from pathlib import Path
import pandas as pd

BLAST_FILE = Path("results/orf/albertii_orfs_vs_ecoli_candidates.tsv")
ANNOTATION_FILE = Path("results/reports/all_ecoli_cds_annotations.tsv")
OUTPUT = Path("results/reports/candidate_homolog_trace.tsv")

TARGETS = {
    "ecoli_Sakai__NP_310248.1",
    "ecoli_Sakai__NP_303550.1",
}

columns = [
    "qseqid", "sseqid", "pident", "length",
    "qlen", "slen", "qstart", "qend",
    "sstart", "send", "evalue", "bitscore"
]

blast = pd.read_csv(BLAST_FILE, sep="\t", names=columns)
annotations = pd.read_csv(ANNOTATION_FILE, sep="\t")

annotation_map = (
    annotations.drop_duplicates("protein_id")
    .set_index("protein_id")
    .to_dict("index")
)

rows = []

for _, hit in blast.iterrows():
    parts = str(hit["qseqid"]).split("|")

    if len(parts) < 2:
        continue

    parent_id = parts[1]

    if parent_id not in TARGETS:
        continue

    query_info = annotation_map.get(parent_id, {})
    hit_info = annotation_map.get(hit["sseqid"], {})

    rows.append({
        "query_candidate": parent_id,
        "query_gene": query_info.get("gene"),
        "query_product": query_info.get("product"),
        "query_seqid": query_info.get("seqid"),
        "query_start": query_info.get("start"),
        "query_end": query_info.get("end"),
        "orf_id": parts[0],
        "hit_candidate": hit["sseqid"],
        "hit_sample": hit_info.get("sample"),
        "hit_gene": hit_info.get("gene"),
        "hit_product": hit_info.get("product"),
        "hit_seqid": hit_info.get("seqid"),
        "hit_start": hit_info.get("start"),
        "hit_end": hit_info.get("end"),
        "hit_strand": hit_info.get("strand"),
        "identity": hit["pident"],
        "alignment_length": hit["length"],
        "query_coverage": round(
            (hit["qend"] - hit["qstart"] + 1) / hit["qlen"] * 100, 2
        ),
        "target_coverage": round(
            (hit["send"] - hit["sstart"] + 1) / hit["slen"] * 100, 2
        ),
        "evalue": hit["evalue"],
        "bitscore": hit["bitscore"],
    })

result = pd.DataFrame(rows)

if result.empty:
    print("No matching records.")
else:
    result = result.sort_values(
        ["query_candidate", "bitscore"],
        ascending=[True, False]
    )

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(OUTPUT, sep="\t", index=False)

    for candidate in sorted(TARGETS):
        subset = result[result["query_candidate"] == candidate]

        print(f"\n=== {candidate} ===")

        if subset.empty:
            print("No hits.")
            continue

        print(
            subset.drop_duplicates("hit_candidate")[
                [
                    "orf_id", "hit_candidate", "hit_sample",
                    "hit_gene", "hit_product", "hit_seqid",
                    "hit_start", "hit_end", "hit_strand",
                    "identity", "alignment_length",
                    "query_coverage", "target_coverage",
                    "evalue", "bitscore"
                ]
            ].head(20).to_string(index=False)
        )

    print("\nUnmapped hit annotations:",
          result["hit_gene"].isna().sum())
    print("Saved:", OUTPUT)
