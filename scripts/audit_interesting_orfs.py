
import pandas as pd
from Bio import SeqIO

candidates = pd.read_csv(
    "results/synteny/candidate_coordinates.tsv", sep="\t"
)

albertii = pd.read_csv(
    "results/synteny/albertii_cds_annotations.tsv", sep="\t"
)

intervals = pd.read_csv(
    "results/synteny/albertii_intervals.fasta.tsv", sep="\t"
)

orfs = pd.read_csv(
    "results/orf/albertii_orfs.tsv", sep="\t"
)

fasta_records = list(
    SeqIO.parse("results/orf/albertii_orfs.faa", "fasta")
)

if len(orfs) != len(fasta_records):
    raise ValueError(
        "Number of ORFs in TSV and FASTA does not match"
    )

orfs["orf_id"] = [record.id for record in fasta_records]

seqids = dict(zip(intervals["candidate_id"], intervals["seqid"]))
orfs["seqid"] = orfs["candidate_id"].map(seqids)

blast = pd.read_csv(
    "results/orf/albertii_orfs_vs_ecoli_candidates.tsv",
    sep="\t",
    names=[
        "orf_id", "hit_candidate", "identity", "alignment_length",
        "orf_length", "candidate_length", "qstart", "qend",
        "evalue", "bitscore"
    ]
)

interesting = [
    "ecoli_REL606__ACT39985.1",
    "ecoli_Sakai__NP_303550.1",
    "ecoli_Sakai__NP_308259.1",
    "ecoli_Sakai__NP_310689.1",
    "ecoli_Sakai__NP_310248.1",
    "ecoli_Sakai__NP_309366.1",
    "ecoli_Sakai__NP_311896.1",
    "ecoli_Sakai__NP_309554.1"
]

print("\n=== E. COLI CANDIDATE ANNOTATIONS ===")

print(
    candidates[candidates["protein_id"].isin(interesting)][
        ["protein_id", "length", "seqid", "start", "end",
         "strand", "gene", "locus_tag", "product"]
    ].to_string(index=False)
)

blast["origin_candidate"] = blast["orf_id"].str.split("|").str[1]

selected = blast[
    blast["origin_candidate"].isin(interesting)
].copy()

best = (
    selected.sort_values("bitscore", ascending=False)
    .drop_duplicates("origin_candidate")
)

print("\n=== BEST BLAST HIT PER INTERESTING LOCUS ===")

print(
    best[
        ["origin_candidate", "orf_id", "hit_candidate",
         "identity", "alignment_length", "orf_length",
         "candidate_length", "evalue", "bitscore"]
    ].to_string(index=False)
)

print("\n=== ALBERTII CDS OVERLAPPING BEST ORFS ===")

audit_rows = []

for _, row in best.iterrows():
    matches = orfs[orfs["orf_id"] == row["orf_id"]]

    if matches.empty:
        print("\nORF not found:", row["orf_id"])
        continue

    o = matches.iloc[0]

    overlaps = albertii[
        (albertii["seqid"] == o["seqid"])
        & (albertii["start"] <= o["genomic_end"])
        & (albertii["end"] >= o["genomic_start"])
    ].copy()

    orf_length_bp = (
        int(o["genomic_end"]) - int(o["genomic_start"]) + 1
    )

    print("\nCandidate:", row["origin_candidate"])
    print("ORF:", row["orf_id"])
    print(
        "ORF coordinates:",
        o["seqid"], o["genomic_start"], o["genomic_end"],
        "strand:", o["strand"],
        "length:", orf_length_bp, "bp"
    )

    if overlaps.empty:
        print("No annotated CDS overlaps this ORF")

        audit_rows.append({
            "origin_candidate": row["origin_candidate"],
            "orf_id": row["orf_id"],
            "orf_seqid": o["seqid"],
            "orf_start": o["genomic_start"],
            "orf_end": o["genomic_end"],
            "orf_strand": o["strand"],
            "orf_length_bp": orf_length_bp,
            "albertii_protein_id": None,
            "albertii_gene": None,
            "albertii_product": None,
            "cds_start": None,
            "cds_end": None,
            "cds_strand": None,
            "overlap_bp": 0,
            "same_strand": None
        })
        continue

    overlaps["overlap_bp"] = (
        overlaps["end"].clip(upper=int(o["genomic_end"]))
        - overlaps["start"].clip(lower=int(o["genomic_start"]))
        + 1
    )

    overlaps = overlaps[overlaps["overlap_bp"] > 0].copy()

    overlaps["same_strand"] = (
        overlaps["strand"] == o["strand"]
    )

    if (overlaps["overlap_bp"] > orf_length_bp).any():
        raise ValueError(
            f"Overlap exceeds ORF length for {row['orf_id']}"
        )

    overlaps = overlaps.sort_values(
        "overlap_bp", ascending=False
    )

    print(
        overlaps[
            ["protein_id", "gene", "product", "start",
             "end", "strand", "overlap_bp", "same_strand"]
        ].to_string(index=False)
    )

    for _, a in overlaps.iterrows():
        audit_rows.append({
            "origin_candidate": row["origin_candidate"],
            "orf_id": row["orf_id"],
            "orf_seqid": o["seqid"],
            "orf_start": o["genomic_start"],
            "orf_end": o["genomic_end"],
            "orf_strand": o["strand"],
            "orf_length_bp": orf_length_bp,
            "albertii_protein_id": a["protein_id"],
            "albertii_gene": a["gene"],
            "albertii_product": a["product"],
            "cds_start": a["start"],
            "cds_end": a["end"],
            "cds_strand": a["strand"],
            "overlap_bp": a["overlap_bp"],
            "same_strand": a["same_strand"]
        })

pd.DataFrame(audit_rows).to_csv(
    "results/orf/interesting_orf_annotation_audit.tsv",
    sep="\t",
    index=False
)

print("\nAudit saved to:")
print("results/orf/interesting_orf_annotation_audit.tsv")

