
import pandas as pd
from Bio.Seq import Seq
from pathlib import Path

out = Path("results/orf")
out.mkdir(parents=True, exist_ok=True)

intervals = pd.read_csv(
    "results/synteny/albertii_intervals.fasta.tsv",
    sep="\t"
)

candidates = pd.read_csv(
    "results/synteny/candidate_coordinates.tsv",
    sep="\t"
)

synteny = pd.read_csv(
    "results/synteny/validated_flanking_synteny.tsv",
    sep="\t"
)

lengths = dict(zip(candidates["protein_id"], candidates["length"]))
statuses = dict(zip(synteny["candidate_id"], synteny["status"]))

starts = {"ATG", "GTG", "TTG"}
stops = {"TAA", "TAG", "TGA"}

MIN_AA = 50

results = []

for _, row in intervals.iterrows():

    candidate = row["candidate_id"]
    dna_seq = str(row["sequence"]).upper()

    genomic_start = int(row["start"])
    genomic_end = int(row["end"])
    interval_strand = row["strand"]

    interval_length = len(dna_seq)
    candidate_length = lengths.get(candidate)

    if interval_strand not in {"+", "-"}:
        raise ValueError(
            f"Unknown interval strand for {candidate}: {interval_strand}"
        )

    sequences = [
        ("+", dna_seq),
        ("-", str(Seq(dna_seq).reverse_complement()))
    ]

    for sequence_strand, dna in sequences:

        for frame in range(3):

            coding = dna[frame:]

            for i in range(0, len(coding) - 2, 3):

                if coding[i:i + 3] not in starts:
                    continue

                stop_position = None

                for j in range(i + 3, len(coding) - 2, 3):

                    if coding[j:j + 3] in stops:
                        stop_position = j
                        break

                if stop_position is None:
                    continue

                orf_dna = coding[i:stop_position + 3]
                aa = str(Seq(orf_dna).translate())

                if len(aa) < MIN_AA:
                    continue

                rel_start = frame + i + 1
                rel_end = frame + stop_position + 3

                if interval_strand == "+":

                    g_start = genomic_start + rel_start - 1
                    g_end = genomic_start + rel_end - 1

                    genomic_strand = sequence_strand

                else:

                    g_start = genomic_end - rel_end + 1
                    g_end = genomic_end - rel_start + 1

                    genomic_strand = (
                        "-" if sequence_strand == "+" else "+"
                    )

                results.append({
                    "candidate_id": candidate,
                    "synteny_status": statuses.get(candidate),
                    "candidate_length_aa": candidate_length,
                    "orf_length_aa": len(aa),
                    "length_ratio": (
                        round(len(aa) / candidate_length, 3)
                        if pd.notna(candidate_length) and candidate_length
                        else None
                    ),
                    "interval_length": interval_length,
                    "interval_strand": interval_strand,
                    "sequence_strand": sequence_strand,
                    "strand": genomic_strand,
                    "frame": frame + 1,
                    "start_codon": orf_dna[:3],
                    "stop_codon": orf_dna[-3:],
                    "genomic_start": g_start,
                    "genomic_end": g_end,
                    "protein_sequence": aa
                })

df = pd.DataFrame(results)

if not df.empty:

    df = df.sort_values(
        ["candidate_id", "orf_length_aa"],
        ascending=[True, False]
    ).reset_index(drop=True)

    df["orf_id"] = [
        f"ORF{i + 1}|{row.candidate_id}|{row.strand}|"
        f"{row.genomic_start}-{row.genomic_end}"
        for i, row in df.iterrows()
    ]

df.to_csv(
    out / "albertii_orfs.tsv",
    sep="\t",
    index=False
)

with open(out / "albertii_orfs.faa", "w") as f:

    for _, row in df.iterrows():

        f.write(f">{row.orf_id}\n")
        f.write(row.protein_sequence + "\n")

print("Intervals scanned:", len(intervals))
print("ORFs found:", len(df))

print(
    "Candidates with ORFs:",
    df["candidate_id"].nunique() if not df.empty else 0
)

print(
    "Candidates without ORFs:",
    len(intervals) - (
        df["candidate_id"].nunique() if not df.empty else 0
    )
)

if not df.empty:

    print("\nORFs by synteny:")
    print(df.groupby("synteny_status").size().to_string())

    print("\nORFs by interval orientation:")
    print(df.groupby("interval_strand").size().to_string())

    print("\nORFs by genomic strand:")
    print(df.groupby("strand").size().to_string())

    print("\nLongest ORF per candidate:")

    best = df.groupby("candidate_id").first()

    print(
        best[
            [
                "synteny_status",
                "candidate_length_aa",
                "orf_length_aa",
                "length_ratio",
                "start_codon",
                "interval_strand",
                "strand"
            ]
        ].to_string()
    )

print("\nSaved:")
print(out / "albertii_orfs.tsv")
print(out / "albertii_orfs.faa")

