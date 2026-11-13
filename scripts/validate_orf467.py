
from Bio import SeqIO
from Bio.Seq import Seq
from pathlib import Path

orf_file = "results/orf/albertii_orfs.faa"
candidate_file = "data/processed/orphan_candidates.faa"

target_orf = "ORF467|ecoli_Sakai__NP_308259.1|-|4007735-4007911"
target_candidate = "ecoli_Sakai__NP_308259.1"

orfs = {
    record.id: record
    for record in SeqIO.parse(orf_file, "fasta")
}

candidates = {
    record.id: record
    for record in SeqIO.parse(candidate_file, "fasta")
}

if target_orf not in orfs:
    raise ValueError(f"ORF not found: {target_orf}")

if target_candidate not in candidates:
    raise ValueError(f"Candidate not found: {target_candidate}")

orf = orfs[target_orf]
candidate = candidates[target_candidate]

print("ORF467")
print("Length:", len(orf.seq), "aa")
print("Sequence:", orf.seq)

print("\nOriginal E. coli candidate")
print("Length:", len(candidate.seq), "aa")
print("Sequence:", candidate.seq)

out = Path("results/orf")
out.mkdir(parents=True, exist_ok=True)

with open(out / "orf467_vs_candidate.faa", "w") as f:
    SeqIO.write([orf, candidate], f, "fasta")

print("\nSaved:")
print("results/orf/orf467_vs_candidate.faa")

