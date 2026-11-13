import os
import csv
from Bio import AlignIO

R = "results/orthofinder_run2/Results_Nov12"
LIST = f"{R}/Orthogroups/Orthogroups_SingleCopyOrthologs.txt"
FOLDER = f"{R}/MultipleSequenceAlignments"
OUT = "results/phylogeny"

species = ["MG1655", "W3110", "REL606", "Sakai", "E_albertii"]

with open(LIST) as f:
    ogs = [x.strip() for x in f if x.strip()]

parts = {s: [] for s in species}
manifest = []

for og in ogs:
    aln = AlignIO.read(f"{FOLDER}/{og}.fa", "fasta")
    records = {}

    for rec in aln:
        matches = [s for s in species if rec.id.startswith(s + "_")]
        if len(matches) != 1:
            raise ValueError(f"Cannot identify species: {rec.id} in {og}")
        records[matches[0]] = str(rec.seq)

    if set(records) != set(species):
        raise ValueError(f"Missing species in {og}: {set(species)-set(records)}")

    length = aln.get_alignment_length()
    if len(set(map(len, records.values()))) != 1:
        raise ValueError(f"Unequal sequence lengths in {og}")

    for s in species:
        parts[s].append(records[s])

    manifest.append((og, length))

path = f"{OUT}/single_copy_supermatrix.faa"
with open(path, "w") as f:
    for s in species:
        f.write(f">{s}\n{''.join(parts[s])}\n")

with open(f"{OUT}/gene_manifest.tsv", "w") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["orthogroup", "alignment_length"])
    w.writerows(manifest)

print("Orthogroups:", len(manifest))
print("Total alignment length:", sum(x[1] for x in manifest))
print("Species:", ", ".join(species))
print("Saved:", path)
