from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

GENOMES = {
    "Sakai": {
        "file": "data/raw/ecoli_Sakai/unpacked/ncbi_dataset/data/GCF_000008865.2",
        "seqid": "NC_002695.2",
        "start": 262678,
        "end": 262902,
    },
    "albertii": {
        "file": "data/raw/escherichia_albertii/unpacked/ncbi_dataset/data/GCF_028622335.1",
        "seqid": "NZ_CP117568.1",
        "start": 4007735,
        "end": 4007911,
    },
}

START_CODONS = {"ATG", "GTG", "TTG"}
STOP_CODONS = {"TAA", "TAG", "TGA"}


for sample, info in GENOMES.items():
    fasta = list(Path(info["file"]).glob("*.fna"))[0]
    genome = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))

    sequence = genome[info["seqid"]].seq[
        info["start"] - 1:info["end"]
    ].reverse_complement()

    dna = str(sequence).upper()
    codons = [dna[i:i + 3] for i in range(0, len(dna), 3)]

    print(f"\n{sample}")
    print(f"Start codon: {codons[0]}")
    print(f"Start codon accepted: {codons[0] in START_CODONS}")
    print(f"Stop codon: {codons[-1]}")
    print(f"Terminal stop accepted: {codons[-1] in STOP_CODONS}")
    print(f"First 15 nt: {dna[:15]}")
    print(f"Last 15 nt: {dna[-15:]}")
