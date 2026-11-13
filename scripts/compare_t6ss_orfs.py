from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
import csv

GENOMES = {
    "Sakai": {
        "file": "data/raw/ecoli_Sakai/unpacked/ncbi_dataset/data/GCF_000008865.2",
        "seqid": "NC_002695.2",
        "start": 262678,
        "end": 262902,
        "strand": "-",
    },
    "albertii": {
        "file": "data/raw/escherichia_albertii/unpacked/ncbi_dataset/data/GCF_028622335.1",
        "seqid": "NZ_CP117568.1",
        "start": 4007735,
        "end": 4007911,
        "strand": "-",
    },
}

OUTPUT = Path("results/synteny/t6ss_orf_comparison.tsv")


def get_sequence(sample, info):
    files = list(Path(info["file"]).glob("*.fna"))

    if not files:
        raise FileNotFoundError(f"No FASTA for {sample}")

    genome = SeqIO.to_dict(SeqIO.parse(files[0], "fasta"))

    if info["seqid"] not in genome:
        raise ValueError(f"Missing sequence {info['seqid']}")

    sequence = genome[info["seqid"]].seq[
        info["start"] - 1:info["end"]
    ]

    if info["strand"] == "-":
        sequence = sequence.reverse_complement()

    return str(sequence).upper()


def main():
    sequences = {}

    for sample, info in GENOMES.items():
        sequence = get_sequence(sample, info)
        sequences[sample] = sequence

        protein = str(Seq(sequence).translate())

        print(f"\n{sample}")
        print(f"DNA length: {len(sequence)}")
        print(f"Protein including stop: {protein}")
        print(f"Internal stops: {protein[:-1].count('*')}")
        print(f"Terminal stop: {protein.endswith('*')}")

    sakai = sequences["Sakai"]
    albertii = sequences["albertii"]

    sakai_protein = str(Seq(sakai).translate()).rstrip("*")
    albertii_protein = str(Seq(albertii).translate()).rstrip("*")

    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -1

    alignment = aligner.align(
        sakai_protein,
        albertii_protein
    )[0]

    aligned_length = sum(
        min(end1 - start1, end2 - start2)
        for (start1, end1), (start2, end2)
        in zip(*alignment.aligned)
    )

    identical = 0

    for (start1, end1), (start2, end2) in zip(
        alignment.aligned[0],
        alignment.aligned[1]
    ):
        for i in range(min(end1 - start1, end2 - start2)):
            if sakai_protein[start1 + i] == albertii_protein[start2 + i]:
                identical += 1

    identity = (
        identical / aligned_length * 100
        if aligned_length else 0
    )

    print("\n=== PROTEIN ALIGNMENT ===")
    print(alignment)
    print(f"Identical positions: {identical}")
    print(f"Aligned positions: {aligned_length}")
    print(f"Identity: {identity:.2f}%")

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)

    with OUTPUT.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "sakai_length_nt",
            "albertii_length_nt",
            "sakai_length_aa",
            "albertii_length_aa",
            "identical_positions",
            "aligned_positions",
            "identity_percent",
        ])
        writer.writerow([
            len(sakai),
            len(albertii),
            len(sakai_protein),
            len(albertii_protein),
            identical,
            aligned_length,
            round(identity, 2),
        ])

    print(f"\nSaved: {OUTPUT}")


if __name__ == "__main__":
    main()
