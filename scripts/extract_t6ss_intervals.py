from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

GENOMES = {
    "Sakai": {
        "seqid": "NC_002695.2",
        "start": 262628,
        "end": 262936,
        "strand": "-",
        "candidate": "NP_308259.1",
    },
    "albertii": {
        "seqid": "NZ_CP117568.1",
        "start": 4007643,
        "end": 4007936,
        "strand": "-",
        "candidate": "ORF467",
    },
}

OUTPUT = Path("results/synteny/t6ss_intervals.fasta")


def find_genome(sample):
    matches = list(Path("data/raw").glob(f"{sample}/**/*.fna"))

    if not matches:
        raise FileNotFoundError(f"No genome FASTA found for {sample}")

    return matches[0]


def main():
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)

    with OUTPUT.open("w") as out:
        for sample, info in GENOMES.items():
            genome_path = find_genome(
                "ecoli_Sakai" if sample == "Sakai"
                else "escherichia_albertii"
            )

            records = SeqIO.to_dict(SeqIO.parse(genome_path, "fasta"))

            if info["seqid"] not in records:
                print(f"Available IDs in {sample}:")
                print(list(records.keys())[:10])
                raise ValueError(f"Sequence {info['seqid']} not found")

            genome = records[info["seqid"]].seq

            sequence = genome[info["start"] - 1:info["end"]]

            if info["strand"] == "-":
                sequence = sequence.reverse_complement()

            print(
                f"{sample}: {info['seqid']} "
                f"{info['start']}-{info['end']} "
                f"length={len(sequence)}"
            )

            out.write(
                f">{sample}|{info['candidate']}|"
                f"{info['seqid']}:{info['start']}-{info['end']}\n"
            )

            sequence_text = str(sequence)
            for i in range(0, len(sequence_text), 70):
                out.write(sequence_text[i:i + 70] + "\n")

    print(f"\nFASTA saved: {OUTPUT}")


if __name__ == "__main__":
    main()
