from Bio import SeqIO
from pathlib import Path

samples = [
    "ecoli_MG1655",
    "ecoli_W3110",
    "ecoli_REL606",
    "ecoli_Sakai",
]

base = Path("data/raw")

for sample in samples:
    sample_dir = base / sample / "unpacked"

    print("\n" + "=" * 55)
    print(f"SAMPLE: {sample}")
    print("=" * 55)

    if not sample_dir.exists():
        print("ERROR: directory not found")
        continue

    genome_files = list(sample_dir.rglob("*.fna"))
    protein_files = list(sample_dir.rglob("*.faa"))
    annotation_files = (
        list(sample_dir.rglob("*.gff"))
        + list(sample_dir.rglob("*.gff3"))
    )

    for label, files in [
        ("GENOME", genome_files),
        ("PROTEINS", protein_files),
    ]:
        print(f"\n{label}:")

        if not files:
            print("File not found")
            continue

        for path in files:
            count = 0
            first_record = None

            for record in SeqIO.parse(path, "fasta"):
                count += 1
                if first_record is None:
                    first_record = record

            print("File:", path.name)
            print("Number of sequences:", count)

            if first_record:
                print("First ID:", first_record.id)
                print("First sequence length:", len(first_record.seq))

    print("\nANNOTATION:")

    if not annotation_files:
        print("File not found")
    else:
        for path in annotation_files:
            with open(path) as f:
                count = sum(
                    1 for line in f
                    if line.strip() and not line.startswith("#")
                )

            print("File:", path.name)
            print("Annotation records:", count)
