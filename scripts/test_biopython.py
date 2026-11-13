from Bio import SeqIO

fasta_path = "tests/test_data/example.fasta"

records = list(SeqIO.parse(fasta_path, "fasta"))

print("Number of sequences:", len(records))

for record in records:
    print("ID:", record.id)
    print("Sequence:", record.seq)
    print("Length:", len(record.seq))
