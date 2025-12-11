import sys
from Bio import SeqIO

if len(sys.argv) != 4:
    print("Usage: python script.py fasta_file1.fasta fasta_file2.fasta")
    sys.exit(1)

fasta_file1 = SeqIO.parse(sys.argv[1], "fasta")
fasta_file2 = SeqIO.parse(sys.argv[2], "fasta")
fasta_name = str(sys.argv[3])
matches = []

seq_record2_id= []
seq_record2_seq = []
for seq_record in fasta_file2:
    seq_record2_id.append(seq_record.id)
    seq_record2_seq.append(seq_record.seq)
j=0
for seq_record1 in fasta_file1:
    flag=True
    for i in range(0,len(seq_record2_seq)):
        if ''.join(["|genome_name:",fasta_name,"|"]) in seq_record2_id[i]:
            if flag:
                if str(seq_record2_seq[i]) in str(seq_record1.seq):
                    matches.append([seq_record1.id,"\t","unique"])
                    flag=False
    if flag:
        matches.append([seq_record1.id,"\t","known"])

    j+=1

print("contig\tfamiliarity")
for header in matches:
    print(''.join(header))
