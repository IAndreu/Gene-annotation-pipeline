import os
import sys
from Bio import SeqIO

db = sys.argv[1]

# Get mean of the proteins family
records = list(SeqIO.parse(db, "fasta"))
mean=0
for i in range(len(records)):
    mean+=len(records[i].seq)
mean=mean/(len(records))

# Get a file with rare sequences (being rare if they are smaller or higher then 60% of the average) and a file with good sequences
rare_sequences=[]
good_sequences=[]
for i in range(len(records)):
    if len(records[i].seq)<mean*0.6 or len(records[i].seq)>mean/0.6:
        rare_sequences.append(records[i])
    else:
        good_sequences.append(records[i])
rare_file = db.replace("_db.fasta","rare_db.fasta")   
good_file = db.replace("_db.fasta","good_db.fasta")       
SeqIO.write(rare_sequences, rare_file, "fasta")
SeqIO.write(good_sequences, good_file, "fasta")