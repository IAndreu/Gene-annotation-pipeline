# How to run: $ python3 run_blast.py query database evalue
import os
import sys
import pandas as pd


############ INPUTS ############
db = sys.argv[1]
query = sys.argv[2]
evalue = sys.argv[3]
# threads = sys.argv[5]



# Generate blast database from fasta file
os.system("makeblastdb -dbtype prot -in %s" %db)

#output name = gene family name + _blast_output (for example OBP_blast_otuput)
output_name=  query.replace("_db.fasta","")+"_blast_otuput"
# use OBP as test:
#OBP_gf = gene_families_db[3]

os.system("blastp -query %s -db %s -out %s -evalue %s -outfmt '6 std qlen slen'"% (query, db, output_name, evalue))
