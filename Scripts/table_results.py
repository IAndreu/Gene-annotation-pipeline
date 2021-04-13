import os
import os.path
from os import path
import sys
import pandas as pd

Genome = sys.argv[1]

# Get family names from 'gene_families.xlsx'
df = pd.read_excel('Data/gene_families.xlsx') 
# Store Gene families:
gene_families= [df['Gene family'][i].replace(" ", "") for i in range(len(df))]
average= [str(df['Average number in ant genomes'][i]) for i in range(len(df))]

with open("table_results.txt", "w") as fp:
    fp.write("Gene family\tstep 1\tBitacora\ttotal\tAverage number in ant genomes\n")
    fp.close()
for i in range(len(gene_families)):
    # check if the output file from that Gene family exists
    result="Data/Genomes/"+Genome+"/gene_families_pipeline/"+gene_families[i]+"/Result/"+Genome+"_genecounts_summary.txt"
    if str(path.isfile(result))=='True':
        f=open(result)
        lines=f.readlines()
        gf, step1, bit, tot = lines[1].split()[:4]
        with open("table_results.txt", "a") as fp:
            fp.write(gf+"\t"+step1+"\t"+bit+"\t"+tot+"\t"+average[i]+"\n")
            fp.close()

