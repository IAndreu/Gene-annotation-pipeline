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
average= [str(df['Average number of genes'][i]) for i in range(len(df))]

with open("Data/Genomes/"+Genome+"/table_results.txt", "w") as fp:
    fp.write("Gene family\tstep 1\tBitacora\tAverage number in ant genomes\n")
    fp.close()
for i in range(len(gene_families)):
    step1 = len(open("Data/Genomes/"+Genome+"/gene_families_pipeline/"+gene_families[i]+"/"+gene_families[i]+"_merged.txt").readlines())
    # check if the output file from that Gene family exists
    result="Data/Genomes/"+Genome+"/gene_families_pipeline/"+gene_families[i]+"/Result/"+Genome+"_genecounts_summary.txt"
    if str(path.isfile(result))=='True':
        f=open(result)
        lines=f.readlines()
        gf, n1, n2, tot = lines[1].split()[:4]
        with open("Data/Genomes/"+Genome+"/table_results.txt", "a") as fp:
            fp.write(gf+"\t"+str(step1)+"\t"+tot+"\t"+average[i]+"\n")
            fp.close()
