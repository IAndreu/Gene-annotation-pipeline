import os
import os.path
from os import path
import sys
import pandas as pd

Genome = sys.argv[1]
xlsx = sys.argv[2]
genome_dir = sys.argv[3]
# Get family names from 'gene_families.xlsx'
df = pd.read_excel(xlsx) 
# Store Gene families:
gene_families= [df['Gene family'][i].replace(" ", "") for i in range(len(df))]
average= [str(df['Average number of genes'][i]) for i in range(len(df))]

with open(genome_dir+"/table_results.txt", "w") as fp:
    fp.write("Gene/Gene family\tStep_1_complete\tStep_1_partial\tStep_1_chimera\tStep_1_total\tNumber of annotated genes identified\tNumber of putative not annotated genes\tTotal number of identified genes (Annotated + Genomic)\tTotal number of identified genes clustering identical sequences\tTotal number of identified genes clustering highly identical sequences and proteins shorter than 30 aa\tAverage number in ant genomes\n")
    fp.close()
for i in range(len(gene_families)):
    merged = genome_dir+"gene_families_pipeline/"+gene_families[i]+"/"+gene_families[i]+"_merged.txt"
    step1 = len(open(merged).readlines())
    complete = open(merged).read().count("complete")
    chimeras = open(merged).read().count("chimera")
    parcial = open(merged).read().count("parcial")
    # check if the output file from that Gene family exists
    result= genome_dir+"gene_families_pipeline/"+gene_families[i]+"/Result/"+Genome+"_genecounts_summary.txt"
    if str(path.isfile(result))=='True':
        f=open(result)
        lines=f.readlines()
        gf, n1, n2, n3, n4, n5 = lines[1].split()
        with open(genome_dir+"/table_results.txt", "a") as fp:
            fp.write(gf+"\t"+str(complete)+"\t"+str(parcial)+"\t"+str(chimeras)+"\t"+str(step1)+"\t"+n1+"\t"+n2+"\t"+n3+"\t"+n4+"\t"+n5+"\t"+average[i]+"\n")
            fp.close()

