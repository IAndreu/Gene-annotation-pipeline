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

os.system("mkdir GFF")

for i in range(len(gene_families)):
    # check if the output file from that Gene family exists
    gff="Data/Genomes/"+Genome+"/gene_families_pipeline/"+gene_families[i]+"/"+gene_families[i]+"_annot_genes.gff3"
    gff_t="Data/Genomes/"+Genome+"/gene_families_pipeline/"+gene_families[i]+"/"+gene_families[i]+"_annot_genes_trimmed.gff3"
    os.system("cp %s GFF" % (gff))
    os.system("cp %s GFF" % (gff_t))
