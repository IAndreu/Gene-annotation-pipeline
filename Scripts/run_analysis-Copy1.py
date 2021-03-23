#############################################################
#                                                           #
#     Run first analysis, using blastp and InterproScan     #
#                                                           #
#############################################################
# How to run: $ python3 run_analysis.py proteome

import os
import sys
import pandas as pd
import numpy as np

genome = sys.argv[1]
genome_directory = sys.argv[2]
# Get family names from 'GAGA_gene_families.xlsx'
df = pd.read_excel('Data/GAGA_gene_families.xlsx') # Location in ERDA: 'GAGA/Gene_annotation/Gene_family_pipeline/GAGA_gene_families.xlsx'

# Store Gene families:
gene_families= [df['Gene family'][i] for i in range(len(df))]

# Get the gene families fasta database file:
gene_families_db = ['Data/Gene_families/'+i+'/'+i+'_db.fasta' for i in gene_families]
# Copy the gene family db in a folder inside the genome
'''for i in range(len(gene_families_db)):
    new = genome_directory+gene_families[i]
    os.system("mkdir -p %s && cp %s %s" % (new, gene_families_db[i], new))
    gene_families_db[i] = new + '/' + gene_families[i]+'_db.fasta'
    '''
os.system("mkdir -p %s && cp %s %s" % ('Data/Genomes/testset/OBP', 'Data/Gene_families/OBP/OBP_db.fasta', 'Data/Genomes/testset/OBP'))

# Get if it is possible to run blast for each gene family
blast= [df['Blast'][i] for i in range(len(df))]

# Get InterPro and Pfam domains
InterPro = [df['InterPro Domain'][i] for i in range(len(df))]
Pfam = [df['Pfam domain'][i] for i in range(len(df))]

# Get Average length of the proteins of each gene family and peptide average length
prot_avg_len = [float(df['Average protein length'][i]) for i in range(len(df))]
pept_avg_len = [float(df['Peptide start length'][i]) for i in range(len(df))]

# Get e-values for each gene family to run blast
evalues = [float(df['E-value'][i]) for i in range(len(df))]


# Run analysis of the gene family database and divide the database into rare and good sequences
'''
for i in range(len(gene_families)):
    os.system("python3 Scripts/db_analysis.py %s" % (gene_families_db[i]))
'''
os.system("python3 Scripts/db_analysis.py %s" % (genome_directory+'OBP/OBP_db.fasta'))

# run blast for every gene family if able
'''
for i in range(len(gene_families)):
    if blast[i]=='Yes':
        os.system("python3 Scripts/run_blast.py %s %s %s" % (genome, gene_families_db[i], evalues[i]))
        name = gene_families_db[i].replace("_db.fasta","_blast_otuput.txt")
        outname = gene_families_db[i].replace("_db.fasta","")
        os.system("perl Scripts/get_blastp_parsed_newv2.pl %s %s %s" % (name, outname, evalues[i]))
    else: # create a blank file called GeneFamilyName_parsed_blast.txt
        with open(gene_families_db[i].replace("_db.fasta","")+'_parsed_blast.txt', 'w') as fp: 
            pass
        
    '''
# EXAMPLE ONLY WITH OBP GENE FAMILY
# run blast
os.system("python3 Scripts/run_blast.py %s %s %s" % (genome, genome_directory+'OBP/OBP_db.fasta', 0.001))
# Parse blast
os.system("perl Scripts/get_blastp_parsed_newv2.pl %s %s %s" % (genome_directory+'OBP/OBP_blast_otuput.txt', genome_directory+'OBP/OBP', 0.001))

# read InterproScan output
tsv = pd.read_csv(genome+'.tsv', sep='\t', header=None, names = list(range(0,14)))

# For each gene family, produce a table with the proteins find in interproscan output if able:
'''for i in range(len(InterPro)):
    if Pfam[i][0]=='P': #or InterPro[i][0]=='I' if there is domain id from InterPro or Pfam, then:
        name= gene_families_db[i].replace("_db.fasta","")+'_parsed_domain.txt'
        table=tsv.loc[tsv[4] == Pfam[i]]
        new_table = table[[0, 6, 7, 2, 3]]
        new_table.to_csv(path_or_buf= name ,header=None, index=None, sep='\t', mode='a')
    else: # otherwise create a blank file
        with open((name, 'w') as fp: 
            pass
        '''

# EXAMPLE ONLY WITH OBP GENE FAMILY
# Get the proteins that have the interproscan identifier:
table=tsv.loc[tsv[4] == 'PF01395'] # tsv[11] == 'IPR006170'
# Get the columns of interest and remove repeated: Protein Name, Start location, Stop location, Sequence length, Signature accession 
new_table = table[[0, 6, 7, 2, 3]]
# Save the table as GeneFamily_parsed_domain.txt
new_table.to_csv(path_or_buf=genome_directory+'OBP/OBP_parsed_domain.txt',header=None, index=None, sep='\t', mode='a')


# Merge blast and interproscan results into a table and retrieve the fasta sequences
'''for i in range(len(gene_families)):
    blast= gene_families_db[i].replace("_db.fasta","")+'_parsed_blast.txt'
    inter=  gene_families_db[i].replace("_db.fasta","")+'_parsed_domain.txt'
    prot_avg= prot_avg_len[i]
    pept_avg= pept_avg_len[i]
    os.system("python3 Scripts/merge_blast_domain.py %s %s %s %s %s" % (genome, blast, inter, prot_avg, pept_avg))
    '''
# EXAMPLE ONLY WITH OBP GENE FAMILY
# run blast
os.system("python3 Scripts/merge_blast_domain.py %s %s %s %s %s" % (genome,genome_directory+'OBP/OBP_parsed_blast.txt', genome_directory+'OBP/OBP_parsed_domain.txt', 146, 40))