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
df = pd.read_excel('Data/GAGA_gene_families.xlsx') 

# Store Gene families:
gene_families= [df['Gene family'][i] for i in range(len(df))]

# Get the gene families fasta database file:
gene_families_db = ['Data/Gene_families/'+i+'/'+i+'_db.fasta' for i in gene_families]
# Copy the gene family db in a folder inside the genome
for i in range(len(gene_families_db)):
    new = genome_directory+gene_families[i]
    os.system("mkdir -p %s && cp %s %s" % (new, gene_families_db[i], new))
    gene_families_db[i] = new + '/' + gene_families[i]+'_db.fasta'


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
for i in range(len(gene_families)):
    os.system("python3 Scripts/db_analysis.py %s" % (gene_families_db[i]))

# run blast for every gene family if able
for i in range(len(gene_families)):
    if blast[i]=='Yes':
        os.system("python3 Scripts/run_blast.py %s %s %s" % (genome, gene_families_db[i], evalues[i]))
        name = gene_families_db[i].replace("_db.fasta","_blast_otuput.txt")
        outname = gene_families_db[i].replace("_db.fasta","")
        os.system("perl Scripts/get_blastp_parsed_newv2.pl %s %s %s" % (name, outname, evalues[i]))
    else:# create a blank file called GeneFamilyName_parsed_blast.txt
        with open(gene_families_db[i].replace("_db.fasta","")+'_parsed_blast.txt', 'w') as fp:
            fp.write('delete')
            fp.write("\t")
            fp.write("1")
            fp.write("\t")
            fp.write("10")
            fp.write("\t")
            fp.write("10")
            fp.write("\t")
            fp.write("blastp")
            fp.write("\n")
            fp.write('delete')
            fp.write("\t")
            fp.write("1")
            fp.write("\t")
            fp.write("10")
            fp.write("\t")
            fp.write("10")
            fp.write("\t")
            fp.write("blastp")
            fp.write("\n")
            fp.close()

# read InterproScan output
tsv = pd.read_csv(genome+'.tsv', sep='\t', header=None, names = list(range(0,14)))

# For each gene family, produce a table with the proteins find in interproscan output if able:
for i in range(len(InterPro)):
    name= gene_families_db[i].replace("_db.fasta","")+'_parsed_domain.txt'
    if str(Pfam[i])[0]=='P': #or InterPro[i][0]=='I' if there is domain id from InterPro or Pfam, then:
        table=tsv.loc[tsv[4] == Pfam[i]]
        new_table = table[[0, 6, 7, 2, 3]]
        new_table.to_csv(path_or_buf= name ,header=None, index=None, sep='\t', mode='a')
    else:
        # otherwise create a blank file
        with open(name, 'w') as fp: 
            fp.write('delete')
            fp.write("\t")
            fp.write("1")
            fp.write("\t")
            fp.write("10")
            fp.write("\t")
            fp.write("10")
            fp.write("\t")
            fp.write("Pfam")
            fp.write("\n")
            fp.write('delete')
            fp.write("\t")
            fp.write("1")
            fp.write("\t")
            fp.write("10")
            fp.write("\t")
            fp.write("10")
            fp.write("\t")
            fp.write("Pfam")
            fp.write("\n")
            fp.close()
                 

# Merge blast and interproscan results into a table and retrieve the fasta sequences
for i in range(len(gene_families)):
    blast= gene_families_db[i].replace("_db.fasta","")+'_parsed_blast.txt'
    inter=  gene_families_db[i].replace("_db.fasta","")+'_parsed_domain.txt'
    prot_avg= prot_avg_len[i]
    pept_avg= pept_avg_len[i]
    os.system("python3 Scripts/merge_blast_domain.py %s %s %s %s %s" % (genome, blast, inter, prot_avg, pept_avg))

# Produce HMM profile of the good+the parsed sequences
for i in range(len(gene_families)):
    fasta = gene_families_db[i].replace("_db.fasta","")+'good_db.fasta'
    alignment = gene_families_db[i].replace("_db.fasta","")+'_alignment.txt'
    profile = gene_families_db[i].replace("_db.fasta","")+'good_db.hmm'
    os.system("mafft --auto %s > %s" % (fasta, alignment))
    os.system("hmmbuild %s  %s" % (profile, alignment))    
