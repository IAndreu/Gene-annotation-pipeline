#############################################################
#                                                           #
#     Run first analysis, using blastp and InterproScan     #
#                                                           #
#############################################################


import os
import sys
import pandas as pd
import numpy as np

path = sys.argv[1]
proteome = sys.argv[2]
genome_directory = sys.argv[3]
gff = sys.argv[4]
genome = sys.argv[5]
run_bitacora = sys.argv[6]
genome_name = sys.argv[7]
threads = sys.argv[8]

# Get family names from 'gene_families.xlsx'
df = pd.read_excel('Data/gene_families.xlsx') 

# Store Gene families:
gene_families= [df['Gene family'][i].replace(" ", "") for i in range(len(df))]

# Get if it is possible to run blast for each gene family
blast= [df['Blast'][i] for i in range(len(df))]

# Get the gene families fasta database file:
gene_families_db = ['Data/Gene_families/'+i+'_db.fasta' for i in gene_families]
# Copy the gene family db in a folder inside the proteome
for i in range(len(gene_families_db)):
    if blast[i]=='Yes':
        new = genome_directory+'/gene_families_pipeline/'+gene_families[i]
        os.system("mkdir -p %s && cp %s %s" % (new, gene_families_db[i], new))
        gene_families_db[i] = new + '/' + gene_families[i]+'_db.fasta'
    else:
        new = genome_directory+'/gene_families_pipeline/'+gene_families[i]
        os.system("mkdir -p %s && echo > %s" % (new, new+'/' + gene_families[i] + '_db.fasta'))
        gene_families_db[i] = new + '/' + gene_families[i]+'_db.fasta'

# Get InterPro and Pfam domains
InterPro = [df['InterPro Domain'][i] for i in range(len(df))]
Pfam = [df['Pfam domain'][i] for i in range(len(df))]

# Get Average length of the proteins of each gene family and peptide average length
prot_avg_len = [float(df['Average protein length'][i]) for i in range(len(df))]
pept_avg_len = [float(df['Peptide start length'][i]) for i in range(len(df))]

# Get e-values for each gene family to run blast
evalues = [float(df['E-value'][i]) for i in range(len(df))]


# run blast for every gene family if able
for i in range(len(gene_families)):
    if blast[i]=='Yes':
        os.system("python3 Scripts/run_blast.py %s %s %s %s" % (proteome, gene_families_db[i], evalues[i], threads))
        name = gene_families_db[i].replace("_db.fasta","_blast_otuput.txt")
        # check if blast output is empty then create blank file
        if os.stat(name).st_size == 0:
            with open(gene_families_db[i].replace("_db.fasta","_parsed_blast.txt"), 'w') as fp:
                fp.write('delete\t1\t10\t10\tblastp\n')
                fp.write('delete\t1\t10\t10\tblastp\n')
                fp.close()
        else:
            outname = gene_families_db[i].replace("_db.fasta","")
            os.system("perl Scripts/get_blastp_parsed_newv2.pl %s %s %s" % (name, outname, evalues[i]))
        with open(gene_families_db[i].replace("_db.fasta","_parsed_blast.txt"), 'a') as fp:
            fp.write('delete\t1\t10\t10\tblastp\n')
            fp.write('delete\t1\t10\t10\tblastp\n')
            fp.close()
    else:# create a blank file called GeneFamilyName_parsed_blast.txt
        with open(gene_families_db[i].replace("_db.fasta","_parsed_blast.txt"), 'w') as fp:
            fp.write('delete\t1\t10\t10\tblastp\n')
            fp.write('delete\t1\t10\t10\tblastp\n')
            fp.close()

# read InterproScan output
tsv = pd.read_csv(proteome+'.tsv', sep='\t', header=None, names = list(range(0,14)))

# For each gene family, produce a table with the proteins find in interproscan output if able:
# Check if there is Pfam domain, if not check if there is InterPro domain. Also, if there are more than one ID then retrieve all of them.
for i in range(len(InterPro)):
    name= gene_families_db[i].replace("_db.fasta","_parsed_domain.txt")
    # Check if there is Pfam domain
    if str(Pfam[i])[0]=='P':
        # Obtain a table with all the pfam domains
        l=[]
        if "&&" not in Pfam[i].split():
            for j in Pfam[i].split():
                if j[0]=='P':
                    l.append(tsv.loc[tsv[4] == j])
            new_table = pd.concat(l)[[0, 6, 7, 2, 3]]
        else:
            for j in Pfam[i].split():
                if j[0]=='P':
                    l.append(tsv.loc[tsv[4] == j].drop_duplicates(subset=[0], keep='first'))
            new_table = pd.concat(l)
            new_table = new_table[new_table.duplicated(subset=[0], keep=False)][[0, 6, 7, 2, 3]]
        if len(new_table)!=0:
            new_table.to_csv(path_or_buf= name ,header=None, index=None, sep='\t', mode='a')
        else:  # if no matches then write a blank file
            with open(name, 'w') as fp: 
                fp.write('delete\t1\t10\t10\tPfam\n')
                fp.write('delete\t1\t10\t10\tPfam\n')
                fp.close()       
    # if there is no Pfam then Interpro:
    elif str(InterPro[i])[0]=='I':
        # Obtain a table with all the domains
        l=[]
        if "&&" not in InterPro[i].split():
            for j in InterPro[i].split():
                if j[0]=='I':
                    l.append(tsv.loc[tsv[11] == j])
            new_table = pd.concat(l)[[0, 6, 7, 2, 3]]
        else:
            for j in InterPro[i].split():
                if j[0]=='I':
                    l.append(tsv.loc[tsv[11] == j].drop_duplicates(subset=[0], keep='first'))
            new_table = pd.concat(l)
            new_table = new_table[new_table.duplicated(subset=[0], keep=False)][[0, 6, 7, 2, 3]]
        if len(new_table)!=0:
            new_table.to_csv(path_or_buf= name ,header=None, index=None, sep='\t', mode='a')
        else:  # if no matches then write a blank file
            with open(name, 'w') as fp: 
                fp.write('delete\t1\t10\t10\tPfam\n')
                fp.write('delete\t1\t10\t10\tPfam\n')
                fp.close() 
    else: # otherwise create a blank file
        with open(name, 'w') as fp:
            fp.write('delete\t1\t10\t10\tPfam\n')
            fp.write('delete\t1\t10\t10\tPfam\n')
            fp.close()
    with open(name, 'a') as fp:
        fp.write('delete\t1\t10\t10\tPfam\n')
        fp.write('delete\t1\t10\t10\tPfam\n')
        fp.close()  

# Merge blast and interproscan results into a table and retrieve the fasta sequences
for i in range(len(gene_families)):
    blast= gene_families_db[i].replace("_db.fasta","_parsed_blast.txt")
    inter=  gene_families_db[i].replace("_db.fasta","_parsed_domain.txt")
    prot_avg= prot_avg_len[i]
    pept_avg= pept_avg_len[i]
    directory = gene_families_db[i].replace(gene_families[i]+"_db.fasta", "Result")
    os.system("mkdir -p %s" % (directory))
    output = directory+"/"+gene_families[i]+'_db.fasta'
    os.system("python3 Scripts/merge_blast_domain.py %s %s %s %s %s %s" % (proteome, blast, inter, prot_avg, pept_avg, output))

# Produce HMM profile of the gene family db (good) + the parsed sequences
for i in range(len(gene_families)):
    fasta = gene_families_db[i].replace(gene_families[i]+"_db.fasta","Result/"+gene_families[i]+"_db.fasta")
    alignment = fasta.replace("_db.fasta","_alignment.txt")
    profile = fasta.replace(".fasta",".hmm")
    os.system("mafft --auto %s > %s" % (fasta, alignment))
    os.system("hmmbuild %s  %s" % (profile, alignment))
    os.system("rm %s" % (alignment))

# Run bitacora
for i in range(len(gene_families)):
    directory = gene_families_db[i].replace(gene_families[i]+"_db.fasta","Result/")
    os.chdir("%s" % (directory))
    os.system("bash %s -q %s -g %s -f %s -p %s -n %s -t %s" % (run_bitacora, directory, genome, gff, proteome, genome_name, threads))

os.chdir("%s" % (path))
os.system("python3 Scripts/table_results.py %s" % (genome_name))
