# Pipeline for gene family annotation across hundreds of genomes

This pipeline automates and standardizes gene family annotation for a number of gene families in dataset of newly generated genomes. This pipeline allows to obtain the most accurate number of gene copies and minimizes methodical biases that would otherwise perturb downstream comparative analyses. _BITACORA_ and _GeMoMa_ are the main tools used for the identification and annotation of gene families in genome assemblies, toghether with a first step that identifies and curate gene models using _Blastp_ and _InterProScan_ based on an input file with information of the gene families to annotate.

### Contents
1. Prerequisites
2. Installation
3. Computational Requirements
4. Usage
    
    4.1 Preparing data
    
    4.2 Run pipeline
    
    4.3 Output
5. Example


## 1. Prerequisites
The requiered dependencies necessary for running the pipeline are:
- **Perl**: Perl is installed by default in most operating systems. See https://learn.perl.org/installing/ for installation instructions.
- **Python**: Download the newest version available from https://www.python.org/downloads/
- **BLAST**: Download blast executables from: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
- **HMMER**: The easiest way to install HMMER in your system is to type one of the following commands in your terminal:
```
brew install hmmer               # OS/X, HomeBrew
port install hmmer               # OS/X, MacPorts
apt install hmmer                # Linux (Ubuntu, Debian...)
dnf install hmmer                # Linux (Fedora)
yum install hmmer                # Linux (older Fedora)
conda install -c bioconda hmmer  # Anaconda
```
Or compile HMMER binaries from the source code: http://hmmer.org/

HMMER and BLAST binaries require to be added to the PATH environment variable. Specify the correct path to bin folders in the master script "runBITACORA_command_line.sh", if necessary.


- **GeMoMa**: By default, _BITACORA_ reconstructs new gene models using the GeMoMa algorithm (Keilwagen et al., 2016; 2018). The GeMoMa jar file (i.e. GeMoMa-1.6.4.jar) must be specified in GEMOMAP variable in "runBITACORA_command_line.sh". GeMoMa is implemented in Java using Jstacs and can be downloaded from: http://www.jstacs.de/index.php/GeMoMa. The latest version of GeMoMa tested in our pipeline is v1.6.4.
```
GEMOMAP=/path/to/GeMoMa.jar  (within runBITACORA.sh script)
```
- **MAFFT**: Download the newest version available from https://mafft.cbrc.jp/alignment/software/

Also, if the TSV output file obtained from InterProScan, of the genome predicted proteins included in the GFF in FASTA format is not provided, then InterProScan must be installed and the script "submit_interpro.sh" executed before running the pipeline.
- **InterProScan**: Download it from https://interproscan-docs.readthedocs.io/en/latest/UserDocs.html#obtaining-a-copy-of-interproscan

## 2. Installation

The pipeline is distributed as a multiplatform shell script (run\_pipeline.sh) that calls several other Python and Perl scripts, which include all functions responsible for performing all pipeline tasks. _BITACORA_ is already installed inside the pipeline, hence, it does not require any installation or compilation step apart from the ones stated in **Prerequisites**.

To run the pipeline edit the master script "run\_pipeline.sh" variables and also the script of _BITACORA_ "runBITACORA\_command\_line.sh" path variables: BLAST and HMMER (if needed), BITACORA script folder and GeMoMa jar file.

## 3. Computational requirements

The pipeline have been tested in UNIX-based platforms (both in Mac OS and Linux operating systems). Multiple threading can be set in blast searches, which is the most time-consuming step, by editing the variable THREADS in "run\_pipeline.sh".

For a typical good quality genome (~2Gb in size and ~10,000 scaffolds) and a 40Gb RAM machine, it is able to analyze a whole genome and a total of 94 gene and gene families in 10 hours using 4 cores.

## 4. Usage
#### 4.1. Preparing data
The input files required to run a full analysis (update the complete path to these files if needed in the master script "run\_pipeline.sh") are the following:
##### In the directory _"Gene-annotation-pipeline/Data/"_
- **Excel file called "gene_families.xlsx" that contains all the gene families names to be annotated together with the information requiered by the pipeline about them.**

Example of the format of this file:

| Function/Classification | Gene family | Blast | InterPro Domain | Pfam domain | Average number of genes |Average protein length | Peptide start length | E-value |
| ----------- | ----------- | ----------- | ----------- | ----------- | ----------- | ----------- | ----------- | ----------- |
| Chemosensory receptors | OR | Yes | IPR004117 | PF02949 | 300-400 | 400 | 40 | 1,00E-5 |
| Chemosensory and others | CD36 | Yes | IPR002159 | PF01130 | 9 | 350 | 40 | 1,00E-5 |
| ... | ... | ... | ... | ... | ... | ... | ... | ... |

 **NOTE:** Not all the fields from the table have to be filled for a gene family to be annotated. "Function/Classification", "InterPro Domain" and "Pfam domain" can be blank cells (specially if there is no domain information about the gene family). "Blast" cell can also be blank if there is no input fasta database of the gene family to run a *blastp*.

##### In the directory _"Gene-annotation-pipeline/Data/Gene_families/"_
- **Files containing the query protein databases (GENEFAMILY-NAME\_db.fasta) in FASTA format, where the “GENEFAMILY-NAME” label is a gene family name from the Excel file. The addition of ”_db” to the database name with its proper extension is mandatory.**

##### In the directory _"Gene-annotation-pipeline/Data/Genomes/"_
Here is where all the genomes that have to be annotated must be placed in a folder with their name and all the required files inside.  Example:

_.../Genomes/Genome-name1/files_

_.../Genomes/Genome-name2/files_

...

_.../Genomes/Genome-nameN/files_

Where the files for each genome are:

- **File with genomic sequences in FASTA format**
- **File with structural annotations in GFF3 format**
- **File with predicted proteins in FASTA format.**
- **File with predicted domains from InterPro in TSV format**

Those files are some of the variables that must be specified in the main script "run_pipeline.sh".
#### 4.2. Run pipeline

Once edited the scripts as indicated in the **Installation** step and prepared the data, the pipeline can be executed with the following command:
```
bash run_pipeline.sh
``` 
#### 4.3. Output

The pipeline generates for each GENOME AND GENE FAMILY, two outputs of interest:

1. The files called **GENE-FAMILY-NAME_parsed_sequences.fasta** and **GENE-FAMILY-NAME_merged.txt**. 

This output files are generated in the first step before running bitacora and are located in the directory _"Gene-annotation-pipeline/Data/Genomes/GENOME-NAME/gene\_families\_pipeline/GENE-FAMILY-NAME"_. The fasta file contains identified and curated sequences to be annotated and the txt file is a summary table of the sequences.

2. The output folder created by _BITACORA_.

It is located in the directory _"Gene-annotation-pipeline/Data/Genomes/GENOME-NAME/gene\_families\_pipeline/GENE-FAMILY-NAME/Result"_. The description of this output can be found in the Documentation of _BITACORA_ https://github.com/molevol-ub/bitacora.git

For a general overview of the results, a file called **table_results.txt** in the main directory is produced and contains the number of genes identified in the first step, after _BITACORA_, and the expected, per family and genome.

#### 5. Example
An example to run the pipeline can be found in Example folder. It consists of two chemosensory-related gene families in insects: Odorant receptors (ORs), and the CD36-SNMP gene family; will be searched in the chromosome 2R of Drosophila melanogaster. The GFF3 and protein files are modified from original annotations, deleting some gene models, to allow that BITACORA can identify novel not-annotated genes.

First, unzip the Example_files.zip to obtain the necessary files. Move the fasta databases files from _"DB"_ folder to _"Gene-annotation-pipeline/Data/Gene_families"_. Then, create a folder named "Drosophila\_melanogaster" (the genome name) in the directory _"Gene-annotation-pipeline/Data/Genomes"_. Finally, move the files from _Files"_ folder to "Gene-annotation-pipeline/Data/Genomes/Drosophila\_melanogaster"_.

To run the example, edit the script "runBITACORA\_command\_line.sh" from _BITACORA_ to add the paths variables and also edit the "run_pipeline.sh" script to add your path to the pipeline and the genome name (in this example "Drosophila\_melanogaster"). Then in the command line:

```
bash run_pipeline.sh
``` 

#### 6. Citation

Joel Vizueta, Alejandro Sánchez-Gracia, Julio Rozas. BITACORA: A comprehensive tool for the identification and annotation of gene families in genome assemblies. Molecular Ecology Resources, 2020. https://doi.org/10.1111/1755-0998.13202

Joel Vizueta, Julio Rozas, Alejandro Sánchez-Gracia; Comparative Genomics Reveals Thousands of Novel Chemosensory Genes and Massive Changes in Chemoreceptor Repertories across Chelicerates, Genome Biology and Evolution, Volume 10, Issue 5, 1 May 2018, Pages 1221–1236, https://doi.org/10.1093/gbe/evy081

Moreover, if you use GeMoMa, please cite:

J. Keilwagen, M. Wenk, J. L. Erickson, M. H. Schattat, J. Grau, and F. Hartung. Using intron position conservation for homology-based gene prediction. Nucleic Acids Research, 2016. https://doi.org/10.1093/nar/gkw092

J. Keilwagen, F. Hartung, M. Paulini, S. O. Twardziok, and J. Grau Combining RNA-seq data and homology-based gene prediction for plants, animals and fungi. BMC Bioinformatics, 2018. https://doi.org/10.1186/s12859-018-2203-5





