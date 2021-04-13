# Pipeline for gene family annotation across hundreds of genomes

This pipeline automates and standardizes gene family annotation for a number of gene families in dataset of newly generated genomes. This pipeline allows to obtain the most accurate number of gene copies and minimizes methodical biases that would otherwise perturb downstream comparative analyses. _BITACORA_ and _GeMoMa_ are the main tools used for the identification and annotation of gene families in genome assemblies, toghether with a first preprocessing step that identifies and curate gene models using _Blastp_ and _InterProScan_ based on an input file with information of the gene families to annotate.

## Prerequisites
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

## Installation

The pipeline is distributed as a multiplatform shell script (run_pipeline.sh) that calls several other Python and Perl scripts, which include all functions responsible of performing all pipeline tasks. Hence, it does not require any installation or compilation step.

To run the pipeline edit the master script "run\_pipeline.sh" variables described in Prerequisites, and also the master script of _BITACORA_ "runBITACORA_command_line.sh" variables.




