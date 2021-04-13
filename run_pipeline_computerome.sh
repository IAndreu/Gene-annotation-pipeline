#!/bin/sh
### Note: No commands may be executed until after the #PBS lines
### Account information
#PBS -W group_list=ku_00039 -A ku_00039
### Job name (comment out the next line to get the name of the script used as the job name)
##PBS -N test
### Output files (comment out the next 2 lines to get the job name used instead)
##PBS -e test.err
##PBS -o test.log
### Only send mail when job is aborted or terminates abnormally
#PBS -m n
### Number of nodes
#PBS -l nodes=1:ppn=4
### Memory
#PBS -l mem=40gb
### Requesting time - format is <days>:<hours>:<minutes>:<seconds> (here 12:00:00, 12 hours)
#PBS -l walltime=24:00:00
### Forward X11 connection (comment out if not needed)
##PBS -X

# Go to the directory from where the job was submitted (initial directory is $HOME)
echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

### Here follows the user commands:
# Define number of processors
NPROCS='wc -l < $PBS_NODEFILE'
echo This job has allocated $NPROCS nodes

# Load all required modules for the job
module load tools
module load ngs
module load anaconda3/4.4.0
module load openjdk/15.0.1
module load perl
module load ncbi-blast/2.2.31+
module load hmmer/3.2.1
module load java/1.8.0 
module load interproscan/5.47-82.0
module load mafft/7.453

# This is where the work is done
# Make sure that this script is not bigger than 64kb ~ 150 lines, otherwise put in seperat script and execute from here

###### EDIT VARIABLES ######
path="/home/projects/ku_00039/people/igngod/Project"
proteome="${path}/Data/Genomes/Mpha/Mpha_annotation_fixed.renamed.representative.pep"
genome_directory="${path}/Data/Genomes/Mpha/"
gff="${path}/Data/Genomes/Mpha/Mpha_annotation_fixed.renamed.gff3"
genome="${path}/Data/Genomes/Mpha/GAGA-0245_Monomorium_pharaonis_chromosome.fasta"
run_bitacora="${path}/bitacora-master/runBITACORA_command_line.sh"
genome_name="Mpha"
num_threads=4
############################

python3 Scripts/run_analysis.py $proteome $genome_directory $gff $genome $run_bitacora $genome_name $num_threads
