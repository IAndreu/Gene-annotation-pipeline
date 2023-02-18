###### EDIT VARIABLES ######

path="/path/to/Gene-annotation-pipeline"
table="/path/to/gene_families.xlsx"
gene_fam_db="/path/to/Gene_families_db"
genome_name="Genome-name"
output_directory="${path}/Data/Genomes/Genome-name"
proteome="${genome_directory}/Proteome-in-fasta-file"
gff="${genome_directory}/gff3-file"
genome="${genome_directory}/Genome-in-fasta-file"
interpro="${genome_directory}/predicted-domains-from-InterPro-in-TSV-format"
run_bitacora="${path}/bitacora-master/runBITACORA_command_line.sh"
num_threads=4
############################

python3 "$path"/Scripts/run_analysis.py $path $table $gene_fam_db $proteome $interpro $gff $genome $run_bitacora $genome_name $output_directory $num_threads
