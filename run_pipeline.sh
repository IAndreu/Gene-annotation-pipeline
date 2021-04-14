###### EDIT VARIABLES ######

path="/path/to/Gene-annotation-pipeline"
genome_name="Genome-name"
genome_directory="${path}/Data/Genomes/Genome-name"
proteome="${genome_directory}/Proteome-in-fasta-file"
gff="${genome_directory}/gff3-file"
genome="${genome_directory}/Genome-in-fasta-file"
run_bitacora="${path}/bitacora-master/runBITACORA_command_line.sh"
num_threads=4

############################

python3 Scripts/run_analysis.py $path $proteome $genome_directory $gff $genome $run_bitacora $genome_name $num_threads
