# GAGA-project
Gene family annotation pipeline:

1. Conduct a first search of gene family members using homology-based methods
2. Identify the characteristic protein domain of the gene family
3. Create a database with the identified and curated sequences
4. Annotate and identify genes in the genome assembly
  Using BITACORA https://github.com/molevol-ub/bitacora
5. Validate the annotations

## Files:

**/data/GAGA_gene_families.xlsx**  Table with gene families information such: Function/Classification, InterPro Domain, Pfam Domain, etc.
**/data/Mpha_annotation_fixed.renamed.representative.pep**  Proteome of Monomorium pharaonis to test the pipeline.
**/data/OBP_db.fasta**  Fasta file with the OBP gene family sequences annotated for other ants and insects.
**/scripts/run_analysis.py**  Python script that runs a first analysis of the pipeline to identify curated sequences, using Blastp and HMMER.
**/scripts/run_blast.py**  Python script that performs a blastp of a given gene_family proteins and a proteome.
**/scripts/parse_blast.py**  Python script that parses the obtained blastp output and generates a table with the proteins of interest.
