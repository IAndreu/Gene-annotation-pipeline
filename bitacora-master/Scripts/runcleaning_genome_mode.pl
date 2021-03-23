#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
my $dirname = dirname(__FILE__);
#########################################################
#####
#
# Clean output folders
#
#####

## Input variables

my $name = $ARGV[0];
my $chemdir = $ARGV[1];

my @chemosensory;
system("ls $chemdir\/*_db.fasta > Query_genes.txt");
open(File, "<", "Query_genes.txt");
while(<File>){
	chomp;
	my $line = $_;
	next if ($line !~ /\S+/);
	my $id = "";
	if ($line =~ /([^\/]+)\_db.fasta/){
		$id = $1;
	} else {
		die "Are you sure you renamed you Query DB to QUERY_db.fasta?\nCannot find QUERY_db.fasta in $line\n";
	}
	push (@chemosensory, $id);
}
close File;

foreach my $chem (@chemosensory){
	system ("mkdir -p $chem\/Intermediate_files");
	system ("rm -rf $chem\/Intermediate_files/hmmer 2>/dev/null");
	
	# Moving Files
	system("mv $chem\/* $chem\/Intermediate_files/ 2>/dev/null");
	system("mv $chem\/Intermediate_files/*bed $chem\/");
	system("mv $chem\/Intermediate_files/*genomic_genes_trimmed* $chem\/ 2>/dev/null");
	system("mv $chem\/Intermediate_files/*genomic_proteins_trimmed* $chem\/ 2>/dev/null");	

	system ("rm -rf $chem\/*_renamed\.fasta $chem\/*\.fasta\.* 2>/dev/null");


}
#system ("rm -rf GeMoMa_version.txt 2>/dev/null");




