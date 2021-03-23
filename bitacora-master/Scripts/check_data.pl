#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use lib dirname (__FILE__);
use Readgff qw (readgff);
my $dirname = dirname(__FILE__);

# Script to check if input data is OK and parseable by BITACORA

# usage: perl Scripts/check_data.pl $gff $genome $proteins query_directory


my ($line, $name, $nameout);
my $gff = $ARGV[0];
my $genome = $ARGV[1];
my $proteome = $ARGV[2];
my $querydir = $ARGV[3];
my $gemoma = $ARGV[4];

print "\n----------------- Checking input data and prerequisites\n";

# Checking GFF

my ($gffgeneref, $gffcdsref, $gffscafcdsref) = &readgff($gff);
my %gffcds = %$gffcdsref; 
my %gffgene = %$gffgeneref;
my %gffscafcds = %$gffscafcdsref;

## Checking if GFF contains overlapped genes (i.e. isoforms or bad-annotations)

my $overlapped = "0";
system("perl $dirname/Tools/get_overlapping_genes_fromgff.pl $gff");
open (Ofile, "<", "$gff\_overlapping_genes.txt");
while (<Ofile>){
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	$overlapped++;
}
close Ofile;
if ($overlapped > 0	){
	print "\nWarning: The input GFF contains coding sequences in overlapping positions (such as isoforms). BITACORA will continue running but take into account that the final annotation will include all isoforms included in the input GFF and protein fasta files (ignore this message if the input protein fasta contains one representative sequence per gene)\nYou can check the overlapping genes in file $gff\_overlapping_genes.txt\n\n";
}

print "GFF parsed correctly\n";


#Checking if IDs from Fasta and GFF match

my $checkgenome = "0";
my $scafnames = "";
open (File , "<", $genome); 
while (<File>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	if ($line =~ />(\S+)/){
		my $gene = $1;
		$scafnames .= "$gene\_\_";
		$checkgenome++;
	}	
}
close File;

if ($checkgenome == 0){
	die "ERROR in $dirname/check_data.pl: Genome file is not found or empty! $genome\n";
}

my $checkprot = "0";
open (File , "<", $proteome); 
while (<File>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);

	if ($line =~ />(\S+)/){
		my $gene = $1;
		$checkprot++;
		if (exists $gffgene{$gene}){
			#OK
			my @subline = split(/\t/, $gffgene{$gene});
			if ($scafnames =~ /$subline[0]\_\_/){
				#OK
			} else {
				die "ERROR in $dirname/check_data.pl: Genome sequence $subline[0] from GFF is not found in the genome fasta file\n";
			}

			if (exists $gffcds{$gene}){
				#OK
			} else {
				die "ERROR in $dirname/check_data.pl: Gene $gene does not contain any CDS sequence in GFF file. Please check your GFF\n";
			}


		} elsif ($gene =~ /(\S+)\-P(\S+)/) {
			my $genenew = "$1"."-R$2"; ## Change P for R if protein is named as transcript

			if (exists $gffgene{$genenew}){
				#OK
				my @subline = split(/\t/, $gffgene{$genenew});
				if ($scafnames =~ /$subline[0]\_\_/){
					#OK
				} else {
					die "ERROR in $dirname/check_data.pl: Genome sequence $subline[0] from GFF is not found in the genome fasta file\n";
				}

				if (exists $gffcds{$genenew}){
					#OK
				} else {
					die "ERROR in $dirname/check_data.pl: Gene $gene does not contain any CDS sequence in GFF file. Please check your GFF\n";
				}
			}
			else {
				die "ERROR in $dirname/check_data.pl: Protein gene $gene nor $genenew is not found in the GFF3\nIf your GFF comes from NCBI, you probably need to reformat the GFF to assign protein ID in mRNA and CDS fields. For that, you can use the following Script and input the formatted GFF in BITACORA\nperl $dirname/Tools/reformat_ncbi_gff.pl $gff\n\nOtherwise, if this error keeps showing, You can run the following script to check which proteins are not annotated in the GFF, or contain different name in mRNA|transcript field in the GFF to explore the data. You can use the outfile_ingff.fasta tu use only annotated sequences in BITACORA\nperl $dirname/Tools/get_proteins_notfound_ingff.pl $gff $proteome\n\nFinally, if all your proteins are named different as in the GFF3, you can codify the protein sequences directly from your GFF3 file with the following script, and use the resulting proteins in BITACORA:\nperl $dirname/gff2fasta_v3.pl $genome $gff $gff\n\n";
			}

		} else {
			die "ERROR in $dirname/check_data.pl: Protein gene $gene is not found in the GFF3\nIf your GFF comes from NCBI, you probably need to reformat the GFF to assign protein ID in mRNA and CDS fields. For that, you can use the following Script and input the formatted GFF in BITACORA\nperl $dirname/Tools/reformat_ncbi_gff.pl $gff\n\nOtherwise, if this error keeps showing, You can run the following script to check which proteins are not annotated in the GFF, or contain different name in mRNA|transcript field in the GFF to explore the data. You can use the outfile_ingff.fasta tu use only annotated sequences in BITACORA\nperl $dirname/Tools/get_proteins_notfound_ingff.pl $gff $proteome\n\nFinally, if all your proteins are named different as in the GFF3, you can codify the protein sequences directly from your GFF3 file with the following script, and use the resulting proteins in BITACORA:\nperl $dirname/gff2fasta_v3.pl $genome $gff $gff\n\n";
		}
	}	
}
close File;

if ($checkprot == 0){
	die "ERROR in $dirname/check_data.pl: Protein file is not found or empty! $proteome\n";
}

print "Genome and protein fasta parsed correctly\n";


# Check if hmmer is installed

system ("hmmsearch -h > testHMMER.out 2> testHMMER.err");
my $hmmerr = 0;
open (File, "<", "testHMMER.err");
while(<File>){
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);	
	$hmmerr++;
}
close File;
my $hmmok = 0;
open (File, "<", "testHMMER.out");
while(<File>){
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);	
	if ($line =~ /domtblout/){
		$hmmok++;
	}
}
close File;

if ($hmmerr == 0 && $hmmok > 0){
	#OK
} else {
	die "ERROR in $dirname/check_data.pl: hmmsearch could not be found\nAre you sure you set the path to HMMER bin correctly and it is properly installed?\n";
}

print "HMMER installed correctly\n";

# Check if blast is installed

system ("blastp -help > testBLAST.out 2> testBLAST.err");
system ("tblastn -help >> testBLAST.out 2>> testBLAST.err");
system ("makeblastdb -help >> testBLAST.out 2>> testBLAST.err");
my $blasterr = 0;
open (File, "<", "testBLAST.err");
while(<File>){
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);	
	$blasterr++;
}
close File;
my $blastok = 0;
open (File, "<", "testBLAST.out");
while(<File>){
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);	
	if ($line =~ /outfmt/){
		$blastok++;
	}
}
close File;

if ($blasterr == 0 && $blastok > 0){
	#OK
} else {
	die "ERROR in $dirname/check_data.pl: BLAST binaries could not be found\nAre you sure you set the path to BLAST bin correctly and it is properly installed?\n";
}

print "BLAST installed correctly\n";


# Check if GeMoMa is installed (if needed)

if ($gemoma =~ /T/){
	my $gpath = $ARGV[5];
	system ("java -jar $gpath CLI > testGemoma.out 2> testGemoma.err");
	my $gemerr = 0;
	my $gemomav = "0";
	my $tempv = "";
	open (File, "<", "testGemoma.err");
	while(<File>){
		chomp;
		$line = $_;
		next if ($line !~ /\S+/);	
		if ($line =~ /Unable/){
	#		$gemerr++;
		}
		if ($line =~ /(GeMoMa\-\S+).jar/){
			$tempv = $1;
			if ($tempv =~ /GeMoMa\-(\d)\.(\d)\.(\d)/){
				$gemomav = "$1"."$2"."$3";
			} elsif ($tempv =~ /GeMoMa\-(\d)\.(\d)/){
				$gemomav = "$1"."$2"."0";
			} elsif ($tempv =~ /GeMoMa\-(\d)/){
				$gemomav = "$1"."0"."0";
			}
		}
	}
	close File;
	my $gemok = 0;
	open (File, "<", "testGemoma.out");
	while(<File>){
		chomp;
		$line = $_;
		next if ($line !~ /\S+/);	
		if ($line =~ /GeMoMa/){
			$gemok++;
		}
	}
	close File;

	if ($gemerr == 0 && $gemok > 0){
		#OK
	} else {
		die "ERROR in $dirname/check_data.pl: GeMoMa could not be found\nAre you sure you set the path to GeMoMa jar file correctly and java is installed in your system?\n";
	}

	if ($gemomav > 100){
		#OK
	} else {
		die "ERROR in $dirname/check_data.pl: GeMoMa version $tempv could not be assessed\nAre you sure you are using a version of GeMoMa tested and working in BITACORA? Check the latest supported version in the readme\n";
	}
	if ($gemomav < 164){
		die "ERROR in $dirname/check_data.pl: An old version of GeMoMa has been detected $tempv\nPlease download a more recent version (at least GeMoMa-1.6.4)\n";
	}	

	print "$tempv is installed correctly\n";

	open (Gemfile, ">", "GeMoMa_version.txt");
	print Gemfile "$gemomav\n";
	close Gemfile;
}



# Check if the query dir contains proper renamed fasta and HMM files

my $nfile = 0;
my $dbfiles = "";
system("ls $querydir\/*_db.fasta > Query_genes.txt");
open(File, "<", "Query_genes.txt");
while(<File>){
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	my $id = "";
	if ($line =~ /([^\/]+)\_db.fasta/){
		$id = $1;
	} else {
		die "Are you sure you renamed you Query DB to QUERY_db.fasta? Avoid using special characters in QUERY as\nCannot find QUERY_db.fasta in $line\n";
	}

	system ("wc -l $line > testdb.out 2> testdb.err");

	my $dbout = 0;
	open (File2, "<", "testdb.out");
	while(<File2>){
		chomp;
		my $line2 = $_;
		next if ($line !~ /\S+/);	
		if ($line2 =~ /(\d+) /){
			$dbout = $1;
		} else {
			die "ERROR in $dirname/check_data.pl: reading testdb.out in $line2\n";
		}
	}
	close File2;

	if ($dbout == 0){
		die "ERROR in $dirname/check_data.pl: $line is empty!\n";		
	}

	# Checking now HMM file

	system ("wc -l $querydir\/$id\_db.hmm > testhmm.out 2> testhmm.err");

	$hmmerr = 0;
	open (File2, "<", "testhmm.err");
	while(<File2>){
		chomp;
		my $line2 = $_;
		next if ($line2 !~ /\S+/);	
		$hmmerr++;
	}
	close File2;
	$hmmok = 0;
	open (File2, "<", "testhmm.out");
	while(<File2>){
		chomp;
		my $line2 = $_;
		next if ($line2 !~ /\S+/);	
		#my @subline = split (/\s/, $line);
		#$hmmok = $subline[0];
		if ($line2 =~ /(\d+) /){
			$hmmok = $1;
		} else {
			die "ERROR in $dirname/check_data.pl: reading testhmm.out in $line2\n";
		}

	}
	close File2;

	if ($hmmerr == 0 && $hmmok > 0){
		#OK
	} elsif ($hmmerr == 0 && $hmmok == 0){
		die "ERROR in $dirname/check_data.pl: $querydir\/$id\_db.hmm is empty!\n";		
	} else {
		die "ERROR in $dirname/check_data.pl: It was not found the expected HMM file: $querydir\/$id\_db.hmm\n";
	}

	$nfile++;
	$dbfiles .= "$id ";

}
close File;

if ($nfile > 0){
	#OK
} else {
	die "ERROR in $dirname/check_data.pl: It was not found any file formatted correctly in $querydir\nFiles must be named as ID_db.fasta and ID_db.hmm as detailed in manual\n";
}

print "Query directory files found and named correctly in $querydir\nFound $nfile query databases to search: $dbfiles\n";

print "Everything looks fine\n----------------- DONE\n\n";

system("rm testHMMER.err testHMMER.out testBLAST.err testBLAST.out testhmm.err testhmm.out testdb.err testdb.out testGemoma.err testGemoma.out");



