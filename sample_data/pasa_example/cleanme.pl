#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;


## we delete all files we don't need in this directory. Be careful in case users try running it somewhere else, outside this dir.
chdir $FindBin::Bin or die "error, cannot cd to $FindBin::Bin";



my @files_to_keep = qw (
cleanme.pl
genome.fasta.gz
pasa_assemblies.fasta.gz
pasa_assemblies.gff3.gz
runMe.sh
pasa_assemblies_described.txt.gz
Makefile
                        );


my %keep = map { + $_ => 1 } @files_to_keep;


foreach my $file (<*>) {
	
	if (! $keep{$file}) {
		print STDERR "-removing file: $file\n";
		unlink($file);
	}
}


`rm -rf ./pasa_assemblies.fasta.transdecoder_dir/`;
`rm -rf ./pasa_assemblies.fasta.transdecoder_dir.__checkpoints`;
`rm -rf ./pasa_assemblies.fasta.transdecoder_dir.__checkpoints_longorfs`;

`rm -rf ./pasa.transdecoder_workdir*`;

exit(0);
