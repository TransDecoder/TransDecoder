#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;


## we delete all files we don't need in this directory. Be careful in case users try running it somewhere else, outside this dir.
chdir $FindBin::Bin or die "error, cannot cd to $FindBin::Bin";



my @files_to_keep = qw (
cleanme.pl
runMe.sh
Trinity.fasta.gz
Makefile
runMe.Docker.sh
genome_alignments.gmap.gff3.gz
                        );


my %keep = map { + $_ => 1 } @files_to_keep;


foreach my $file (<*>) {
	
	if (! $keep{$file}) {
		print STDERR "-removing file: $file\n";
		unlink($file);
	}
}


`rm -rf ./Trinity.fasta.transdecoder_dir/`;
`rm -rf ./Trinity.fasta.transdecoder_dir.__checkpoints`;
`rm -rf ./Trinity.fasta.transdecoder_dir.__checkpoints_longorfs`;
`sudo rm -rf ./transdecoder.docker_outdir`;


exit(0);
