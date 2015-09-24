#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "\n\n\tusage: $0 transcripts.fasta.transdecoder.gff3\n\n\n";

my $gff3_file = $ARGV[0] or die $usage;


main: {

    open (my $fh, $gff3_file) or die $usage;
    while (<$fh>) {
        unless (/\w/) { next; }
        my @x = split(/\t/);
        my $type = $x[2];
        if ($type eq 'mRNA') {
            my $info = $x[8];
            if ($info =~ /complete/) {
                $info =~ /ID=([^;]+);/ or die "Error, cannot parse ID from $info";
                my $isoform_id = $1;
                print "$isoform_id\n";
            }
        }
    }
    close $fh;

    exit(0);
}

