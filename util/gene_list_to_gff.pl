#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use GFF3_utils2;
use Carp;

$|++;

my $usage = "\n\nusage: $0 gene_ID_list_file genes.gff3\n\n";

my $gene_ID_list = $ARGV[0] or die $usage;
my $gff3_file = $ARGV[1] or die $usage;

my $gene_obj_indexer_href = {};

&GFF3_utils2::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer_href);

open (my $fh, $gene_ID_list) or die $!;
while (<$fh>) {
    unless (/\w/) { next;}
    chomp;
	my ($gene_id, @rest) = split (/\t/);
    $gene_id =~ s/\s+//;
    
    my $gene_obj = $gene_obj_indexer_href->{$gene_id} or confess "Error, didn't find gene obj for $gene_id in $gff3_file";
    
    print $gene_obj->to_GFF3_format(source => "transdecoder") . "\n";
    
}
close $fh;


exit(0);




