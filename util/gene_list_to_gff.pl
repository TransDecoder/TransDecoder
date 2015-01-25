#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use Gene_obj_indexer;
use GFF3_utils;
use Carp;

$|++;

my $usage = "\n\nusage: $0 gene_ID_list_file gene_db.inx_file  \n\nNote, you must first run index_gff3_files.pl\n\n";

my $gene_ID_list = $ARGV[0] or die $usage;
my $inx_file = $ARGV[1] or die $usage;

my $gene_obj_indexer = new Gene_obj_indexer( { "use" => "$inx_file" } );

open (my $fh, $gene_ID_list) or die $!;
while (<$fh>) {
    unless (/\w/) { next;}
    chomp;
	my ($gene_id, $com_name) = split (/\t/);
    $gene_id =~ s/\s+//;
    
    my $gene_obj = $gene_obj_indexer->get_gene($gene_id);
	
	if (defined ($com_name) && $com_name =~ /\w/) {
		$gene_obj->{com_name} = $com_name;
	}

    print $gene_obj->to_GFF3_format(source => "transdecoder") . "\n";
    
}
close $fh;


exit(0);




