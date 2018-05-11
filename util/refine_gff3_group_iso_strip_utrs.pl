#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use GFF3_utils2;
use Carp;

my $usage = "usage: $0  transdecoder.gff3\n\n";

my $gff3_file = $ARGV[0] or die $usage;


main: {

    my $indexer = {};
    my $asmbl_id_to_gene_list_href = &GFF3_utils2::index_GFF3_gene_objs($gff3_file, $indexer);

    
    
    foreach my $asmbl_id (keys %$asmbl_id_to_gene_list_href) {

        my @gene_ids = @{$asmbl_id_to_gene_list_href->{$asmbl_id}};

        my $primary_gene_obj;
        
        foreach my $gene_id (@gene_ids) {

            my $gene_obj = $indexer->{$gene_id};
            $gene_obj->trim_UTRs();

            if ($primary_gene_obj) {
                $primary_gene_obj->add_isoform($gene_obj);
            }
            else {
                $primary_gene_obj = $gene_obj;
            }
        }
    
        $primary_gene_obj->refine_gene_object();
        print $primary_gene_obj->to_GFF3_format() . "\n";
        
        
    }
    

    exit(0);
}


