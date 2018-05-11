#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use GFF3_utils2;
use Carp;

my $usage = "\n\nusage: $0 gff3_file\n\n";

my $gff3_file = $ARGV[0] or die $usage;


main: {

    my $gene_obj_indexer_href = {};
    
    my $asmbl_id_to_gene_list_href = &GFF3_utils2::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer_href);

    my %seen;

    
    foreach my $asmbl_id (sort keys %$asmbl_id_to_gene_list_href) {
        
        my @gene_ids = @{$asmbl_id_to_gene_list_href->{$asmbl_id}};
        
        #print "ASMBL: $asmbl_id, gene_ids: @gene_ids\n";
        my @gene_entries;
        
        foreach my $gene_id (@gene_ids) {
            
            my $gene_obj = $gene_obj_indexer_href->{$gene_id};
            
            my $cds_token = &get_CDS_token($gene_obj);
            if (! $seen{$cds_token}) {
                
                print $gene_obj->to_GFF3_format(source => "transdecoder") . "\n";
                
                $seen{$cds_token} = 1;
                
            }
            else {
                print STDERR "-ignoring entry $cds_token, already represented by another transcript\n";
            }
            
        }
    }
    
    
    exit(0);

}



####
sub get_CDS_token {
    my ($gene_obj) = @_;

    my $cds_text = $gene_obj->{asmbl_id};
    
    my @exons = $gene_obj->get_exons();
    foreach my $exon (@exons) {
        if (my $cds_obj = $exon->get_CDS_obj()) {
            my ($cds_end5, $cds_end3) = $cds_obj->get_coords();
            $cds_text .= ":$cds_end5-$cds_end3";
        }
    }

    return($cds_text);
}
    
