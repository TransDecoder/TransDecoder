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

my $usage = "\n\nusage: $0 gff3_file [ gff3_file, ... ]\n\n"
    . "\tgenes are indexed by model feat_name.\n\n";


my @gff3_files = @ARGV;
unless (@gff3_files) { die $usage; }

my $index_file = "gene_structures.inx";
if (scalar @gff3_files == 1) {
    $index_file = $gff3_files[0] . ".inx";
}

my $gene_obj_indexer = new Gene_obj_indexer( { "create" => $index_file } );

my %seen; #track gene_ids already visited

foreach my $gff3_file (@gff3_files) {

#    print STDERR "// indexing $gff3_file\n";
    
    ## associate gene identifiers with contig id's.
    my $temp_inx = "tmp.inx";
    my $temp_gene_indexer = new Gene_obj_indexer( { "create" => $temp_inx } );
    
    my $asmbl_id_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($gff3_file, $temp_gene_indexer);
    
    
    foreach my $asmbl_id (sort keys %$asmbl_id_to_gene_list_href) {
        
        my @gene_ids = @{$asmbl_id_to_gene_list_href->{$asmbl_id}};
        
        #print "ASMBL: $asmbl_id, gene_ids: @gene_ids\n";
        
        foreach my $gene_id (@gene_ids) {
            
            if ($seen{$gene_id}) {
                croak "Error, already stored gene_id: [$gene_id], not allowed to have the same gene id multiple GFF3 files.\n";
            }
            $seen{$gene_id} = 1;
            
            my $gene_obj_ref = $temp_gene_indexer->get_gene($gene_id);
            
            foreach my $gene_obj ($gene_obj_ref, $gene_obj_ref->get_additional_isoforms()) {
                
                $gene_obj->delete_isoforms(); # unbundle the model object!
                
                my $model_id = $gene_obj->{Model_feat_name};
                $gene_obj_indexer->store_gene($model_id, $gene_obj);
                print STDERR "\r Indexed $model_id ";
            }
        }
        
    }
    print STDERR "\n";

    unlink ($temp_inx);
}


exit(0);

