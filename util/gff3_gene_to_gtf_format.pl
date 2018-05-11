#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use Fasta_retriever;
use GFF3_utils2;
use Carp;


my $usage = "usage: $0 file.GFF3 genome.fasta > file.GTF\n\n";

my $gff3_file = $ARGV[0] or die $usage;
my $genome_fasta = $ARGV[1] or die $usage;


main: {

    my $fasta_retriever = new Fasta_retriever($genome_fasta);
     
    my $inx_file = "$gff3_file.tmp.inx";
    my $gene_obj_indexer = {};
    
    my $asmbl_id_to_gene_list_href = &GFF3_utils2::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer);

    foreach my $asmbl_id (sort keys %$asmbl_id_to_gene_list_href) {
        
        ## get the genome sequence
        my $genome_seq = $fasta_retriever->get_seq($asmbl_id); 

        my @gene_ids = @{$asmbl_id_to_gene_list_href->{$asmbl_id}};
        
        #print "ASMBL: $asmbl_id, gene_ids: @gene_ids\n";
        
        foreach my $gene_id (@gene_ids) {
            
            ## note models of isoforms are bundled into the same gene object.
            my $gene_obj_ref = $gene_obj_indexer->{$gene_id} or confess "Error, no gene_obj for $gene_id in $gff3_file";
            
            foreach my $gene_obj ($gene_obj_ref, $gene_obj_ref->get_additional_isoforms()) {
                
                $gene_obj->delete_isoforms(); # unbundle the model object!
                
                my $gtf_text = "";

                eval {
                    $gtf_text = $gene_obj->to_GTF_format(\$genome_seq);
                };
                if ($@) {
                    # do it in pseudogene mode
                    $gene_obj->{is_pseudogene} = 1;
                    $gtf_text = $gene_obj->to_GTF_format(\$genome_seq);
                }
                print "$gtf_text\n";
            }
        }
    }

    unlink $inx_file;
    
    exit(0);
}

