#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use GFF3_utils2;
use Carp;

$|++;

my $usage = "\n\nusage: $0 gff3_file\n\n";

my $gff3_file = $ARGV[0] or die $usage;


main: {

    my $gene_obj_indexer_href = {};
    
    my $asmbl_id_to_gene_list_href = &GFF3_utils2::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer_href);
    
    foreach my $asmbl_id (sort keys %$asmbl_id_to_gene_list_href) {
        
        my @gene_ids = @{$asmbl_id_to_gene_list_href->{$asmbl_id}};
        
        #print "ASMBL: $asmbl_id, gene_ids: @gene_ids\n";
        my @gene_entries;
        
        foreach my $gene_id (@gene_ids) {
            
            my $gene_obj_ref = $gene_obj_indexer_href->{$gene_id};
            
            my ($lend, $rend) = sort {$a<=>$b} $gene_obj_ref->get_model_span();
            
            my $struct = { gene_obj => $gene_obj_ref,
                           lend => $lend,
                           rend => $rend,
                           length => $rend - $lend + 1,
            };
            
            push (@gene_entries, $struct);
            
            
        }
        
        @gene_entries = reverse sort {$a->{length}<=>$b->{length}} @gene_entries;
        
        my @largest_orfs = shift @gene_entries;
        
        while (@gene_entries) {
            my $next_gene = shift @gene_entries;
            
            my ($next_lend, $next_rend) = ($next_gene->{lend}, $next_gene->{rend});
            
            
            my $found_eclipsed = 0;
            
            foreach my $gene (@largest_orfs) {
                
                my ($lend, $rend) = ($gene->{lend}, $gene->{rend});
                
                if ($next_lend > $lend && $next_rend < $rend) {
                    ## eclipsed
                    my $model_feat_name = $gene->{gene_obj}->{Model_feat_name};
                    print STDERR "\tECLIPSE: $model_feat_name is eclipsed by longer ORF, removing it.\n";
                    $found_eclipsed = 1;
                    last;
                }
            }
            
            unless ($found_eclipsed) {
                push (@largest_orfs, $next_gene);
            }
        }
        

        foreach my $struct (@largest_orfs) {
            my $gene_obj = $struct->{gene_obj};
            
            print $gene_obj->to_GFF3_format(source => "transdecoder") . "\n";
        }
        
        
        
        
    }
    
    
    exit(0);

}
