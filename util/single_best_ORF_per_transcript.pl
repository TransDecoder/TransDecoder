#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use Gene_obj_indexer;
use GFF3_utils;
use Carp;
use Data::Dumper;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);



my $usage = <<__EOUSAGE__;

################################################################################
#
# Required:
#
#  --gff3_file <string>          gff3 file w/ predicted ORFs
#
# Optional:
#
#  --blast_hits <string>         blast hits file
#
#  --pfam_hits <string>          pfam hits file
#
################################################################################


__EOUSAGE__

    ;

my $gff3_file;
my $blast_hits_file;
my $pfam_hits_file;


&GetOptions('gff3_file=s' => \$gff3_file,
            'blast_hits=s' => \$blast_hits_file,
            'pfam_hits=s' => \$pfam_hits_file,
    );

unless ($gff3_file) {
    die $usage;
}


main: {

    my %blast_hits;
    if ($blast_hits_file) {
        %blast_hits = &parse_blastp_hits_file($blast_hits_file);
    }

    my %pfam_hits;
    if ($pfam_hits_file) {
        %pfam_hits = &parse_pfam_hits_file($pfam_hits_file);
    }
        
    my $gene_obj_indexer_href = {};
    
    my $asmbl_id_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer_href);
    
    foreach my $asmbl_id (sort keys %$asmbl_id_to_gene_list_href) {
        
        my @gene_ids = @{$asmbl_id_to_gene_list_href->{$asmbl_id}};
        
        #print "ASMBL: $asmbl_id, gene_ids: @gene_ids\n";
        my @gene_entries;
        
        foreach my $gene_id (@gene_ids) {
            
            my $gene_obj_ref = $gene_obj_indexer_href->{$gene_id};

            my $model_id = $gene_obj_ref->{Model_feat_name};

            my $homology_count = 0;
            if ($blast_hits{$model_id}) {
                $homology_count++;
            }
            if ($pfam_hits{$model_id}) {
                $homology_count++;
            }
            
            
            my $struct = { gene_obj => $gene_obj_ref,
                           length => $gene_obj_ref->get_CDS_length(),
                           homology_count => $homology_count,
            };
            
            push (@gene_entries, $struct);
            
            
        }
        
        @gene_entries = sort {
        
            $b->{homology_count} <=> $a->{homology_count}
            ||
            $b->{length} <=> $a->{length}
            
        } @gene_entries;
        
        if (scalar @gene_entries > 1) {
            print STDERR "ORFs prioritized as:\n";
            foreach my $entry (@gene_entries) {
                print STDERR "\t" . join("\t", $entry->{gene_obj}->{Model_feat_name}, 
                                         "homology_count: " . $entry->{homology_count},
                                         "len: " . $entry->{length}) . "\n";
            }
            
        }
        
        my $best_gene_entry = shift @gene_entries;

        my $gene_obj = $best_gene_entry->{gene_obj};
            
        print $gene_obj->to_GFF3_format(source => "transdecoder") . "\n";
        
        foreach my $remaining_gene (@gene_entries) {
            print STDERR "-discarding " . $remaining_gene->{gene_obj}->{Model_feat_name} . " as not single best orf on trans\n";
        }
        
    }
    
    
    exit(0);

}



## note borrowed code below from TransDecoder.predict and should consolidate.   ## //FIXME

sub parse_pfam_hits_file {
    my ($pfam_hits_file) = @_;

    my %has_pfam_hit;

    if (! -e $pfam_hits_file) {
        die "Error, cannot find pfam hits file: $pfam_hits_file";
    }

    print "PFAM output found and processing...\n";
    # capture those proteins having pfam hits
    open (my $fh, $pfam_hits_file) or die "Error, cannot open file: $pfam_hits_file";
    while (my $ln=<$fh>) {
        next if $ln=~/^\#/;
        my @x = split(/\s+/,$ln);
        next unless $x[3];  # domtbl
        my $orf_acc = $x[3];
        $has_pfam_hit{$orf_acc} = 1;
    }
    close $fh;


    return(%has_pfam_hit);
}

####
sub parse_blastp_hits_file {
    my ($blastp_file) = @_;

    unless (-e $blastp_file) {
        die "Error, cannot find file $blastp_file";
    }

    my %blastp_hits;

    open (my $fh, $blastp_file) or die "Error, cannot open file $blastp_file";
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $id = $x[0];

        $blastp_hits{$id} = 1;
    }
    close $fh;

    return(%blastp_hits);
}


