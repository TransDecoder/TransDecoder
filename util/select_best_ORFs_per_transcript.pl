#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use GFF3_utils2;
use Carp;
use Data::Dumper;
use DelimParser;
use List::Util qw (max);
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
require "overlapping_nucs.ph";

srand(1234);

my $min_length_auto_accept = 1000000; # infinity effectively

my $usage = <<__EOUSAGE__;

################################################################################
#
# Required:
#
#  --gff3_file <string>          gff3 file w/ predicted ORFs
#
#  --cds_scores <string>         cds scores file
#
# Optional:
#
#  --min_length_auto_accept <int>  min length (nt) of orf to automatically accept
#
#  --blast_hits <string>         blast hits file
#
#  --pfam_hits <string>          pfam hits file
#
#  --single_best_orf          retain single longest orf that meets selection criteria
#
#  --verbose                  show orf rankings in stderr
#
#       Initial ORF selection requirement:
#
#            blast_hit |  pfam_hit | (frame_score[0] > 0 and frame_score[0] > max(frame_score[1], frame_score[2]) )
#
#       Subsequent ORF prioritization scheme:
#
#            homology_count, frame_score[0], orf_length
# 
################################################################################


__EOUSAGE__

    ;

my $gff3_file;
my $blast_hits_file;
my $pfam_hits_file;
my $cds_scores_file;
my $SINGLE_BEST_ORF_FLAG = 0;

my $MAX_PCT_OVERLAP = 10;
my $help_flag;
my $verbose_flag = 0;

&GetOptions( # required
             'gff3_file=s' => \$gff3_file,
             'cds_scores=s' => \$cds_scores_file,
             
             # optional

             'min_length_auto_accept=i' => \$min_length_auto_accept,
             
             'blast_hits=s' => \$blast_hits_file,
             'pfam_hits=s' => \$pfam_hits_file,
             
             'single_best_orf' => \$SINGLE_BEST_ORF_FLAG,

             'help|h' => \$help_flag,
             'verbose|v' => \$verbose_flag,
    );

if (@ARGV) {
    die "Error, dont understand params: @ARGV";
}

if ($help_flag) {
    die $usage;
}

unless ($gff3_file && $cds_scores_file) {
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

    my %cds_scores;
    if ($cds_scores_file) {
      %cds_scores = &parse_CDS_scores_file ($cds_scores_file);
    }
        
    my $gene_obj_indexer_href = {};
    
    my $asmbl_id_to_gene_list_href = &GFF3_utils2::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer_href);
    
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

            my $orf_length = $gene_obj_ref->get_CDS_length();
            my $cds_scores_aref = $cds_scores{$model_id};


            ## apply selection criteria
            unless ($homology_count 
                    ||
                    $orf_length >= $min_length_auto_accept
                    ||
                    ($cds_scores_aref->[0] > 0 &&  $cds_scores_aref->[0] > max($cds_scores_aref->[1], $cds_scores_aref->[2]))
                ) {
                # doesn't meet selection criteria. skipping it.
                next;
            }
            
            # orf accepted as preliminary candidate
            
            my $struct = { gene_obj => $gene_obj_ref,
                           length => $orf_length, 
                           homology_count => $homology_count,
                           cds_scores => $cds_scores_aref,
            };
            
            push (@gene_entries, $struct);
            
            
        }

        unless (@gene_entries) {
            # skip contig.
            next;
        }
        
        ## prioritize orfs for further screening.
        
        @gene_entries = sort {
            $b->{homology_count} <=> $a->{homology_count}
            ||
            $b->{cds_scores}->[0] <=> $a->{cds_scores}->[0]
            ||
            $b->{length} <=> $a->{length}
            
        } @gene_entries;


        if ($verbose_flag) {
            print STDERR "ORFs prioritized as:\n";
            foreach my $entry (@gene_entries) {
                print STDERR "\t" . join("\t", $entry->{gene_obj}->{Model_feat_name}, 
                                         "homology_count: " . $entry->{homology_count},
                                         "cds_score: " . $entry->{cds_scores}[0],
                                         "len: " . $entry->{length}) . "\n";
            }
        }
                
        if ($SINGLE_BEST_ORF_FLAG) {
            ## re-sort and take the longest one w/ homology info:
            @gene_entries = sort {
                $b->{homology_count} <=> $a->{homology_count}
                ||
                $b->{length} <=> $a->{length}
            } @gene_entries;
            
            @gene_entries = ($gene_entries[0]);

            if ($verbose_flag) {
                print STDERR "- ** selecting best orf, prioritized by homology then length:\n";
                my $entry = $gene_entries[0];
                print STDERR "\t" . join("\t", $entry->{gene_obj}->{Model_feat_name}, 
                                         "homology_count: " . $entry->{homology_count},
                                         "cds_score: " . $entry->{cds_scores}[0],
                                         "len: " . $entry->{length}) . "\n";
            }
        }
        else {
            @gene_entries = &remove_overlapping_orfs(@gene_entries);
        }
        
        foreach my $gene_entry (@gene_entries) {
            
            my $gene_obj = $gene_entry->{gene_obj};

            my $model_id = $gene_obj->{Model_feat_name};

            my $com_name = $gene_obj->{com_name};

            my $adj_comname = $com_name . ",score=" . $gene_entry->{cds_scores}[0];
            
            if (my $blast_hits = $blast_hits{$model_id}) {
                $adj_comname .= $blast_hits;
            }
            if (my $pfam_hits = $pfam_hits{$model_id}) {
                $adj_comname .= $pfam_hits;
            }


            $gene_obj->{com_name} =  $adj_comname;            
            
            print $gene_obj->to_GFF3_format(source => "transdecoder") . "\n";
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

    print STDERR "PFAM output found ($pfam_hits_file) and processing...\n";
    # capture those proteins having pfam hits
    open (my $fh, $pfam_hits_file) or die "Error, cannot open file: $pfam_hits_file";
    while (my $ln=<$fh>) {
        next if $ln=~/^\#/;
        my @x = split(/\s+/,$ln);

        # hmmscan format:
        #0       Fe-ADH_2
        #1       PF13685.5
        #2       250
        #3       CUFF.50.1.p1
        #4       -

        # hmmsearch format:
        #0       CUFF.50.1.p1
        #1       -
        #2       423
        #3       Fe-ADH_2
        #4       PF13685.5
        #5       250

        
        next unless $x[3];  # domtbl

        if ($x[1] =~ /^PF\d+/) {
            # hmmscan formatting 
            my $orf_acc = $x[3]; # CUFF.50.1.p1
            my $pfam_hit = $x[0]; # Fe-ADH_2
            my $pfam_acc = $x[1]; # PF13685.5
            my $domain_evalue = $x[12];
            $has_pfam_hit{$orf_acc} .= "," . join("|", $pfam_hit, $pfam_acc, $domain_evalue);
        }
        elsif ($x[4] =~ /^PF\d+/) {
            # hmmsearch formatting:
            my $orf_acc = $x[0]; # CUFF.50.1.p1
            my $pfam_hit = $x[3]; # Fe-ADH_2
            my $pfam_acc = $x[4]; # PF13685.5
            my $domain_evalue = $x[12];
            $has_pfam_hit{$orf_acc} .= "," . join("|", $pfam_hit, $pfam_acc, $domain_evalue);
        }
        else {
            confess "Error, cannot decipher pfam hit formatting: $ln ";
        }
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
    print STDERR "blastp output found ($blastp_file) and processing...\n";
    my %blastp_hits;

    open (my $fh, $blastp_file) or die "Error, cannot open file $blastp_file";
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $id = $x[0];
        my $hit_acc = $x[1];
        my $per_id = $x[2];
        my $Evalue = $x[10];
        
        $blastp_hits{$id} .= "," . join("|", $hit_acc, $per_id, $Evalue);
    }
    close $fh;

    return(%blastp_hits);
}

sub parse_CDS_scores_file {
    my ($cds_scores_file) = @_;
    
    unless (-e $cds_scores_file) {
        die "Error, cannot find file $cds_scores_file";
    }
    
    my %cds_scores;
    
    open (my $fh, $cds_scores_file) or die "Error, cannot open file $cds_scores_file";
    my $tabreader = new DelimParser::Reader($fh, "\t");
    while (my $row = $tabreader->get_row()) {
        
        my ($acc, @scores) = ($row->{'#acc'}, 
                              $row->{'score_1'}, $row->{'score_2'}, $row->{'score_3'},
                              $row->{'score_4'}, $row->{'score_5'}, $row->{'score_6'});

        # make sure we have score vals
        foreach my $score (@scores) {
            unless (defined $score) {
                confess "Error, missing a score entry for: " . Dumper($row);
            }
        }
                
        $cds_scores{$acc} = [@scores];
    }
    close $fh;
    
    return(%cds_scores);
}



####
sub remove_overlapping_orfs {
    my (@gene_entries) = @_;

    my @selected_entries;
    
    foreach my $gene_entry (@gene_entries) {
        if ($gene_entry->{homology_count} ||  ! &has_sufficient_overlap($gene_entry, \@selected_entries)) {
            push (@selected_entries, $gene_entry);
        }
    }
    
    return(@selected_entries);
}

####
sub has_sufficient_overlap {
    my ($gene_entry, $other_entries_aref) = @_;


    my $gene_obj = $gene_entry->{gene_obj};
    my ($lend, $rend) = sort {$a<=>$b} $gene_obj->get_model_span();
    my $gene_len = $rend - $lend + 1; 
    
    foreach my $other_entry (@$other_entries_aref) {
        my ($other_lend, $other_rend) = sort {$a<=>$b} $other_entry->{gene_obj}->get_model_span();
        
        my $other_len = $other_rend - $other_lend + 1;

        if (&coordsets_overlap([$lend, $rend], [$other_lend, $other_rend])) {

            my $overlap_len = &nucs_in_common($lend, $rend, $other_lend, $other_rend);

            if ($overlap_len / $gene_len * 100 > $MAX_PCT_OVERLAP
                ||
                $overlap_len / $other_len * 100 > $MAX_PCT_OVERLAP) {

                return(1);
            }
        }
    }

    return(0); # insufficient overlap detected
}

    
