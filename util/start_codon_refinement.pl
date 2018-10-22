#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use GFF3_utils2;
use Fasta_reader;
use Nuc_translator;
use Carp;
use Data::Dumper;
use List::Util qw (max);
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use PWM;
use File::Basename;


my $atg_pwm_pos = 20;
my $adj_dist = 30;
my $adj_pct = 15;

my $usage = <<__EOUSAGE__;

#####################################################################################
#
#  --transcripts <string>   transcripts.fasta file targeted by transdecoder
#
#  --gff3_file <string>     target gff3 file
#
#  optional:
#
#  --workdir <string>       TransDecoder working dir (default: (--transcripts val) + ".transdecoder_dir")
#
#  --adj_dist <int>         distance allowed for start adjustment (default: $adj_dist)
# 
#  --adj_pct <int>          pecentage of orf length for examining start adjustment (default: $adj_pct)
#
#  --atg_pos <int>          atg index position in pwm (default: $atg_pwm_pos)
#
#  --debug                  verbose
#
#######################################################################################



__EOUSAGE__

    ;


my $transcripts_file;
my $gff3_file;
my $DEBUG = 0;
my $workdir;

&GetOptions('transcripts=s' => \$transcripts_file,
            'gff3_file=s' => \$gff3_file,
            'adj_dist=i' => \$adj_dist,
            'atg_pos=i'=> \$atg_pwm_pos,
            'atg_pct=i' => \$adj_pct,
            'debug' => \$DEBUG,
            'workdir=s' => \$workdir,
    );


unless ($transcripts_file && $gff3_file) {
    die $usage;
}

if ($adj_pct > 30 || $adj_pct < 0) {
    die "Error, --adj_pct is out of range...  must be between 0 and 30 ";
}

main: {

    # should be in the output directory
    
    # get transdecoder start train files
    my $transdecoder_dir = $workdir;
    unless ($workdir) {
        $workdir = basename($transcripts_file) . ".transdecoder_dir";
    }
    unless (-d $transdecoder_dir) {
        die "Error, cannot locate transdecoder working directory as: $transdecoder_dir";
    }
    

    my $pwm_plus = "${transdecoder_dir}/start_refinement.enhanced.+.pwm";
    my $pwm_minus = "${transdecoder_dir}/start_refinement.-.pwm";
    
    my $pwm_plus_obj = new PWM();
    $pwm_plus_obj->load_pwm_from_file($pwm_plus);

    my $pwm_minus_obj = new PWM();
    $pwm_minus_obj->load_pwm_from_file($pwm_minus);
        
    my $roc_file = "${transdecoder_dir}/start_refinement.enhanced.feature.scores.roc";
    my $auc_file = "${transdecoder_dir}/start_refinement.enhanced.feature.scores.roc.auc";

    my ($pwm_range, $min_threshold) = &parse_range_and_thresholds($auc_file, $roc_file);
    print STDERR "-best pwm: $pwm_range, with score thresh: $min_threshold\n" if $DEBUG;
    
    my ($pwm_range_left, $pwm_range_right) = split(",", $pwm_range); # extent around the atg
    my $pwm_range_aref = [$atg_pwm_pos - $pwm_range_left -1, $atg_pwm_pos + 2 + $pwm_range_right -1]; # zero based
    
    my $start_scores_log_file = "${transdecoder_dir}/start_refinement.alt_start_scores";
    open(my $ofh_start_scores, ">$start_scores_log_file") or die "Error, cannot write to $start_scores_log_file";
    
    print STDERR "-reading transcripts: $transcripts_file\n" if $DEBUG;
    my $fasta_reader = new Fasta_reader($transcripts_file);

    my %seqs = $fasta_reader->retrieve_all_seqs_hash();

    print STDERR "-parsing orf annotations: $gff3_file\n" if $DEBUG;
    my $gene_obj_indexer_href = {};

    my $asmbl_id_to_gene_list_href = &GFF3_utils2::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer_href);

    my $num_starts_revised = 0;
    
    foreach my $transcript_acc (sort keys %$asmbl_id_to_gene_list_href) {

        print STDERR "-processing: $transcript_acc\n" if $DEBUG;
        
        my @gene_ids = @{$asmbl_id_to_gene_list_href->{$transcript_acc}};
        

        my $transcript_seq = uc $seqs{$transcript_acc};

        foreach my $gene_id (@gene_ids) {

            my $gene_obj = $gene_obj_indexer_href->{$gene_id};
            
            my $revised_start_flag = &refine_start_codon_position($transcript_acc, $gene_id,
                                                                  $gene_obj, $pwm_plus_obj, $pwm_minus_obj, 
                                                                  $pwm_range_aref, $min_threshold, $transcript_seq,
                                                                  $ofh_start_scores);

            if ($revised_start_flag) {
                $num_starts_revised++;

                ## update naming convention based on now having an updated start codon.
                if ($gene_obj->{com_name} =~ /internal/) {
                    $gene_obj->{com_name} =~ s/internal/3prime_partial/;
                }
                elsif ($gene_obj->{com_name} =~ /5prime_partial/) {
                    $gene_obj->{com_name} =~ s/5prime_partial/complete/;
                }
            } 
            print $gene_obj->to_GFF3_format(source => "transdecoder") . "\n";
            

        }
    }

    close $ofh_start_scores;
    
    print STDERR "-number of revised start positions: $num_starts_revised\n";

    exit(0);
    
}

####
sub parse_range_and_thresholds {
    my ($auc_file, $roc_file) = @_;

    ## get the most accurate search range:
    my $best_range = "";
    my $best_range_score = 0;
    {
        open(my $fh, $auc_file) or die "Error, cannot open file: $auc_file";
        while (<$fh>) {
            chomp;
            my ($range, $score) = split(/\t/);
            if ($score > $best_range_score) {
                $best_range_score = $score;
                $best_range = $range;
            }
        }
        close $fh;
    }

    # get most accurate threshold
    my $best_min_threshold = undef;
    my $best_F1 = 0;
    
    {
        open(my $fh, $roc_file) or die "Error, cannot open file $roc_file";
        my $header = <$fh>;
        while (<$fh>) {
            chomp;
            my @x = split(/\t/);
            my $cat = $x[0];

            unless ($cat eq $best_range) { next; }
            
            my $thresh = $x[1];
            my $F1 = $x[8];
            if ($F1 > $best_F1) {
                $best_F1 = $F1;
                $best_min_threshold = $thresh;
            }
        }
        close $fh;
    }
    
    return($best_range, $best_min_threshold);

}


####
sub refine_start_codon_position {
    my ($transcript_acc, $gene_id,
        $gene_obj, $pwm_plus_obj, $pwm_minus_obj, 
        $pwm_range_aref, $min_threshold, $transcript_seq, $ofh_start_scores) = @_;
    
    my $revised_start_flag = 0;
    
    my $orient = $gene_obj->get_orientation();
    
    my ($lend, $rend) = sort {$a<=>$b} $gene_obj->get_model_span();
    my $orf_len = $rend - $lend + 1;
    if ($orf_len % 3 != 0) {
        die "Error, $orf_len is not mod 3 " . $gene_obj->toString();
    }
    
    my $orig_start_pos = $lend;
    
    my $start_pos = $lend;
    if ($orient eq '-') {
        $transcript_seq = &reverse_complement($transcript_seq);
        $start_pos = length($transcript_seq) - $rend + 1;
    }


    # only work on 5' partials
    my $start_index = $start_pos - 1; # zero based

    if (substr($transcript_seq, $start_index, 3) eq "ATG") {
        return(0);
    }
    
    my @alt_starts;


    
    my $max_search_pos = max($start_index + $adj_dist, $start_index + int($adj_pct * $orf_len / 100));
    

    while ($transcript_seq =~ /(ATG)/g) {
        my $pos = $-[0];
        if ($pos > $max_search_pos) { last; } # too far
        if ($pos > $start_index 
            && 
            ($pos - $start_index) % 3  == 0) { # in frame start

            push (@alt_starts, $pos);
        }
    }

    unless (@alt_starts) {
        return($revised_start_flag);
    }
    
    my $pwm_len = $pwm_plus_obj->get_pwm_length();
    
        
    my $best_alt_start = undef;
    my $best_alt_start_score = undef;

    my @alt_start_scores;
    
    foreach my $alt_start (@alt_starts) {
        my $feature_seq_start = $alt_start - $atg_pwm_pos;
        if ($feature_seq_start > 0) {
            my $feature_seq = substr($transcript_seq, $feature_seq_start, $pwm_len);
            unless (length($feature_seq) == $pwm_len) { next; }
            
            my $alt_start_score = $pwm_plus_obj->score_plus_minus_pwm($feature_seq, $pwm_minus_obj,
                                                                      pwm_range => $pwm_range_aref);

            if ($alt_start_score eq "NA") { next; }
     
            $alt_start_score = sprintf("%.3f", $alt_start_score);
            
            my $short_feature_seq = &translate_sequence(substr($transcript_seq, $alt_start, 15), 1);
                        
            push (@alt_start_scores, "${alt_start}_${short_feature_seq}_${alt_start_score}");
            
            if ($alt_start_score >= $min_threshold

                &&
                (  (! defined $best_alt_start_score) || $alt_start_score > $best_alt_start_score) ) {
                
                $best_alt_start = $alt_start;
                $best_alt_start_score = $alt_start_score;
            }
        }
    }
    
    
    if ($best_alt_start) {
        $best_alt_start++; # make 1-based coord
        my $new_start = $best_alt_start;
        if ($orient eq '-') {
            $new_start = length($transcript_seq) - $best_alt_start + 1;
        }
        
        my ($exon_obj) = $gene_obj->get_exons();
        my $cds_obj = $exon_obj->get_CDS_obj();
        $cds_obj->{end5} = $new_start;
        
        $gene_obj->refine_gene_object();

        $revised_start_flag = 1;
        
        print STDERR "# refined start codon: $orig_start_pos -> $new_start\n" if $DEBUG;
        print "# refined start codon: $orig_start_pos -> $new_start (score: $best_alt_start_score)\n";
    }
    

    if (@alt_start_scores) {
        unshift(@alt_start_scores, $transcript_acc, $gene_id);
        print $ofh_start_scores join("\t", @alt_start_scores) . "\n";
    }
    
    return($revised_start_flag);
    
}

