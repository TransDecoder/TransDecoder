#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Fasta_reader;
use Nuc_translator;
use PWM;

my $pwm_left = 20;
my $pwm_right = 10;


my $usage = <<__EOUSAGE__;

#####################################################################
#
#  --transcripts <string>     target transcripts fasta file
#
#  --selected_orfs <string>   longest_orfs.cds.top_longest_5000.nr80
#
#  --out_prefix <string>      output prefix for pwm and feature sequences
#
## Optional
#
#  --pwm_left <int>           default: $pwm_left
#
#  --pwm_right <int>          default: $pwm_right  
#
#
#####################################################################


__EOUSAGE__

    ;


my $help_flag;
my $transcripts_file;
my $selected_orfs_file;
my $out_prefix = undef;


&GetOptions ( 'h' => \$help_flag, 
              'transcripts=s' => \$transcripts_file,
              'selected_orfs=s' => \$selected_orfs_file,
              'pwm_left=i' => \$pwm_left,
              'pwm_right=i' => \$pwm_right,
              'out_prefix=s' => \$out_prefix,
    
    );

if ($help_flag) {
    die $usage;
}

unless ($transcripts_file && $selected_orfs_file && $out_prefix) {
    die $usage;
}

my $pwm_length = $pwm_left + 3 + $pwm_right;

main: {

    my $fasta_reader = new Fasta_reader($transcripts_file);

    my %seqs = $fasta_reader->retrieve_all_seqs_hash();

    my @starts = &parse_starts($selected_orfs_file);



    my $pwm_plus = new PWM();
    my $pwm_minus = new PWM();
    
    
    my $features_plus_file = "$out_prefix.+.features";
    open (my $ofh_features_plus, ">$features_plus_file") or die "Error, cannot write to $features_plus_file";

    my $features_minus_file = "$out_prefix.-.features";
    open (my $ofh_features_minus, ">$features_minus_file") or die "Error, cannot write to $features_minus_file";
        
    
    foreach my $start_info (@starts) {
        my ($transcript_acc, $start_coord, $orient) = @$start_info;
        
        my $transcript_seq = uc $seqs{$transcript_acc};

        if ($orient eq '-') {
            # convert info to (+) strand
            $transcript_seq = &reverse_complement($transcript_seq);
            $start_coord = length($transcript_seq) - $start_coord + 1;
        }

        my $start_codon = substr($transcript_seq, $start_coord - 1, 3);
        if ($start_codon ne "ATG") { next; }
                
        
        my @trans_chars = split(//, uc $transcript_seq);

        my $feature_seq = &get_feature_seq(\@trans_chars, $start_coord - 1);

        if ($feature_seq) {
            
            print $ofh_features_plus join("\t", $feature_seq, $transcript_acc) . "\n";
            
            $pwm_plus->add_feature_seq_to_pwm($feature_seq);
        }
        
        ######################################
        ## contribute rest to the negative set.
        $transcript_seq = substr($transcript_seq, $start_coord + 1);
        
        &add_all_starts_to_pwm($pwm_minus, $transcript_seq, $ofh_features_minus);
        
    }

    close $ofh_features_plus;
    close $ofh_features_minus;
    
    # build and write plus model
    $pwm_plus->build_pwm();
    
    my $pwm_plus_file = "$out_prefix.+.pwm";
    $pwm_plus->write_pwm_file($pwm_plus_file);

    # build and write minus model
    $pwm_minus->build_pwm();

    my $pwm_minus_file = "$out_prefix.-.pwm";
    $pwm_minus->write_pwm_file($pwm_minus_file);

        
    exit(0);
    

}


####
sub get_feature_seq {
    my ($seq_aref, $start_coord) = @_;
    
    my $begin_pwm = $start_coord - $pwm_left;
    my $end_pwm = $begin_pwm + $pwm_length - 1;
        
    if ($begin_pwm < 0 || $end_pwm > $#$seq_aref) {
        # out of bounds on sequence
        return(undef);
    }
    
    
    my $feature_seq = "";
    for (my $i = $begin_pwm; $i <= $end_pwm; $i++) {
        my $char = $seq_aref->[$i];
        $feature_seq .= $char;
    }
    
    return($feature_seq);
}



####
sub add_all_starts_to_pwm {
    my ($pwm_obj, $sequence, $ofh) = @_;

    my @seqarray = split(//, $sequence);
    
    my @candidates;
    
    while ($sequence =~ /(ATG)/g) {
        
        my $pos_start = $-[0];
        push (@candidates, $pos_start);
    }

    foreach my $start_pos (@candidates) {
        my $feature_seq = &get_feature_seq(\@seqarray, $start_pos);
        
        if ($feature_seq) {
            $pwm_obj->add_feature_seq_to_pwm($feature_seq);
            print $ofh "$feature_seq\n";
        }
    }
    
    return;
}
    


####
sub parse_starts {
    my ($selected_orfs_file) = @_;

    my @starts;

    open(my $fh, $selected_orfs_file) or die "Error, cannot open file $selected_orfs_file";
    while (<$fh>) {
        chomp;
        if (/^>/) {
            if (/(\S+):(\d+)-(\d+)\(([+-])\)/) {
                my $transcript = $1;
                my $prime5 = $2;
                my $prime3= $3;
                my $orient = $4;

                
                push (@starts, [$transcript, $prime5, $orient]);
            }
            else {
                die "Error, couldn't parse transcript coordinate info from $_";
            }

            
        }
        
    }
    close $fh;
    
    return(@starts);
}


####
sub convert_counts_to_prob_vals {
    my (@vals) = @_;

    my $index = shift @vals;
    
    my $sum = 0;
    foreach my $val (@vals) {
        $sum += $val;
    }

    my @probs = ($index);
    foreach my $val (@vals) {
        my $prob = sprintf("%.3f", $val / $sum);
        
        push (@probs, $prob);
    }

    return(@probs);
}
