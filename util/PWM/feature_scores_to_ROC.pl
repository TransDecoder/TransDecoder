#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 feature.scores [num_ROC_threshold_intervals=10]\n\n";

my $feature_score_file = $ARGV[0] or die $usage;

my $num_intervals = $ARGV[1] || 10;

main: {

    my ($min_val, $max_val);

    my %cat_type_to_scores;

    print STDERR "-parsing scores\n";
    open(my $fh, $feature_score_file) or die $!;
    while (<$fh>) {
        chomp;
        my ($category, $pos_or_neg, $score) = split(/\t/);

        if ($score eq 'NA') { next; }
        
        push (@{$cat_type_to_scores{$category}}, [$pos_or_neg, $score]);
    
        if (! defined $min_val) {
            $min_val = $score;
            $max_val = $score;
        }
        else {
            if ($score < $min_val) {
                $min_val = $score;
            }
            if ($score > $max_val) {
                $max_val = $score;
            }
        }
        

    }
    close $fh;


    my $delta = ($max_val - $min_val) / $num_intervals;

    print join("\t", "cat", "thresh", "TP", "TN", "FP", "FN", "TPR", "FPR", "F1") . "\n";
    
    foreach my $cat_type (keys %cat_type_to_scores) {

        my @type_n_score = @{$cat_type_to_scores{$cat_type}};


        for (my $i = $min_val; $i < $max_val; $i += $delta) {

            my ($TP, $TN, $FP, $FN) = &score_cat_type(\@type_n_score, $i);


            my $TPR = $TP / ($TP + $FN);
            my $FPR = $FP / ($FP + $TN);
            my $F1 = (2 * $TP) / (2 * $TP + $FP + $FN);
            
            print join("\t", $cat_type, $i, $TP, $TN, $FP, $FN, $TPR, $FPR, $F1) . "\n";
        }
    }


    exit(0);
}


####
sub score_cat_type {
    my ($type_n_score_aref, $min_val) = @_;

    my ($TP, $TN, $FP, $FN) = (0,0,0,0);

    foreach my $type_n_score (@$type_n_score_aref) {
        my ($type, $score) = @$type_n_score;

        if ($score eq "NA") { next; }
        
        if ($type eq 'pos') {
            if ($score >= $min_val) {
                $TP++;
            }
            else {
                $FN++;
            }
        }
        elsif ($type eq 'neg') {
            if ($score >= $min_val) {
                $FP++;
            }
            else {
                $TN++;
            }
        }
        else {
            die "Error, don't understand type: $type";
        }
    }


    return($TP, $TN, $FP, $FN);
}

