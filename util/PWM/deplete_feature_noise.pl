#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use PWM;
use List::Util qw(shuffle);

my $pct_features_retain = 30;


my $usage = <<__EOUSAGE__;


############################################################################
#
#  --features_plus <string>      features plus file
#
#  --pwm_minus <string>          pwm minus file
#
#  --out_prefix <string>         output prefix (for pwm and feature files)
#
#  optional:
#
#  --pct_features_retain <int>   percent of features to retain in optimal set
#                                (default: $pct_features_retain)
#
##############################################################################


__EOUSAGE__

    ;


my $help_flag;

my $features_plus_file;
my $pwm_minus_file;
my $out_prefix;

&GetOptions ( 'h' => \$help_flag,
              'features_plus=s' => \$features_plus_file,
              'pct_features_retain=i' => \$pct_features_retain,
              'pwm_minus=s' => \$pwm_minus_file,
              'out_prefix=s' => \$out_prefix,
    );


if ($help_flag) {
    die $usage;
}


unless ($features_plus_file && $pwm_minus_file && $out_prefix) {
    die $usage;
}

if ($pct_features_retain < 1 || $pct_features_retain > 99) {
    die "Error, need to set --pct_features_retain between 1 and 99 percent, use an integer value.";
}


main: {

    my $pwm_minus = new PWM();
    $pwm_minus->load_pwm_from_file($pwm_minus_file);
    
    
    my @features = &parse_feature_seqs($features_plus_file);
    
    @features = shuffle(@features);
    
    
    my $num_features = scalar(@features);
    my $num_incorporate = int($num_features * $pct_features_retain/100);

    print STDERR "num features: $num_features\tnum_incorporate: $num_incorporate\n";


    my $pwm_plus = new PWM();

    my @init_features;
    for (1..$num_incorporate) {
        my $feature_seq = pop @features;
        push (@init_features, $feature_seq);
        $pwm_plus->add_feature_seq_to_pwm($feature_seq);
    }
    $pwm_plus->build_pwm();
    
    
    ## score features in the initial motif
    my @scored_features;
    foreach my $feature (@init_features) {
        my $score = $pwm_plus->score_plus_minus_pwm($feature, $pwm_minus);
        if ($score ne "NA") {
            push (@scored_features, 
                  { score => $score,
                    seq => $feature,
                  } );
        }
    }
    
    # sort by score ascendingly
    @scored_features = sort {$a->{score} <=> $b->{score}} @scored_features;
    
    # examine rest of the features
    # if one scores better than an existing entry, replace it.

    my $num_feature_swaps = 0;
    
    for my $feature (@features) {
        my $score = $pwm_plus->score_plus_minus_pwm($feature, $pwm_minus);
        
        my $worst_score = $scored_features[0]->{score};

        if ($score ne 'NA' && $score > $worst_score) {
            
            # purge the worst feature
            my $purge_feature = shift @scored_features;
                        
            # add better feature
            push(@scored_features, { score => $score,
                                     seq => $feature,
                 } );

    
            # adjust the pwm
            $pwm_plus->remove_feature_seq_from_pwm($purge_feature->{seq});
            $pwm_plus->add_feature_seq_to_pwm($feature);
            $pwm_plus->build_pwm();
        

            print STDERR "-feature swap of score: $score instead of $worst_score\n";
            $num_feature_swaps++;
            
            ## rescore all currently integrated features:
            foreach my $scored_feature (@scored_features) {
                $scored_feature->{score} = $pwm_plus->score_plus_minus_pwm($scored_feature->{seq}, $pwm_minus);
            }
            
            @scored_features = sort {$a->{score}<=> $b->{score}} @scored_features;

        }
    }

    print STDERR "-num feature swaps: $num_feature_swaps\n";

    ## prune non-pos scoring features from model
    my $num_features_removed = 0;
    foreach my $scored_feature (@scored_features) {
        if ($scored_feature->{score} <= 0) {
            $pwm_plus->remove_feature_seq_from_pwm($scored_feature->{seq});
            $num_features_removed++;
        }
    }

    if ($num_features_removed > 0) {
        $pwm_plus->build_pwm();
        print STDERR "-removed $num_features_removed / " . scalar(@scored_features) . " from PWM based on low scores\n";
    }
    
    my $out_pwm_file = "$out_prefix.+.pwm";
    $pwm_plus->write_pwm_file($out_pwm_file);
    
    my $out_features_file = "$out_prefix.+.features";
    {
        open(my $ofh, ">$out_features_file") or die "Error, cannot write to file: $out_features_file";
        foreach my $feature (@scored_features) {
            print $ofh $feature->{seq} . "\n";
        }
    }
        
                
    exit(0);
}

####
sub parse_feature_seqs {
    my ($file) = @_;

    my @features;
    
    open(my $fh, $file) or die "Error, cannot open file: $file";
    while (<$fh>) {
        chomp;
        my ($feature_seq, @rest) = split(/\s+/);
        
        push (@features, $feature_seq);
    }
    close $fh;
    
    return (@features);
}
