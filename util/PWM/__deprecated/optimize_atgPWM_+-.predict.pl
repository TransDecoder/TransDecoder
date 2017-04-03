#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Fasta_reader;
use Nuc_translator;
use PWM;
use Process_cmd;


my $num_rounds = 5;
my $fraction_train = 0.75;
my $atg_position = 10;


my $usage = <<__EOUSAGE__;

#####################################################################
#
#  --features <string>        features file from initial pwm construction
#
#  --atg_pwm_minus <string>   negative start codon PWM model
#
#  optional:
#
#  --atg_position <int>       position of ATG in the PWM (default: $atg_position)  (zero-coord sys)
#
#  --num_rounds <int>         number of rounds to iterate training / evaluations (default: $num_rounds)
#
#  --fraction_train <float>   fraction of features to use for training (rest used for eval) (default: $fraction_train)
#
#####################################################################


__EOUSAGE__

    ;


my $help_flag;
my $features_file;
my $atg_pwm_minus_file;



&GetOptions ( 'h' => \$help_flag,
              'features=s' => \$features_file,
              'atg_pwm_minus=s' => \$atg_pwm_minus_file,
              'atg_position=i' => \$atg_position,
              'num_rounds=i' => \$num_rounds,
              'fraction_train=f' => \$fraction_train,
    );

if ($help_flag) {
    die $usage;
}

unless ($features_file && $features_file && $atg_pwm_minus_file) {
    die $usage;
}

unless ($fraction_train > 0 && $fraction_train < 1) {
    die "Error, --fraction_train must be between 0 and 1 ";
}

main: {

    my $pwm_minus = new PWM();
    $pwm_minus->load_pwm_from_file($atg_pwm_minus_file);
    
    my $pwm_length = $pwm_minus->get_pwm_length();
    my $pwm_upstream_max = $atg_position; # positions 0..$atg_position-1
    my $pwm_downstream_max = $pwm_length - ($atg_position +1 + 2); # positions available after the ATG

    my @up_down_combos;
    for (my $i = 1; $i <= $pwm_upstream_max; $i++) {
        for (my $j = 1; $j <= $pwm_downstream_max; $j++) {
            push (@up_down_combos, "$i,$j");
        }
    }
        

    my %motif_range_to_vals;
    
    open(my $fh, $features_file) or die "Error, cannot open file $features_file";
    open(my $ofh, ">$features_file.+-.pwm_range_scores") or die "Error, cannot write to $features_file.+-.pwm_range_scores";
    
    print $ofh join("\t", "", @up_down_combos) . "\n";
    
    my $counter = 0;

    my @feature_seqs;

    {
        while (<$fh>) {
            chomp;
            my ($feature_seq, $trans_name) = split(/\s+/);
            push (@feature_seqs, $feature_seq);
        }
        close $fh;
    }

    ## do 5 rounds
    use List::Util qw(shuffle);
    my $num_features = scalar(@feature_seqs);
    for my $round (1..$num_rounds) {

        # shuffle sequences, select 75% for training and the other half for testing
        my @features_to_use = shuffle(@feature_seqs);

        my $num_features_train = int($fraction_train * $num_features);
        my @train_features = @feature_seqs[0..$num_features_train-1];
        my @test_features = @feature_seqs[$num_features_train..$#feature_seqs]; # rest of them

        my $pwm_plus = &build_pwm(@train_features);
        
        foreach my $feature_seq (@test_features) {
      
            my @row = ("$feature_seq.$round." . ++$counter);
            
            foreach my $up_down (@up_down_combos) {
                my ($up, $down) = split(/,/, $up_down);
                my $range_left = $atg_position - $up;
                my $range_down = $atg_position + 2 + $down; # the -1's are to convert from 1-based to 0-based coord sys
                my $score = $pwm_plus->score_plus_minus_pwm($feature_seq, $pwm_minus, 
                                                            pwm_range => [$range_left, $range_down]);

                my $local_pwm_len = $down + $up;
                
                if ($score ne "NA") {
                    $score = sprintf("%.3f", $score/$local_pwm_len);  # normalizing for length of pwm
                } 
                push (@row, $score); 
            }
            print $ofh join("\t", @row) . "\n";
            
        }
    }
    close $fh;
    close $ofh;
    
    ## examine distribution:
    {
        my $Rscript = "$features_file.+-.pwm_range_scores.Rscript";
        open (my $ofh, ">$Rscript") or die "Error, cannot write to $Rscript";
        print $ofh "data = read.table(\"$features_file.+-.pwm_range_scores\")\n"
                 . "pdf(\"$features_file.+-.pwm_range_scores.boxplot.pdf\")\n"
                 . "boxplot(data, outline=F, las=2, cex.axis=0.4, main='log_odds/pwm_length, order by comparison')\n"
                 . "abline(h=0, col='blue')\n"
                 . "cd = apply(data, 2, median)\n"
                 . "boxplot(data[,rev(order(cd))], outline=F, las=2, cex.axis=0.4, main='log_odds/pwm_length, order by median desc')\n"
                 . "abline(h=0, col='blue')\n"
                 . "dev.off()\n";

        close $ofh;

        &process_cmd("Rscript --vanilla $Rscript");
    }
    
    
    exit(0);
}


####
sub parse_base_freqs {
    my ($base_freqs_file) = @_;

    my %base_freqs;
    
    open(my $fh, $base_freqs_file) or die "Error, cannot open file $base_freqs_file";
    while (<$fh>) {
        chomp;
        my ($base, $count, $relfreq) = split(/\t/);
        
        $base_freqs{$base} = $relfreq;
    }
    close $fh;

    return(%base_freqs);
}


####
sub build_pwm {
    my (@features) = @_;

    my $pwm = new PWM();

    foreach my $feature (@features) {
        $pwm->add_feature_seq_to_pwm($feature);
    }

    $pwm->build_pwm();

    return($pwm);
}


