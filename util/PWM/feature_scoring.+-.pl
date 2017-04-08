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
use Process_cmd;


my $num_rounds = 5;
my $fraction_train = 0.75;
my $atg_position = 20;
my $max_feature_select = 1000;

my $usage = <<__EOUSAGE__;


#####################################################################
#
#  --features_plus <string>        features plus file
#
#  --features_minus <string>       features minus file
#
#  optional:
#
#  --atg_position <int>       position of ATG in the PWM (default: $atg_position)  (zero-coord sys)
#
#  --num_rounds <int>         number of rounds to iterate training / evaluations (default: $num_rounds)
#
#  --fraction_train <float>   fraction of features to use for training (rest used for eval) (default: $fraction_train)
#
#  --max_feature_select <int>  max number of features to test on (default: $max_feature_select)
#
#####################################################################


__EOUSAGE__

    ;


my $help_flag;
my $features_plus_file;
my $features_minus_file;
my $atg_pwm_minus_file;


&GetOptions ( 'h' => \$help_flag,
              'features_plus=s' => \$features_plus_file,
              'features_minus=s' => \$features_minus_file,
              
              'atg_position=i' => \$atg_position,
              'num_rounds=i' => \$num_rounds,
              'fraction_train=f' => \$fraction_train,
    );

if ($help_flag) {
    die $usage;
}

unless ($features_plus_file && $features_minus_file) {
    die $usage;
}

unless ($fraction_train > 0 && $fraction_train < 1) {
    die "Error, --fraction_train must be between 0 and 1 ";
}

main: {


    my @features_plus = &parse_features_file($features_plus_file);
    my @features_minus = &parse_features_file($features_minus_file);
        
    my $pwm_length = length($features_plus[0]);
    my $pwm_upstream_max = $atg_position; # positions 0..$atg_position-1
    my $pwm_downstream_max = $pwm_length - ($atg_position +1 + 2); # positions available after the ATG

    my @up_down_combos;
    for (my $i = 1; $i <= $pwm_upstream_max; $i++) {
        for (my $j = 1; $j <= $pwm_downstream_max; $j++) {
            push (@up_down_combos, "$i,$j");
        }
    }
        

    my %motif_range_to_vals;

    my $counter = 0;

    my @feature_seqs;


    ## do 5 rounds
    use List::Util qw(shuffle);
    my $num_features = scalar(@feature_seqs);
    for my $round (1..$num_rounds) {

        print STDERR "-round: $round\n";
        
        my ($plus_train_seqs_aref, $plus_test_seqs_aref) = &sample_features(@features_plus);
        my ($minus_train_seqs_aref, $minus_test_seqs_aref) = &sample_features(@features_minus);
        
        my $pwm_plus = &build_pwm(@$plus_train_seqs_aref);
        my $pwm_minus = &build_pwm(@$minus_train_seqs_aref);
                
        &score_features($plus_test_seqs_aref, $pwm_plus, $pwm_minus, \@up_down_combos, 'pos');
        &score_features($minus_test_seqs_aref, $pwm_plus, $pwm_minus, \@up_down_combos, 'neg');
        
    }



=remove
    
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

=cut    
    
    exit(0);
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

####
sub parse_features_file {
    my ($features_file) = @_;

    my @feature_seqs;
    open(my $fh, $features_file) or die "Error, cannot open file: $features_file";
    while (<$fh>) {
        chomp;
        my ($feature_seq, @rest) = split(/\s+/);
        push (@feature_seqs, $feature_seq);
    }
    close $fh;
    
    return(@feature_seqs);
}
    

####
sub sample_features {
    my @feature_seqs = @_;
    
    @feature_seqs = shuffle(@feature_seqs);
    my $num_features = scalar(@feature_seqs);
    
    my $num_features_train = int($fraction_train * $num_features);
    my @train_features = @feature_seqs[0..$num_features_train-1];
    my @test_features = @feature_seqs[$num_features_train..$#feature_seqs]; # rest of them

    return(\@train_features, \@test_features);
}


####
sub score_features {
    my ($features_aref, $pwm_plus, $pwm_minus, $up_down_combos_aref, $feature_set_type) = @_;
    
    foreach my $up_down (@$up_down_combos_aref) {
        #print STDERR "\t$up_down\n";
        my ($up, $down) = split(/,/, $up_down);
        my $range_left = $atg_position - $up;
        my $range_down = $atg_position + 2 + $down; # the -1's are to convert from 1-based to 0-based coord sys
        
        my $counter = 0;
        foreach my $feature_seq (@$features_aref) {

            $counter++;
            if ($counter > $max_feature_select) { last; }
            
            my $score = $pwm_plus->score_plus_minus_pwm($feature_seq, $pwm_minus, 
                                                        pwm_range => [$range_left, $range_down]);
            
            my $local_pwm_len = $down + $up;
            
            if ($score ne "NA") {
                $score = sprintf("%.3f", $score/$local_pwm_len);  # normalizing for length of pwm
            } 
            print join("\t", $up_down, $feature_set_type, $score) . "\n";
        }
                
    }
    
    return;
}
