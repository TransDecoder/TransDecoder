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

my $usage = <<__EOUSAGE__;

#####################################################################
#
#  --features <string>        features file from initial pwm construction
#
#  --atg_pwm <string>         start codon PWM file
#
#  --base_freqs <string>      base freqs dat file
#
#  --atg_position <int>       position of ATG in the PWM (default: 10)
#
#####################################################################


__EOUSAGE__

    ;


my $help_flag;
my $features_file;
my $atg_pwm_file;
my $base_freqs_file;
my $atg_position = 10;

&GetOptions ( 'h' => \$help_flag,
              'features=s' => \$features_file,
              'atg_pwm=s' => \$atg_pwm_file,
              'base_freqs=s' => \$base_freqs_file,
              'atg_position=i' => \$atg_position,
    );

if ($help_flag) {
    die $usage;
}

unless ($features_file && $atg_pwm_file && $base_freqs_file) {
    die $usage;
}



main: {

    my $pwm = new PWM();
    $pwm->load_pwm_from_file($atg_pwm_file);

    my $pwm_length = $pwm->get_pwm_length();
    my $pwm_upstream_max = $atg_position + 1 - 1;
    my $pwm_downstream_max = $pwm_length - ($atg_position +1 + 2);

    my @up_down_combos;
    for (my $i = 1; $i <= $pwm_upstream_max; $i++) {
        for (my $j = 1; $j <= $pwm_downstream_max; $j++) {
            push (@up_down_combos, "$i,$j");
        }
    }
        
    my @mask = ($atg_position, $atg_position+1, $atg_position+2);
    
    my %base_freqs = &parse_base_freqs($base_freqs_file);
    
    my %motif_range_to_vals;
    
    open(my $fh, $features_file) or die "Error, cannot open file $features_file";
    open(my $ofh, ">$features_file.pwm_range_scores") or die "Error, cannot write to $features_file.pwm_range_scores";

    print $ofh join("\t", "", @up_down_combos) . "\n";
    
    my $counter = 0;
    while (<$fh>) {
        chomp;
        my ($feature_seq, $trans_name) = split(/\s+/);
        
        my @row = ("$feature_seq." . ++$counter);

        foreach my $up_down (@up_down_combos) {
            my ($up, $down) = split(/,/, $up_down);
            my $range_left = $atg_position - $up;
            my $range_down = $atg_position + 2 + $down; # the -1's are to convert from 1-based to 0-based coord sys
            my $score = $pwm->score_pwm($feature_seq, \%base_freqs, 
                                        mask => \@mask,
                                        pwm_range => [$range_left, $range_down]);

            if ($score ne "NA") {
                $score = sprintf("%.3f", $score);
            } 
            push (@row, $score); 
        }
        print $ofh join("\t", @row) . "\n";
        
    }
    close $fh;
    close $ofh;

    ## examine distribution:
    {
        my $Rscript = "$features_file.pwm_range_scores.Rscript";
        open (my $ofh, ">$Rscript") or die "Error, cannot write to $Rscript";
        print $ofh "data = read.table(\"$features_file.pwm_range_scores\")\n"
                 . "pdf(\"$features_file.pwm_range_scores.boxplot.pdf\")\n"
                 . "boxplot(data, outline=F, las=2)\n"
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


    
