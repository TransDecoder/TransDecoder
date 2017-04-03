#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use PWM;

my $usage = <<__EOUSAGE__;

#####################################################################
#
#  --features <string>        features file
#
#  --pwm_out <string>         pwm output file
#
#####################################################################


__EOUSAGE__

    ;


my $help_flag;

my $features_file;
my $pwm_output_file;

&GetOptions ( 'h' => \$help_flag,
              'pwm_out=s' => \$pwm_output_file,
              'features=s' => \$features_file,
    );

if ($help_flag) {
    die $usage;
}

unless ($pwm_output_file && $features_file) {
    die $usage;
}


main: {

    my $pwm = new PWM();

    open(my $fh, $features_file) or die "Error, cannot open file $features_file";
    while (<$fh>) {
        chomp;
        my ($feature_seq, @rest) = split(/\s+/);
        $pwm->add_feature_seq_to_pwm($feature_seq);
    }
    close $fh;

    $pwm->build_pwm();

    $pwm->write_pwm_file($pwm_output_file);
        

    exit(0);
}

