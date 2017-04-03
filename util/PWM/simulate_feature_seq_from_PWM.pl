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
#  --pwm <string>             pwm file
#
#  --num_features <int>       number of features to simulate
#
#####################################################################


__EOUSAGE__

    ;


my $help_flag;
my $pwm_file;
my $num_features;

&GetOptions ( 'h' => \$help_flag,
              'pwm=s' => \$pwm_file,
              'num_features=i' => \$num_features,
    );

if ($help_flag) {
    die $usage;
}

unless ($pwm_file && $num_features) {
    die $usage;
}


main: {

    my $pwm = new PWM();
    $pwm->load_pwm_from_file($pwm_file);

    for (1..$num_features) {
        my $feature_seq = $pwm->simulate_feature();
        print "$feature_seq\n";
    }

    exit(0);
}

