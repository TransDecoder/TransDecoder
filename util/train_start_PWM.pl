#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use FindBin;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use lib ("$FindBin::Bin/../PerlLib");
use Pipeliner;

my $utildir = "$FindBin::Bin/PWM";

my $help_flag;


my $pwm_left = 20;
my $pwm_right = 10;

my $usage = <<__EOUSAGE__;


####################################################################
#
#  --transcripts <string>     target transcripts fasta file
#
#  --selected_orfs <string>   longest_orfs.cds.top_longest_5000.nr80
#
#  Optional:
#
#  --pwm_left <int>           default: $pwm_left
#
#  --pwm_right <int>          default: $pwm_right
#
#
#####################################################################

__EOUSAGE__

    ;


my $transcripts_file;
my $selected_orfs_file;

&GetOptions ( 'h' => \$help_flag,
              'transcripts=s' => \$transcripts_file,
              'selected_orfs=s' => \$selected_orfs_file,
              'pwm_left=i' => \$pwm_left,
              'pwm_right=i' => \$pwm_right,
    );



unless ($transcripts_file && $selected_orfs_file) {
    die $usage;
}

main: {

    my $checkpoints_dir = "__checkpoints";
    if (! -d $checkpoints_dir) {
        mkdir($checkpoints_dir) or die $!;
    }
    
    my $pipeliner = new Pipeliner(-verbose => 2);

    my $cmd = "$utildir/build_atgPWM_+-.pl "
        . " --transcripts $transcripts_file "
        . " --selected_orfs $selected_orfs_file "
        . " --out_prefix atg"
        . " --pwm_left $pwm_left --pwm_right $pwm_right ";
    
    $pipeliner->add_commands(new Command($cmd, "$checkpoints_dir/built_init.ok"));

    $cmd = "$utildir/feature_scoring.+-.pl "
        . " --features_plus atg.+.features "
        . " --features_minus atg.-.features "
        . " --atg_position $pwm_left "
        . " > atg.feature.scores";
    
    $pipeliner->add_commands(new Command($cmd, "$checkpoints_dir/score_features.ok"));

    $cmd = "$utildir/feature_scores_to_ROC.pl atg.feature.scores > atg.feature.scores.roc";
    $pipeliner->add_commands(new Command($cmd, "$checkpoints_dir/roc_features.ok"));

    $cmd = "$utildir/plot_ROC.Rscript atg.feature.scores.roc";
    $pipeliner->add_commands(new Command($cmd, "$checkpoints_dir/rocplot.ok"));
    
    $cmd = "$utildir/compute_AUC.pl atg.feature.scores.roc";
    $pipeliner->add_commands(new Command($cmd, "$checkpoints_dir/aucplot.ok"));

    
    
    ## motif enhancer
    $cmd = "$utildir/deplete_feature_noise.pl "
        . " --features_plus atg.+.features "
        . " --pwm_minus atg.-.pwm "
        . " --out_prefix enhanced";
    
    $pipeliner->add_commands(new Command($cmd, "$checkpoints_dir/enhance.ok"));

    
    $cmd = "$utildir/feature_scoring.+-.pl "
        . " --features_plus enhanced.+.features "
        . " --features_minus atg.-.features "
        . " --atg_position $pwm_left "
        . " > enhanced.feature.scores";
    
    $pipeliner->add_commands(new Command($cmd, "$checkpoints_dir/enhanced_score_features.ok"));

    
    $cmd = "$utildir/feature_scores_to_ROC.pl enhanced.feature.scores > enhanced.feature.scores.roc";
    $pipeliner->add_commands(new Command($cmd, "$checkpoints_dir/enhanced_roc_features.ok"));

    $cmd = "$utildir/plot_ROC.Rscript enhanced.feature.scores.roc";
    $pipeliner->add_commands(new Command($cmd, "$checkpoints_dir/enhanced_rocplot.ok"));
    
    $cmd = "$utildir/compute_AUC.pl enhanced.feature.scores.roc";
    $pipeliner->add_commands(new Command($cmd, "$checkpoints_dir/enhanced_aucplot.ok"));
    

    
    $pipeliner->run();

    exit(0);
    
}

    
