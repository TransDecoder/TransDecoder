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

# Training needs this many start-site features (training ORFs with a complete ATG window).
# Below it, start refinement is skipped: start_codon_refinement.pl passes the ORFs through.
my $MIN_TRAINING_FEATURES = 50;

my $usage = <<__EOUSAGE__;


####################################################################
#
#  --transcripts <string>     target transcripts fasta file
#
#  --selected_orfs <string>   longest_orfs.cds.top_longest_5000.nr80
#
#  --out_prefix <string>      output prefix
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
my $out_prefix;

&GetOptions ( 'help|h' => \$help_flag,
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

main: {

    my $checkpoints_dir = "${out_prefix}_" . time() . "_checkpoints";
    if (! -d $checkpoints_dir) {
        mkdir($checkpoints_dir) or die $!;
    }
    
    my $pipeliner = new Pipeliner(-verbose => 2, -checkpoint_dir => $checkpoints_dir);
    
    my $cmd = "$utildir/build_atgPWM_+-.pl "
        . " --transcripts $transcripts_file "
        . " --selected_orfs $selected_orfs_file "
        . " --out_prefix $out_prefix"
        . " --pwm_left $pwm_left --pwm_right $pwm_right ";
    
    $pipeliner->add_commands(new Command($cmd, "built_init.ok"));
    $pipeliner->run();

    my $skip_file = "$out_prefix.SKIPPED";
    unlink($skip_file) if -e $skip_file;
    my $num_features = 0;
    open(my $fh, "$out_prefix.+.features") or die "Error, cannot open $out_prefix.+.features";
    while (<$fh>) { $num_features++ if /\S/; }
    close $fh;
    if ($num_features < $MIN_TRAINING_FEATURES) {
        my $msg = "Only $num_features start-site training features (minimum $MIN_TRAINING_FEATURES); "
            . "skipping start codon refinement. This happens when training ORFs start at or near the "
            . "transcript 5' end (e.g. transcripts lacking 5' UTRs).";
        print STDERR "WARNING: $msg\n";
        open(my $ofh, ">$skip_file") or die "Error, cannot write $skip_file";
        print $ofh "$msg\n";
        close $ofh;
        exit(0);
    }

    $cmd = "$utildir/feature_scoring.+-.pl "
        . " --features_plus $out_prefix.+.features "
        . " --features_minus $out_prefix.-.features "
        . " --atg_position $pwm_left "
        . " > $out_prefix.feature.scores";
    
    $pipeliner->add_commands(new Command($cmd, "score_features.ok"));

    $cmd = "$utildir/feature_scores_to_ROC.pl $out_prefix.feature.scores > $out_prefix.feature.scores.roc";
    $pipeliner->add_commands(new Command($cmd, "roc_features.ok"));

    $cmd = "$utildir/plot_ROC.Rscript $out_prefix.feature.scores.roc || :";
    $pipeliner->add_commands(new Command($cmd, "rocplot.ok"));
    
    $cmd = "$utildir/compute_AUC.pl $out_prefix.feature.scores.roc";
    $pipeliner->add_commands(new Command($cmd, "aucplot.ok"));

    $cmd = "$utildir/make_seqLogo.Rscript $out_prefix.+.pwm || :";
    $pipeliner->add_commands(new Command($cmd, "seqlogo.+.ok"));

    $cmd = "$utildir/make_seqLogo.Rscript $out_prefix.-.pwm || :";
    $pipeliner->add_commands(new Command($cmd, "seqlogo.-.ok"));

        
    #################
    ## motif enhancer

    
    $cmd = "$utildir/deplete_feature_noise.pl "
        . " --features_plus $out_prefix.+.features "
        . " --pwm_minus $out_prefix.-.pwm "
        . " --out_prefix $out_prefix.enhanced";
    
    $pipeliner->add_commands(new Command($cmd, "enhance.ok"));

    
    $cmd = "$utildir/feature_scoring.+-.pl "
        . " --features_plus $out_prefix.enhanced.+.features "
        . " --features_minus $out_prefix.-.features "
        . " --atg_position $pwm_left "
        . " > $out_prefix.enhanced.feature.scores";
    
    $pipeliner->add_commands(new Command($cmd, "enhanced_score_features.ok"));

    
    $cmd = "$utildir/feature_scores_to_ROC.pl $out_prefix.enhanced.feature.scores > $out_prefix.enhanced.feature.scores.roc";
    $pipeliner->add_commands(new Command($cmd, "enhanced_roc_features.ok"));

    $cmd = "$utildir/plot_ROC.Rscript $out_prefix.enhanced.feature.scores.roc || :";
    $pipeliner->add_commands(new Command($cmd, "enhanced_rocplot.ok"));
    
    $cmd = "$utildir/compute_AUC.pl $out_prefix.enhanced.feature.scores.roc";
    $pipeliner->add_commands(new Command($cmd, "enhanced_aucplot.ok"));
    
    $cmd = "$utildir/make_seqLogo.Rscript $out_prefix.enhanced.+.pwm || :";
    $pipeliner->add_commands(new Command($cmd, "seqlogo.enhanced.+.ok"));
        
    
    $pipeliner->run();

    exit(0);
    
}

    
