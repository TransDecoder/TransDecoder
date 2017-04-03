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

my $usage = <<__EOUSAGE__;

#####################################################################
#
#  --transcripts <string>     target transcripts fasta file
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
my $transcripts_file;
my $atg_pwm_file;
my $base_freqs_file;
my $atg_position = 10;

&GetOptions ( 'h' => \$help_flag,
              'transcripts=s' => \$transcripts_file,
              'atg_pwm=s' => \$atg_pwm_file,
              'base_freqs=s' => \$base_freqs_file,
              'atg_position=i' => \$atg_position,
    );

if ($help_flag) {
    die $usage;
}

unless ($transcripts_file && $atg_pwm_file && $base_freqs_file) {
    die $usage;
}



main: {

    my $pwm = new PWM();
    $pwm->load_pwm_from_file($atg_pwm_file);
    
    my %base_freqs = &parse_base_freqs($base_freqs_file);

    my $fasta_reader = new Fasta_reader($transcripts_file);
    while (my $seq_obj = $fasta_reader->next()) {
        my $seq_acc = $seq_obj->get_accession();
        my $sequence = uc $seq_obj->get_sequence();
        
        &score_seq_using_pwm($seq_acc, $sequence, $pwm, \%base_freqs);
        
    }

    exit(0);
}

sub score_seq_using_pwm {
    my ($seq_acc, $sequence, $pwm_obj, $base_freqs_href) = @_;

        
    my @candidates;
    
    while ($sequence =~ /(ATG)/g) {
        
        my $pos_start = $-[0];
        push (@candidates, $pos_start);
    }
    
    my $seq_len = length($sequence);
    my $pwm_length = $pwm_obj->get_pwm_length();

    my @mask = ($atg_position, $atg_position+1, $atg_position+2);
    my $target_prefix_length = $atg_position + 1 - 1;

    
    foreach my $start_pos (@candidates) {

        my $target_seq = substr($sequence, $start_pos - $target_prefix_length, $pwm_length);

        unless (length($target_seq) == $pwm_length) { next; } # need to skip it
        
        my $motif_score = $pwm_obj->score_pwm($target_seq, $base_freqs_href, 'mask' => \@mask);
        
        print "$seq_acc\t$start_pos\t$motif_score\n";
    }
    print "\n";

    return;
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


    
