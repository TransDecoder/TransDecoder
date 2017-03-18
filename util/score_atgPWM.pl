#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Fasta_reader;
use Nuc_translator;

my $usage = <<__EOUSAGE__;

#####################################################################
#
#  --transcripts <string>     target transcripts fasta file
#
#  --atg_pwm <string>         start codon PWM file
#
#  --base_freqs <string>      base freqs dat file
#
#####################################################################


__EOUSAGE__

    ;


my $help_flag;
my $transcripts_file;
my $atg_pwm_file;
my $base_freqs_file;


&GetOptions ( 'h' => \$help_flag,
              'transcripts=s' => \$transcripts_file,
              'atg_pwm=s' => \$atg_pwm_file,
              'base_freqs=s' => \$base_freqs_file
    );

if ($help_flag) {
    die $usage;
}

unless ($transcripts_file && $atg_pwm_file && $base_freqs_file) {
    die $usage;
}



main: {

    my %pwm = &get_pwm($atg_pwm_file);

    
    my %base_freqs = &parse_base_freqs($base_freqs_file);

    my $fasta_reader = new Fasta_reader($transcripts_file);
    while (my $seq_obj = $fasta_reader->next()) {
        my $seq_acc = $seq_obj->get_accession();
        my $sequence = uc $seq_obj->get_sequence();

        &score_seq_using_pwm($seq_acc, $sequence, \%pwm, \%base_freqs);
                
        
    }

    exit(0);
}

sub score_seq_using_pwm {
    my ($seq_acc, $sequence, $pwm_href, $base_freqs_href) = @_;

    my $pwm_lend = $pwm_href->{'pwm_left'};
    my $pwm_rend = $pwm_href->{'pwm_right'};
    my $pwm_len = $pwm_lend + 3 + $pwm_rend;

    
    my @candidates;
    
    while ($sequence =~ /(ATG)/g) {
        
        my $pos_start = $-[0];
        push (@candidates, $pos_start);
    }
    
    my $seq_len = length($sequence);
    
    my @seq_chars = split(//, $sequence);
    
    
    foreach my $start_pos (@candidates) {
        
        my $motif_score = 0;
        
        my $begin = $start_pos - $pwm_lend;
        
        if ($begin < 0 || $begin + $pwm_len >= $seq_len) { next; }
        
        my $pwm_pos = 0;
        for (my $i = $begin; $i < $begin + $pwm_len; $i++) {
            
            if ($pwm_pos > $pwm_lend && $pwm_pos <= $pwm_lend + 3) { next; } # skip the ATG itself
            
            my $char = $seq_chars[$i];
            my $prob = $pwm_href->{pwm}->{$pwm_pos}->{$char};
            my $prob_rand = $base_freqs_href->{$char};
            
            my $loglikelihood = log($prob/$prob_rand);
            $motif_score += $loglikelihood;
            
            $pwm_pos++;
        }
        
        print "$seq_acc\t$start_pos\t$motif_score\n";
    }
    print "\n";
}
    

####
sub get_pwm {
    my ($pwm_file) = @_;

    open(my $fh, $pwm_file) or die "Error, cannot open file: $pwm_file";
    my $header = <$fh>;
    chomp $header;
    my @x = split(/\s+/, $header);
    my $pwm_left = $x[1];
    my $pwm_right = $x[3];

    my %pwm;
    $pwm{'pwm_left'} = $pwm_left;
    $pwm{'pwm_right'} = $pwm_right;

    my $seq_header = <$fh>;
    my $i = 0;
    while(<$fh>) {
        chomp;
        my ($relpos, $G, $A, $T, $C) = split(/\t/);
        $pwm{'pwm'}->{$i} = { G => $G,
                              A => $A,
                              T => $T,
                              C => $C };
        $i++;
    }

    return(%pwm);
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


    
