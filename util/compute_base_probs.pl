#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Fasta_reader;
use Nuc_translator;

my $usage = "usage: $0 transcripts_file [top_strand_only]\n\n";

my $transcripts_file = $ARGV[0] or die $usage;
my $top_strand_only_flag = $ARGV[1] || 0;

main: {

    my %base_counter;

    my $fasta_reader = new Fasta_reader($transcripts_file);
    while (my $seq_obj = $fasta_reader->next()) {

        my $sequence = uc $seq_obj->get_sequence();

        &count_bases($sequence, \%base_counter);

        unless ($top_strand_only_flag) {
            $sequence = &reverse_complement($sequence);
            &count_bases($sequence, \%base_counter);
        }
        
    }

    
    my $sum = 0;
    foreach my $count (values %base_counter) {
        $sum += $count;
    }

    foreach my $base (sort keys %base_counter) {
        
        my $count = $base_counter{$base};

        my $ratio = $count/$sum;
        
        print join("\t", $base, $count, sprintf("%.3f", $ratio)) . "\n";
    }
    

    exit(0);
}

####
sub count_bases {
    my ($sequence, $base_counter_href) = @_;

    my @chars = split(//, $sequence);
    
    foreach my $char (@chars) {
        if ($char =~ /[GATC]/) {
            $base_counter_href->{$char}++;
        }
    }
    
    return;
}

