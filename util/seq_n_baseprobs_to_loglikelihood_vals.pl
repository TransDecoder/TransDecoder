#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Fasta_reader;
use Nuc_translator;

## hexamer stats
my %framed_kmers;

my %background_base_probs;


my $usage = "usage: $0 targetCDSs base_probs.dat\n\n";

my $target_CDS = $ARGV[0] or die $usage;
my $base_probs_dat_file = $ARGV[1] or die $usage;
my $debug = $ARGV[2];

main: {

	&parse_targetCDSs($target_CDS);

	&parse_background($base_probs_dat_file);

	&report_logliklihood_ratios();
	
	exit(0);
}





####
sub report_logliklihood_ratios {


	print join("\t", "#framed_kmer", "kmer_count", "kminus1_prefix_count", "loglikelihood") . "\n";
    
	## Markov-chain based probabilities
    
	foreach my $framed_kmer (sort keys %framed_kmers) {
		my ($kmer, $frame) = split (/-/, $framed_kmer);
	
        if ($kmer =~ /[^GATC]/) {
            ## ignoring hexamers containing non-GATC bases
            next;
        }
        my $kmer_length = length($kmer);
        my $framed_kmer_count = $framed_kmers{"${kmer}-${frame}"};

        my $kminus1mer_count = undef; # set below
        my $kminus1mer_frame = $frame - 1;
        if ($kminus1mer_frame < 0) {
            $kminus1mer_frame = 2;
        }
        
        if ($kmer_length > 1) {
            my $kminus1mer = substr($kmer, 0, $kmer_length-1);
            $kminus1mer_count = $framed_kmers{"${kminus1mer}-${kminus1mer_frame}"} || 0;
        }
        else {
            $kminus1mer_count = $framed_kmers{"FRAME-${kminus1mer_frame}"} || 0;
        }
         
        
		my $markov_prob_framed = ($framed_kmer_count + 1) / ($kminus1mer_count + 4); # adding pseudocounts, 1 per kmer possibility
        
		my $last_base = substr($kmer, -1);
		my $background_prob = $background_base_probs{$last_base} or die "Error, no background probability set for base: $last_base of kmer $framed_kmer";;
        
        my $logliklihood = log($markov_prob_framed / $background_prob);
        
        print "$framed_kmer\t$framed_kmer_count\t$kminus1mer_count\t$logliklihood\n";
        
    }
    
	return;
}


	
####
sub parse_targetCDSs {
	my ($seqFile) = @_;

	my $fasta_reader = new Fasta_reader($seqFile);
	
	while (my $seq_obj = $fasta_reader->next()) {
		
		my $accession = $seq_obj->get_accession();
		print STDERR "\r     Target: processing $accession           " if $debug;
		
		my $sequence = uc $seq_obj->get_sequence();

		my $seq_len = length($sequence);

        for my $markov_order (0..5) {
            
            for (my $i = $markov_order; $i < $seq_len; $i++) {
                my $frame = $i % 3;

                ## avoid stop codons!
                if ($i == $seq_len - 2 - 1 # last triplet
                    &&
                    $frame == 0 # first position of a codon
                    ) {

                    my $codon = substr($sequence, $i, 3);
                    # stops: UAA, UAG, UGA
                    if ($codon =~ /^(TAA|TAG|TGA)$/) {
                        last;
                    }
                }
                
                if ($markov_order == 0) {
                    # include counts of number of framed positions, needed for Markov chain computes later.
                    $framed_kmers{"FRAME-${frame}"}++;
                }
                
                my $kmer = substr($sequence, $i-$markov_order, $markov_order+1);
                $framed_kmers{"${kmer}-${frame}"}++;
                
                
            }
        }
    }

	print "\r     CDS base frequency processing complete.           \n" if $debug;
    
	return;
}

#### 
sub parse_background {
	my ($base_probs_dat_file) = @_;

    open (my $fh, $base_probs_dat_file) or die "Error, cannot open file $base_probs_dat_file";
    while (<$fh>) {
        chomp;
        my ($base, $count, $ratio) = split(/\t/);
        $background_base_probs{$base} = $ratio;
    }
    close $fh;
        	
	return;
}

