#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Fasta_reader;
use Nuc_translator;
use List::Util qw(min);

my $usage = "usage: $0 CDS hexamerScores\n\n";


my $MAX_REPEAT_COUNT = 3;

my $cds_file = $ARGV[0] or die $usage;
my $kmer_scores_file = $ARGV[1] or die $usage;


my %scores = &parse_kmer_scores($kmer_scores_file);



main: {

    print join("\t", "#acc", "Markov_order", "seq_length", 
               "score_1", "score_2", "score_3", "score_4", "score_5", "score_6") . "\n";
    
	my $fasta_reader = new Fasta_reader($cds_file);
	while (my $seq_obj = $fasta_reader->next()) {
		&score_seq($seq_obj);
    }
	
	exit(0);

}

####
sub score_seq {
	my ($seq_obj) = @_;
        
    my $accession = $seq_obj->get_accession();
    my $sequence = uc $seq_obj->get_sequence();
    
    for my $markov_order (5) {
        
		my $score1 = &score_CDS_via_Markov($sequence, $markov_order);				
		my $score2 = &score_CDS_via_Markov(substr($sequence, 1), $markov_order);
		my $score3 = &score_CDS_via_Markov(substr($sequence, 2), $markov_order);

		my $rev_seq = &reverse_complement($sequence, $markov_order);
		
		my $score4 = &score_CDS_via_Markov($rev_seq, $markov_order);
		my $score5 = &score_CDS_via_Markov(substr($rev_seq, 1), $markov_order);
		my $score6 = &score_CDS_via_Markov(substr($rev_seq, 2), $markov_order);
		
		printf("$accession\t$markov_order\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",
			   length($sequence), 
               $score1, $score2, $score3,
			   $score4, $score5, $score6);
                
    }
    
    return;
}


####
sub score_CDS_via_Markov {
    my ($sequence, $markov_order) = @_;
    
	my $seq_length = length($sequence);
	
	if ($seq_length < $markov_order + 1) {
		return(0);
	}

    my %repeat_counter;
    
    my $score = 0;
    for (my $i = 0; $i < $seq_length; $i++) {
        my $frame = $i % 3;
        my $markov_use = min($i, $markov_order);
        my $kmer = substr($sequence, $i-$markov_use, $markov_use + 1);

        ## avoid stop codons!
        if ($i == $seq_length - 2 - 1 # last triplet
            &&
            $frame == 0 # first position of a codon
            ) {
            
            my $codon = substr($sequence, $i, 3);
            # stops: UAA, UAG, UGA
            if ($codon =~ /^(TAA|TAG|TGA)$/) {
                last;
            }
        }
        

        
        my $framed_kmer = "${kmer}-${frame}";


        $repeat_counter{$framed_kmer}++;
        
        #if ($repeat_counter{$framed_kmer} > $MAX_REPEAT_COUNT) { next; }
        

        my $loglikelihood = $scores{$framed_kmer} || 0;


        #print "$i\t$framed_kmer\t$loglikelihood\n";

        
        
        $score += $loglikelihood;
    }
    
    #print "Score: $score\n";
    

    return($score);

}


####
sub parse_kmer_scores {
	my ($kmer_scores_file) = @_;

	my %scores;
	open (my $fh, $kmer_scores_file) or die "Error, cannot open $kmer_scores_file";
	while (<$fh>) {
		chomp;
		if (/^\#/) { next; } # skip header and comments
        my ($token, $count, $count_kmerminus1, $score) = split (/\t/);
		$scores{$token} = $score;
	}
	close $fh;

	return (%scores);
}

