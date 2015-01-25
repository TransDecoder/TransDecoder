#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Fasta_reader;
use Nuc_translator;

my $usage = "usage: $0 CDS hexamerScores\n\n";

my $cds_file = $ARGV[0] or die $usage;
my $hexamer_scores_file = $ARGV[1] or die $usage;


my %scores = &parse_hexamer_scores($hexamer_scores_file);

main: {
		
	my $fasta_reader = new Fasta_reader($cds_file);
	while (my $seq_obj = $fasta_reader->next()) {
		
		my $accession = $seq_obj->get_accession();
		my $sequence = uc $seq_obj->get_sequence();

		my $score1 = &score_CDS($sequence);				
		my $score2 = &score_CDS(substr($sequence, 1));
		my $score3 = &score_CDS(substr($sequence, 2));

		my $rev_seq = &reverse_complement($sequence);
		
		my $score4 = &score_CDS($rev_seq);
		my $score5 = &score_CDS(substr($rev_seq, 1));
		my $score6 = &score_CDS(substr($rev_seq, 2));
		
		printf("$accession\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",
			   length($sequence), 
               $score1, $score2, $score3,
			   $score4, $score5, $score6);
		
	}
	
	exit(0);

}

####
sub score_CDS {
	my ($sequence) = @_;
	
	my $seq_length = length($sequence);
	
	if ($seq_length < 5) {
		return(0);
	}

	## init score to first pentamer
	my $pentamer = substr($sequence, 0, 5);
	my $framed_pentamer = "${pentamer}-0";
	my $score = $scores{$framed_pentamer} || 0;
	
	for (my $i = 5; $i <= $seq_length - 6; $i++) {
		my $hexamer = substr($sequence, $i, 6);
		my $frame = $i % 3;
		my $framed_hexamer = "${hexamer}-${frame}";
		my $hex_score = $scores{$framed_hexamer} || 0;
		$score += $hex_score;
	}
	
	return ($score);
}


####
sub parse_hexamer_scores {
	my ($hexamer_scores_file) = @_;

	my %scores;
	open (my $fh, $hexamer_scores_file) or die "Error, cannot open $hexamer_scores_file";
	while (<$fh>) {
		chomp;
		my ($token, $score) = split (/\t/);
		$scores{$token} = $score;
	}
	close $fh;

	return (%scores);
}

