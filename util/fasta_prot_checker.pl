#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Fasta_reader;


my $usage = "usage: $0 proteins.fa\n";

my $fasta_file = $ARGV[0];
unless ($fasta_file) {
    die $usage;
}

my $fasta_reader = new Fasta_reader($fasta_file);


## make intervening stop codons fatal:

my $seen_intervening_stop = 0;

while (my $seq_obj = $fasta_reader->next()) {
	
	my $error;
	
    my $header = $seq_obj->get_header();
    my $sequence = $seq_obj->get_sequence();
	    
	unless ($sequence =~ /\*$/) {
		$error .= "\tNo stop codon.";
    }
	
	chop $sequence; # get rid of any trailing stop
	my $numstops = 0;
	while ($sequence =~ /\*/g) {
		$numstops++;
    }
	if ($numstops || !($sequence)) {
		$error .= "\t\*$numstops";
        $seen_intervening_stop++;
    } 
	
    unless ($sequence =~ /^m/i) {
		$error .= "\tDoesn't start with M.";
    }
    
    if ($error && $numstops) {
		print "$header\tERRORS: $error\n";
	}
    
}

exit($seen_intervening_stop);

