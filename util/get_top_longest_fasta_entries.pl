#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Fasta_reader;

my $usage = "usage: $0 file.fasta numTopLongest\n\n";

my $file = $ARGV[0] or die $usage;
my $num_longest = $ARGV[1] or die $usage;

main: {

	my @entries;

	my $fasta_reader = new Fasta_reader($file);
	while (my $seq_obj = $fasta_reader->next()) {

		my $seq = $seq_obj->get_sequence();
		my $len = length($seq);
		
		push (@entries, [$seq_obj, $len]);

	}

	@entries = reverse sort {$a->[1]<=>$b->[1]} @entries;
	
	my $counter = 0;
	foreach my $entry (@entries) {
		
		my ($seq_obj, $num) = @$entry;
		
		print $seq_obj->get_FASTA_format();
		
		$counter++;

		if ($counter >= $num_longest) {
			last;
		}

	}
		
	exit(0);
}

