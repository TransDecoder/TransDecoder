#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Fasta_reader;
use Nuc_translator;

## hexamer stats
my %framed_hexamers;
my %background_hexamers;


## pentamer stats
my %framed_pentamers;
my %background_base_probs;
my %framed_all_pentamer_counts;

my $usage = "usage: $0 targetCDSs base_probs.dat\n\n";

my $target_CDS = $ARGV[0] or die $usage;
my $base_probs_dat_file = $ARGV[1] or die $usage;
my $debug = $ARGV[2];

main: {

	&parse_targetCDSs($target_CDS);

	&parse_background($base_probs_dat_file);

	&add_pseudocounts();

	&report_logliklihood_ratios();
	
	exit(0);
}





####
sub report_logliklihood_ratios {
	

	## Markov-based probabilities (5th order markov chain):
	
	foreach my $framed_hexamer (sort keys %framed_hexamers) {
		my ($hexamer, $frame) = split (/-/, $framed_hexamer);
	
        if ($hexamer =~ /[^GATC]/) {
            ## ignoring hexamers containing non-GATC bases
            next;
        }
        
		my $pentamer = substr($hexamer, 0, 5);

		my $framed_hexamer_count = $framed_hexamers{$framed_hexamer};
		my $framed_pentamer_count = $framed_pentamers{"${pentamer}-${frame}"};

		my $markov_prob_framed = $framed_hexamer_count / $framed_pentamer_count;

		my $last_base = substr($hexamer, 5, 1);
		my $background_prob = $background_base_probs{$last_base} or die "Error, no background probability set for base: $last_base of hexamer $hexamer";;
        
        my $logliklihood = log($markov_prob_framed / $background_prob);
            
        print "$framed_hexamer\t$logliklihood\n";
        
    }
    


	## The Initialization Matrix based on framed pentamer frequencies.

	foreach my $framed_pentamer (sort keys %framed_pentamers) {
		
        if ($framed_pentamer =~ /[^GATC]/) { 
            next;
        }
        
		my ($pentamer, $frame) = split (/-/, $framed_pentamer);

		my $frame_counts = $framed_all_pentamer_counts{$frame};
		my $framed_pentamer_counts = $framed_pentamers{$framed_pentamer};

		my $prob_framed_pentamer = $framed_pentamer_counts / $frame_counts;

		## now background
		my @bases = split(//, $pentamer);
        my $prob_background_pentamer = 1;
        foreach my $base (@bases) {
            $prob_background_pentamer *= $background_base_probs{$base};
        }
        
		my $logliklihood = log($prob_framed_pentamer / $prob_background_pentamer);

		print "$framed_pentamer\t$logliklihood\n";

	}

	return;
}

####
sub add_pseudocounts {
	
	foreach my $framed_hexamer (keys %framed_hexamers) {
		my ($hexamer, $frame) = split (/-/, $framed_hexamer);
		
		my $pentamer = substr($hexamer, 0, 5);
		
		$framed_hexamers{$framed_hexamer}++;
		$framed_pentamers{"${pentamer}-${frame}"}++;
		$framed_all_pentamer_counts{$frame}++;

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

		for (my $i = 0; $i <= $seq_len - 5; $i++) {
			my $frame = $i % 3;
			my $pentamer = substr($sequence, $i, 5);
			$framed_pentamers{"${pentamer}-${frame}"}++;
			$framed_all_pentamer_counts{$frame}++;
			
			if ($i <= $seq_len - 6) { 
				# got a hexamer
				my $hexamer = substr($sequence, $i, 6);
				$framed_hexamers{"${hexamer}-${frame}"}++;
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

