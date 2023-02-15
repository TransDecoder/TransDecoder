#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Fasta_reader;

my $usage = "\n\n\tusage: $0 all_long_orfs.pep\n\n";

my $pep_file = $ARGV[0] or die $usage;

 main: {

     my %contig_to_longest_prot;
     
     my $fasta_reader = new Fasta_reader($pep_file);
     while (my $seq_obj = $fasta_reader->next()) {
         my $header = $seq_obj->get_header();
         my $sequence = $seq_obj->get_sequence();
         my $len = length($sequence);

         if ($header =~ / (\S+):\d+-\d+\([+-]\)/) {
             my $contig = $1;
             &store_entry(\%contig_to_longest_prot, $seq_obj, $len, $contig);
         }
         else {
             die "Error, cannot parse header info: $header ";
         }
     }
     
     foreach my $struct (values %contig_to_longest_prot) {
         my $fa = $struct->{seq_obj}->get_FASTA_format();

         print $fa;
     }

     exit(0);
}

####
sub store_entry {
    my ($contig_to_longest_prot_href, $seq_obj, $len, $contig) = @_;

    my $struct = $contig_to_longest_prot_href->{$contig};
    
    if ( (! $struct) || $struct->{length} < $len) {
        
        $struct = { seq_obj => $seq_obj,
                    length => $len,
        };

        $contig_to_longest_prot_href->{$contig} = $struct;
    }

    return;
}


    
