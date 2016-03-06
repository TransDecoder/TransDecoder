#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Fasta_reader;

my $usage = "\n\n\tusage: $0 file.transdecoder.pep\n\n";

my $pep_file = $ARGV[0] or die $usage;

 main: {

     my %contig_to_longest_prot;
     
     my $fasta_reader = new Fasta_reader($pep_file);
     while (my $seq_obj = $fasta_reader->next()) {
         my $header = $seq_obj->get_header();

         #  len:365 (+) asmbl_104:2-1096(+)
         if ($header =~ /len:(\d+) \([+-]\) (\S+):\d+-\d+\([+-]\)/) {
             my $len = $1;
             my $contig = $2;


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


    
