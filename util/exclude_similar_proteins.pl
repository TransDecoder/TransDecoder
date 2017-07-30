#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Carp;
use Nuc_translator;
use Fasta_reader;

my $usage = "\n\nusage: $0 fasta_file [DEBUG]\n\n";

my $fasta_file = $ARGV[0] or die $usage;
my $DEBUG = $ARGV[1] || 0;

my $WMER_SIZE = 5;

my %CAPTURED_WMERS;
my $MAX_PER_ID = 80;
my $MIN_UNIQUE_WMERS_PER_LEN = 0;

main: {

    
    my $fasta_reader = new Fasta_reader($fasta_file);
    
    my @seq_data;
    
    while (my $seq_obj = $fasta_reader->next()) {

        my $acc = $seq_obj->get_accession();
        my $sequence = $seq_obj->get_sequence();
        
        push (@seq_data, { acc => $acc,
                           seq => $sequence,
                           seqlen => length($sequence),
                           seq_obj => $seq_obj,
              } );
    }

    @seq_data = reverse sort {$a->{seqlen}<=>$b->{seqlen}} @seq_data;


    my $count_retained = 0;
    
    foreach my $seq_struct (@seq_data) {

        my $acc = $seq_struct->{acc};
        my $sequence = $seq_struct->{seq};

        my $seq_obj = $seq_struct->{seq_obj};
        
        my $prot = translate_sequence($sequence, 1);

        # get spaced wmers for this target sequence
        my %wmers;
        for (my $i = 0; $i <= length($prot) - $WMER_SIZE; $i+=$WMER_SIZE) {
            my $wmer = substr($prot, $i, $WMER_SIZE);
            $wmers{$wmer} = 1;
        }
        
        # see if wmer profile has been already seen
        my %prots_seen;
        foreach my $wmer (keys %wmers) {
            if (exists $CAPTURED_WMERS{$wmer}) {
                foreach my $prot (@{$CAPTURED_WMERS{$wmer}}) {
                    $prots_seen{$prot}++;
                }
            }
        }

        # determine best match
        my $num_unique_wmers = scalar(keys %wmers);
        my $pct_unique_wmers = sprintf("%.2f", $num_unique_wmers / (length($prot) / $WMER_SIZE) * 100);

        if ($pct_unique_wmers < 30) {
            print STDERR "-skipping $acc as likely low complexity sequence\n";
            next;
        }
        
        my $proxy_per_id = 0;
        if (%prots_seen) {
            my @counts_seen = reverse sort {$prots_seen{$a}<=>$prots_seen{$b}} keys %prots_seen;
            my $max_seen_acc = shift @counts_seen;
            my $max_seen_count = $prots_seen{$max_seen_acc};
            $proxy_per_id = sprintf("%.2f", $max_seen_count / $num_unique_wmers * 100);
            print STDERR "-analyze($acc):  max_seen: $max_seen_count in $max_seen_acc, proxy_per_id: $proxy_per_id\n" if $DEBUG;
        }
        
        if ($proxy_per_id <= $MAX_PER_ID) {
            # ok, keep it.
            
            # recompute wmers as overlapping, unspaced
            my %wmers;
            for (my $i = 0; $i <= length($prot) - $WMER_SIZE; $i++) {
                my $wmer = substr($prot, $i, $WMER_SIZE);
                $wmers{$wmer} = 1;
            }
            # store wmers for this protein
            for my $wmer (keys %wmers) {
                push(@{$CAPTURED_WMERS{$wmer}}, $acc);
            }
            
            print $seq_obj->get_FASTA_format();
            $count_retained++;
        }
        else {
            print STDERR "-skipping training candidate: $acc, not unique enough\n";
        }
    }

    print STDERR "\n\t-redundancy-minimized set includes $count_retained / " . scalar(@seq_data) . " = " . sprintf("%.2f", ($count_retained/scalar(@seq_data) * 100)) . "%\n\n";


    exit(0);
    
}

