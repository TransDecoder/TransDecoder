#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use Gene_obj_indexer;
use GFF3_utils;
use Carp;
use Fasta_retriever;
use Nuc_translator;

$|++;

my $usage = "\n\nusage: $0 gene_db.inx_file transcripts.fasta train_accs.list pre_length post_length\n\n";

my $inx_file = $ARGV[0] or die $usage;
my $transcripts_fasta = $ARGV[1] or die $usage;
my $train_accs_file = $ARGV[2] or die $usage;
my $pre_length = $ARGV[3] or die $usage;
my $post_length = $ARGV[4] or die $usage;


main: {

    my $gene_obj_indexer = new Gene_obj_indexer( { "use" => "$inx_file" } );

    my $fasta_retriever = new Fasta_retriever($transcripts_fasta);

    my @train_accs = `cat $train_accs_file`;
    chomp @train_accs;
    
    my @nuc_frequency;
    
    my $seq_counter = 0;

    foreach my $gene_id (@train_accs) {

        my $gene_obj = $gene_obj_indexer->get_gene($gene_id);
        
        my $com_name = $gene_obj->{com_name};

        unless ($com_name =~ /type:complete/ || $com_name =~ /type:3prime_partial/) { next; } # need start codon


        my $orient = $gene_obj->get_orientation();
        my ($lend, $rend) = sort {$a<=>$b} $gene_obj->get_model_span();
        my $contig_id = $gene_obj->{asmbl_id};
        
        my $trans_seq = uc $fasta_retriever->get_seq($contig_id);

        my $seq_length = length($trans_seq);
        
        if ($orient eq '-') {
            $trans_seq = &reverse_complement($trans_seq);
            $lend = $seq_length - $rend + 1;
        }

        my $atg = substr($trans_seq, $lend-1, 3);
        unless ($atg eq "ATG") {
            die "Error, missing ATG start: $atg ";
        }
        
        my $begin_profile = $lend - 1 - $pre_length;
        my $profile_length = $pre_length + 3 + $post_length;

        if ($begin_profile >= 1) {
            my $start_codon_context_substring = substr($trans_seq, $begin_profile, $profile_length);
            #print "$start_codon_context_substring\n";
            &add_to_profile(\@nuc_frequency, $start_codon_context_substring);
            $seq_counter++;
        }
        
        
        
    }
 

    &write_PSSM($pre_length, $seq_counter, \@nuc_frequency);
    
   
    exit(0);
}


####
sub add_to_profile {
    my ($nuc_freq_aref, $context_string) = @_;

    my @chars = split(//, $context_string);

    for (my $i = 0; $i <= $#chars; $i++) {
        my $char = $chars[$i];
        $nuc_freq_aref->[$i]->{$char}++;
    }
    
    return;
}


####
sub write_PSSM {
    my ($pre_length, $num_seqs, $nuc_freq_aref) = @_;
    
    my @nucs = qw(G A T C);
    
    my $atg_pos = $pre_length + 1;
    print "#ATG=$atg_pos\n";
    print "#" . join("\t", @nucs) . "\n";
    
    for (my $i = 0; $i < $#$nuc_freq_aref; $i++) {

        my %char_counts = %{$nuc_freq_aref->[$i]};
        
        my @rel_freqs;
        
        foreach my $nuc (@nucs) {
            my $count = $char_counts{$nuc} || 0;

            my $relative_freq = sprintf("%.3f", $count/$num_seqs);

            push (@rel_freqs, $relative_freq);
        }

        print join("\t", @rel_freqs) . "\n";
    }

    return;
}
