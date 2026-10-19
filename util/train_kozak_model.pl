#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Fasta_reader;
use Nuc_translator;
use Kozak_model;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);

# Train the start-codon (Kozak) context model used by start_codon_refinement.pl: log-odds of
# bases at -6..-1 and +4..+5 around the ATG, training ORF starts vs in-frame ATGs inside the
# training ORFs.

my $MIN_TRAINING_STARTS = 50;

my $usage = <<__EOUSAGE__;

#####################################################################################
#
#  --transcripts <string>       transcripts.fasta file targeted by transdecoder
#
#  --training_orfs <string>     training ORF CDS fasta (headers include acc:lend-rend(strand))
#
#  --out_prefix <string>        writes:
#                                 <prefix>.kozak.tsv         log-odds model (read by start_codon_refinement.pl)
#                                 <prefix>.kozak.counts.tsv  base counts at start and internal ATG sites
#                                 <prefix>.kozak.logo.svg    sequence logo
#                               or, with fewer than $MIN_TRAINING_STARTS training ORFs starting with ATG,
#                                 <prefix>.kozak.SKIPPED
#
#######################################################################################

__EOUSAGE__

    ;

my ($transcripts_file, $training_orfs_file, $out_prefix);

&GetOptions('transcripts=s' => \$transcripts_file,
            'training_orfs=s' => \$training_orfs_file,
            'out_prefix=s' => \$out_prefix,
    );

unless ($transcripts_file && $training_orfs_file && $out_prefix) {
    die $usage;
}


main: {

    my %seqs = new Fasta_reader($transcripts_file)->retrieve_all_seqs_hash();

    my $pos_counts = &Kozak_model::new_counts();
    my $neg_counts = &Kozak_model::new_counts();
    my $num_starts = 0;

    my $fasta_reader = new Fasta_reader($training_orfs_file);
    while (my $seq_obj = $fasta_reader->next()) {
        unless ($seq_obj->get_header() =~ /(\S+):(\d+)-(\d+)\(([+-])\)/) { next; }
        my ($acc, $lend, $rend, $strand) = ($1, $2, $3, $4);
        unless (exists $seqs{$acc}) { next; }
        ($lend, $rend) = sort {$a<=>$b} ($lend, $rend);

        my $seq = uc $seqs{$acc};
        if ($strand eq '-') {
            $seq = &reverse_complement($seq);
            ($lend, $rend) = (length($seq) - $rend + 1, length($seq) - $lend + 1);
        }
        my ($start, $end) = ($lend - 1, $rend);
        unless (substr($seq, $start, 3) eq "ATG") { next; }
        $num_starts++;

        &Kozak_model::add_site($pos_counts, $seq, $start);
        for (my $i = $start + 3; $i < $end - 3; $i += 3) {
            if (substr($seq, $i, 3) eq "ATG") {
                &Kozak_model::add_site($neg_counts, $seq, $i);
            }
        }
    }

    unlink("$out_prefix.kozak.SKIPPED", "$out_prefix.kozak.tsv");
    if ($num_starts < $MIN_TRAINING_STARTS) {
        my $msg = "only $num_starts training ORFs start with ATG (minimum $MIN_TRAINING_STARTS); "
            . "start refinement will use frame contrast only";
        print STDERR "WARNING: $msg\n";
        open(my $ofh, ">$out_prefix.kozak.SKIPPED") or die "Error, cannot write $out_prefix.kozak.SKIPPED";
        print $ofh "$msg\n";
        close $ofh;
        exit(0);
    }

    my $model = &Kozak_model::log_odds_model($pos_counts, $neg_counts);
    &Kozak_model::write_model($model, $num_starts, "$out_prefix.kozak.tsv");
    &Kozak_model::write_counts($pos_counts, $neg_counts, "$out_prefix.kozak.counts.tsv");
    &Kozak_model::write_logo_svg($pos_counts, $model, $num_starts, "$out_prefix.kozak.logo.svg");

    print STDERR "-Kozak model from $num_starts training starts: $out_prefix.kozak.tsv, $out_prefix.kozak.logo.svg\n";

    exit(0);
}
