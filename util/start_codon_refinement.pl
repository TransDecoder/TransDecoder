#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use GFF3_utils2;
use Fasta_reader;
use Nuc_translator;
use Carp;
use List::Util qw (min max);
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use File::Basename;
use Kozak_model;
use Frame_contrast;

# An ORF that runs to the 5' end of its transcript (no start codon) is either truncated
# (5' partial) or extends into the 5' UTR. Two scores are compared for the segment before its
# first in-frame ATG, both in bases' worth of typical coding signal (mu, the training ORFs'
# median per-base frame contrast):
#
#  coding score = frame contrast of the segment / mu
#                 (coding Markov score in frame minus the better of the scores shifted by one and
#                 two bases, as in ORF selection: about 1 per base for coding sequence, 0 for UTR)
#  UTR score    = contrast_alpha * segment_length + kozak_weight * kozak_bits
#                 (kozak_bits: log2-odds of the bases at -6..-1 and +4..+5 around the ATG, A of
#                 ATG = +1, trained on the training ORFs' starts against in-frame ATGs inside them;
#                 train_kozak_model.pl, --kozak_model)
#
# The ORF is trimmed to the ATG when the UTR score is higher; otherwise it stays 5' partial.
# Each bit of Kozak score lets the ATG overrule kozak_weight bases of coding signal.

my $contrast_alpha = 0.4;
my $kozak_weight = 17;

my $usage = <<__EOUSAGE__;

#####################################################################################
#
#  --transcripts <string>       transcripts.fasta file targeted by transdecoder
#
#  --gff3_file <string>         target gff3 file
#
#  --training_orfs <string>     training ORF CDS fasta (headers include acc:lend-rend(strand))
#
#  --hexamer_scores <string>    Markov model log-likelihood file (hexamer.scores)
#
#  optional:
#
#  --kozak_model <string>       Kozak log2-odds model from train_kozak_model.pl (if absent, coding score only)
#
#  --workdir <string>           TransDecoder working dir (default: (--transcripts val) + ".transdecoder_dir")
#
#  --contrast_alpha <float>     UTR score per base of the segment before the first ATG: the fraction
#                               of typical coding signal it may show and still be called UTR
#                               (default: $contrast_alpha)
#
#  --kozak_weight <float>       UTR score per bit of the ATG's Kozak score, in bases of coding
#                               signal (default: $kozak_weight)
#
#  --genetic_code|G <string>    genetic code (default: Universal)
#
#  --debug                      verbose
#
#######################################################################################

__EOUSAGE__

    ;


my $transcripts_file;
my $gff3_file;
my $training_orfs_file;
my $hexamer_scores_file;
my $kozak_model_file;
my $genetic_code = "Universal";
my $DEBUG = 0;
my $workdir;

&GetOptions('transcripts=s' => \$transcripts_file,
            'gff3_file=s' => \$gff3_file,
            'training_orfs=s' => \$training_orfs_file,
            'hexamer_scores=s' => \$hexamer_scores_file,
            'kozak_model=s' => \$kozak_model_file,
            'contrast_alpha=f' => \$contrast_alpha,
            'kozak_weight=f' => \$kozak_weight,
            'genetic_code|G=s' => \$genetic_code,
            'debug' => \$DEBUG,
            'workdir=s' => \$workdir,
    );


unless ($transcripts_file && $gff3_file && $training_orfs_file && $hexamer_scores_file) {
    die $usage;
}

&Nuc_translator::use_specified_genetic_code($genetic_code);
my $FRAME_CONTRAST = new Frame_contrast($hexamer_scores_file, [ &Nuc_translator::get_stop_codons() ]);


main: {

    unless ($workdir) {
        $workdir = basename($transcripts_file) . ".transdecoder_dir";
    }
    unless (-d $workdir) {
        die "Error, cannot locate transdecoder working directory as: $workdir";
    }

    print STDERR "-reading transcripts: $transcripts_file\n" if $DEBUG;
    my $fasta_reader = new Fasta_reader($transcripts_file);
    my %seqs = $fasta_reader->retrieve_all_seqs_hash();

    my $coding_contrast = $FRAME_CONTRAST->median_training_contrast($training_orfs_file);
    unless (defined $coding_contrast) {
        print STDERR "WARNING: no training ORFs for start refinement; ORFs emitted unchanged\n";
        open(my $fh, $gff3_file) or die "Error, cannot open $gff3_file";
        print while (<$fh>);
        close $fh;
        exit(0);
    }
    printf STDERR "-training ORFs: median frame contrast %.4f per base; open 5' ends trimmed to the first ATG when %.2f x segment length + %g x Kozak bits exceeds the segment's coding score\n",
        $coding_contrast, $contrast_alpha, $kozak_weight;

    my $kozak_model;
    if ($kozak_model_file && -s $kozak_model_file) {
        $kozak_model = &Kozak_model::read_model($kozak_model_file);
    }
    else {
        print STDERR "-no Kozak model; start refinement uses frame contrast only\n";
    }

    print STDERR "-parsing orf annotations: $gff3_file\n" if $DEBUG;
    my $gene_obj_indexer_href = {};
    my $asmbl_id_to_gene_list_href = &GFF3_utils2::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer_href);

    my $num_starts_revised = 0;

    foreach my $transcript_acc (sort keys %$asmbl_id_to_gene_list_href) {

        my $transcript_seq = uc $seqs{$transcript_acc};

        foreach my $gene_id (@{$asmbl_id_to_gene_list_href->{$transcript_acc}}) {

            my $gene_obj = $gene_obj_indexer_href->{$gene_id};

            if (&refine_start_codon_position($gene_obj, $transcript_seq, $coding_contrast, $kozak_model)) {
                $num_starts_revised++;

                ## update naming convention based on now having a start codon.
                if ($gene_obj->{com_name} =~ /internal/) {
                    $gene_obj->{com_name} =~ s/internal/3prime_partial/;
                }
                elsif ($gene_obj->{com_name} =~ /5prime_partial/) {
                    $gene_obj->{com_name} =~ s/5prime_partial/complete/;
                }
            }
            print $gene_obj->to_GFF3_format(source => "transdecoder") . "\n";
        }
    }

    print STDERR "-number of revised start positions: $num_starts_revised\n";

    exit(0);
}


####
sub refine_start_codon_position {
    my ($gene_obj, $transcript_seq, $coding_contrast, $kozak_model) = @_;

    my $orient = $gene_obj->get_orientation();
    my ($lend, $rend) = sort {$a<=>$b} $gene_obj->get_model_span();
    my $orf_len = $rend - $lend + 1;

    my $start_index = $lend - 1; # zero based, oriented
    if ($orient eq '-') {
        $transcript_seq = &reverse_complement($transcript_seq);
        $start_index = length($transcript_seq) - $rend;
    }
    my $orf_seq = substr($transcript_seq, $start_index, $orf_len);

    # only ORFs open at their 5' end
    if (substr($orf_seq, 0, 3) eq "ATG") {
        return(0);
    }

    my $atg;
    for (my $i = 3; $i < $orf_len - 3; $i += 3) {
        if (substr($orf_seq, $i, 3) eq "ATG") {
            $atg = $i;
            last;
        }
    }
    unless (defined $atg) {
        return(0);
    }

    # both scores in bases' worth of typical coding signal
    my $coding_score = $FRAME_CONTRAST->frame_contrast(substr($orf_seq, 0, $atg)) / $coding_contrast;
    my $kozak_score = ($kozak_model) ? &Kozak_model::score($kozak_model, $transcript_seq, $start_index + $atg) : 0;
    my $utr_score = $contrast_alpha * $atg + $kozak_weight * $kozak_score;

    if ($coding_score >= $utr_score) {
        return(0); # upstream segment looks coding: 5' partial
    }

    my $new_start = $start_index + $atg + 1; # 1-based, oriented
    if ($orient eq '-') {
        $new_start = length($transcript_seq) - $new_start + 1;
    }

    my ($exon_obj) = $gene_obj->get_exons();
    my $cds_obj = $exon_obj->get_CDS_obj();
    $cds_obj->{end5} = $new_start;
    $gene_obj->refine_gene_object();

    printf STDERR "# refined start codon: %s -> %d (coding score %.1f, UTR score %.1f, kozak %.2f bits)\n",
        $gene_obj->{Model_feat_name}, $new_start, $coding_score, $utr_score, $kozak_score if $DEBUG;

    return(1);
}


