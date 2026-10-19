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

# An ORF that runs to the 5' end of its transcript (no start codon) is either truncated
# (5' partial) or extends into the 5' UTR. Two signals decide, and the ORF is trimmed to start
# at its first in-frame ATG unless the segment before that ATG looks coding:
#
#  - frame contrast of the segment: its coding Markov score in frame minus the better of the
#    scores shifted by one and two bases (as in ORF selection). Coding sequence scores well
#    above its shifted frames; UTR does not.
#  - Kozak context of the ATG: log-odds of the bases at -6..-1 and +4..+5 (A of ATG = +1),
#    trained on the training ORFs' starts against in-frame ATGs inside them
#    (train_kozak_model.pl; --kozak_model).
#
#  trim if  contrast < (contrast_alpha + kozak_weight * kozak_log_odds) * mu * segment_length
#
# where mu is the median per-base frame contrast of the training ORFs.

my $contrast_alpha = 0.4;
my $kozak_weight = 0.15;

my $MARKOV_ORDER = 5;

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
#  --kozak_model <string>       Kozak log-odds model from train_kozak_model.pl (if absent, frame contrast only)
#
#  --workdir <string>           TransDecoder working dir (default: (--transcripts val) + ".transdecoder_dir")
#
#  --contrast_alpha <float>     fraction of the training ORFs' per-base frame contrast the segment
#                               before the first ATG must reach to stay 5' partial (default: $contrast_alpha)
#
#  --kozak_weight <float>       weight of the ATG's Kozak log-odds in that cutoff (default: $kozak_weight)
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
my %STOP_CODONS = map { $_ => 1 } &Nuc_translator::get_stop_codons();

my %MARKOV_SCORES = &parse_kmer_scores($hexamer_scores_file);


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

    my $coding_contrast = &median_training_contrast($training_orfs_file);
    unless (defined $coding_contrast) {
        print STDERR "WARNING: no training ORFs for start refinement; ORFs emitted unchanged\n";
        open(my $fh, $gff3_file) or die "Error, cannot open $gff3_file";
        print while (<$fh>);
        close $fh;
        exit(0);
    }
    printf STDERR "-training ORFs: median frame contrast %.4f per base; open 5' ends kept when the upstream segment reaches %.2f of it\n",
        $coding_contrast, $contrast_alpha;

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

    my $contrast = &frame_contrast(substr($orf_seq, 0, $atg));
    my $kozak_score = ($kozak_model) ? &Kozak_model::score($kozak_model, $transcript_seq, $start_index + $atg) : 0;
    my $cutoff = $contrast_alpha + $kozak_weight * $kozak_score;

    if ($contrast >= $cutoff * $coding_contrast * $atg) {
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

    printf STDERR "# refined start codon: %s -> %d (frame contrast %.3f of training, kozak %.2f)\n",
        $gene_obj->{Model_feat_name}, $new_start, $contrast / ($coding_contrast * $atg), $kozak_score if $DEBUG;

    return(1);
}


####
sub frame_contrast {
    my ($seq) = @_;

    # in-frame Markov score minus the better of the two shifted-frame scores
    return(&score_CDS_via_Markov($seq)
           - max(&score_CDS_via_Markov(substr($seq, 1)), &score_CDS_via_Markov(substr($seq, 2))));
}


####
sub score_CDS_via_Markov {
    my ($sequence) = @_;

    # as in score_CDS_likelihood_all_6_frames.pl
    my $seq_length = length($sequence);
    if ($seq_length < $MARKOV_ORDER + 1) {
        return(0);
    }

    my $score = 0;
    for (my $i = 0; $i < $seq_length; $i++) {
        my $frame = $i % 3;
        my $markov_use = min($i, $MARKOV_ORDER);
        my $kmer = substr($sequence, $i - $markov_use, $markov_use + 1);

        if ($i == $seq_length - 3 && $frame == 0 && $STOP_CODONS{substr($sequence, $i, 3)}) {
            last; # terminal stop codon
        }
        $score += $MARKOV_SCORES{"${kmer}-${frame}"} || 0;
    }
    return($score);
}


####
sub parse_kmer_scores {
    my ($kmer_scores_file) = @_;

    my %scores;
    open (my $fh, $kmer_scores_file) or die "Error, cannot open $kmer_scores_file";
    while (<$fh>) {
        chomp;
        if (/^\#/) { next; }
        my ($token, $count, $count_kmerminus1, $score) = split (/\t/);
        $scores{$token} = $score;
    }
    close $fh;

    return (%scores);
}


####
sub parse_training_orfs {
    my ($training_orfs_file) = @_;

    # returns list of [transcript_acc, lend, rend, strand, cds_seq]
    my @orfs;
    my $fasta_reader = new Fasta_reader($training_orfs_file);
    while (my $seq_obj = $fasta_reader->next()) {
        my $header = $seq_obj->get_header();
        my $cds_seq = uc $seq_obj->get_sequence();
        if ($header =~ /(\S+):(\d+)-(\d+)\(([+-])\)/) {
            my ($acc, $lend, $rend, $strand) = ($1, $2, $3, $4);
            ($lend, $rend) = sort {$a<=>$b} ($lend, $rend);
            push (@orfs, [$acc, $lend, $rend, $strand, $cds_seq]);
        }
        else {
            push (@orfs, [undef, undef, undef, undef, $cds_seq]);
        }
    }
    return(@orfs);
}


####
sub median_training_contrast {
    my ($training_orfs_file) = @_;

    my @per_base;
    foreach my $orf (&parse_training_orfs($training_orfs_file)) {
        my $cds_seq = $orf->[4];
        if (length($cds_seq) >= 30) {
            push (@per_base, &frame_contrast($cds_seq) / length($cds_seq));
        }
    }
    unless (@per_base) {
        return(undef);
    }
    @per_base = sort {$a<=>$b} @per_base;
    my $mid = int(scalar(@per_base) / 2);
    return( (scalar(@per_base) % 2) ? $per_base[$mid] : ($per_base[$mid-1] + $per_base[$mid]) / 2 );
}
