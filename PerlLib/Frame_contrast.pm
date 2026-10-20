package Frame_contrast;

# Frame contrast of a coding Markov model: the in-frame score of a sequence minus the better of
# its scores shifted by one and two bases (as in ORF selection), and the median per-base contrast
# of the training ORFs. Used by start codon refinement and by the short-ORF test in ORF selection.
# Same computation as PyTransDecoder (predict.py: _frame_contrast, _coding_contrast).

use strict;
use warnings;
use Carp;
use List::Util qw(min max);
use Fasta_reader;

my $MARKOV_ORDER = 5;


####
sub new {
    my ($packagename, $hexamer_scores_file, $stop_codons_aref) = @_;

    my $self = { scores => {}, stops => { map { $_ => 1 } @$stop_codons_aref } };
    open (my $fh, $hexamer_scores_file) or confess "Error, cannot open $hexamer_scores_file";
    while (<$fh>) {
        chomp;
        if (/^\#/) { next; }
        my ($token, $count, $count_kmerminus1, $score) = split (/\t/);
        $self->{scores}->{$token} = $score;
    }
    close $fh;

    bless($self, $packagename);
    return($self);
}


####
sub score_CDS_via_Markov {
    my ($self, $sequence) = @_;

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

        if ($i == $seq_length - 3 && $frame == 0 && $self->{stops}->{substr($sequence, $i, 3)}) {
            last; # terminal stop codon
        }
        $score += $self->{scores}->{"${kmer}-${frame}"} || 0;
    }
    return($score);
}


####
sub frame_contrast {
    my ($self, $seq) = @_;

    # in-frame Markov score minus the better of the two shifted-frame scores
    return($self->score_CDS_via_Markov($seq)
           - max($self->score_CDS_via_Markov(substr($seq, 1)), $self->score_CDS_via_Markov(substr($seq, 2))));
}


####
sub median_training_contrast {
    my ($self, $training_orfs_file) = @_;

    # median per-base frame contrast of the training ORFs (undef if none)
    my @per_base;
    my $fasta_reader = new Fasta_reader($training_orfs_file);
    while (my $seq_obj = $fasta_reader->next()) {
        my $cds_seq = uc $seq_obj->get_sequence();
        if (length($cds_seq) >= 30) {
            push (@per_base, $self->frame_contrast($cds_seq) / length($cds_seq));
        }
    }
    unless (@per_base) {
        return(undef);
    }
    @per_base = sort {$a<=>$b} @per_base;
    my $mid = int(scalar(@per_base) / 2);
    return( (scalar(@per_base) % 2) ? $per_base[$mid] : ($per_base[$mid-1] + $per_base[$mid]) / 2 );
}

1;
