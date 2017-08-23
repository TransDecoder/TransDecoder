#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 hexamer.scores\n\n";

my $hexamer_file = $ARGV[0] or die $usage;


open(my $fh, $hexamer_file) or die "Error, cannot open file $hexamer_file";
my $header = <$fh>;

my @hex_scores;
while (<$fh>) {
    chomp;
    my ($framed_kmer, $kmer_count, $prefix_count, $loglikelihood) = split(/\t/);
    my $struct = { framed_kmer => $framed_kmer,
                   kmer_count => $kmer_count,
                   prefix_count => $prefix_count,
                   loglikelihood => $loglikelihood,
    };

    push (@hex_scores, $struct);
}

@hex_scores = sort {$a->{loglikelihood} <=> $b->{loglikelihood}} @hex_scores;

# adjust extremes
# define IQR
my $left_iqr = $hex_scores[int($#hex_scores * 0.25)]->{loglikelihood};
my $right_iqr = $hex_scores[int($#hex_scores * 0.75)]->{loglikelihood};
my $center_iqr = ($left_iqr + $right_iqr)/2;
my $iqr_1pt5 = 1.5 * ($right_iqr - $left_iqr);
my $lower_1pt5_iqr = $left_iqr - $iqr_1pt5;
my $upper_1pt5_iqr = $right_iqr + $iqr_1pt5;

print STDERR "left_iqr: $left_iqr\n"
    . "right_iqr: $right_iqr\n"
    . "lower_1pt5_iqr = $lower_1pt5_iqr\n"
    . "upper_1pt5_iqr = $upper_1pt5_iqr\n";


# adjust vals, reassign vals to outliers
for my $struct (@hex_scores) {
    if ($struct->{loglikelihood} < $lower_1pt5_iqr) {
        $struct->{loglikelihood} = $lower_1pt5_iqr;
    }
    elsif ($struct->{loglikelihood} > $upper_1pt5_iqr) {
        $struct->{loglikelihood} = $upper_1pt5_iqr;
    }
}


# resort by kmer name
@hex_scores = sort {$a->{framed_kmer} cmp $b->{framed_kmer} } @hex_scores;

# output results

print $header;
foreach my $struct (@hex_scores) {
    print join("\t", $struct->{framed_kmer},
               $struct->{kmer_count},
               $struct->{prefix_count},
               $struct->{loglikelihood}) . "\n";
}


exit(0);


                   
