package Decoys;

# Decoy ORFs from shuffled transcripts, for target-decoy FDR estimates of candidate ORFs.
#
# Each transcript is shuffled preserving its dinucleotide counts (Altschul & Erickson 1985; the random
# Eulerian path is drawn with Wilson's algorithm for the last-exit tree), so decoys keep the
# transcript's length and base composition. Within strata of ORF length (log-spaced bins), transcript
# length (log2 bins) and ORF type (complete vs partial), a target ORF with frame contrast c has
#     FDR(c) = decoys with contrast >= c / (shuffles x targets with contrast >= c)
# and its q-value is the lowest FDR of any cutoff that accepts it. Strata without decoys give q = 0.
#
# Same computation, random numbers and output as PyTransDecoder (core/prng.py, core/decoys.py):
# xoshiro128** on 32-bit words seeded by FNV-1a hashes, with all arithmetic below 2^53.

use strict;
use warnings;
use Carp;

my $M32 = 0xFFFFFFFF;
our $LENGTH_BIN_RATIO = 1.1;
our $LENGTH_BIN_START = 50;  # codons
our $TX_BIN_START = 500;     # nt


####
sub fnv1a32 {
    my ($text) = @_;
    my $h = 0x811C9DC5;
    foreach my $b (unpack("C*", $text)) {
        $h = (($h ^ $b) * 0x01000193) & $M32;
    }
    return($h);
}

####
sub _rotl {
    my ($x, $k) = @_;
    return((($x << $k) | ($x >> (32 - $k))) & $M32);
}

####
sub new_rng {
    my ($seed) = @_;
    my @s = map { &fnv1a32("$seed#$_") } (0..3);
    unless ($s[0] || $s[1] || $s[2] || $s[3]) {
        $s[0] = 1;
    }
    return(\@s);
}

####
sub next32 {
    my ($s) = @_;
    my $result = (&_rotl(($s->[1] * 5) & $M32, 7) * 9) & $M32;
    my $t = ($s->[1] << 9) & $M32;
    $s->[2] ^= $s->[0];
    $s->[3] ^= $s->[1];
    $s->[1] ^= $s->[2];
    $s->[0] ^= $s->[3];
    $s->[2] ^= $t;
    $s->[3] = &_rotl($s->[3], 11);
    return($result);
}

####
sub below {
    # uniform integer in [0, n) by rejection
    my ($s, $n) = @_;
    my $limit = int(4294967296 / $n) * $n;
    while (1) {
        my $r = &next32($s);
        if ($r < $limit) {
            return($r % $n);
        }
    }
}

####
sub shuffle_in_place {
    # Fisher-Yates
    my ($s, $aref) = @_;
    for (my $i = $#$aref; $i > 0; $i--) {
        my $j = &below($s, $i + 1);
        @$aref[$i, $j] = @$aref[$j, $i];
    }
}

####
sub dinucleotide_shuffle {
    # random sequence with the same first base, last base and dinucleotide counts as $seq
    my ($seq, $rng) = @_;
    my $n = length($seq);
    if ($n < 3) {
        return($seq);
    }
    my @bases = split(//, $seq);
    my %edges;
    for (my $i = 0; $i < $n - 1; $i++) {
        push (@{$edges{$bases[$i]}}, $bases[$i+1]);
    }
    my $end = $bases[-1];
    # random last-exit tree rooted at the final base (Wilson's loop-erased random walks)
    my %in_tree = ($end => 1);
    my (%nxt, %last);
    foreach my $v (sort keys %edges) {
        my $u = $v;
        while (! $in_tree{$u}) {
            my $out = $edges{$u};
            $nxt{$u} = $out->[ &below($rng, scalar(@$out)) ];
            $u = $nxt{$u};
        }
        $u = $v;
        while (! $in_tree{$u}) {
            $in_tree{$u} = 1;
            $last{$u} = $nxt{$u};
            $u = $nxt{$u};
        }
    }
    my %order;
    foreach my $v (sort keys %edges) {
        my @out = @{$edges{$v}};
        if (exists $last{$v}) {
            # remove the first occurrence of the last-exit edge
            for (my $i = 0; $i <= $#out; $i++) {
                if ($out[$i] eq $last{$v}) {
                    splice(@out, $i, 1);
                    last;
                }
            }
        }
        &shuffle_in_place($rng, \@out);
        if (exists $last{$v}) {
            push (@out, $last{$v});
        }
        $order{$v} = \@out;
    }
    my %pos;
    my @res = ($bases[0]);
    my $u = $bases[0];
    for (1 .. $n - 1) {
        my $w = $order{$u}->[ $pos{$u}++ ];
        push (@res, $w);
        $u = $w;
    }
    return(join("", @res));
}

####
sub length_bin {
    # largest k with codons >= START * RATIO**k (powers by repeated multiplication); -1 below START
    my ($codons) = @_;
    if ($codons < $LENGTH_BIN_START) {
        return(-1);
    }
    my $k = 0;
    my $edge = $LENGTH_BIN_START * $LENGTH_BIN_RATIO;
    while ($codons >= $edge) {
        $k++;
        $edge *= $LENGTH_BIN_RATIO;
    }
    return($k);
}

####
sub tx_bin {
    # floor(log2(tx_len / 500)) + 1, 0 below 500 nt: the bit length of int(tx_len / 500)
    my ($tx_len) = @_;
    if ($tx_len < $TX_BIN_START) {
        return(0);
    }
    my $x = int($tx_len / $TX_BIN_START);
    my $bits = 0;
    while ($x > 0) {
        $bits++;
        $x >>= 1;
    }
    return($bits);
}

####
sub stratum {
    # stratum key as written by PyTransDecoder: "<len bin>,<tx bin>,<True|False>"
    my ($codons, $tx_len, $complete) = @_;
    return(join(",", &length_bin($codons), &tx_bin($tx_len), ($complete ? "True" : "False")));
}

####
sub qvalues {
    # $targets_href: orf_id => [stratum, contrast]; $decoys_aref: [ [stratum, contrast], ... ]
    my ($targets_href, $decoys_aref, $n_shuffles) = @_;
    my (%tby, %dby);
    foreach my $orf_id (keys %$targets_href) {
        my ($k, $c) = @{$targets_href->{$orf_id}};
        push (@{$tby{$k}}, $c);
    }
    foreach my $d (@$decoys_aref) {
        push (@{$dby{$d->[0]}}, $d->[1]);
    }
    my %q_of;
    foreach my $k (keys %tby) {
        my @tc = sort {$a <=> $b} @{$tby{$k}};
        my @dc = sort {$a <=> $b} @{$dby{$k} || []};
        my $nd = scalar(@dc);
        my $nt = scalar(@tc);
        my @q;
        my $run_min;
        my $j = 0;  # decoys below the current contrast (bisect_left)
        for (my $i = 0; $i < $nt; $i++) {
            my $c = $tc[$i];
            while ($j < $nd && $dc[$j] < $c) { $j++; }
            my $fdr = ($nd - $j) / $n_shuffles / ($nt - $i);
            $run_min = $fdr if (! defined($run_min) || $fdr < $run_min);
            push (@q, $run_min);
        }
        $q_of{$k} = [\@tc, \@q];
    }
    my %out;
    foreach my $orf_id (keys %$targets_href) {
        my ($k, $c) = @{$targets_href->{$orf_id}};
        my ($tc, $q) = @{$q_of{$k}};
        # first index with tc >= c (bisect_left)
        my ($lo, $hi) = (0, scalar(@$tc));
        while ($lo < $hi) {
            my $mid = int(($lo + $hi) / 2);
            if ($tc->[$mid] < $c) { $lo = $mid + 1; } else { $hi = $mid; }
        }
        $out{$orf_id} = $q->[$lo];
    }
    return(%out);
}

####
sub format_g5 {
    # Python's format(x, '.5g')
    my ($x) = @_;
    return(sprintf("%.5g", $x));
}

1;
