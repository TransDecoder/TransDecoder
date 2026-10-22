package Kozak_model;

# Start-codon (Kozak) context model used by start codon refinement.
#
# Positions -6..-1 and +4..+5 around an ATG (A of ATG = +1; the ATG itself is not scored).
# Log2-odds (bits, +1 pseudocounts) of training ORF start contexts vs in-frame ATGs inside the
# training ORFs: +1 bit = twice as start-like as an average internal ATG. Same model and file formats as PyTransDecoder (core/kozak.py).

use strict;
use warnings;
use Carp;
use List::Util qw(max);

our $UPSTREAM = 6;
our $DOWNSTREAM = 2;
our @BASES = qw(A C G T);


####
sub position_labels {
    return( (map { "-$_" } reverse(1..$UPSTREAM)), (map { "+$_" } (4..(3 + $DOWNSTREAM))) );
}


####
sub window {
    my ($seq, $atg_index) = @_;

    # bases at -6..-1 and +4..+5 around the ATG at zero-based $atg_index (undef outside the sequence)
    my @positions = (($atg_index - $UPSTREAM) .. ($atg_index - 1),
                     ($atg_index + 3) .. ($atg_index + 2 + $DOWNSTREAM));
    my @bases;
    foreach my $j (@positions) {
        my $base = ($j >= 0 && $j < length($seq)) ? substr($seq, $j, 1) : undef;
        push (@bases, (defined($base) && $base =~ /^[ACGT]$/) ? $base : undef);
    }
    return(@bases);
}


####
sub new_counts {
    return( [ map { {A=>0, C=>0, G=>0, T=>0} } (1..($UPSTREAM + $DOWNSTREAM)) ] );
}


####
sub add_site {
    my ($counts_aref, $seq, $atg_index) = @_;

    my @window = &window($seq, $atg_index);
    for my $k (0..$#window) {
        if (defined $window[$k]) {
            $counts_aref->[$k]->{$window[$k]}++;
        }
    }
}


####
sub log_odds_model {
    my ($pos_counts, $neg_counts) = @_;

    my @model;
    for my $k (0..$#$pos_counts) {
        my ($pt, $nt) = (4, 4);
        foreach my $b (@BASES) {
            $pt += $pos_counts->[$k]->{$b};
            $nt += $neg_counts->[$k]->{$b};
        }
        my %log_odds;
        foreach my $b (@BASES) {
            $log_odds{$b} = log( (($pos_counts->[$k]->{$b} + 1) / $pt) / (($neg_counts->[$k]->{$b} + 1) / $nt) ) / log(2);
        }
        push (@model, \%log_odds);
    }
    return(\@model);
}


####
sub write_model {
    my ($model, $num_starts, $filename) = @_;

    my @labels = &position_labels();
    open(my $ofh, ">$filename") or confess "Error, cannot write to $filename";
    print $ofh "#position\tA\tC\tG\tT  (log2-odds in bits, start vs internal in-frame ATG; $num_starts training starts)\n";
    for my $k (0..$#$model) {
        print $ofh join("\t", $labels[$k], map { sprintf("%.3f", $model->[$k]->{$_}) } @BASES) . "\n";
    }
    close $ofh;
}


####
sub read_model {
    my ($filename) = @_;

    my @model;
    open(my $fh, $filename) or confess "Error, cannot open $filename";
    while (<$fh>) {
        if (/^\#/) { next; }
        chomp;
        my ($label, @vals) = split(/\t/);
        unless (scalar(@vals) == 4) { next; }
        my %log_odds;
        @log_odds{@BASES} = @vals;
        push (@model, \%log_odds);
    }
    close $fh;
    unless (scalar(@model) == $UPSTREAM + $DOWNSTREAM) {
        confess "Error, $filename does not have " . ($UPSTREAM + $DOWNSTREAM) . " positions";
    }
    return(\@model);
}


####
sub write_counts {
    my ($pos_counts, $neg_counts, $filename) = @_;

    my @labels = &position_labels();
    open(my $ofh, ">$filename") or confess "Error, cannot write to $filename";
    print $ofh join("\t", "#position", (map { "start_$_" } @BASES), (map { "internal_$_" } @BASES)) . "\n";
    for my $k (0..$#$pos_counts) {
        print $ofh join("\t", $labels[$k], (map { $pos_counts->[$k]->{$_} } @BASES),
                        (map { $neg_counts->[$k]->{$_} } @BASES)) . "\n";
    }
    close $ofh;
}


####
sub score {
    my ($model, $seq, $atg_index) = @_;

    my @window = &window($seq, $atg_index);
    my $score = 0;
    for my $k (0..$#window) {
        if (defined $window[$k]) {
            $score += $model->[$k]->{$window[$k]};
        }
    }
    return($score);
}


####
sub write_logo_svg {
    my ($pos_counts, $model, $num_starts, $filename) = @_;

    # Two panels over -6..-1, ATG, +4..+5 (the ATG itself is labeled, not drawn):
    #  top:    information content (bits) of the training start contexts
    #  bottom: log2-odds vs internal in-frame ATGs, positive above and negative below the axis
    # Both y scales fit the data (at least 0.5 bits / 1 log2 unit).

    my %color = (A => '#109648', C => '#255C99', G => '#F7B32B', T => '#D62839');
    my @labels = &position_labels();
    my @cols = (@labels[0..$UPSTREAM-1], 'A', 'T', 'G', @labels[$UPSTREAM..$#labels]);
    my $ncol = scalar(@cols);
    my ($colw, $left, $top_h, $bot_h, $gap) = (40, 56, 150, 220, 56);

    # per-column stacks
    my (@ic_stack, @up_stack, @dn_stack);
    my ($max_ic, $max_lo) = (0.5, 1);
    for my $k (0..$UPSTREAM + $DOWNSTREAM - 1) {
        my $n = 0; $n += $pos_counts->[$k]->{$_} for @BASES;
        my %f = map { $_ => ($pos_counts->[$k]->{$_} + 0.5) / ($n + 2) } @BASES;
        my $entropy = 0; $entropy -= $f{$_} * log($f{$_}) / log(2) for @BASES;
        my $ic = 2 - $entropy;
        $ic_stack[$k] = [ map { [$_, $f{$_} * $ic] } sort { $f{$a} <=> $f{$b} } @BASES ];
        $max_ic = $ic if $ic > $max_ic;

        my %lo = map { $_ => $model->[$k]->{$_} } @BASES;
        $up_stack[$k] = [ map { [$_, $lo{$_}] } sort { $lo{$a} <=> $lo{$b} } grep { $lo{$_} > 0 } @BASES ];
        $dn_stack[$k] = [ map { [$_, -$lo{$_}] } sort { $lo{$b} <=> $lo{$a} } grep { $lo{$_} < 0 } @BASES ];
        my ($u, $d) = (0, 0);
        $u += $_->[1] for @{$up_stack[$k]};
        $d += $_->[1] for @{$dn_stack[$k]};
        $max_lo = max($max_lo, $u, $d);
    }
    my $ic_top = &_nice_ceiling($max_ic);
    my $lo_top = &_nice_ceiling($max_lo);

    my $top_y0 = 40 + $top_h;                 # baseline of the top panel
    my $bot_axis = $top_y0 + $gap + $bot_h / 2;
    my $height = $top_y0 + $gap + $bot_h + 40;
    my $width = $left + $ncol * $colw + 20;
    my $ic_px = $top_h / $ic_top;
    my $lo_px = ($bot_h / 2) / $lo_top;

    my @svg;
    push (@svg, qq(<svg xmlns="http://www.w3.org/2000/svg" width="$width" height="$height" font-family="Helvetica,Arial,sans-serif">));
    push (@svg, qq(<rect width="100%" height="100%" fill="white"/>));
    push (@svg, qq(<text x="$left" y="22" font-size="14">Start-codon context of $num_starts training ORF starts (bits)</text>));
    push (@svg, sprintf(qq(<text x="%d" y="%d" font-size="14">Log2-odds vs in-frame ATGs inside the training ORFs</text>), $left, $top_y0 + $gap - 14));

    # bold Helvetica/Arial glyphs at 100px: cap height ~72px; advance widths below
    my %glyph_w = (A => 72, C => 72, G => 78, T => 61);
    my $letter = sub {
        my ($b, $x, $y_base, $h) = @_;
        if ($h < 0.5) { return; }
        push (@svg, sprintf(qq(<text x="0" y="0" font-size="100" font-weight="bold" fill="%s" transform="translate(%.1f,%.1f) scale(%.4f,%.4f)">%s</text>),
                            $color{$b}, $x + 2, $y_base, ($colw - 4) / $glyph_w{$b}, $h / 72, $b));
    };

    for my $c (0..$ncol-1) {
        my $x = $left + $c * $colw;
        my $is_atg = ($c >= $UPSTREAM && $c < $UPSTREAM + 3);
        push (@svg, sprintf(qq(<text x="%.1f" y="%.1f" font-size="%d" text-anchor="middle"%s>%s</text>),
                            $x + $colw / 2, $height - 14, ($is_atg ? 16 : 12), ($is_atg ? ' font-weight="bold" fill="#888888"' : ''), $cols[$c]));
        if ($is_atg) { next; }
        my $k = ($c < $UPSTREAM) ? $c : $c - 3;

        my $y = $top_y0;
        foreach my $e (@{$ic_stack[$k]}) {
            my $h = $e->[1] * $ic_px;
            $letter->($e->[0], $x, $y, $h);
            $y -= $h;
        }
        my $yup = $bot_axis;
        foreach my $e (@{$up_stack[$k]}) {
            my $h = $e->[1] * $lo_px;
            $letter->($e->[0], $x, $yup, $h);
            $yup -= $h;
        }
        my $ydn = $bot_axis;
        foreach my $e (@{$dn_stack[$k]}) {
            my $h = $e->[1] * $lo_px;
            $ydn += $h;
            $letter->($e->[0], $x, $ydn, $h);
        }
    }

    # axes and ticks
    my $xr = $left + $ncol * $colw;
    push (@svg, qq(<line x1="$left" y1="$top_y0" x2="$xr" y2="$top_y0" stroke="black"/>));
    push (@svg, qq(<line x1="$left" y1="$bot_axis" x2="$xr" y2="$bot_axis" stroke="black"/>));
    push (@svg, sprintf(qq(<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="black"/>), $left - 4, $top_y0 - $top_h, $left - 4, $top_y0));
    push (@svg, sprintf(qq(<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="black"/>), $left - 4, $bot_axis - $bot_h / 2, $left - 4, $bot_axis + $bot_h / 2));
    foreach my $t ([0, $top_y0], [$ic_top / 2, $top_y0 - $top_h / 2], [$ic_top, $top_y0 - $top_h],
                   ["+$lo_top", $bot_axis - $bot_h / 2], [0, $bot_axis], ["-$lo_top", $bot_axis + $bot_h / 2]) {
        push (@svg, sprintf(qq(<text x="%d" y="%.1f" font-size="11" text-anchor="end">%s</text>), $left - 8, $t->[1] + 4, $t->[0]));
    }
    push (@svg, "</svg>");

    open(my $ofh, ">$filename") or confess "Error, cannot write to $filename";
    print $ofh join("\n", @svg) . "\n";
    close $ofh;
}


####
sub _nice_ceiling {
    my ($v) = @_;
    foreach my $step (0.25, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 8, 10) {
        if ($v <= $step) { return($step); }
    }
    return(int($v) + 1);
}

1;
