#!/usr/bin/env perl
# Decoy shuffles and strata: golden values shared with PyTransDecoder (tests/test_decoys.py).
use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Decoys;

my $fail = 0;
sub check { my ($ok, $what) = @_; print(($ok ? "ok" : "FAIL"), "\t$what\n"); $fail++ unless $ok; }

&check(&Decoys::fnv1a32("") == 0x811C9DC5 && &Decoys::fnv1a32("a") == 0xE40C292C, "fnv1a32");
my $g = &Decoys::new_rng("42:x:0");
&check(join(",", map { &Decoys::next32($g) } 1..5) eq "330423170,325826645,3944918289,989894671,3051725679", "xoshiro128**");
&check(&Decoys::dinucleotide_shuffle("ATGCGCGTATTAGCCGATAGGCTAGCTAACG", &Decoys::new_rng("42:x:0"))
       eq "ATTATAATAGCGACGCCTAGGCGCGTGCTAG", "dinucleotide shuffle");
&check(&Decoys::stratum(49, 400, 1) eq "-1,0,True" && &Decoys::stratum(50, 500, 0) eq "0,1,False"
       && &Decoys::stratum(56, 999, 1) eq "1,1,True" && &Decoys::stratum(100, 4000, 1) eq "7,4,True", "strata");
my %t = map { ("t$_" => ["s", $_ + 1]) } (0..9);
my %q = &Decoys::qvalues(\%t, [["s", 1.5], ["s", 2.5]], 1);
&check($q{t9} == 0 && $q{t2} == 0 && $q{t1} == 1/9 && $q{t0} == 2/10, "q-values");
exit($fail ? 1 : 0);
