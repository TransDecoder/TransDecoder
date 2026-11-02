#!/usr/bin/env perl

# Target-decoy q-values of candidate ORFs (TransDecoder.Predict ORF selection).
#
# Decoy mode: each transcript is dinucleotide-shuffled (--shuffles copies; 0 = enough for ~20 Mb of decoy
# sequence, at most 50), decoy ORFs are found by TransDecoder.LongOrfs with the settings of the
# original run (longorfs.params.json) and scored with the run's Markov model. Eligible ORFs (in-frame
# score > 0 and frame contrast > 0, contrast = score_1 - max(score_2, score_3)) get q-values within
# strata of ORF length, transcript length and ORF type (PerlLib/Decoys.pm). Writes
# <workdir>/decoy_qvalues.tsv, identical to PyTransDecoder's for the same inputs and seed.
#
# Preset mode (--preset): no decoys; an eligible ORF gets q = 0 when its frame contrast reaches the
# preset's threshold for its stratum at --max_qvalue, else q = 1 (written to <workdir>/preset_qvalues.tsv).

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Carp;
use List::Util qw(max min);
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use Fasta_reader;
use Decoys;
use Presets;

my $DECOY_TARGET_NT = 20_000_000;
my $DECOY_MAX_SHUFFLES = 50;

my $usage = <<__EOUSAGE__;

usage: $0 --transcripts <fasta> --workdir <transdecoder_dir> --cds_scores <file> --hexamer_scores <file>
          [--shuffles 0] [--seed 42]
       $0 --transcripts <fasta> --workdir <transdecoder_dir> --cds_scores <file> --preset <dir> --max_qvalue <q>

__EOUSAGE__
    ;

my ($transcripts_file, $workdir, $cds_scores_file, $hexamer_scores_file, $preset_dir, $max_qvalue);
my $n_shuffles = 0;
my $seed = 42;

&GetOptions('transcripts=s' => \$transcripts_file,
            'workdir=s' => \$workdir,
            'cds_scores=s' => \$cds_scores_file,
            'hexamer_scores=s' => \$hexamer_scores_file,
            'shuffles=i' => \$n_shuffles,
            'seed=i' => \$seed,
            'preset=s' => \$preset_dir,
            'max_qvalue=f' => \$max_qvalue,
    );

if (@ARGV) {
    die "Error, don't understand options: @ARGV";
}
unless ($transcripts_file && $workdir && $cds_scores_file && ($hexamer_scores_file || ($preset_dir && defined $max_qvalue))) {
    die $usage;
}

main: {

    my %transcripts;
    {
        my $fasta_reader = new Fasta_reader($transcripts_file);
        while (my $seq_obj = $fasta_reader->next()) {
            $transcripts{ $seq_obj->get_accession() } = uc $seq_obj->get_sequence();
        }
    }

    my %targets = &target_strata("$workdir/longest_orfs.cds", $cds_scores_file, \%transcripts);

    if ($preset_dir) {
        my $preset = new Presets($preset_dir);
        my $level = $preset->level($max_qvalue);
        my $outfile = "$workdir/preset_qvalues.tsv";
        open (my $ofh, ">$outfile") or confess "Error, cannot write $outfile";
        print $ofh "#orf_id\tstratum\tcontrast\tq\n";
        my $n_pass = 0;
        foreach my $orf_id (sort keys %targets) {
            my ($k, $c) = @{$targets{$orf_id}};
            my $q = ($c >= $preset->threshold($k, $max_qvalue)) ? 0 : 1;
            $n_pass++ unless $q;
            printf $ofh "%s\t%s\t%.2f\t%s\n", $orf_id, $k, $c, $q;
        }
        close $ofh;
        printf STDERR "-preset thresholds at q %g (stored level %g): %d of %d eligible candidates pass\n",
            $max_qvalue, $level, $n_pass, scalar(keys %targets);
        exit(0);
    }

    my %params = &read_longorfs_params($workdir);

    my $total_nt = 0;
    foreach my $seq (values %transcripts) {
        $total_nt += length($seq);
    }
    unless ($n_shuffles) {
        my $per = int(($DECOY_TARGET_NT + max($total_nt, 1) - 1) / max($total_nt, 1));  # ceiling
        $n_shuffles = max(1, min($DECOY_MAX_SHUFFLES, $per));
    }

    ## decoy transcripts
    my $decoy_dir = "$workdir/decoys";
    unless (-d $decoy_dir) {
        mkdir($decoy_dir) or confess "Error, cannot mkdir $decoy_dir";
    }
    my $decoy_fasta = "$decoy_dir/decoys.fa";
    my %decoy_tx_len;
    open (my $dfh, ">$decoy_fasta") or confess "Error, cannot write $decoy_fasta";
    foreach my $tid (sort keys %transcripts) {
        my $seq = $transcripts{$tid};
        for my $k (0 .. $n_shuffles - 1) {
            my $shuf = &Decoys::dinucleotide_shuffle($seq, &Decoys::new_rng("$seed:$tid:$k"));
            my $decoy_id = "decoy${k}_$tid";
            $decoy_tx_len{$decoy_id} = length($shuf);
            print $dfh ">$decoy_id\n$shuf\n";
        }
    }
    close $dfh;

    ## decoy ORFs: the LongOrfs finder with the run's settings
    my $cmd = "$FindBin::Bin/TransDecoder.LongOrfs -t $decoy_fasta -O $decoy_dir -m $params{min_protein_length} "
        . " -G $params{genetic_code} "
        . ($params{strand_specific} ? " -S " : "")
        . ($params{complete_orfs_only} ? " --complete_orfs_only " : "")
        . " > $decoy_dir/longorfs.log 2>&1";
    &process_cmd($cmd);
    my $decoy_cds = "$decoy_dir/decoys.fa.transdecoder_dir/longest_orfs.cds";
    my $decoy_scores = "$decoy_cds.scores";
    &process_cmd("$FindBin::Bin/score_CDS_likelihood_all_6_frames.pl $decoy_cds $hexamer_scores_file "
                 . " --genetic_code $params{genetic_code} > $decoy_scores");

    my %decoys_eligible = &target_strata($decoy_cds, $decoy_scores, undef, \%decoy_tx_len);
    my $n_decoy_orfs = &count_fasta_entries($decoy_cds);
    my @decoys = values %decoys_eligible;

    my %q = &Decoys::qvalues(\%targets, \@decoys, $n_shuffles);

    my $outfile = "$workdir/decoy_qvalues.tsv";
    open (my $ofh, ">$outfile") or confess "Error, cannot write $outfile";
    print $ofh "#orf_id\tstratum\tcontrast\tq\n";
    foreach my $orf_id (sort keys %q) {
        my ($k, $c) = @{$targets{$orf_id}};
        printf $ofh "%s\t%s\t%.2f\t%s\n", $orf_id, $k, $c, &Decoys::format_g5($q{$orf_id});
    }
    close $ofh;

    printf STDERR "-decoys: %d shuffle(s) per transcript (seed %d), %d decoy ORFs (%d eligible); %d eligible candidates (%s)\n",
        $n_shuffles, $seed, $n_decoy_orfs, scalar(@decoys), scalar(keys %targets), $outfile;

    exit(0);
}


####
sub target_strata {
    # (stratum, frame contrast) of each eligible ORF: in-frame score > 0 and contrast > 0
    my ($cds_file, $scores_file, $transcripts_href, $tx_len_href) = @_;

    my (%type, %tid, %codons);
    {
        my $fasta_reader = new Fasta_reader($cds_file);
        while (my $seq_obj = $fasta_reader->next()) {
            my $acc = $seq_obj->get_accession();
            my $header = $seq_obj->get_header();
            my @x = split(/\s+/, $header);
            # header: <orf_id> type:<type> <transcript>:<lend>-<rend>(<strand>)
            my ($t) = $x[1] =~ /^type:(\S+)/ or confess "Error, cannot parse ORF type from: $header";
            my $loc = $x[2];
            $loc =~ s/:[^:]*$//;
            $type{$acc} = $t;
            $tid{$acc} = $loc;
            $codons{$acc} = int(length($seq_obj->get_sequence()) / 3);
        }
    }

    my %out;
    open (my $fh, $scores_file) or confess "Error, cannot open $scores_file";
    while (<$fh>) {
        chomp;
        next if /^\#/;
        my ($acc, $order, $len, @s) = split(/\t/);
        next unless exists $type{$acc};
        my $c = $s[0] - max($s[1], $s[2]);
        if ($s[0] > 0 && $c > 0) {
            my $tx_len = ($transcripts_href)
                ? length($transcripts_href->{ $tid{$acc} } // "")
                : ($tx_len_href->{ $tid{$acc} } // 0);
            $out{$acc} = [ &Decoys::stratum($codons{$acc}, $tx_len, $type{$acc} eq 'complete'), $c ];
        }
    }
    close $fh;

    return(%out);
}

####
sub read_longorfs_params {
    my ($workdir) = @_;
    my $params_file = "$workdir/longorfs.params.json";
    my %params;
    if (-s $params_file) {
        open (my $fh, $params_file) or confess "Error, cannot open $params_file";
        my $json = join("", <$fh>);
        close $fh;
        ($params{min_protein_length}) = $json =~ /"min_protein_length":\s*(\d+)/;
        ($params{genetic_code}) = $json =~ /"genetic_code":\s*"([^"]+)"/;
        ($params{strand_specific}) = $json =~ /"strand_specific":\s*(true|false)/;
        ($params{complete_orfs_only}) = $json =~ /"complete_orfs_only":\s*(true|false)/;
        foreach my $k (keys %params) {
            unless (defined $params{$k}) {
                confess "Error, cannot parse $k from $params_file";
            }
            $params{$k} = 1 if $params{$k} eq 'true';
            $params{$k} = 0 if $params{$k} eq 'false';
        }
    }
    else {
        # LongOrfs from before the params file: infer its settings from its outputs
        my ($min_aa, $minus, $all_complete) = (1e9, 0, 1);
        open (my $fh, "$workdir/longest_orfs.pep") or confess "Error, cannot open $workdir/longest_orfs.pep";
        my ($gc, $pep) = ("universal", "");
        while (<$fh>) {
            chomp;
            if (/^>/) {
                if ($pep ne "") { (my $p = $pep) =~ s/\*$//; $min_aa = min($min_aa, length($p)); }
                $pep = "";
                $minus = 1 if /\(-\)\s*$/;
                $all_complete = 0 unless / type:complete /;
                $gc = $1 if / gc:(\S+)/;
            }
            else {
                $pep .= $_;
            }
        }
        if ($pep ne "") { (my $p = $pep) =~ s/\*$//; $min_aa = min($min_aa, length($p)); }
        close $fh;
        %params = (min_protein_length => $min_aa, genetic_code => $gc,
                   strand_specific => ($minus ? 0 : 1), complete_orfs_only => $all_complete);
        print STDERR "-warning: no longorfs.params.json; decoy ORF settings inferred: "
            . join(", ", map { "$_=$params{$_}" } sort keys %params) . "\n";
    }
    # LongOrfs' default name for the standard code
    $params{genetic_code} = "Universal" if lc($params{genetic_code}) eq 'universal';
    return(%params);
}

####
sub count_fasta_entries {
    my ($file) = @_;
    my $n = 0;
    open (my $fh, $file) or confess "Error, cannot open $file";
    while (<$fh>) {
        $n++ if /^>/;
    }
    close $fh;
    return($n);
}

####
sub process_cmd {
    my ($cmd) = @_;
    print STDERR "CMD: $cmd\n";
    my $ret = system($cmd);
    if ($ret) {
        confess "Error, cmd: $cmd died with ret $ret";
    }
}
