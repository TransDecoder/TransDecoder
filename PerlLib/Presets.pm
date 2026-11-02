package Presets;

# Presets: coding model, start model and ORF acceptance thresholds learned from a large reference
# transcriptome, for inputs too small (or too partial) to learn from. Same files and lookups as
# PyTransDecoder (core/presets.py); bundled presets are in PerlLib/presets/<name>.
#
# A preset directory holds
#   preset.json            name, source, genetic code, GC %, long-ORF auto-accept length (nt or null),
#                          decoy strata, q levels
#   hexamer.scores.gz      the reference run's Markov coding model
#   training_orfs.cds.gz   the reference run's training ORFs
#   kozak.tsv              the reference run's start-context model (optional)
#   thresholds.tsv.gz      per stratum and q level, the smallest frame contrast accepted
# (uncompressed files are read too).

use strict;
use warnings;
use Carp;
use FindBin;
use File::Basename;

my $INF = 9**9**9;
our $BUNDLED_DIR = dirname(__FILE__) . "/presets";


####
sub bundled_presets {
    my @names;
    if (opendir(my $dh, $BUNDLED_DIR)) {
        @names = sort grep { -e "$BUNDLED_DIR/$_/preset.json" } readdir($dh);
        closedir $dh;
    }
    return(@names);
}

####
sub resolve {
    # a preset directory or a bundled preset name
    my ($name_or_dir) = @_;
    if (-e "$name_or_dir/preset.json") {
        return($name_or_dir);
    }
    if (-e "$BUNDLED_DIR/$name_or_dir/preset.json") {
        return("$BUNDLED_DIR/$name_or_dir");
    }
    my $avail = join(", ", &bundled_presets()) || "none installed";
    confess "Error, no preset '$name_or_dir': not a preset directory or a bundled preset ($avail)";
}

####
sub find_file {
    # <dir>/<name>.gz or <dir>/<name>
    my ($dir, $name) = @_;
    foreach my $f ("$dir/$name.gz", "$dir/$name") {
        return($f) if -e $f;
    }
    return(undef);
}

####
sub open_file {
    my ($file) = @_;
    my $fh;
    if ($file =~ /\.gz$/) {
        open ($fh, "-|", "gzip", "-dc", $file) or confess "Error, cannot read $file";
    }
    else {
        open ($fh, $file) or confess "Error, cannot open $file";
    }
    return($fh);
}

####
sub copy_plain {
    # copy a preset file, decompressing .gz
    my ($src, $dest) = @_;
    my $in = &open_file($src);
    open (my $out, ">$dest") or confess "Error, cannot write $dest";
    while (<$in>) {
        print $out $_;
    }
    close $in;
    close $out;
}

####
sub new {
    my ($packagename, $dir) = @_;
    $dir = &resolve($dir);
    open (my $fh, "$dir/preset.json") or confess "Error, cannot open $dir/preset.json";
    my $json = join("", <$fh>);
    close $fh;
    my %meta;
    foreach my $key (qw(name source genetic_code strata)) {
        ($meta{$key}) = $json =~ /"$key":\s*"([^"]*)"/;
    }
    ($meta{gc_pct}) = $json =~ /"gc_pct":\s*([\d.]+)/;
    my ($auto) = $json =~ /"auto_accept_nt":\s*([\d.]+|null)/;
    $meta{auto_accept_nt} = (defined($auto) && $auto ne 'null') ? $auto : undef;

    my $self = { dir => $dir,
                 meta => \%meta,
                 hexamer_scores => &find_file($dir, "hexamer.scores"),
                 training_orfs => &find_file($dir, "training_orfs.cds"),
                 kozak => (-e "$dir/kozak.tsv") ? "$dir/kozak.tsv" : undef,
                 thr => {},          # level => { stratum => min contrast }
                 order => {},        # level => [ strata in file order ]
    };
    my $tfile = &find_file($dir, "thresholds.tsv") or confess "Error, no thresholds.tsv in $dir";
    my $tfh = &open_file($tfile);
    while (<$tfh>) {
        chomp;
        next if /^\#/;
        my ($stratum, $q, $c, $n) = split(/\t/);
        $q += 0;
        $self->{thr}->{$q}->{$stratum} = ($c eq 'inf') ? $INF : $c + 0;
        push (@{$self->{order}->{$q}}, $stratum);
    }
    close $tfh;
    bless($self, $packagename);
    return($self);
}

####
sub auto_accept_nt {
    my ($self) = @_;
    my $a = $self->{meta}->{auto_accept_nt};
    return(defined($a) ? $a : 1000000000);
}

####
sub level {
    # largest stored q level <= q (the smallest level if q is below all of them)
    my ($self, $q) = @_;
    my @levels = sort {$a <=> $b} keys %{$self->{thr}};
    my @ok = grep { $_ <= $q + 1e-12 } @levels;
    return(@ok ? $ok[-1] : $levels[0]);
}

####
sub threshold {
    # minimum frame contrast at q for a stratum ("lenbin,txbin,True|False"); strata absent from the
    # reference use the nearest ORF-length bin with the same remaining keys, then any transcript-length
    # bin of the same type (first in file order among equals, as in PyTransDecoder)
    my ($self, $stratum, $q) = @_;
    my $level = $self->level($q);
    my $table = $self->{thr}->{$level};
    if (exists $table->{$stratum}) {
        return($table->{$stratum});
    }
    my @s = split(/,/, $stratum);
    my @all = @{$self->{order}->{$level}};
    my $rest = join(",", @s[1..$#s]);
    my @same = grep { my @k = split(/,/, $_); join(",", @k[1..$#k]) eq $rest } @all;
    if (! @same && scalar(@s) == 3) {
        @same = grep { my @k = split(/,/, $_); $k[2] eq $s[2] } @all;
    }
    if (! @same) {
        @same = @all;
    }
    my ($best, @best_key);
    foreach my $k (@same) {
        my @kk = split(/,/, $k);
        my @key = (abs($kk[0] - $s[0]), (scalar(@kk) > 1 ? abs($kk[1] - $s[1]) : 0));
        if (! defined($best) || $key[0] < $best_key[0] || ($key[0] == $best_key[0] && $key[1] < $best_key[1])) {
            $best = $k;
            @best_key = @key;
        }
    }
    return($table->{$best});
}

1;
