#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use Fasta_reader;
use GFF3_utils2;
use Carp;
use Nuc_translator;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);


my $usage = <<__EOUSAGE__;

####################################################
#
# Required:
#
#  --gff3 <string>          gff3 file
#
#  --fasta <string>         fasta file corresponding to gff3 file
#
##
#  Optional:
#
#  --seqType <string>        prot|CDS|cDNA|gene,  default=prot
#
#  --genetic_code  <string>   universal (default)
#                             Euplotes, Tetrahymena, Candida
#                             Acetabularia, Mitochondrial-Canonical
#                             Mitochondrial-Vertebrates, Mitochondrial-Arthropods
#                             Mitochondrial-Echinoderms, Mitochondrial-Molluscs
#                             Mitochondrial-Ascidians, Mitochondrial-Nematodes
#                             Mitochondrial-Platyhelminths,Mitochondrial-Yeasts
#                             Mitochondrial-Euascomycetes, Mitochondrial-Protozoans
#
#  --batch_size <int>        number of GFF3 sequence IDs to process per batch
#                             in streaming mode (default: 5000)
#
#  --tmpdir <string>         directory for temporary GFF3 batch files.  If not
#                             supplied, a temporary directory is created and
#                             cleaned up automatically.  User-supplied tmpdirs
#                             retain their generated batch subdirectory.
#
#  --legacy_memory_mode      use the original all-in-memory execution path
#                             (loads all FASTA records and indexes the full GFF3)
#
###################################################


__EOUSAGE__

    ;


my $gff3_file;
my $fasta_db;
my $seq_type = 'prot';
my $genetic_code = '';
my $batch_size = 5000;
my $tmpdir = "";
my $legacy_memory_mode = 0;

&GetOptions ( 'gff3=s' => \$gff3_file,
              'fasta=s' => \$fasta_db,
              'seqType=s' => \$seq_type,
              'genetic_code=s' => \$genetic_code,
              'batch_size=i' => \$batch_size,
              'tmpdir=s' => \$tmpdir,
              'legacy_memory_mode' => \$legacy_memory_mode,
    );

unless ($gff3_file && $fasta_db) {
    die $usage;
}

unless ($seq_type =~ /^(prot|CDS|cDNA|gene)$/) {
    die "Error, don't understand sequence type [$seq_type]\n\n$usage";
}

if ($genetic_code) {
    &Nuc_translator::use_specified_genetic_code($genetic_code);
}

if ($legacy_memory_mode) {
    run_legacy_memory_mode($gff3_file, $fasta_db, $seq_type);
}
else {
    run_streaming_mode($gff3_file, $fasta_db, $seq_type, $batch_size, $tmpdir);
}

exit(0);


sub run_legacy_memory_mode {
    my ($gff3_file, $fasta_db, $seq_type) = @_;

    ## read genome
    my $fasta_reader = new Fasta_reader($fasta_db);
    my %genome = $fasta_reader->retrieve_all_seqs_hash();


    my $gene_obj_indexer_href = {};

    ## associate gene identifiers with contig id's.
    my $contig_to_gene_list_href = &GFF3_utils2::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer_href);

    foreach my $asmbl_id (sort keys %$contig_to_gene_list_href) {

        my $genome_seq = $genome{$asmbl_id} or die "Error, no sequence for $asmbl_id";

        my @gene_ids = @{$contig_to_gene_list_href->{$asmbl_id}};

        foreach my $gene_id (@gene_ids) {
            my $gene_obj_ref = $gene_obj_indexer_href->{$gene_id};
            emit_gene_sequences($gene_obj_ref, \$genome_seq, $asmbl_id, $gene_id, $seq_type);
        }
    }
}


sub run_streaming_mode {
    my ($gff3_file, $fasta_db, $seq_type, $batch_size, $tmpdir) = @_;

    require File::Spec;
    require File::Temp;

    my $cleanup_tmpdir = 1;
    if ($tmpdir) {
        if (! -d $tmpdir) {
            die "Error: tmpdir does not exist: $tmpdir\n";
        }
        $cleanup_tmpdir = 0;
        $tmpdir = File::Temp::tempdir("gff3_to_proteins_XXXXXX", DIR => $tmpdir, CLEANUP => $cleanup_tmpdir);
    }
    else {
        $tmpdir = File::Temp::tempdir("gff3_to_proteins_XXXXXX", TMPDIR => 1, CLEANUP => $cleanup_tmpdir);
    }

    my @asmbl_ids = collect_asmbl_ids_from_gff3($gff3_file);
    my %wanted_asmbl = map { $_ => 1 } @asmbl_ids;

    my $fasta_offsets_href = build_fasta_offset_index($fasta_db, \%wanted_asmbl);

    foreach my $asmbl_id (@asmbl_ids) {
        unless (defined $fasta_offsets_href->{$asmbl_id}) {
            die "Error: no FASTA entry found for GFF3 sequence id: $asmbl_id\n";
        }
    }

    my $batch_files_aref = split_gff3_into_batches(
        $gff3_file,
        \@asmbl_ids,
        $batch_size,
        $tmpdir,
    );

    open(my $fasta_fh, "<", $fasta_db) or die "Error: cannot open FASTA $fasta_db";
    binmode($fasta_fh);

    foreach my $batch_gff3 (@$batch_files_aref) {

        next unless -s $batch_gff3;

        my $gene_obj_indexer_href = {};

        my $contig_to_gene_list_href = &GFF3_utils2::index_GFF3_gene_objs($batch_gff3, $gene_obj_indexer_href);

        foreach my $asmbl_id (sort keys %$contig_to_gene_list_href) {

            my $offset = $fasta_offsets_href->{$asmbl_id};
            defined $offset or die "Error: no FASTA offset for $asmbl_id\n";

            my $genome_seq = fetch_fasta_seq_by_offset($fasta_fh, $offset, $asmbl_id);

            my @gene_ids = @{$contig_to_gene_list_href->{$asmbl_id}};

            foreach my $gene_id (@gene_ids) {

                my $gene_obj_ref = $gene_obj_indexer_href->{$gene_id};

                emit_gene_sequences(
                    $gene_obj_ref,
                    \$genome_seq,
                    $asmbl_id,
                    $gene_id,
                    $seq_type,
                );

                delete $gene_obj_indexer_href->{$gene_id};
            }

            delete $contig_to_gene_list_href->{$asmbl_id};
            undef $genome_seq;
        }
    }

    close $fasta_fh;
}


sub collect_asmbl_ids_from_gff3 {
    my ($gff3_file) = @_;

    my %seen;

    open(my $fh, "<", $gff3_file) or die "Error: cannot open $gff3_file";

    while (my $line = <$fh>) {
        last if $line =~ /^##FASTA/;
        next if $line =~ /^\s*#/;
        next unless $line =~ /\S/;

        my @x = split(/\t/, $line);
        next unless scalar(@x) >= 9;

        my $asmbl_id = $x[0];

        $seen{$asmbl_id} = 1;
    }

    close $fh;

    return sort keys %seen;
}


sub build_fasta_offset_index {
    my ($fasta_file, $wanted_href) = @_;

    my %offsets;

    open(my $fh, "<", $fasta_file) or die "Error: cannot open FASTA $fasta_file";
    binmode($fh);

    while (1) {
        my $offset = tell($fh);
        my $line = <$fh>;
        last unless defined $line;

        if ($line =~ /^>(\S+)/) {
            my $acc = $1;
            if (! $wanted_href || exists $wanted_href->{$acc}) {
                $offsets{$acc} = $offset;
            }
        }
    }

    close $fh;

    return \%offsets;
}


sub fetch_fasta_seq_by_offset {
    my ($fh, $offset, $expected_acc) = @_;

    seek($fh, $offset, 0) or die "Error: seek failed at FASTA offset $offset";

    my $header = <$fh>;
    defined $header or die "Error: no FASTA header at offset $offset";

    $header =~ /^>(\S+)/
        or die "Error: expected FASTA header at offset $offset, got: $header";

    my $acc = $1;

    if (defined $expected_acc && $acc ne $expected_acc) {
        die "Error: FASTA index mismatch. Expected $expected_acc but found $acc at offset $offset";
    }

    my $seq = "";

    while (my $line = <$fh>) {
        last if $line =~ /^>/;
        chomp $line;
        $line =~ s/\s+//g;
        $seq .= $line;
    }

    return $seq;
}


sub split_gff3_into_batches {
    my ($gff3_file, $asmbl_ids_aref, $batch_size, $tmpdir) = @_;

    die "Error: batch_size must be positive\n" unless $batch_size && $batch_size > 0;

    my %asmbl_to_batch;
    my @batch_files;

    for (my $i = 0; $i < scalar(@$asmbl_ids_aref); $i++) {
        my $batch_no = int($i / $batch_size);
        my $asmbl_id = $asmbl_ids_aref->[$i];

        $asmbl_to_batch{$asmbl_id} = $batch_no;

        if (! defined $batch_files[$batch_no]) {
            $batch_files[$batch_no] = File::Spec->catfile($tmpdir, "batch_$batch_no.gff3");
        }
    }

    open(my $in_fh, "<", $gff3_file) or die "Error: cannot open $gff3_file";

    my %fh_cache;
    my @fh_order;
    my $max_open_fhs = 64;

    while (my $line = <$in_fh>) {
        last if $line =~ /^##FASTA/;
        next if $line =~ /^\s*#/;
        next unless $line =~ /\S/;

        my @x = split(/\t/, $line);
        next unless scalar(@x) >= 9;

        my $asmbl_id = $x[0];
        my $batch_no = $asmbl_to_batch{$asmbl_id};

        next unless defined $batch_no;

        my $out_fh = $fh_cache{$batch_no};

        if (! $out_fh) {
            if (scalar(keys %fh_cache) >= $max_open_fhs) {
                my $old_batch = shift @fh_order;
                close $fh_cache{$old_batch};
                delete $fh_cache{$old_batch};
            }

            my $batch_file = $batch_files[$batch_no];
            my $new_file = ! -e $batch_file;

            open($out_fh, ">>", $batch_file)
                or die "Error: cannot write $batch_file";

            if ($new_file) {
                print $out_fh "##gff-version 3\n";
            }

            $fh_cache{$batch_no} = $out_fh;
            push @fh_order, $batch_no;
        }

        print $out_fh $line;
    }

    close $in_fh;

    foreach my $batch_no (keys %fh_cache) {
        close $fh_cache{$batch_no};
    }

    return \@batch_files;
}


sub emit_gene_sequences {
    my ($gene_obj_ref, $genome_seq_sref, $asmbl_id, $gene_id, $seq_type) = @_;

    my %params;
    if ($seq_type eq "gene") {
        $params{unspliced_transcript} = 1;
    }

    $gene_obj_ref->create_all_sequence_types($genome_seq_sref, %params);

    my $counter = 0;
    foreach my $isoform ($gene_obj_ref, $gene_obj_ref->get_additional_isoforms()) {

        $counter++;

        my $orientation = $isoform->get_orientation();
        my ($model_lend, $model_rend) = sort {$a<=>$b} $isoform->get_model_span();
        my ($gene_lend, $gene_rend) = sort {$a<=>$b} $isoform->get_gene_span();

        my $isoform_id = $isoform->{Model_feat_name};

        my $seq = "";

        if ($seq_type eq "prot") {
            $seq = $isoform->get_protein_sequence();
        }
        elsif ($seq_type eq "CDS") {
            $seq = $isoform->get_CDS_sequence();
        }
        elsif ($seq_type eq "cDNA") {
            $seq = $isoform->get_cDNA_sequence();
        }
        elsif ($seq_type eq "gene" && $counter == 1) {
            $seq = $isoform->get_gene_sequence();
        }

        unless ($seq) {
            print STDERR "-warning, no $seq_type sequence for $isoform_id\n";
            next;
        }

        my $seqlen = length($seq);
        if ($seq =~ /\*$/) {
            $seqlen -= 1; # dont count stop codon
        }

        $seq =~ s/(\S{60})/$1\n/g; # make fasta format
        chomp $seq;

        my $com_name = $isoform->{com_name} || "";

        if ($com_name eq $isoform_id) { $com_name = ""; } # no sense in repeating it

        my $locus = $isoform->{pub_locus};
        my $model_locus = $isoform->{model_pub_locus};

        my $locus_string = "";
        if ($model_locus) {
            $locus_string .= $model_locus;
        }
        if ($locus) {
            $locus_string .= " $locus";
        }
        if ($locus_string) {
            $locus_string .= " "; # add spacer
        }

        #if ($seq_type eq 'prot' || $seq_type eq 'CDS') {  # this was a bad idea, just use the original id.
        #    $isoform_id = "cds.$isoform_id";
        #}

        print ">$isoform_id $gene_id $locus_string $com_name len:$seqlen $asmbl_id:$model_lend-$model_rend($orientation)\n$seq\n";
    }
}
