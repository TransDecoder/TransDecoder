#!/usr/bin/env perl

=pod

=head1 NAME

extract_FL_subset.pl

=head1 USAGE

Required:

 --transcripts, -t <string>                            transcripts.fasta

 --gff3, -g <string>                            transcripts.fasta.transdecoder.gff3

=cut


use strict;
use warnings;
use FindBin;
use Pod::Usage;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use Data::Dumper;
use File::Basename;

use lib ("$FindBin::RealBin/../PerlLib");
use Gene_obj;


my $UTIL_DIR = "$FindBin::RealBin";

my $help;

my $transcripts_file;
my $gff3_file;

&GetOptions( 'transcripts|t=s' => \$transcripts_file,
             'gff3|g=s' => \$gff3_file,

             'h' => \$help,
    );



pod2usage(-verbose => 2, -output => \*STDERR) if ($help || ! ($transcripts_file && $gff3_file));

if (@ARGV) {
    die "Error, don't understand options: @ARGV";
}


main: {
    
    my $output_prefix = basename($transcripts_file);

    # get the FL entries:
    my $FL_accs_file = "$output_prefix.complete_only";
    my $cmd = "$UTIL_DIR/get_FL_accs.pl $gff3_file > $FL_accs_file";
    &process_cmd($cmd);

    # index the current gff file:
    $cmd = "$UTIL_DIR/index_gff3_files_by_isoform.pl $gff3_file";
    &process_cmd($cmd);
        
    # retrieve the best entries:
    my $complete_gff3_filename = $FL_accs_file . ".gff3";
    $cmd = "$UTIL_DIR/gene_list_to_gff.pl $FL_accs_file $gff3_file.inx > $FL_accs_file.gff3";
    &process_cmd($cmd);
    
    
    ##############################
    ## Generate the final outputs.
    ##############################
    
    {
                        
        ## write final outputs:
        $gff3_file = "$FL_accs_file.gff3";
        

        ## make a BED file for viewing in IGV
        my $bed_file = $gff3_file;
        $bed_file =~ s/\.gff3$/\.bed/;
        $cmd = "$UTIL_DIR/gff3_file_to_bed.pl $gff3_file > $bed_file";
        &process_cmd($cmd);
        
    
        # make a peptide file:
        my $best_pep_file = $gff3_file;
        $best_pep_file =~ s/\.gff3$/\.pep/;
        $cmd = "$UTIL_DIR/gff3_file_to_proteins.pl --gff3 $gff3_file --fasta $transcripts_file --seqType prot > $best_pep_file";
        &process_cmd($cmd);
        
        

        # make a CDS file:
        my $best_cds_file = $best_pep_file;
        $best_cds_file =~ s/\.pep$/\.cds/;
        $cmd = "$UTIL_DIR/gff3_file_to_proteins.pl --gff3 $gff3_file --fasta $transcripts_file --seqType CDS > $best_cds_file";
        &process_cmd($cmd);
        
        # make a CDS file:
        my $best_cdna_file = $best_pep_file;
        $best_cdna_file =~ s/\.pep$/\.mRNA/;
        $cmd = "$UTIL_DIR/gff3_file_to_proteins.pl --gff3 $gff3_file --fasta $transcripts_file --seqType cDNA > $best_cdna_file";
        &process_cmd($cmd);
    
    }
    
    print STDERR "Done. See output files $FL_accs_file.\*\n\n\n";
    
    
    
    exit(0);
}


####
sub process_cmd {
	my ($cmd) = @_;

	print "CMD: $cmd\n";
	my $ret = system($cmd);

	if ($ret) {
		die "Error, cmd: $cmd died with ret $ret";
	}
	
	return;

}

