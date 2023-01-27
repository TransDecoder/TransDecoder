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
###################################################


__EOUSAGE__

    ;


my $gff3_file;
my $fasta_db;
my $seq_type = 'prot';
my $genetic_code = '';

&GetOptions ( 'gff3=s' => \$gff3_file,
              'fasta=s' => \$fasta_db,
              'seqType=s' => \$seq_type,
              'genetic_code=s' => \$genetic_code,
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
		
        my %params;
        if ($seq_type eq "gene") {
            $params{unspliced_transcript} = 1;
        }
        
        $gene_obj_ref->create_all_sequence_types(\$genome_seq, %params);
        
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
}


exit(0);


