#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use Fasta_reader;
use GFF3_utils;
use Carp;
use Nuc_translator;

my $usage = "\n\nusage: $0 gff3_file genome_db [prot|CDS|cDNA|gene,default=prot] [flank=0]\n\n";

my $gff3_file = $ARGV[0] or die $usage;
my $fasta_db = $ARGV[1] or die $usage;
my $seq_type = $ARGV[2] || "prot";
my $flank = $ARGV[3] || 0;

my ($upstream_flank, $downstream_flank) = (0,0);

if ($flank) {
	if ($flank =~ /:/) {
		($upstream_flank, $downstream_flank) = split (/:/, $flank);
	}
	else {
		($upstream_flank, $downstream_flank) = ($flank, $flank);
	}
}

if ($upstream_flank < 0 || $downstream_flank < 0) {
	die $usage;
}



unless ($seq_type =~ /^(prot|CDS|cDNA|gene)$/) {
    die "Error, don't understand sequence type [$seq_type]\n\n$usage";
}


## read genome
my $fasta_reader = new Fasta_reader($fasta_db);
my %genome = $fasta_reader->retrieve_all_seqs_hash();


my $gene_obj_indexer_href = {};

## associate gene identifiers with contig id's.
my $contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer_href);

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
				if ($upstream_flank || $downstream_flank) {
					$seq = &add_flank($seq, $upstream_flank, $downstream_flank, $model_lend, $model_rend, $orientation, \$genome_seq);
				}
			}
            elsif ($seq_type eq "cDNA") {
                $seq = $isoform->get_cDNA_sequence();
				if ($upstream_flank || $downstream_flank) {
					$seq = &add_flank($seq, $upstream_flank, $downstream_flank, $gene_lend, $gene_rend, $orientation, \$genome_seq);
				}
			}
            elsif ($seq_type eq "gene" && $counter == 1) {
                $seq = $isoform->get_gene_sequence();
				if ($upstream_flank || $downstream_flank) {
					$seq = &add_flank($seq, $upstream_flank, $downstream_flank, $gene_lend, $gene_rend, $orientation, \$genome_seq);
				}
			}
            
            unless ($seq) {
                print STDERR "-warning, no $seq_type sequence for $isoform_id\n";
                next;
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
            
            print ">$isoform_id $gene_id $locus_string $com_name $asmbl_id:$model_lend-$model_rend($orientation)\n$seq\n";
        }
    }
}


exit(0);


####
sub add_flank {
	my ($seq, $upstream_flank, $downstream_flank, $lend, $rend, $orientation, $genome_seq_ref) = @_;
	
	my $far_left = ($orientation eq '+') ? $lend - $upstream_flank : $lend - $downstream_flank;
	
	if ($far_left < 1) { $far_left = 1; }
	
	my $flank_right = ($orientation eq '+') ? $downstream_flank : $upstream_flank;

	my $left_seq = substr($$genome_seq_ref, $far_left - 1, $lend - $far_left);

	my $right_seq = substr($$genome_seq_ref, $rend, $flank_right);
	
	if ($orientation eq '+') {
		return (lc($left_seq) . uc($seq) . lc($right_seq));
	}
	else {
		return (lc(&reverse_complement($right_seq)) . uc($seq) . lc(&reverse_complement($left_seq)));
	}
}


