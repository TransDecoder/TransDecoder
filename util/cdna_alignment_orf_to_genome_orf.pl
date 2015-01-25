#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/PerlLib");
use Gene_obj;
use GFF3_utils;
use Data::Dumper;
use Fasta_reader;

my $usage = "\nusage: $0 cdna_orfs.genes.gff3 cdna_genome.alignments.gff3 cdna.fasta\n\n";

my $cdna_orfs_gff3 = $ARGV[0] or die $usage;
my $cdna_genome_gff3 = $ARGV[1] or die $usage;
my $cdna_fasta = $ARGV[2] or die $usage;


my $WARNING_COUNT = 0; # count those orfs that appear to be on strand opposite from the transcribed strand.

main: {

    my %cdna_seq_lengths = &parse_cdna_seq_lengths($cdna_fasta);
    
    my %orf_counter;

    my %cdna_acc_to_transcript_structure = &parse_transcript_alignment_info($cdna_genome_gff3);

    ## parse ORFs on cDNAs

    my $gene_obj_indexer_href = {};
    ## associate gene identifiers with contig id's.
    my $contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($cdna_orfs_gff3, $gene_obj_indexer_href);

    foreach my $asmbl_id (sort keys %$contig_to_gene_list_href) {
    
        my @gene_ids = @{$contig_to_gene_list_href->{$asmbl_id}};
    
        foreach my $gene_id (@gene_ids) {
            my $gene_obj_ref = $gene_obj_indexer_href->{$gene_id};
            
            my $asmbl_id = $gene_obj_ref->{asmbl_id};
            
            
            ## pasa stuff
            
            if ($asmbl_id =~ /(S\d+)_(asmbl_\d+)/) { 
                
                my $subcluster = $1;
                $asmbl_id = $2;
            }
            
            my $transcript_struct = $cdna_acc_to_transcript_structure{$asmbl_id} or die "Error, no cdna struct for $asmbl_id";
            
            #print Dumper($transcript_struct) . $gene_obj_ref->toString();

            eval {
                my $new_orf_gene = &place_orf_in_cdna_alignment_context($transcript_struct, $gene_obj_ref, \%cdna_seq_lengths);
                
                if ($new_orf_gene) {
                    
                    my $orf_count = $orf_counter{$asmbl_id}++;
                    $new_orf_gene->{asmbl_id} = $transcript_struct->{contig};
                    #$new_orf_gene->{TU_feat_name} = "t.$asmbl_id.$orf_count";
                    #$new_orf_gene->{Model_feat_name} = "m.$asmbl_id.$orf_count";
                    $new_orf_gene->{com_name} = "ORF";
                    
                    $new_orf_gene->{TU_feat_name} = $gene_id;
                    $new_orf_gene->{Model_feat_name} = $gene_obj_ref->{Model_feat_name};
                                        

                    print $new_orf_gene->to_GFF3_format(source => "transdecoder") . "\n";
                    
                }
            };

            if ($@) {
                
                print STDERR "Error occurred.\n";
                
                print STDERR Dumper($transcript_struct);
                print STDERR $gene_obj_ref->toString();
                print STDERR "$@";
                die;
            }
            
        }
    }

    exit(0);

}

####
sub parse_transcript_alignment_info {
    my ($cdna_align_gff3) = @_;

    my %cdna_alignments;

    open (my $fh, $cdna_align_gff3) or die "Error, cannot open file $cdna_align_gff3";
    while (<$fh>) {
        unless (/\w/) { next; }

        my @x = split(/\t/);
        my $contig = $x[0];
        my $lend = $x[3];
        my $rend = $x[4];
        my $orient = $x[6];
        my $info = $x[8];

        $info =~ /Target=(\S+)/ or die "Error, cannot parse ID from $info";
        my $asmbl = $1;
        
        if (my $struct = $cdna_alignments{$asmbl}) {
            push (@{$struct->{coords}}, [$lend, $rend]);
        }
        else {
            # first time
            my $struct = { asmbl => $asmbl,
                           contig => $contig,
                           
                           coords => [ 
                                       [$lend, $rend]
                                     ],
                           orient => $orient,
                                   };

            $cdna_alignments{$asmbl} = $struct;
        }

    }

    close $fh;

    return(%cdna_alignments);
}


####
sub place_orf_in_cdna_alignment_context {
    my ($transcript_struct, $orf_gene_obj, $cdna_seq_lengths_href) = @_;

    my $trans_seq_length = $cdna_seq_lengths_href->{ $transcript_struct->{asmbl} } or die "Error, no length for " . Dumper($transcript_struct);
    


    ## unwrap the gene
    my @cds_coords;
    my $orf_orient = $orf_gene_obj->get_orientation();
        
    foreach my $exon ($orf_gene_obj->get_exons()) {
        
        if (my $cds_exon = $exon->get_CDS_obj()) {

            my ($lend, $rend) = sort {$a<=>$b} $cds_exon->get_coords();
            push (@cds_coords, [$lend, $rend]);
        }
    }

    @cds_coords = sort {$a->[0]<=>$b->[0]} @cds_coords;

    my $cds_span_lend = $cds_coords[0]->[0];
    my $cds_span_rend = $cds_coords[$#cds_coords]->[1];
    
    if ($cds_span_rend > $trans_seq_length) {
        $cds_span_rend = $trans_seq_length;
    }
    
    
    my @exon_coords = @{$transcript_struct->{coords}};
    @exon_coords = sort {$a->[0]<=>$b->[0]} @exon_coords;
    my $trans_orient = $transcript_struct->{orient};

    ## examine each potential context of orf in alignment.
    
    my ($cds_genome_lend, $cds_genome_rend);
    my $transcribed_orient;

    if ($orf_orient eq '+') {


        if ($trans_orient eq '+') { 

            $cds_genome_lend = &from_cdna_lend($cds_span_lend, \@exon_coords);
            $cds_genome_rend = &from_cdna_lend($cds_span_rend, \@exon_coords);
            $transcribed_orient = '+';

        }
    
        elsif ($trans_orient eq '-') {
            
            $cds_genome_lend = &from_cdna_rend($cds_span_rend, \@exon_coords);
            $cds_genome_rend = &from_cdna_rend($cds_span_lend, \@exon_coords);
            $transcribed_orient = '-';

        }
        
    }
    
    else {
        ## orf orient is '-'
        if (scalar(@exon_coords) > 1) {
            # any correct ORF should be in the '+' orientation here.... must be a false positive orf or transcript structure is wrong
            $WARNING_COUNT++;
            print STDERR "Warning [$WARNING_COUNT], shouldn't have a minus-strand ORF on a spliced transcript structure. Skipping entry.\n";
            
            return undef;
        }
        
        if ($trans_orient eq '+') { 
            
            
            $cds_genome_lend = &from_cdna_lend($cds_span_rend, \@exon_coords);
            $cds_genome_rend = &from_cdna_lend($cds_span_lend, \@exon_coords);
            $transcribed_orient = '-';
            
        }
        
        elsif ($trans_orient eq '-') {
            
            $cds_genome_lend = &from_cdna_rend($cds_span_rend, \@exon_coords);
            $cds_genome_rend = &from_cdna_rend($cds_span_lend, \@exon_coords);
            $transcribed_orient = '+';
        }
        

        
    }
    
    my $new_gene_obj = new Gene_obj();
    $new_gene_obj->build_gene_obj_exons_n_cds_range(\@exon_coords, $cds_genome_lend, $cds_genome_rend, $transcribed_orient);
    
    return ($new_gene_obj);
}


####
sub from_cdna_lend {
    my ($pt, $coords_aref) = @_;

    my $lend_accrue = 0;

    my @coords = sort {$a->[0]<=>$b->[0]} @$coords_aref;

    foreach my $coordset (@coords) {
        my ($lend, $rend) = @$coordset;

        my $seg_len = $rend - $lend + 1;

        my $rend_accrue = $lend_accrue + $seg_len;
        $lend_accrue++;

        
        if ($pt >= $lend_accrue && $pt <= $rend_accrue) {
            
            my $pos = $lend + ($pt - $lend_accrue);
            return($pos);
        }
        
        $lend_accrue = $rend_accrue;
    }
    

    die "Error, couldn't localize pt $pt within coordsets: " . Dumper($coords_aref);

    return;
}

####
sub from_cdna_rend {
    my ($pt, $coords_aref) = @_;
    
    my $lend_accrue = 0;
    
    my @coords = reverse sort {$a->[0]<=>$b->[0]} @$coords_aref;
    
    foreach my $coordset (@coords) {
        my ($lend, $rend) = @$coordset;
        
        my $seg_len = $rend - $lend + 1;
        
        my $rend_accrue = $lend_accrue + $seg_len;
        $lend_accrue++;
        
        
        if ($pt >= $lend_accrue && $pt <= $rend_accrue) {
            
            my $pos = $rend - ($pt - $lend_accrue);
            return($pos);
        }
                
        $lend_accrue = $rend_accrue;
    }
    
    
    die "Error, couldn't localize pt $pt within coordsets: " . Dumper($coords_aref);
    
    return;
}

####
sub parse_cdna_seq_lengths {
    my ($fasta_file) = @_;

    my %seq_lengths;

    my $fasta_reader = new Fasta_reader($fasta_file);
    while (my $seq_obj = $fasta_reader->next()) {

        my $acc = $seq_obj->get_accession();
        
        my $asmbl = $acc;
        
        if ($acc =~ /(asmbl_\d+)/) {
            # pasa stuff
            $asmbl = $1;
        }
        
        my $sequence = $seq_obj->get_sequence();

        $seq_lengths{$asmbl} = length($sequence);
    }
    
    return(%seq_lengths);
}

