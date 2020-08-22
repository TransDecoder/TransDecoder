#!/usr/bin/env perl

use strict;
use warnings;
use Carp;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use GFF3_utils2;
use Data::Dumper;
use Fasta_reader;

my $usage = "\nusage: $0 cdna_orfs.genes.gff3 cdna_genome.alignments.gff3 cdna.fasta [DEBUG]\n\n";

my $cdna_orfs_gff3 = $ARGV[0] or die $usage;
my $cdna_genome_gff3 = $ARGV[1] or die $usage;
my $cdna_fasta = $ARGV[2] or die $usage;

# ensure can find the inputs
foreach my $file ($cdna_orfs_gff3, $cdna_genome_gff3, $cdna_fasta) {
    unless (-s $file) {
        die "Error, cannot locate file: $file";
    }
}


my $WARNING_COUNT = 0; # count those orfs that appear to be on strand opposite from the transcribed strand.

my $DEBUG = $ARGV[3];

main: {

    print STDERR "-parsing cdna lengths\n" if $DEBUG;
    my %cdna_seq_lengths = &parse_cdna_seq_lengths($cdna_fasta);
    
    my %orf_counter;

    print STDERR "-parse transcript alignment info\n" if $DEBUG;
    my %cdna_acc_to_transcript_structure = &parse_transcript_alignment_info($cdna_genome_gff3);

    ## parse ORFs on cDNAs

    print STDERR "-index $cdna_orfs_gff3\n" if $DEBUG;
    my $gene_obj_indexer_href = {};
    ## associate gene identifiers with contig id's.
    my $contig_to_gene_list_href = &GFF3_utils2::index_GFF3_gene_objs($cdna_orfs_gff3, $gene_obj_indexer_href);


    my %isolated_gene_id_to_new_genes;

    my $count_propagated = 0;
    my $total = 0;
    
    foreach my $asmbl_id (sort keys %$contig_to_gene_list_href) {
        ## $asmbl_id is the actual Transcript identifier from which ORFs were predicted.
        
        my @gene_ids = @{$contig_to_gene_list_href->{$asmbl_id}};
    
        foreach my $gene_id (@gene_ids) { # gene identifiers as given by transdecoder on the transcript sequences
            my $gene_obj_ref = $gene_obj_indexer_href->{$gene_id};

            $total++;
            
            my $transcript_struct = $cdna_acc_to_transcript_structure{$asmbl_id};
            unless ($transcript_struct) {
                print STDERR "-WARNING, $asmbl_id has no genome representation... skipping\n";
                next;
            }
            
            #print "$gene_obj_ref\n";

            
            eval {
                my $new_orf_gene = &place_orf_in_cdna_alignment_context($transcript_struct, $gene_obj_ref, \%cdna_seq_lengths);
                
                if ($new_orf_gene) {
                    
                    my $orf_count = $orf_counter{$asmbl_id}++;
                    $new_orf_gene->{asmbl_id} = $transcript_struct->{contig};
                    #$new_orf_gene->{TU_feat_name} = "t.$asmbl_id.$orf_count";
                    #$new_orf_gene->{Model_feat_name} = "m.$asmbl_id.$orf_count";
                    $new_orf_gene->{com_name} = $gene_obj_ref->{com_name}; #"ORF";
                    
                    my $use_gene_id = $transcript_struct->{gene_id};
                    unless ($use_gene_id) {
                        ## extract the orig gene id from the incoming gene obj
                        my ($isolated_gene_id, @rest) = split(/::/, $gene_id);
                        $use_gene_id = $isolated_gene_id;
                    }

                    ## some annoying gene annotations are showing up with same ids but different locations or strands.
                    ## going to address that by including the chr and orient info in the gene identifier.
                    my $chr = $new_orf_gene->{asmbl_id};
                    my $orient = $new_orf_gene->get_orientation(); # note, my have been reoriented
                    $use_gene_id = "${use_gene_id}^$chr^$orient";
                                        
                    $new_orf_gene->{TU_feat_name} = $use_gene_id;
                    $new_orf_gene->{Model_feat_name} = $gene_obj_ref->{Model_feat_name};

                    print STDERR "New orf: " . $new_orf_gene->toString() if $DEBUG;
                    push (@{$isolated_gene_id_to_new_genes{$use_gene_id}}, $new_orf_gene);

                    $count_propagated++;
                }
            };

            if ($@) {
                
                print STDERR "Error, $gene_id couldn't be fully propagated:\n";
                
                print STDERR Dumper($transcript_struct);
                print STDERR $gene_obj_ref->toString();
                print STDERR "$@";
            }
            
        }
    }
    
    ## output results:
    foreach my $gene_id (sort keys %isolated_gene_id_to_new_genes) {

        my @gene_objs = @{$isolated_gene_id_to_new_genes{$gene_id}};

        foreach my $gene_obj (@gene_objs) {
            $gene_obj->set_CDS_phases_from_init_phase(0);
        }

        my $parent_gene_obj = shift @gene_objs;
        foreach my $gene_obj (@gene_objs) {
            $parent_gene_obj->add_isoform($gene_obj);
        }

        $parent_gene_obj->refine_gene_object();
        
        print $parent_gene_obj->to_GFF3_format(source => "transdecoder") . "\n";
    }


    print STDERR "\n\n\tDone.  $count_propagated / $total transcript orfs could be propagated to the genome\n\n";
    
    exit(0);

}

####
sub parse_transcript_alignment_info {
    my ($cdna_align_gff3) = @_;

    my %cdna_alignments;

    open (my $fh, $cdna_align_gff3) or die "Error, cannot open file $cdna_align_gff3";
    while (<$fh>) {
        if (/^\#/) { next; }
        unless (/\w/) { next; }

        my $line = $_;
        
        my @x = split(/\t/);
        my $contig = $x[0];
        my $lend = $x[3];
        my $rend = $x[4];
        my $orient = $x[6];
        my $info = $x[8];

        $info =~ /Target=(\S+)/ or die "Error, cannot parse ID from info [$info] of line $line";
        my $asmbl = $1;

        my $trans_id = "";
        my $gene_id = "";
           
        ## trick for retaining gene and trans identifier information from cufflinks data 
        if ($asmbl =~ /^GENE\^(\S+),TRANS\^(\S+)/) {
            $gene_id = $1;
            $trans_id = $2;

            $asmbl = $2;
        }
        
        
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
                               trans_id => $trans_id,
                               gene_id => $gene_id,
                               
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

    #print $orf_gene_obj->toString();
    
    my $trans_seq_length = $cdna_seq_lengths_href->{ $transcript_struct->{asmbl} } or confess "Error, no length for " . Dumper($transcript_struct) . " Please be sure to use a cDNA fasta file and not a genome fasta file for your commandline parameter.";
    


    ## unwrap the gene
    my @cds_coords;
    my $orf_orient = $orf_gene_obj->get_orientation();
        
    foreach my $exon ($orf_gene_obj->get_exons()) {
        
        if (my $cds_exon = $exon->get_CDS_obj()) {

            my ($lend, $rend) = sort {$a<=>$b} $cds_exon->get_coords();
            push (@cds_coords, [$lend, $rend]);
        }
    }

    #print Dumper(\@cds_coords);
    
    @cds_coords = sort {$a->[0]<=>$b->[0]} @cds_coords;

    my $cds_span_lend = $cds_coords[0]->[0];
    my $cds_span_rend = $cds_coords[$#cds_coords]->[1];

    unless ($cds_span_lend && $cds_span_rend) {
        die "Error, missing cds span for gene: " . $orf_gene_obj->toString();
    }
    
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


            
            ##   ===========>......==========>..........===========>   exons
            ##          ---->......---------->..........---->          CDS

            $cds_genome_lend = &from_cdna_genome_lend($cds_span_lend, \@exon_coords);
            $cds_genome_rend = &from_cdna_genome_lend($cds_span_rend, \@exon_coords);
            $transcribed_orient = '+';

        }
        
        elsif ($trans_orient eq '-') {


            
            ##   <===========......<==========..........<===========   exons
            ##          <----......<----------..........<----          CDS

                        
            $cds_genome_lend = &from_cdna_genome_lend($trans_seq_length - $cds_span_rend + 1, \@exon_coords);
            $cds_genome_rend = &from_cdna_genome_lend($trans_seq_length - $cds_span_lend + 1, \@exon_coords);
            $transcribed_orient = '-';
            
        }
        
    }
    
    else {

        
        ####################
        ## orf orient is '-' w/ respect to transcript orientation
        ####################
        
        if (scalar(@exon_coords) > 1) {
            # any correct ORF should be in the '+' orientation here.... must be a false positive orf or transcript structure is wrong
            $WARNING_COUNT++;
            print STDERR "-Warning [$WARNING_COUNT], shouldn't have a minus-strand ORF on a plus-strand spliced transcript structure. Skipping entry $orf_gene_obj->{Model_feat_name}.\n";
            
            return undef;
        }
        
        
        if ($trans_orient eq '+') { 



            ## broken
            
            ##   ===========>......==========>..........===========>   exons
            ##          <----......<----------..........<----          CDS

            $WARNING_COUNT++;
            print STDERR "-Warning [$WARNING_COUNT], orf found on opposite strand of single-exon transcript. Reorienting transcribed direction accordingly to match orf orient on genome (-).\n";

            print STDERR "CDS coords: $cds_span_lend, $cds_span_rend\n";
            
            $cds_genome_lend = &from_cdna_genome_lend($cds_span_lend, \@exon_coords);
            $cds_genome_rend = &from_cdna_genome_lend($cds_span_rend, \@exon_coords);
            
            
            $transcribed_orient = '-';
            
        }
        
        elsif ($trans_orient eq '-') {

                        
            ##   <===========......<==========..........<===========   exons
            ##          ---->......---------->..........---->          CDS


            $WARNING_COUNT++;
            print STDERR "-Warning [$WARNING_COUNT], orf found on opposite strand of single-exon transcript. Reorienting transcribed direction accordingly to match orf orient on genome (+).\n";
                        
            $cds_genome_lend = &from_cdna_genome_lend($trans_seq_length - $cds_span_rend +1, \@exon_coords);
            $cds_genome_rend = &from_cdna_genome_lend($trans_seq_length - $cds_span_lend +1, \@exon_coords);
            $transcribed_orient = '+';
        }
        
    }
    
    my $new_gene_obj = new Gene_obj();
    $new_gene_obj->build_gene_obj_exons_n_cds_range(\@exon_coords, $cds_genome_lend, $cds_genome_rend, $transcribed_orient);
    
    return ($new_gene_obj);
}


####
sub from_cdna_genome_lend {
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
    

    die "Error, couldn't localize pt [$pt] within coordsets: " . Dumper($coords_aref);

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
        
        #if ($acc =~ /(asmbl_\d+)/) {
            # pasa stuff
        #    $asmbl = $1;
        #}
        
        my $sequence = $seq_obj->get_sequence();

        $seq_lengths{$asmbl} = length($sequence);
    }
    
    return(%seq_lengths);
}

