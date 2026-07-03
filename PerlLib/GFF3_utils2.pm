#!/usr/local/bin/perl

package main;
our $SEE;


package GFF3_utils2; # removes DB_File requirement, no use of Gene_obj_indexer module.

use strict;
use warnings;
use Gene_obj;
use Carp;
use URI::Escape;
use Data::Dumper;
use Scalar::Util qw(refaddr);

my %CONTIG_ITERATOR_STATE;


####
sub index_GFF3_gene_objs {
    
    my ($gff_filename, $gene_obj_indexer, $contig_id) = @_;
    # contig_id is optional.
    
    unless (ref $gene_obj_indexer eq 'HASH') {
        confess "Error, \$gene_obj_indexer must be a hashref";
    }
    
    ## note can use either a gene_obj_indexer or a hash reference.
    
    my %asmbl_id_to_gene_id_list;

    open (my $fh, $gff_filename) or die $!;

    my $parse_context = _new_gff3_parse_context();
    
    my $counter = 0;
    # print STDERR "\n-parsing file $gff_filename\n";
    while (my $line = <$fh>) {
        my ($line_text, $fields_aref) = _parse_gff3_line_fields($line);
        next unless $fields_aref;
        _parse_gff3_feature($line_text, $fields_aref, $parse_context, $contig_id);
    }
    close $fh;
    
    ## 
    # print STDERR "\n-caching genes.\n";
    foreach my $asmbl_id (sort keys %{$parse_context->{gene_coords}}) {
        my $gene_ids_aref = _store_gene_objs_for_asmbl($asmbl_id, $parse_context, $gene_obj_indexer);
        $asmbl_id_to_gene_id_list{$asmbl_id} = $gene_ids_aref if @$gene_ids_aref;
    }
    #print STDERR "\n";
    return (\%asmbl_id_to_gene_id_list);
}

####
# Scalar context returns [$asmbl_id, \@gene_ids] or undef at EOF.
# List context returns ($asmbl_id, \@gene_ids) or an empty list at EOF.
sub next_contig_gene_objs {

    my ($fh, $gene_obj_indexer) = @_;

    unless (ref $gene_obj_indexer eq 'HASH') {
        confess "Error, \$gene_obj_indexer must be a hashref";
    }

    my $fh_key = ref($fh) ? refaddr($fh) : "$fh";
    my $iterator_state = $CONTIG_ITERATOR_STATE{$fh_key} ||= {
        pending_record => undef,
        seen_asmbl_ids => {},
    };

    my $parse_context = _new_gff3_parse_context();
    my $current_asmbl_id;

    while (1) {
        my $record = delete $iterator_state->{pending_record};

        unless ($record) {
            my $line = <$fh>;

            unless (defined $line) {
                delete $CONTIG_ITERATOR_STATE{$fh_key};

                if (defined $current_asmbl_id) {
                    my $gene_ids_aref = _store_gene_objs_for_asmbl($current_asmbl_id, $parse_context, $gene_obj_indexer);
                    return wantarray ? ($current_asmbl_id, $gene_ids_aref) : [$current_asmbl_id, $gene_ids_aref];
                }

                return wantarray ? () : undef;
            }

            my ($line_text, $fields_aref) = _parse_gff3_line_fields($line);
            next unless $fields_aref;

            $record = {
                line_text => $line_text,
                fields_aref => $fields_aref,
            };
        }

        my $asmbl_id = $record->{fields_aref}->[0];

        unless (defined $current_asmbl_id) {
            if ($iterator_state->{seen_asmbl_ids}->{$asmbl_id}) {
                confess "Error, GFF3 asmbl_id $asmbl_id appears in non-contiguous blocks";
            }
            $current_asmbl_id = $asmbl_id;
        }

        if ($asmbl_id ne $current_asmbl_id) {
            $iterator_state->{pending_record} = $record;
            $iterator_state->{seen_asmbl_ids}->{$current_asmbl_id} = 1;

            my $gene_ids_aref = _store_gene_objs_for_asmbl($current_asmbl_id, $parse_context, $gene_obj_indexer);
            return wantarray ? ($current_asmbl_id, $gene_ids_aref) : [$current_asmbl_id, $gene_ids_aref];
        }

        _parse_gff3_feature($record->{line_text}, $record->{fields_aref}, $parse_context);
    }
}


sub _new_gff3_parse_context {

    return {
        gene_coords => {},
        transcript_to_gene => {},
        cds_phases => {},
        gene_names => {},
        loci => {},
        gene_id_to_source_type => {},
        source_tracker => {},
    };
}


sub _parse_gff3_line_fields {

    my ($line) = @_;

    chomp $line;

    unless ($line =~ /\w/) { return; } # empty line

    if ($line =~ /^\#/) { return; } # comment entry in gff3

    my @x = split (/\t/, $line);

    unless (scalar @x >= 9) {
        print STDERR "-ignoring line $line\n";
        return;
    }

    return ($line, \@x);
}


sub _parse_gff3_feature {

    my ($line, $x_aref, $parse_context, $contig_id) = @_;

    my ($asmbl_id, $source, $feat_type, $lend, $rend, $orient, $cds_phase, $gene_info) = ($x_aref->[0], $x_aref->[1], $x_aref->[2], $x_aref->[3], $x_aref->[4], $x_aref->[6], $x_aref->[7], $x_aref->[8]);

    if ($contig_id && $asmbl_id ne $contig_id) { return; }

    unless ($feat_type) { die "Error, $line, no feat_type: line\[$line\]"; }

    unless ($feat_type =~ /^(gene|mRNA|CDS|exon)$/) { return;} ## these are the only fields I care about right now.

    $gene_info = uri_unescape($gene_info);

    $gene_info =~ /ID=([^;\s]+);?/;
    my $id = $1 or die "Error, couldn't get the id field $line";

    my $source_tracker_href = $parse_context->{source_tracker};
    if (exists $source_tracker_href->{$id} && $source_tracker_href->{$id} ne $source) {
        confess "Error, gene ID $id is given source $source when previously encountered with source $source_tracker_href->{$id} ";
    }

    if ($feat_type eq 'gene') {
        my $gene_name = "";
        if ($gene_info =~ /Name=\"?([^\;\"]+)\"?/) {
            $gene_name = $1;
        }
        else {
            $gene_name = "";
        }

        if ($gene_info =~ /Note=\"?([^\;\"]+)\"?/) {
            $gene_name .= " $1";
        }

        $parse_context->{gene_names}->{$id} = $gene_name;

    }

    if ($gene_info =~ /Alias=([^;]+)/) {
        my $locus = $1;
        $parse_context->{loci}->{$id} = $locus;
    }


    if ($feat_type eq 'gene') { return;} ## beyond this pt, gene is not needed.

    $gene_info =~ /Parent=([^;\s]+);?/;
    my $parent = $1 or die "Error, couldn't get the parent info $line";

    # print "id: $id, parent: $parent\n";

    if ($feat_type eq 'mRNA') {
        ## just get the identifier info
        $parse_context->{transcript_to_gene}->{$id} = $parent;
        return;
    }

    my $transcript_id = $parent;
    my $gene_id = $parse_context->{transcript_to_gene}->{$transcript_id};
    unless (defined $gene_id) {
        print STDERR "Error, no gene feature found for $transcript_id.... ignoring feature.\n";
        return;
    }


    $parse_context->{gene_id_to_source_type}->{$gene_id} = $source;

    my ($end5, $end3) = ($orient eq '+') ? ($lend, $rend) : ($rend, $lend);

    $parse_context->{gene_coords}->{$asmbl_id}->{$gene_id}->{$transcript_id}->{$feat_type}->{$end5} = $end3;
    # print "$asmbl_id, $gene_id, $transcript_id, $feat_type, $end5, $end3\n";

    if ($cds_phase =~ /^\d+$/) {
        $parse_context->{cds_phases}->{$gene_id}->{$transcript_id}->{$end5} = $cds_phase;
    }
}


sub _store_gene_objs_for_asmbl {

    my ($asmbl_id, $parse_context, $gene_obj_indexer) = @_;

    my $genes_href = $parse_context->{gene_coords}->{$asmbl_id} || {};
    my @gene_ids;

    foreach my $gene_id (keys %$genes_href) {
        #print STDERR "\r-indexing [$gene_id]  ";
        my $transcripts_href = $genes_href->{$gene_id};

        my @gene_objs;

        foreach my $transcript_id (keys %$transcripts_href) {

            my $cds_coords_href = $transcripts_href->{$transcript_id}->{CDS} || {}; # could be a noncoding transcript w/ no CDS
            my $exon_coords_href = $transcripts_href->{$transcript_id}->{exon};

            unless (ref $exon_coords_href) {
                print STDERR Dumper ($transcripts_href);
                die "Error, missing exon coords for $transcript_id, $gene_id\n";
            }

            my $gene_obj = new Gene_obj();


            if (scalar (keys %$cds_coords_href) == 1) {

                ## could be that only the cds span was provided.
                ## break it up across the exon segments

                my ($cds_lend, $cds_rend) = sort {$a<=>$b} %$cds_coords_href;
                my @exon_coords;
                my $orient;
                foreach my $exon_end5 (keys %$exon_coords_href) {
                    my $exon_end3 = $exon_coords_href->{$exon_end5};
                    push (@exon_coords, [$exon_end5, $exon_end3]);
                    if ($exon_end5 < $exon_end3) {
                        $orient = '+';
                    }
                    elsif ($exon_end5 > $exon_end3) {
                        $orient = '-';
                    }
                }

                $gene_obj->build_gene_obj_exons_n_cds_range(\@exon_coords, $cds_lend, $cds_rend, $orient);
            }
            else {

                ## cds and exons specified separately

                $gene_obj->populate_gene_obj($cds_coords_href, $exon_coords_href);
            }

            $gene_obj->{Model_feat_name} = $transcript_id;
            $gene_obj->{TU_feat_name} = $gene_id;
            $gene_obj->{asmbl_id} = $asmbl_id;

            if (my $gene_locus = $parse_context->{loci}->{$gene_id}) {
                $gene_obj->{pub_locus} = $gene_locus;
            }
            if (my $transcript_locus = $parse_context->{loci}->{$transcript_id}) {
                $gene_obj->{model_pub_locus} = $transcript_locus;
            }


            $gene_obj->{com_name} = $parse_context->{gene_names}->{$gene_id} || $transcript_id;

            $gene_obj->{source} = $parse_context->{gene_id_to_source_type}->{$gene_id};

            ## set CDS phase info if available from the gff
            my $cds_phases_href = $parse_context->{cds_phases}->{$gene_id}->{$transcript_id};
            if (ref $cds_phases_href) {
                ## set the cds phases
                my @exons = $gene_obj->get_exons();
                foreach my $exon (@exons) {
                    if (my $cds = $exon->get_CDS_obj()) {
                        my ($end5, $end3) = $cds->get_coords();
                        my $phase = $cds_phases_href->{$end5};
                        unless ($phase =~ /\d+/) {
                            confess "Error, should have phase set for cds $gene_id $transcript_id $end5, but I do not. ";
                        }
                        $cds->set_phase($phase);
                    }
                }
            }

            push (@gene_objs, $gene_obj);
        }

        ## want single gene that includes all alt splice variants here
        my $template_gene_obj = shift @gene_objs;
        foreach my $other_gene_obj (@gene_objs) {
            $template_gene_obj->add_isoform($other_gene_obj);
        }

        $template_gene_obj->refine_gene_object();


        $gene_obj_indexer->{$gene_id} = $template_gene_obj;

        print "GFF3_utils: stored $gene_id\n" if $SEE;

        # add to gene list for asmbl_id
        push (@gene_ids, $gene_id);
    }

    return \@gene_ids;
}


1; #EOM
