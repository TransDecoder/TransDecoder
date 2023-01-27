#!/bin/bash

if [ ! -e Trinity.fasta ]; then
    gunzip -c Trinity.fasta.gz > Trinity.fasta
    gunzip -c genome_alignments.gmap.gff3.gz > genome_alignments.gmap.gff3
fi

../../TransDecoder.LongOrfs -t Trinity.fasta $*

../../TransDecoder.Predict -t Trinity.fasta 

# gmap was used to align the Trinity.fasta transcripts to the genome,
# using the gmap '-f 3' output formatting parameter, generating file 'genome_alignments.gmap.gff3'

../../util/cdna_alignment_orf_to_genome_orf.pl  Trinity.fasta.transdecoder.gff3 genome_alignments.gmap.gff3 Trinity.fasta  > Trinity.fasta.transdecoder.genome.gff3


../../util/fasta_prot_checker.pl Trinity.fasta.transdecoder.pep


exit 0
