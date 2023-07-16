#!/bin/bash

if [ ! -e Trinity.fasta ]; then
    gunzip -c Trinity.fasta.gz > Trinity.fasta
    gunzip -c genome_alignments.gmap.gff3.gz > genome_alignments.gmap.gff3
fi

docker run --rm  -v `pwd`:/data trinityrnaseq/transdecoder:latest TransDecoder.LongOrfs -t /data/Trinity.fasta -O /data

docker run --rm  -v `pwd`:/data trinityrnaseq/transdecoder:latest TransDecoder.Predict -t /data/Trinity.fasta -O /data 

# gmap was used to align the Trinity.fasta transcripts to the genome,
# using the gmap '-f 3' output formatting parameter, generating file 'genome_alignments.gmap.gff3'

docker run --rm  -v `pwd`:/data trinityrnaseq/transdecoder:latest bash -c "/usr/local/bin/util/cdna_alignment_orf_to_genome_orf.pl  /data/Trinity.fasta.transdecoder.gff3 /data/genome_alignments.gmap.gff3 /data/Trinity.fasta  > /data/Trinity.fasta.transdecoder.genome.gff3"


docker run --rm  -v `pwd`:/data trinityrnaseq/transdecoder:latest /usr/local/bin/util/fasta_prot_checker.pl /data/Trinity.fasta.transdecoder.pep


exit 0
