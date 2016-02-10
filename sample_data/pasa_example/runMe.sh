#!/bin/bash

if [ ! -e genome.fasta ]; then
    gunzip -c genome.fasta.gz > genome.fasta
fi

if [ ! -e pasa_assemblies.fasta ]; then
    gunzip -c pasa_assemblies.fasta.gz > pasa_assemblies.fasta
fi

if [ ! -e pasa_assemblies.gff3 ]; then
    gunzip -c pasa_assemblies.gff3.gz > pasa_assemblies.gff3
fi

../../TransDecoder.LongOrfs -t pasa_assemblies.fasta

../../TransDecoder.Predict -t pasa_assemblies.fasta

../../util/cdna_alignment_orf_to_genome_orf.pl  pasa_assemblies.fasta.transdecoder.gff3 pasa_assemblies.gff3 pasa_assemblies.fasta | tee pasa_assemblies.fasta.transdecoder.genome.gff3

exit 0
