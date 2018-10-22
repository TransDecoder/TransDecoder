#!/bin/bash -ve

if [ ! -e genome.fasta ]; then
    gunzip -c genome.fasta.gz > genome.fasta
fi

if [ ! -e pasa_assemblies.fasta ]; then
    gunzip -c pasa_assemblies.fasta.gz > pasa_assemblies.fasta
fi

if [ ! -e pasa_assemblies.gff3 ]; then
    gunzip -c pasa_assemblies.gff3.gz > pasa_assemblies.gff3
fi

if [ ! -e pasa_assemblies_described.txt ]; then
    gunzip -c pasa_assemblies_described.txt.gz > pasa_assemblies_described.txt
fi


# get the gene-to-transcript relationships
cut -f2,3 pasa_assemblies_described.txt > pasa.gene_trans_map.txt

../../TransDecoder.LongOrfs -t pasa_assemblies.fasta --gene_trans_map pasa.gene_trans_map.txt -O pasa.transdecoder_workdir



../../TransDecoder.Predict -t pasa_assemblies.fasta $ARGS -O pasa.transdecoder_workdir

../../util/cdna_alignment_orf_to_genome_orf.pl  pasa_assemblies.fasta.transdecoder.gff3 pasa_assemblies.gff3 pasa_assemblies.fasta  >  pasa_assemblies.fasta.transdecoder.genome.gff3


../../util/fasta_prot_checker.pl pasa_assemblies.fasta.transdecoder.pep


echo "Done.  See pasa_assemblies.fasta.transdecoder.\*"


exit 0
