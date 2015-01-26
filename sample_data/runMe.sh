#!/bin/bash -e

if [ -e test.genome.fasta.gz ] && [ ! -e test.genome.fasta ]; then
    gunzip -c test.genome.fasta.gz > test.genome.fasta
fi

if [ -e test.tophat.sam.gz ] && [ ! -e test.tophat.sam ]; then
    gunzip -c test.tophat.sam.gz > test.tophat.sam
fi

if [ -e transcripts.gtf.gz ] && [ ! -e transcripts.gtf ]; then
    gunzip -c transcripts.gtf.gz > transcripts.gtf
fi

if [ -e blastp.results.outfmt6.gz ] && [ ! -e blastp.results.outfmt6 ]; then
    gunzip -c blastp.results.outfmt6.gz > blastp.results.outfmt6
fi

if [ -e pfam.domtblout.gz ] && [ ! -e pfam.domtblout ]; then
    gunzip -c pfam.domtblout.gz > pfam.domtblout
fi



## generate alignment gff3 formatted output
../util/cufflinks_gtf_to_alignment_gff3.pl transcripts.gtf > transcripts.gff3

## generate transcripts fasta file
../util/cufflinks_gtf_genome_to_cdna_fasta.pl transcripts.gtf test.genome.fasta > transcripts.fasta 

## Extract the long ORFs
../TransDecoder.LongOrfs -t transcripts.fasta


## Predict likely ORFs
if [ $1 ]; then
    ## use pfam and blast results:
    ../TransDecoder.Predict  -t transcripts.fasta --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.results.outfmt6 -v
else
    # just coding metrics
    ../TransDecoder.Predict -t transcripts.fasta 
fi

## convert to genome coordinates
../util/cdna_alignment_orf_to_genome_orf.pl transcripts.fasta.transdecoder.gff3 transcripts.gff3 transcripts.fasta > transcripts.fasta.transdecoder.genome.gff3


## make bed files for viewing with GenomeView

# covert cufflinks gtf to bed
../util/cufflinks_gtf_to_bed.pl transcripts.gtf > transcripts.bed

# convert the genome-based gene-gff3 file to bed
../util/gff3_file_to_bed.pl transcripts.fasta.transdecoder.genome.gff3 > transcripts.fasta.transdecoder.genome.bed

echo
echo
echo Done!  Coding region genome annotations provided as: best_candidates.eclipsed_orfs_removed.genome.gff3
echo
echo 

exit 0
