#!/bin/bash

set -ex

# convert to gff3 alignment format
../../util/gtf_to_alignment_gff3.pl supertranscripts.gtf > supertranscripts.gff3

# extract the transcript sequences
../../util/gtf_genome_to_cdna_fasta.pl supertranscripts.gtf supertranscripts.fasta > transcripts.fasta

# run TransDecoder
../../TransDecoder.LongOrfs -t transcripts.fasta

cmd="../../TransDecoder.Predict -t transcripts.fasta"
if [ $1 ]; then
    cmd="$cmd --no_refine_starts"
fi

eval $cmd

# propagate predicted ORFs to the supertranscript gff3 file:
../../util/cdna_alignment_orf_to_genome_orf.pl transcripts.fasta.transdecoder.gff3 supertranscripts.gff3 transcripts.fasta > supertranscripts.wOrfs.gff3

# make bed file for visualization
../../util/gff3_file_to_bed.pl supertranscripts.wOrfs.gff3 > supertranscripts.wOrfs.bed

# convert to gtf as companion to the gff3 file (some other tools prefer gtf format)
../../util/gff3_gene_to_gtf_format.pl supertranscripts.wOrfs.gff3 supertranscripts.fasta  > supertranscripts.wOrfs.gtf


echo Done
