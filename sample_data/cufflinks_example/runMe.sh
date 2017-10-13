#!/bin/bash

set -ev

export PERL_HASH_SEED=0

if [ ! -e test.genome.fasta ]; then
    gunzip -c test.genome.fasta.gz > test.genome.fasta
fi


if [ ! -e transcripts.gtf ]; then
    gunzip -c transcripts.gtf.gz > transcripts.gtf
fi

if [ ! -e blastp.outfmt6 ]; then
    gunzip -c blastp.outfmt6.gz > blastp.outfmt6
fi

if [ ! -e pfam.domtblout ]; then
    gunzip -c pfam.domtblout.gz > pfam.domtblout
fi


## generate alignment gff3 formatted output
../../util/gtf_to_alignment_gff3.pl transcripts.gtf > transcripts.gff3

## generate transcripts fasta file
../../util/gtf_genome_to_cdna_fasta.pl transcripts.gtf test.genome.fasta > transcripts.fasta 

## Extract the long ORFs
../../TransDecoder.LongOrfs -t transcripts.fasta

cmd=""
## Predict likely ORFs
if [ 1 ]; then   # always doing this now.

    # this is how I would have run blast and pfam but I'm using precomputed results for ease of demonstration.
    #BLASTDB=/seq/RNASEQ/DBs/TRINOTATE_RESOURCES/TRINOTATE_V3/uniprot_sprot.pep
    #PFAMDB=/seq/RNASEQ/DBs/TRINOTATE_RESOURCES/TRINOTATE_V3/Pfam-A.hmm
    #
    ## run blast
    #blastp -query transcripts.fasta.transdecoder_dir/longest_orfs.pep -db $BLASTDB -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > blastp.outfmt6
    #
    ## run pfam
    #hmmscan --domtblout pfam.domtblout $PFAMDB transcripts.fasta.transdecoder_dir/longest_orfs.pep > pfam.log
    
    
    ## use pfam and blast results:
    cmd="../../TransDecoder.Predict  -t transcripts.fasta --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6   -v"
else
    # just coding metrics
    cmd="../../TransDecoder.Predict -t transcripts.fasta"
fi

eval $cmd $ARGS


## convert to genome coordinates
../../util/cdna_alignment_orf_to_genome_orf.pl transcripts.fasta.transdecoder.gff3 transcripts.gff3 transcripts.fasta > transcripts.fasta.transdecoder.genome.gff3


## make bed files for viewing with GenomeView

# covert cufflinks gtf to bed
../../util/gtf_to_bed.pl transcripts.gtf > transcripts.bed

# convert the genome-based gene-gff3 file to bed
../../util/gff3_file_to_bed.pl transcripts.fasta.transdecoder.genome.gff3 > transcripts.fasta.transdecoder.genome.bed


# ensure no fatal problems w/ pep file
../../util/fasta_prot_checker.pl transcripts.fasta.transdecoder.pep

# Done!  Coding region genome annotations provided as: transcripts.fasta.transdecoder.genome.\*


exit 0
