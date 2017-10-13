#!/bin/bash -ve

export PERL_HASH_SEED=0



## generate alignment gff3 formatted output
../../util/gtf_to_alignment_gff3.pl stringtie_merged.gtf > stringtie_merged.gff3


## generate transcripts fasta file
# not including the genome here... too big, but here's how you'd do it.
#../../util/gtf_genome_to_cdna_fasta.pl stringtie_merged.gtf  genome.fasta > stringtie_merged.transcripts.fasta

## Extract the long ORFs
../../TransDecoder.LongOrfs -t stringtie_merged.transcripts.fasta -S


## Predict likely ORFs

../../TransDecoder.Predict -t stringtie_merged.transcripts.fasta $ARGS


## convert to genome coordinates
../../util/cdna_alignment_orf_to_genome_orf.pl stringtie_merged.transcripts.fasta.transdecoder.gff3 stringtie_merged.gff3 stringtie_merged.transcripts.fasta > stringtie_merged.transcripts.fasta.transdecoder.genome.gff3   


## make bed files for viewing with GenomeView

# covert cufflinks gtf to bed
../../util/gtf_to_bed.pl stringtie_merged.gtf > stringtie_merged.bed 

# convert the genome-based gene-gff3 file to bed
../../util/gff3_file_to_bed.pl stringtie_merged.transcripts.fasta.transdecoder.genome.gff3  > stringtie_merged.transcripts.fasta.transdecoder.genome.bed  


../../util/fasta_prot_checker.pl stringtie_merged.transcripts.fasta.transdecoder.pep

# Done!  Coding region genome annotations provided as: transcripts.fasta.transdecoder.genome.\*


exit 0
