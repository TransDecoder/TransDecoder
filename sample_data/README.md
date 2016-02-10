# TransDecoder examples

Here are three examples of how TransDecoder can be applied.

The simplest is in 'simple_transcriptome_target/', where TransDecoder is applied to a small Trinity.fasta file containing some assembled transcript sequences.  Here there is no reference genome to compare to.

Other examples where a reference genome is available include starting with Cufflinks output or PASA output.

The 'cufflinks_example/' contains a small genome sequence and cufflinks transcript annotations in GTF file format.  Here, the transcript sequences are extracted, TransDecoder is executed to find likely coding regions, and the coding region annotations are generated to reflect the coding regions on the genome sequence.

Similarly, in the 'pasa_example/', we start with a set of PASA-reconstructed transcript sequences, the transcript structure annotation on the genome, and a small genome sequence.  Likely coding regions are identified in the PASA assemblies, and a genome annotation file is generated to describe these coding regions in the genomic context.

