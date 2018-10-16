


all: incl_blast_n_pfam simplest single_best best_with_homology

simplest:
	../util/select_best_ORFs_per_transcript.pl --gff3_file longest_orfs.gff3 --cds_scores longest_orfs.cds.scores > test.simplest
	diff test.simplest __test.simplest.expected


single_best:
	../util/select_best_ORFs_per_transcript.pl --gff3_file longest_orfs.gff3 --cds_scores longest_orfs.cds.scores --single_best_orf > test.single_best
	diff test.single_best __test.single_best.expected

incl_blast_n_pfam:
	../util/select_best_ORFs_per_transcript.pl --gff3_file longest_orfs.gff3 --cds_scores longest_orfs.cds.scores  --pfam_hits pfam.domtblout  --blast_hits blastp.outfmt6 > test.wBlastNPfam
	diff test.wBlastNPfam __test.wBlastNPfam.expected


best_with_homology:
	../util/select_best_ORFs_per_transcript.pl --gff3_file longest_orfs.gff3 --cds_scores longest_orfs.cds.scores  --pfam_hits pfam.domtblout  --blast_hits blastp.outfmt6 --single_best_orf > test.best_w_homology
	diff test.best_w_homology __test.best_w_homology.expected





clean:
	rm -f ./test.*



# for checking the trans count dist:
#
#  for file in test.*; do cat $file  | perl -lane 'if ($F[2] eq "CDS") { print;}' | print.pl 0 | token_counter.pl  > __$file.expected.hist; echo done $file; done

