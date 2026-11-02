# Bundled presets

Coding model, start (Kozak) model and q-value thresholds learned from a reference transcriptome, for
inputs too small (fewer than ~1,000 transcripts) or too partial to learn from: `--preset <name>`.

| name | built from | genetic code | GC % | long-ORF auto-accept (nt) |
|---|---|---|---|---|
| mammal | *Homo sapiens* (GENCODE coding + lncRNA) | Universal | 49.8 | 726 |
| insect | *Drosophila melanogaster* | Universal | 50.8 | 804 |
| plant | *Arabidopsis thaliana* | Universal | 42.1 | 663 |
| fungus | *Schizosaccharomyces pombe* | Universal | 37.8 | 615 |
| apicomplexan | *Plasmodium falciparum* | Universal | 22.0 | 519 |
| ciliate | *Tetrahymena thermophila* | Tetrahymena | 25.9 | none |
| green_alga | *Chlamydomonas reinhardtii* | Universal | 68.2 | 1347 |

Each reference holds all annotated complete coding transcripts plus coding fragments of them (and lncRNAs
for human), as used in the TransDecoder paper benchmark; presets built from complete transcripts and
fragments keep sensitivity on partial transcripts. Choose the closest clade: transfer between mammals is
nearly lossless, across distant clades sensitivity drops. For other species, build one from a Predict run
on a related transcriptome with `pytransdecoder build-preset` (TransDecoder Perl reads the same files).
The same presets ship with PyTransDecoder in `pytransdecoder/presets`.
