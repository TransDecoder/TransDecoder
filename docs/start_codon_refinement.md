# Start-codon refinement: complete or 5' partial?

This page explains, step by step, how TransDecoder decides where an ORF starts when the ORF runs
off the 5' end of its transcript. The rule is the same in Perl TransDecoder (v7.0.0 and later)
and PyTransDecoder (v0.2.0 and later).

| | Perl TransDecoder.Predict | PyTransDecoder predict |
|---|---|---|
| per-base UTR allowance (alpha) | `--start_contrast_alpha` (0.4) | `--start-contrast-alpha` (0.4) |
| bases per bit of Kozak score | `--start_kozak_weight` (17) | `--start-kozak-weight` (17) |
| turn refinement off | `--no_refine_starts` | `--no-refine-starts` |

## The problem

LongOrfs reports an ORF as **5' partial** when there is no in-frame stop codon upstream of it. The
ORF then begins at the very first codon of the transcript, not at an ATG. Two different
situations look like this:

1. **The transcript is truncated.** The real start codon lies further upstream, beyond the end of
   the assembled sequence. The ORF should stay 5' partial.
2. **The transcript is complete, but its 5' UTR happens to contain no in-frame stop codon.** The
   ORF runs back through the UTR to the transcript's end. The real start is an ATG further
   downstream, and the ORF should be trimmed to it and called complete.

TransDecoder looks at the ORF's **first in-frame ATG** and at the **upstream segment**: the stretch
between the transcript's 5' end and that ATG. It then asks: does the upstream segment look like
part of the protein (case 1), or like UTR in front of a real start codon (case 2)?

## The unit: "bases' worth of typical coding signal"

Both sides of the comparison are measured in one unit, so they can be compared directly.

TransDecoder trains a coding model on the transcriptome's own longest ORFs. For any stretch of
sequence it computes a **frame contrast**: how much better the stretch scores as coding in the
ORF's reading frame than in the better of the two shifted frames. Real coding sequence scores
clearly higher in its own frame. UTR does not.

The typical (median) frame contrast per base of the training ORFs is one unit: **one base's worth
of typical coding signal**. Measured this way:

- a truly coding segment of 90 bases scores about 90;
- a 90-base UTR scores about 0, or below.

Why divide by the typical value? The raw contrast per base differs between organisms, by about
two-fold across the eight organisms in the TransDecoder benchmark, because codon bias is stronger
in some genomes than others. Dividing by each transcriptome's own typical value puts every
organism on the same scale, so one set of default parameters works everywhere. The typical value
is estimated automatically from each run's training ORFs.

## The Kozak score, in bits

The sequence around a real start codon (the Kozak context) differs from the sequence around ATGs
that sit inside coding regions. TransDecoder learns this from the training ORFs: it compares the
bases at positions -6 to -1 and +4, +5 (the A of the ATG is +1) at their start codons with the
same positions around the ATGs inside them.

The Kozak score of an ATG is in **bits** (log2-odds):

- **0 bits**: the context looks like an average ATG inside a coding region;
- **+1 bit**: twice as start-like; **+2 bits**: four times;
- **-1 bit**: half as start-like.

On the TransDecoder benchmark, real start codons score a median of about +0.5 to +1.3 bits, and
ATGs inside coding regions about -0.5 to -1.6 bits. The signal is real but weak, and the two
groups overlap a good deal.

The trained model is written to the work directory as `start_refinement.kozak.tsv` (bits per base
and position), `start_refinement.kozak.counts.tsv`, and a sequence logo,
`start_refinement.kozak.logo.svg`.

## The rule: two competing scores

**Coding score**: how coding the upstream segment looks.

    coding score = frame contrast of the upstream segment / typical coding signal per base

**UTR score**: the case for "this is UTR, and the ATG is the real start".

    UTR score = 0.4 x upstream length  +  17 x Kozak bits

The UTR score has two parts:

- **Baseline, 0.4 x upstream length** (alpha = 0.4). A truly coding segment would score about
  1.0 x its length, a UTR about 0. The UTR side gets 0.4 per base by default: a segment must look
  more than 40% as coding as typical protein-coding sequence before TransDecoder treats it as
  part of the protein.
- **Start bonus, 17 x Kozak bits.** A favourable start context adds to the UTR side and a poor one
  subtracts from it. In other words, each bit of Kozak score lets the ATG overrule 17 bases' worth
  of coding signal.

**Whichever score is higher wins:**

- **UTR score higher**: the ORF is trimmed to start at the ATG and is relabelled complete (or 3'
  partial if it has no stop codon).
- **Coding score higher, or equal**: the ORF stays 5' partial.

## Worked examples

A 90-base upstream segment, so the baseline is 0.4 x 90 = 36:

| coding score | Kozak score | UTR score | winner |
|---|---|---|---|
| 45 | +1 bit | 36 + 17 = 53 | UTR: **trim to the ATG** |
| 45 | 0 bits | 36 + 0 = 36 | coding: **stay 5' partial** |
| 45 | -1 bit | 36 - 17 = 19 | coding: **stay 5' partial** |
| 10 | -1 bit | 36 - 17 = 19 | UTR: **trim to the ATG** |

The last row shows that a weak start context cannot save a segment that clearly does not look
coding.

How much each part matters depends on the length of the upstream segment:

- **Short segments (under about 60-100 bases)** carry little coding signal, and the baseline is
  small, so the Kozak bonus (typically -35 to +35) usually decides.
- **Long segments (over about 300 bases)** have a baseline of 120 or more, and the coding signal of
  the segment usually decides.

In PyTransDecoder's GFF3 output, every trimmed ORF records the two scores:
`start_revised=true;upstream_coding_score=...;upstream_utr_score=...;kozak_score=...` (Kozak in
bits).

## Details

- Only the first in-frame ATG is considered. Choosing among later ATGs by their start context did
  not improve accuracy on the benchmark.
- If fewer than 50 training ORFs begin with ATG, no Kozak model is trained and the Kozak score is
  0: the decision uses the coding score alone.
- Refinement is skipped when there are no training ORFs to estimate the typical coding signal, or
  when it is turned off.
- The procedure is deterministic: the same input always gives the same result.

## Choosing the parameters

The defaults (0.4 and 17) were evaluated on a benchmark of eight organisms and on human GM12878
assemblies. With complete transcripts and transcript fragments weighted equally, they were the
best of 78 combinations, and a leave-one-organism-out test chose them in 7 of 8 cases. Accuracy
changes little for a Kozak weight between about 7 and 35 bases per bit.

`--start_contrast_alpha` (Perl) / `--start-contrast-alpha` (Python) trades two kinds of error:

- **Higher alpha (0.6-1.0)** trims more ORFs. Suits **full-length transcripts**, such as long-read
  assemblies, where nearly every open 5' end is really UTR.
- **Lower alpha (0.1-0.2)** keeps more ORFs 5' partial. Suits **fragmentary assemblies**, such as
  de novo assemblies of short reads, where many open 5' ends are real truncations.

Releases up to v6.0.0 instead chose among in-frame ATGs with position weight matrices trained on
randomly sampled start sites; their results varied with the random seed.
