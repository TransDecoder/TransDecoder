#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
workdir="$(mktemp -d "${TMPDIR:-/tmp}/gff3_to_proteins_streaming.XXXXXX")"
trap 'rm -rf "$workdir"' EXIT

cat > "$workdir/test.fa" <<'FASTA'
>tx1 first transcript
ATGGCCGCTTAA
>tx2 second transcript
ATGAAACCCGGGTAA
>tx3 third transcript
ATGTTTCCCTAG
FASTA

cat > "$workdir/test.gff3" <<'GFF3'
##gff-version 3
tx1	TD	gene	1	12	.	+	.	ID=gene1;Name=gene one
tx1	TD	mRNA	1	12	.	+	.	ID=mrna1;Parent=gene1
tx1	TD	exon	1	12	.	+	.	ID=exon1;Parent=mrna1
tx1	TD	CDS	1	12	.	+	0	ID=cds1;Parent=mrna1
tx2	TD	gene	1	15	.	+	.	ID=gene2;Name=gene two
tx2	TD	mRNA	1	15	.	+	.	ID=mrna2;Parent=gene2
tx2	TD	exon	1	15	.	+	.	ID=exon2;Parent=mrna2
tx2	TD	CDS	1	15	.	+	0	ID=cds2;Parent=mrna2
tx3	TD	gene	1	12	.	+	.	ID=gene3;Name=gene three
tx3	TD	mRNA	1	12	.	+	.	ID=mrna3;Parent=gene3
tx3	TD	exon	1	12	.	+	.	ID=exon3;Parent=mrna3
tx3	TD	CDS	1	12	.	+	0	ID=cds3;Parent=mrna3
GFF3

for seq_type in prot CDS cDNA gene; do
    "$repo_dir/util/gff3_file_to_proteins.pl" \
        --gff3 "$workdir/test.gff3" \
        --fasta "$workdir/test.fa" \
        --seqType "$seq_type" \
        > "$workdir/default.${seq_type}"

    for batch_size in 1 2 5000; do
        "$repo_dir/util/gff3_file_to_proteins.pl" \
            --batch_size "$batch_size" \
            --gff3 "$workdir/test.gff3" \
            --fasta "$workdir/test.fa" \
            --seqType "$seq_type" \
            > "$workdir/streaming.${seq_type}.${batch_size}"

        diff -u "$workdir/default.${seq_type}" "$workdir/streaming.${seq_type}.${batch_size}"
    done
done

if "$repo_dir/util/gff3_file_to_proteins.pl" \
    --legacy_memory_mode \
    --gff3 "$workdir/test.gff3" \
    --fasta "$workdir/test.fa" \
    > "$workdir/legacy_mode.out" 2> "$workdir/legacy_mode.err"; then
    echo "Expected removed --legacy_memory_mode option to fail" >&2
    exit 1
fi

if ! grep -q "Unknown option: legacy_memory_mode" "$workdir/legacy_mode.err"; then
    echo "Removed --legacy_memory_mode option was not rejected as expected" >&2
    cat "$workdir/legacy_mode.err" >&2
    exit 1
fi

cat > "$workdir/missing.fa" <<'FASTA'
>tx1 first transcript
ATGGCCGCTTAA
FASTA

if "$repo_dir/util/gff3_file_to_proteins.pl" \
    --gff3 "$workdir/test.gff3" \
    --fasta "$workdir/missing.fa" \
    > "$workdir/missing.out" 2> "$workdir/missing.err"; then
    echo "Expected missing FASTA entry check to fail" >&2
    exit 1
fi

if ! grep -q "no FASTA entry found for GFF3 sequence id: tx2" "$workdir/missing.err"; then
    echo "Missing FASTA entry error was not reported as expected" >&2
    cat "$workdir/missing.err" >&2
    exit 1
fi

echo "gff3_file_to_proteins streaming tests passed"
