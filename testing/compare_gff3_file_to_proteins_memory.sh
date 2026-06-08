#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    cat >&2 <<'USAGE'
Usage: testing/compare_gff3_file_to_proteins_memory.sh <num_records> <seq_len>

Generates synthetic transcript FASTA/GFF3 data and compares maximum resident
set size for legacy all-in-memory mode and default streaming mode. Requires
/usr/bin/time with the -v option.
USAGE
    exit 2
fi

num_records="$1"
seq_len="$2"
repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
workdir="$(mktemp -d "${TMPDIR:-/tmp}/gff3_to_proteins_memory.XXXXXX")"
trap 'rm -rf "$workdir"' EXIT

python3 - "$num_records" "$seq_len" "$workdir/test.fa" "$workdir/test.gff3" <<'PY'
import sys
num_records = int(sys.argv[1])
seq_len = int(sys.argv[2])
fasta = sys.argv[3]
gff3 = sys.argv[4]
if seq_len < 6:
    raise SystemExit("seq_len must be at least 6")
seq = "ATG" + ("A" * (seq_len - 6)) + "TAA"
with open(fasta, "w") as fa, open(gff3, "w") as gf:
    gf.write("##gff-version 3\n")
    for i in range(1, num_records + 1):
        tid = f"tx{i}"
        gid = f"gene{i}"
        mid = f"mrna{i}"
        fa.write(f">{tid}\n{seq}\n")
        gf.write(f"{tid}\tTD\tgene\t1\t{seq_len}\t.\t+\t.\tID={gid}\n")
        gf.write(f"{tid}\tTD\tmRNA\t1\t{seq_len}\t.\t+\t.\tID={mid};Parent={gid}\n")
        gf.write(f"{tid}\tTD\texon\t1\t{seq_len}\t.\t+\t.\tID=exon{i};Parent={mid}\n")
        gf.write(f"{tid}\tTD\tCDS\t1\t{seq_len}\t.\t+\t0\tID=cds{i};Parent={mid}\n")
PY

/usr/bin/time -v "$repo_dir/util/gff3_file_to_proteins.pl" \
    --legacy_memory_mode \
    --gff3 "$workdir/test.gff3" \
    --fasta "$workdir/test.fa" \
    > "$workdir/legacy.pep" 2> "$workdir/legacy.time"

/usr/bin/time -v "$repo_dir/util/gff3_file_to_proteins.pl" \
    --batch_size 1000 \
    --gff3 "$workdir/test.gff3" \
    --fasta "$workdir/test.fa" \
    > "$workdir/streaming.pep" 2> "$workdir/streaming.time"

diff -u "$workdir/legacy.pep" "$workdir/streaming.pep"

echo "Legacy mode memory:"
grep "Maximum resident set size" "$workdir/legacy.time"
echo "Streaming mode memory:"
grep "Maximum resident set size" "$workdir/streaming.time"
