#!/usr/bin/env bash
set -euo pipefail

# --- locate project root (…/h9_annotation) from this script path ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# --- allow overrides: bash fix_telogator2_ref.sh [in.fa] [out.fa] ---
REF_IN="${1:-${ROOT_DIR}/2_data/2.2_processed/H9_T2T_v0.1_dip.telogator2.ref.fa}"
REF_OUT="${2:-${ROOT_DIR}/2_data/2.2_processed/H9_T2T_v0.1_dip.telogator2.ref.fixed.fa}"

if [[ ! -s "$REF_IN" ]]; then
  echo "[fix-ref] ERROR: missing input fasta: $REF_IN" >&2
  exit 1
fi

python3 - "$REF_IN" "$REF_OUT" <<'PY'
import re, sys

inf, outf = sys.argv[1], sys.argv[2]

# Accept sample names with underscores; we'll sanitize them for output (no underscores before _chrNp).
pat = re.compile(r"^(.*)_chr([0-9]+|X|Y)_hap([12])([pq])$")

def die(msg: str) -> None:
    raise SystemExit(f"[fix-ref] {msg}")

with open(inf) as f, open(outf, "w") as g:
    for line in f:
        if line.startswith(">"):
            name = line[1:].strip().split()[0]
            m = pat.match(name)
            if not m:
                die(f"unexpected header: {name}")
            sample, chrom, hap, arm = m.groups()
            sample_safe = sample.replace("_", "-")
            new = f"{sample_safe}-hap{hap}_chr{chrom}{arm}"
            g.write(f">{new}\n")
        else:
            g.write(line)

# Try to index (optional)
try:
    import pysam
    pysam.faidx(outf)
    print("[fix-ref] indexed with pysam:", outf + ".fai")
except Exception as e:
    print("[fix-ref] note: no pysam faidx (ok). You can run: samtools faidx", outf)
print("[fix-ref] wrote:", outf)
PY

# --- sanity checks ---
echo "[fix-ref] head:"
grep -E '^>' "$REF_OUT" | head

echo -n "[fix-ref] count chr arms (expect 92 for hap1+hap2, chr1-22+X, p+q): "
grep -E '^>.*_chr([0-9]+|X|Y)[pq]$' "$REF_OUT" | wc -l
