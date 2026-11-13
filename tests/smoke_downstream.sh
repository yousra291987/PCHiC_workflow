#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Downstream smoke test — executes the parts of the pipeline a dry-run cannot.
#
# `snakemake -n` only validates the workflow graph; it never RUNS the R/Perl
# scripts, so runtime/format bugs inside them (wrong argument counts, wrong
# .mat columns feeding juicer, malformed PRE files) sail through CI. This test
# runs the real scripts on tiny fixtures — no HiCUP, bowtie2 or juicer_tools
# needed (juicer is replaced by tests/stub_java.sh) — and asserts the outputs
# have the exact shape downstream tools require.
#
# It guards the concrete regressions we have hit:
#   * BamToFragments.pl must emit the 5-column .mat (Fragment1 ProbeF1 Fragment2
#     ProbeF2 N) with ProbeF1/ProbeF2 as 0/1 bait flags — not empty — so the
#     capture-efficiency step can compute P0/P1/P2 instead of all zeros.
#   * 03_Generate_hic_and_TADs.r must select (Fragment1, Fragment2, N) so every
#     PRE row has 11 columns.
#   * Bam2Juicer.pl must never write a short/malformed PRE row.
#   * align2ibed.pl must accept its 8 args and produce an .ibed.
#
# Usage: bash tests/smoke_downstream.sh
# ---------------------------------------------------------------------------
set -euo pipefail

here="$(cd "$(dirname "$0")" && pwd)"
root="$(cd "$here/.." && pwd)"
fx="$here/fixtures"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

pass() { printf '  OK   %s\n' "$1"; }
fail() { printf '  FAIL %s\n' "$1"; [ "${2:-}" ] && { echo "----"; cat "$2"; echo "----"; }; exit 1; }

echo "== 1) BamToFragments.pl -> .mat =="
# BamToFragments.pl writes the 5-column table: Fragment1 ProbeF1 Fragment2 ProbeF2 N
mat="$work/toy.mat"
perl "$root/scripts/BamToFragments.pl" \
  "$fx/toy.bedpe" "$fx/toy.baitmap" "$fx/toy.rmap" "$mat" 1 2 4 5 2>/dev/null
cols=$(head -1 "$mat" | awk -F'\t' '{print NF}')
[ "$cols" -eq 5 ] || fail ".mat header has $cols columns, expected 5" "$mat"
[ "$(wc -l < "$mat")" -ge 2 ] || fail ".mat produced no interaction rows" "$mat"
head -1 "$mat" | grep -q 'Fragment1' || fail ".mat missing expected Fragment1 header" "$mat"
# ProbeF1/ProbeF2 (cols 2 & 4) must be 0/1 bait flags, not empty: the
# capture-efficiency step reads them and needs real flags (regression: they were
# taken from a non-existent 5th rmap column -> empty -> P0=P1=P2=0, empty plot).
badflag=$(tail -n +2 "$mat" | awk -F'\t' '($2!="0"&&$2!="1")||($4!="0"&&$4!="1"){c++} END{print c+0}')
[ "$badflag" -eq 0 ] || fail "ProbeF1/ProbeF2 are not all 0/1 bait flags" "$mat"
tail -n +2 "$mat" | awk -F'\t' '$2==1||$4==1{f++} END{exit !(f>0)}' \
  || fail ".mat has no fragment flagged as a bait (baitmap lookup failed)" "$mat"
pass ".mat is the 5-column table with 0/1 bait flags and interactions"

echo "== 2) 03_Generate_hic_and_TADs.r (real script, fake java) -> PRE + .hic =="
pre="$work/toy_pre.txt"
hic="$work/toy.hic"
Rscript "$root/scripts/03_Generate_hic_and_TADs.r" \
  --dir "$work" --sample toy \
  --mat "$mat" --rmap "$fx/toy.rmap" \
  --pre "$pre" --hic "$hic" \
  --java "$here/stub_java.sh" --juicer_tools "$work/dummy.jar" \
  --perl perl --bam2juicer "$root/scripts/Bam2Juicer.pl" \
  --genome_id mm10 --mem_gb 1 --q 1 >/dev/null
[ -s "$pre" ] || fail "PRE file not created or empty"
bad=$(awk 'NF!=11{c++} END{print c+0}' "$pre")
[ "$bad" -eq 0 ] || fail "$bad PRE row(s) do not have 11 columns" "$pre"
[ -s "$hic" ] || fail ".hic placeholder not created"
pass "every PRE row has 11 columns and .hic was produced"

echo "== 3) Bam2Juicer.pl with a CRLF (Windows) rmap -> still 11 columns =="
# Regression for the real-world failure: design files saved with CRLF line
# endings left a trailing \r on the rmap's id column, so Bam2Juicer.pl could
# not match fragment ids and emitted short, malformed PRE rows.
crlf_rmap="$work/toy_crlf.rmap"
sed 's/$/\r/' "$fx/toy.rmap" > "$crlf_rmap"
trip="$work/trip.tmp"
tail -n +2 "$mat" | awk -F'\t' 'BEGIN{OFS="\t"}{print $1,$3,$5}' > "$trip"
crlf_pre="$work/pre_crlf.txt"
perl "$root/scripts/Bam2Juicer.pl" "$trip" "$crlf_rmap" "$crlf_pre" 2>/dev/null
[ -s "$crlf_pre" ] || fail "Bam2Juicer produced no PRE rows from a CRLF rmap"
bad=$(awk -F'\t' 'NF!=11{c++} END{print c+0}' "$crlf_pre")
[ "$bad" -eq 0 ] || fail "$bad PRE row(s) from a CRLF rmap do not have 11 columns" "$crlf_pre"
pass "CRLF design files still yield well-formed 11-column PRE rows"

echo "== 4) align2ibed.pl -> .ibed =="
ibed="$work/toy.ibed"
perl "$root/scripts/align2ibed.pl" \
  "$fx/toy.bedpe" "$fx/toy.baitmap" "$fx/toy.rmap" "$ibed" 1 2 4 5 >/dev/null
[ -s "$ibed" ] || fail ".ibed not created or empty"
pass ".ibed was produced"

echo "== 5) BamToFragments.pl tolerates chr-naming differences (bed '1' vs rmap 'chr1') =="
# Regression: a real run on a non-'chr' assembly (e.g. chicken ggal7) had the
# rmap keyed as 'chr1' while the BEDPE used '1', so every coordtofrag() lookup
# missed -> an empty .mat and millions of "uninitialized value" warnings
# (a multi-GB log). Chromosome names must match after stripping a 'chr' prefix.
nochr_bedpe="$work/toy_nochr.bedpe"
sed -E 's/^chr//; s/\tchr/\t/g' "$fx/toy.bedpe" > "$nochr_bedpe"
mm_mat="$work/toy_mm.mat"
warn="$work/bamtofrag.warn"
perl "$root/scripts/BamToFragments.pl" \
  "$nochr_bedpe" "$fx/toy.baitmap" "$fx/toy.rmap" "$mm_mat" 1 2 4 5 2>"$warn" >/dev/null
[ "$(wc -l < "$mm_mat")" -ge 2 ] || fail "chr-mismatch: .mat has no interactions (lookup failed)" "$mm_mat"
if grep -q "uninitialized" "$warn"; then fail "chr-mismatch: coordtofrag still floods warnings" "$warn"; fi
pass "chr-naming difference handled, with no per-read warnings"

echo "== 6) 00_CaptureHiC_CaptureEfficiency_Control.r -> non-empty P0/P1/P2 =="
# Regression: with empty ProbeF1/ProbeF2 flags the stat was Total=N, P0=P1=P2=0
# and the CaptureEfficiency.pdf barplot was blank. P0+P1+P2 must equal Total and
# at least one captured (P1/P2) category must be > 0.
cedir="$work/ce"
mkdir -p "$cedir"
stat="$cedir/toy_CaptureEfficiency.stat"
Rscript "$root/scripts/00_CaptureHiC_CaptureEfficiency_Control.r" \
  --dir "$cedir" --sample toy --mat "$mat" --stats "$stat" >/dev/null 2>&1
[ -s "$stat" ] || fail "capture-efficiency stat not created"
tot=$(awk -F': ' '/^Total:/{print $2}' "$stat")
p0=$(awk -F': ' '/^P0:/{print $2}' "$stat")
p1=$(awk -F': ' '/^P1:/{print $2}' "$stat")
p2=$(awk -F': ' '/^P2:/{print $2}' "$stat")
[ "$((p0 + p1 + p2))" -eq "$tot" ] || fail "P0+P1+P2 ($((p0+p1+p2))) != Total ($tot)" "$stat"
[ "$((p1 + p2))" -gt 0 ] || fail "capture efficiency is 0 (P1=P2=0): bait flags not read" "$stat"
[ -s "$cedir/toy_CaptureEfficiency.pdf" ] || fail "CaptureEfficiency.pdf not created"
pass "capture efficiency computes real P0/P1/P2 and a non-empty plot"

echo
echo "ALL DOWNSTREAM SMOKE CHECKS PASSED"
