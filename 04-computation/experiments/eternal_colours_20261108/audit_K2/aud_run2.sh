#!/bin/bash
# audit K2 of THM-4930 statement 5 / HYP-9520 (2026-11-08): reproduce the auditor's computations.
# Needs audit K1's work products (run ../audit_K1/aud_run.sh first; WORK1 = its work directory with cls_n.txt,
# census_vals_n.txt, anti_n_red.txt, anti_10_red.txt and o11_full.txt).
# Wall time on the loaded mac-mini: validation ~4 min; production (orders 12-13, H <= 1124) ~11 min on 2 niced workers.
set -e
cd "$(dirname "$0")"
WORK1=${WORK1:-../audit_K1/work}; WORK=${WORK:-./work}; mkdir -p "$WORK"
PKT=..
cc -O3 -march=native -DCHECK -o "$WORK/aud_sc13_chk" aud_sc13.c
cc -O3 -march=native -o "$WORK/aud_witness16" aud_witness16.c
# 1. validation with no H bound: must equal the auditor's complete reduced censuses (count for count)
for NO in 5 7 9 11; do
  NE=$((NO-1)); nice -n 10 "$WORK/aud_sc13_chk" $NO 200000 > "$WORK/v13_$NO.txt"
  grep '^E ' "$WORK/v13_$NO.txt" | cut -d' ' -f2- | cmp - "$WORK1/anti_${NE}_red.txt"
  if [ $NO -lt 11 ]; then grep '^O ' "$WORK/v13_$NO.txt" | cut -d' ' -f2- | cmp - "$WORK1/anti_${NO}_red.txt";
  else grep '^O ' "$WORK/v13_$NO.txt" | cut -d' ' -f2- | cmp - "$WORK1/o11_full.txt"; fi
done
# pruned validation at order 11 (LIM 1125)
nice -n 10 "$WORK/aud_sc13_chk" 11 1125 > "$WORK/v13_11_lim.txt"
grep '^O ' "$WORK/v13_11_lim.txt" | cut -d' ' -f2- | cmp - <(awk '$1<1125' "$WORK1/o11_full.txt")
# 2. production: orders 12 and 13 with H <= 1124, two niced workers over interleaved quarters of the base codes
Q=262144
( nice -n 10 "$WORK/aud_sc13_chk" 13 1125 0 $Q > "$WORK/p13_0.txt"; nice -n 10 "$WORK/aud_sc13_chk" 13 1125 $((2*Q)) $((3*Q)) > "$WORK/p13_2.txt" ) &
( nice -n 10 "$WORK/aud_sc13_chk" 13 1125 $Q $((2*Q)) > "$WORK/p13_1.txt"; nice -n 10 "$WORK/aud_sc13_chk" 13 1125 $((3*Q)) $((4*Q)) > "$WORK/p13_3.txt" ) &
wait
python3 - "$WORK" << 'EOF'
import sys, collections
W = sys.argv[1]; cnt = {"E": collections.Counter(), "O": collections.Counter()}; sc = {"E": collections.Counter(), "O": collections.Counter()}
for k in range(4):
    for line in open(f"{W}/p13_{k}.txt"):
        t, h, c, s = line.split(); cnt[t][int(h)] += int(c); sc[t][int(h)] += int(s)
with open(f"{W}/o1213_lt1125.txt", "w") as fo:
    for t in ("E", "O"):
        for h in sorted(cnt[t]): fo.write(f"{t} {h} {cnt[t][h]} {sc[t][h]}\n")
EOF
# 3. lists, witnesses, sandwiches
python3 aud_lists2.py "$WORK1" "$WORK" "$PKT"
"$WORK/aud_witness16" < "$PKT/sc_walk.out"
"$WORK/aud_witness16" < "$PKT/sc_witness.out"
python3 aud_sandwich.py
