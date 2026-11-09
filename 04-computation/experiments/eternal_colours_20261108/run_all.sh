#!/bin/sh
# Reproduce THM-4930's data (mac-mini-2026-11-08-firstrealizers). Large intermediate files go to $WORK (default ./work).
# Total about 30 minutes on one machine (order 10: ~18 min; order 11: four parallel workers).
set -e
cd "$(dirname "$0")"
WORK=${WORK:-./work}; mkdir -p "$WORK"
PKT=../../../05-knowledge/results/first_creators_genealogy_20261107
cc -O2 -o "$WORK/classgen" "$PKT/classgen.c"
( cd "$WORK" && ./classgen 9 cls )                                   # census, orders 1..9 (THM-4922's generator)
for n in 6 8; do cc -O3 -DN=$n -o "$WORK/sc_even$n" sc_even10.c; "$WORK/sc_even$n" > "$WORK/sc_even$n.out"; done
cc -O3 -DN=10 -o "$WORK/sc_even10" sc_even10.c; "$WORK/sc_even10" > sc_even10.out
cp sc_even10.out "$WORK/"
cc -O3 -DCHECK -o "$WORK/sc_odd11_check" sc_odd11_fast.c; "$WORK/sc_odd11_check" 625 > /dev/null
cc -O3 -o "$WORK/sc_odd11_fast" sc_odd11_fast.c
Q=8388608
for k in 0 1 2 3; do "$WORK/sc_odd11_fast" 625 $((k*Q)) $(((k+1)*Q)) > "$WORK/sc11_part$k.out" & done; wait
cat "$WORK"/sc11_part*.out | awk '{s[$2]=(s[$2]||$3)} END {for (h in s) print "E", h, s[h]}' | sort -k2 -n > sc_odd11.out
cp sc_odd11.out "$WORK/"
cc -O2 -o "$WORK/sc_witness" sc_witness.c; "$WORK/sc_witness" 3 13 4 2001 > sc_witness.out
cp sc_witness.out "$WORK/"
python3 sc_holes.py "$WORK" > sc_holes.out
python3 refine_dichotomy.py "$WORK" > refine_dichotomy.out
echo "all done"
