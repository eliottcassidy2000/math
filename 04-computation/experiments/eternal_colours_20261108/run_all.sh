#!/bin/sh
# Reproduce THM-4930's data (mac-mini-2026-11-08-firstrealizers). Large intermediate files go to $WORK (default ./work).
# Core steps take about 45 minutes on one machine; the optional cross-checks (OPTIONAL=1) take several hours.
set -e
cd "$(dirname "$0")"
WORK=${WORK:-./work}; mkdir -p "$WORK"
PKT=../../../05-knowledge/results/first_creators_genealogy_20261107
cc -O2 -o "$WORK/classgen" "$PKT/classgen.c"
( cd "$WORK" && ./classgen 9 cls )                                   # census, orders 1..9 (THM-4922's generator)
for n in 6 8; do cc -O3 -DN=$n -o "$WORK/sc_even$n" sc_even10.c; "$WORK/sc_even$n" > "$WORK/sc_even$n.out"; done
cc -O3 -DN=10 -o "$WORK/sc_even10" sc_even10.c; "$WORK/sc_even10" > sc_even10.out                   # order 10, ~18 min
cc -O3 -DCHECK -o "$WORK/sc_odd11_check" sc_odd11_fast.c; "$WORK/sc_odd11_check" 625 > /dev/null
cc -O3 -o "$WORK/sc_odd11_fast" sc_odd11_fast.c
Q=8388608
for k in 0 1 2 3; do "$WORK/sc_odd11_fast" 625 $((k*Q)) $(((k+1)*Q)) > "$WORK/sc11_part$k.out" & done; wait
cat "$WORK"/sc11_part*.out | awk '{s[$2]=(s[$2]||$3)} END {for (h in s) print "E", h, s[h]}' | sort -k2 -n > sc_odd11.out
# sc_odd13.c: validation builds (orders 6-9 against the census; orders 10, 11 complete), then orders 12, 13 below 1125
for m in 2 3; do cc -O3 -DCHECK -DM=$m -o "$WORK/sc13_m$m" sc_odd13.c; "$WORK/sc13_m$m" 4000 > "$WORK/v13_m$m.out"; done
cc -O3 -DCHECK -DM=4 -o "$WORK/sc13_m4" sc_odd13.c; "$WORK/sc13_m4" 95095 > "$WORK/v13_m4.out"        # ~4 min
awk '$2<=2001' "$WORK/v13_m4.out" > sc_full1011_le2001.out
cc -O3 -DCHECK -DM=5 -o "$WORK/sc13" sc_odd13.c
"$WORK/sc13" 1124 0 524288 > "$WORK/s13_0.out" & "$WORK/sc13" 1124 524288 1048576 > "$WORK/s13_1.out" & wait   # ~10 min
cat "$WORK"/s13_[01].out | awk '{k=$1" "$2; s[k]=(s[k]||$3)} END {for (k in s) print k, s[k]}' | sort -k1,1 -k2,2n > sc_odd13_1124.out
python3 validate_sc_odd13.py "$WORK" > sc_odd13_validation.out                                # full comparisons
cc -O2 -o "$WORK/sc_witness" sc_witness.c; "$WORK/sc_witness" 3 13 4 2001 > sc_witness.out
for f in sc_even10.out sc_odd11.out sc_full1011_le2001.out sc_odd13_1124.out sc_witness.out; do cp $f "$WORK/"; done
[ -f sc_walk.out ] && cp sc_walk.out "$WORK/"
[ -f sc_even12_377.out ] && cp sc_even12_377.out "$WORK/"
[ -d audit_K1 ] && ln -sfn "$(pwd)/audit_K1" "$WORK/audit_K1"
python3 check_order11.py "$WORK" > check_order11.out
python3 sc_holes.py "$WORK" > sc_holes.out
python3 refine_dichotomy.py "$WORK" > refine_dichotomy.out
if [ "${OPTIONAL:-0}" = 1 ]; then
  # order 12 up to 377 without the reduction (all 2^25 labelled bases; ~2-3 h on three workers)
  cc -O3 -DM=3 -o "$WORK/sc_pair8" sc_even_pair.c; "$WORK/sc_pair8" 100000 > "$WORK/sc_pair8.out"
  cc -O3 -DM=5 -o "$WORK/sc_pair12" sc_even_pair.c; T=33554432; P=$((T/3))
  "$WORK/sc_pair12" 377 0 $P > "$WORK/p12_0.out" & "$WORK/sc_pair12" 377 $P $((2*P)) > "$WORK/p12_1.out" & "$WORK/sc_pair12" 377 $((2*P)) $T > "$WORK/p12_2.out" & wait
  cat "$WORK"/p12_*.out | awk '{s[$2]=(s[$2]||$3)} END {for (h in s) print "P", h, s[h]}' | sort -k2 -n > sc_even12_377.out
  # random-walk witnesses at orders 9 (validation) and 12-15 (heuristic; only the witnesses count)
  cc -O3 -o "$WORK/sc_walk" sc_walk.c
  : > sc_walk_runs.out
  for a in "9 1 3357 200000 1" "13 625 720 600000 11" "13 625 900 800000 12" "13 625 1200 1000000 13" "13 625 2001 1500000 14" \
           "12 375 1300 1500000 21" "14 1125 1300 600000 22" "15 1875 2001 500000 23"; do
    echo "# run: $a" >> sc_walk_runs.out; "$WORK/sc_walk" $a >> sc_walk_runs.out; done
  grep '^W ' sc_walk_runs.out | awk '$3 >= 12' | sort -k2,2n -k3,3n | awk '!seen[$2]++' > sc_walk.out
  python3 check_witness_rev.py . sc_walk_runs.out > check_walk.out
  python3 walk_summary.py "$WORK" > walk_summary.out
  cp sc_even12_377.out sc_walk.out "$WORK/"; python3 sc_holes.py "$WORK" > sc_holes.out   # with the optional data
fi
echo "all done"
