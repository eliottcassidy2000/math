#!/bin/sh
# audit K1 of THM-4930 / HYP-9520 (2026-11-08): reproduce the auditor's computations.
# Large intermediate files go to $WORK (default ./work; the census files cls_n.txt are ~30 MB and are NOT committed).
# Wall time on the loaded 8-core mac-mini: about 40 min (order 11, three workers) + 10 min (order 12, one worker).
set -e
cd "$(dirname "$0")"
WORK=${WORK:-./work}; mkdir -p "$WORK"
PKT=..                                                                    # the packet directory (read only)
GEN=../../../../05-knowledge/results/first_creators_genealogy_20261107   # THM-4922's census generator
cc -O2 -o "$WORK/classgen" "$GEN/classgen.c"
( cd "$WORK" && ./classgen 9 cls )
cc -O3 -march=native -o "$WORK/aud_census" aud_census.c
cc -O3 -march=native -DCHECKANTI -o "$WORK/aud_sc_chk" aud_sc.c
cc -O3 -march=native -o "$WORK/aud_sc" aud_sc.c
cc -O3 -march=native -DCHECK -o "$WORK/aud_sc12_chk" aud_sc12.c
cc -O3 -march=native -o "$WORK/aud_witness" aud_witness.c
cc -O3 -march=native -DCHECK -o "$WORK/aud_scpair_chk" aud_scpair.c
cc -O3 -march=native -o "$WORK/aud_find11" aud_find11.c
# 1. census re-verification, orders 1..9 (H, c3, strong, SC, |Aut| per class; orbit-stabilizer)
for n in 1 2 3 4 5 6 7 8 9; do "$WORK/aud_census" $n "$WORK/cls_$n.txt" vals > "$WORK/census_vals_$n.txt"; done
# 2. anti-invariant model, orders 1..9, full and reduced, with the anti-invariance assertion on every tournament
for n in 1 2 3 4 5 6 7 8 9; do for m in full red; do "$WORK/aud_sc_chk" $n $m > "$WORK/anti_${n}_$m.txt"; done; done
python3 aud_compare_le9.py "$WORK"
# 3. order 10 (reduced model, 2^20 codes): compare with the packet and with THM-4922's audit class census
"$WORK/aud_sc" 10 red > "$WORK/anti_10_red.txt"
python3 aud_compare_o10.py "$WORK" "$PKT" "$GEN/audit/aud_o10_byH.txt"
# 4. order 11: ALL self-converse 11-tournaments (reduced model, 2^25 codes, no H bound, no pruning), three workers
A=11184811; B=22369622; T=33554432
"$WORK/aud_sc" 11 red 0 $A > "$WORK/o11_p0.txt" & "$WORK/aud_sc" 11 red $A $B > "$WORK/o11_p1.txt" &
"$WORK/aud_sc" 11 red $B $T > "$WORK/o11_p2.txt" & wait
python3 aud_merge.py "$WORK"/o11_p0.txt "$WORK"/o11_p1.txt "$WORK"/o11_p2.txt > "$WORK/o11_full.txt"
# 5. order 12, H < 625 (supplementary; decides 377), and the validation of the pair-extension method:
#    with no H bound it must reproduce the reduced enumeration exactly at n = 4, 6, 8, 10, and with LIM 625 at n = 10
"$WORK/aud_sc12_chk" 625 0 1048576 > "$WORK/o12_lt625.txt"
for n in 4 6 8 10; do "$WORK/aud_scpair_chk" $n 200000 > "$WORK/pair_$n.txt"; cmp "$WORK/pair_$n.txt" "$WORK/anti_${n}_red.txt"; done
"$WORK/aud_scpair_chk" 10 625 > "$WORK/pair_10_lt625.txt"; awk '$1<625' "$WORK/anti_10_red.txt" | cmp - "$WORK/pair_10_lt625.txt"
"$WORK/aud_scpair_chk" 12 625 > "$WORK/pair12_lt625.txt"; cmp "$WORK/pair12_lt625.txt" "$WORK/o12_lt625.txt"
# explicit order-11 witnesses for 1551 and 1627 (HYP-9520's list)
"$WORK/aud_find11" 1551 1627
# 6. witnesses, lists, refinement table
"$WORK/aud_witness" < "$PKT/sc_witness.out"
python3 aud_lists.py "$WORK" "$PKT" --no12     # the packet's setting
python3 aud_lists.py "$WORK" "$PKT"            # with order 12 below 625
python3 aud_275.py 34 209
python3 aud_refine.py "$WORK"
python3 "$PKT/refine_dichotomy.py" "$WORK"
