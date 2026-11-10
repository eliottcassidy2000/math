#!/bin/sh
# AUDIT 2: rerun the independent checks of audit2/AUDIT2.md.  Binaries and scratch parts go to $W (default
# /tmp/claude-audit2-ensemble); text outputs are rewritten next to this script.  Needs cc, python3 and nauty's gentourng.
# Usage: sh run_all.sh [stage ...]     stages: a b c d e f   (default: all, in this order)
#   a  small orders, about 2 minutes single core: DP-free brute force; vertex cubes of all labelled tournaments m <= 6
#      and of all classes m <= 8 (direct counts); box norms N <= 8; alternating matrices m <= 7; affine spans; L4; null model
#   b  vertex cubes of all 2^21 labelled tournaments of order 7 (direct); alternating matrices m = 8      (P processes)
#   c  vertex cubes of all 191,536 classes of order 9, every out-set counted directly                      (P processes)
#   d  vertex cubes of all 9,733,056 classes of order 10 (insertion formula + one direct control per class) (P processes)
#   e  box norms / deletion pairs from the classes of orders 9 and 10 (route 2 on all placements; route 1 on a
#      random 1/64 of the placements at order 9 and on a slice of order 10)                                (P processes)
#   f  summaries
# With P = 6 on a shared, heavily loaded 8-core machine the stages took about 1, 5, 10, 9, 12 and 0 minutes
# (stage e: about 41 CPU-minutes, 7 minutes on 6 free cores).
set -e
A=$(cd "$(dirname "$0")" && pwd); W=${W:-/tmp/claude-audit2-ensemble}; G=${GENTOURNG:-/opt/homebrew/bin/gentourng}; P=${P:-6}
mkdir -p "$W"; cd "$A"
[ $# -eq 0 ] && set -- a b c d e f
for p in x1_cube x2_box x5_altrank x6_affine; do [ -x "$W/$p" ] && [ "$W/$p" -nt "$p.c" ] || cc -O2 -o "$W/$p" "$p.c"; done
R=$(seq 0 $((P - 1)))
for stage in "$@"; do case $stage in
a)
  python3 -B x3_bruteforce_small.py > x3_bruteforce_small.out
  (for m in 2 3 4 5 6; do "$W/x1_cube" $m labelled 0 1 pgd; done) > x1_raw_labelled_2_6.out
  (for m in 2 3 4 5 6 7 8; do $G -q $m | "$W/x1_cube" $m direct pgd; done) > x1_raw_direct_2_8.out
  (for n in 2 3 4 5 6 7; do "$W/x2_box" $n labelled; done) > x2_raw_labelled_2_7.out
  (for n in 2 3 4 5 6 7 8; do $G -q $n | "$W/x2_box" $n classes 1; done) > x2_raw_classes_2_8.out
  (for m in 2 3 4 5 6 7; do "$W/x5_altrank" $m; done) > "$W/x5_2_7.raw"
  (for m in 2 3 4 5 6 7; do "$W/x6_affine" $m all; done; for m in 8 9 10 11; do "$W/x6_affine" $m random 400 20261110; done) > x6_affine.out
  python3 -B x7_deletion_rank.py 4 5 6 7 > x7_deletion_rank.out
  python3 -B x8_nullmodel.py > x8_nullmodel.out
  python3 -B x10_image.py > x10_image.out ;;
b)
  for r in $R; do ( "$W/x1_cube" 7 labelled $r $P pgd > "$W/x1_lab7_$r.raw" ) & done; wait
  cat $(for r in $R; do echo "$W/x1_lab7_$r.raw"; done) > x1_raw_labelled_7.out
  for r in $R; do ( "$W/x5_altrank" 8 $r $P > "$W/x5_8_$r.raw" ) & done; wait
  cat "$W/x5_2_7.raw" $(for r in $R; do echo "$W/x5_8_$r.raw"; done) > x5_raw_altrank_2_8.out ;;
c)
  for r in $R; do ( $G -q 9 $r/$P | "$W/x1_cube" 9 direct pgd > "$W/x1_d9_$r.raw" ) & done; wait
  cat $(for r in $R; do echo "$W/x1_d9_$r.raw"; done) > x1_raw_direct_9.out ;;
d)
  for r in $R; do ( $G -q 10 $r/$P | "$W/x1_cube" 10 cut 1 d > "$W/x1_c10_$r.raw" ) & done; wait
  cat $(for r in $R; do echo "$W/x1_c10_$r.raw"; done) > x1_raw_cut_10.out ;;
e)
  for r in $R; do ( $G -q 9 $r/$P | "$W/x2_box" 9 classes 2 > "$W/x2_9_$r.raw" ) & done; wait
  cat $(for r in $R; do echo "$W/x2_9_$r.raw"; done) > x2_raw_classes_9.out
  for r in $R; do ( $G -q 10 $r/$P | "$W/x2_box" 10 classes 0 > "$W/x2_10_$r.raw" ) & done; wait
  cat $(for r in $R; do echo "$W/x2_10_$r.raw"; done) > x2_raw_classes_10.out
  # route 1 (explicit cycle counts) against route 2 on a random 1/64 of the placements, for a slice of the order-10 classes
  ( $G -q 10 7/3000 | "$W/x2_box" 10 classes 2 ) > x2_raw_classes_10_slice_check.out ;;
f)
  python3 -B x1_merge.py x1_raw_labelled_2_6.out x1_raw_labelled_7.out x1_raw_direct_2_8.out x1_raw_direct_9.out x1_raw_cut_10.out > x1_summary.out
  python3 -B x2_merge.py x2_raw_labelled_2_7.out x2_raw_classes_2_8.out x2_raw_classes_9.out x2_raw_classes_10.out > x2_summary.out
  python3 -B x5_oddedge_check.py x5_raw_altrank_2_8.out --packet=../out/oddedge_model.out x1_raw_labelled_2_6.out x1_raw_labelled_7.out x1_raw_direct_2_8.out > x5_oddedge_check.out
  python3 -B x9_compare_author.py ../out x1_raw_direct_2_8.out x1_raw_direct_9.out x1_raw_cut_10.out > x9_compare_author.out ;;
esac; done
echo "done: $*"
