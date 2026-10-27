#!/bin/bash
# sequential follow-on runs (one process at a time)
A=/private/tmp/claude-501/-Users-e-Documents-GitHub-math/e3293618-0c68-4b6f-9092-a857de278ee6/scratchpad/audit_thm4812
cd $A/vt
{
  echo "== order-9 regular census"; gentourng -q -d4 -D4 9 | ./vt_checks census 9
  echo "== order-7 regular census"; gentourng -q -d3 -D3 7 | ./vt_checks census 7
  echo "== order-5 regular census"; gentourng -q -d2 -D2 5 | ./vt_checks census 5
  echo "== Cayley"; for g in "3 1" "5 1" "7 1" "9 1" "3 3" "11 1" "13 1" "15 1" "5 3" "17 1"; do ./vt_checks cayley $g; done
  echo "== affine15"; ./vt_checks affine15 A; ./vt_checks affine15 B
} > vt_checks.out 2>&1
cd $A/window
for j in jobs_g_unique/*.job; do ./window_enum $j; done > window_enum_g.out 2>&1
cd $A/vt && cc -O3 -march=native -o tc_census tc_census.c && for n in 3 4 5 6 7 8 9 10; do gentourng -c -q $n | ./tc_census $n 61; done > tc_census.out 2>&1
cd $A/breakit
python3 random_groups.py 20261027 400 4 > rand_tours.txt 2> rand_groups.log
./tcbatch < rand_tours.txt > tcbatch.out 2>&1
for s in 3 4; do gentourng -c -q $s | ./subst_tc $s 13 30 2; done > subst.out 2>&1; gentourng -c -q 5 | ./subst_tc 5 13 30 1 >> subst.out 2>&1
gentourng -c -q 6 | ./subst_tc 6 13 30 1 >> subst.out 2>&1
echo PHASE2 DONE
