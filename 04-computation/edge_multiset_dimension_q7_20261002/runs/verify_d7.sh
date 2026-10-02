#!/bin/bash
# Verification-mode samples (-v) for d = 7 (engine 0, the production engine): every leaf is re-tested from
# scratch (keys recomputed from the definition), every hot-pair kill mask is compared with a direct
# recomputation, and every kill is checked against the true collision status.  Any discrepancy aborts (FATAL).
# usage: runs/verify_d7.sh big1|big2|small      (big1, big2 ~ 75M leaves each, ~9 min; small ~ a few minutes)
cd "$(dirname "$0")/.."
N="nice -n 10"; X=./src/edimsearch; R=data/q6
case "$1" in
  big1)  $N $X 7 12 6 $R/reps_n6_a6.bin 1777 1778 -e 0 -v ;;               # k=12, a=6, b=6, delta=0
  big2)  $N $X 7 13 7 $R/reps_n6_a7.bin 12345 12346 -e 0 -v ;;             # k=13, a=7, b=6, delta=1
  small) $N $X 7 12 7 $R/reps_n6_a7.bin 9000 9002 -e 0 -v                  # k=12, a=7, b=5, delta=2
         $N $X 7 12 8 $R/reps_n6_a8.bin 70000 70020 -e 0 -v                # k=12, a=8, b=4, delta=4
         $N $X 7 13 8 $R/reps_n6_a8.bin 60000 60010 -e 0 -v                # k=13, a=8, b=5, delta=3
         $N $X 7 13 9 $R/reps_n6_a9.bin 300000 300200 -e 0 -v              # k=13, a=9, b=4, delta=5
         $N $X 7 13 10 $R/reps_n6_a10.bin 0 3000 -e 0 -v ;;                # k=13, a=10, b=3, delta=7
esac
echo "VERIFY-DONE $1 rc=$?"
