#!/bin/bash
# Simulated annealing for Q7, sizes 18 and 17 (2 processes at a time, nice), 2 runs x SEC seconds per size:
# one with the default schedule (2M moves per restart), one with 4x longer restarts (8M moves per restart).
# (Calibration runs on sizes where resolving sets are known to exist, Q6 k=15 and Q7 k=19, were made by
#  runs/chain2.sh: logs/anneal_q6_k15_calib.out, logs/anneal_q7_k19_calib.out.)
# usage: runs/anneal.sh [SEC]
cd "$(dirname "$0")/.."
SEC=${1:-900}; N="nice -n 10"; A=./src/anneal
$N $A 7 18 $SEC 11 > logs/anneal_q7_k18_s11.out 2>&1 &
$N $A 7 18 $SEC 12 8000000 > logs/anneal_q7_k18_s12_long.out 2>&1 &
wait
$N $A 7 17 $SEC 21 > logs/anneal_q7_k17_s21.out 2>&1 &
$N $A 7 17 $SEC 22 8000000 > logs/anneal_q7_k17_s22_long.out 2>&1 &
wait
echo "ANNEAL-DONE $(date)" >> logs/anneal_done.out
