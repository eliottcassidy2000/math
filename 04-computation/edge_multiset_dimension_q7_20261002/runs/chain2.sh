#!/bin/bash
# wait for d=7 k<=12 to finish; then 2-minute annealing calibration (2 processes); then the gated post-12
# checks (post12.sh launches k=13 only if all checks pass)
cd "$(dirname "$0")/.."
until grep -qE "ALLDONE 1..12|FAILED" logs/d7_e0.out; do sleep 20; done
if grep -q "ALLDONE 1..12" logs/d7_e0.out && ! grep -q FAILED logs/d7_e0.out; then
  nice -n 10 ./src/anneal 6 15 60 1 > logs/anneal_q6_k15_calib.out 2>&1 &
  nice -n 10 ./src/anneal 7 19 120 1 > logs/anneal_q7_k19_calib.out 2>&1 &
  wait
  runs/post12.sh
else echo "k<=12 run failed; nothing launched" > logs/post12.out; fi
