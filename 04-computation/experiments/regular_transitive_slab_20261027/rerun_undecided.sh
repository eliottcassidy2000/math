#!/bin/sh
# rerun every UNDECIDED regular instance of a regular_geodesic_n<n>_*.out file with a larger budget and several seeds
# usage: ./rerun_undecided.sh n infile budget cands
n=$1; f=$2; budget=$3; cands=$4
grep '^UNDECIDED' "$f" | sed 's/.*result -1)://' | while read rows; do
  ok=0
  for seed in 1 2 3 4 5 6; do
    out=$(ROWS="$rows" ./paley_geodesic_search $n 3 1 $seed $budget $cands | grep 'mode=3')
    case "$out" in *"replayed 1"*) echo "RESOLVED seed $seed: $rows :: $out"; ok=1; break;; esac
  done
  [ $ok = 0 ] && echo "STILL UNDECIDED: $rows"
done
