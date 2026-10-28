#!/bin/sh
# Compare reproduction outputs with the committed .out files: CRLF normalised, timing lines and [N s] suffixes removed.
A="$(cd "$(dirname "$0")" && pwd)"
EXP=/c/Users/Eliott/AppData/Local/Temp/math-wt-s66-eisenstein/04-computation/experiments
cd "$A"
norm() { tr -d '\r' < "$1" | grep -v -E '^\[[0-9]+ s\]$' | sed -E 's/ +\[[0-9]+ s\]$//'; }
for pair in "repro_e4.out:$EXP/eisenstein_four_fourth_powers_20261027.out" \
            "repro_two_grids.out:$EXP/two_grids_class_group_20261027.out" \
            "repro_obstructions.out:$EXP/eisenstein_fourth_power_obstructions_20261027.out"; do
  mine=${pair%%:*}; theirs=${pair#*:}
  norm "$mine" > _a.tmp; norm "$theirs" > _b.tmp
  la=$(wc -l < _a.tmp); lb=$(wc -l < _b.tmp)
  if cmp -s _a.tmp _b.tmp; then echo "IDENTICAL (CRLF and timing ignored; $la lines): $mine vs $(basename "$theirs")";
  else echo "DIFFERENT: $mine ($la lines) vs $(basename "$theirs") ($lb lines)"; diff _a.tmp _b.tmp | head -20; fi
done
rm -f _a.tmp _b.tmp
