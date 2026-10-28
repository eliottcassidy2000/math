#!/bin/sh
# Item 6: reproduce the three committed programs as documented and compare with the committed .out files (timing lines
# ignored).  The author's sources are compiled/run from the read-only worktree; outputs go to this audit folder.
A="$(cd "$(dirname "$0")" && pwd)"
WT=/c/Users/Eliott/AppData/Local/Temp/math-wt-s66-eisenstein
EXP=$WT/04-computation/experiments
cd "$A"
gcc -O2 "$EXP/eisenstein_four_fourth_powers_20261027.c" -o author_e4.exe || exit 1
: > repro_e4.out
for args in "all 20 0 50" "pat 30 22 1 9" "box 30 22" "box 30 38" "box 100 6" "box 150 6" "all 60 51 200" "all 120 51 200"; do
  echo "### ./a.out $args" >> repro_e4.out
  t0=$(date +%s)
  ./author_e4.exe $args >> repro_e4.out
  t1=$(date +%s)
  echo "[$((t1 - t0)) s]" >> repro_e4.out
done
python "$EXP/two_grids_class_group_20261027.py" > repro_two_grids.out 2>&1
python "$EXP/eisenstein_fourth_power_obstructions_20261027.py" > repro_obstructions.out 2>&1
strip() { grep -v -E '^\[[0-9]+ s\]$' "$1" | sed -E 's/  \[[0-9]+ s\]$//'; }
for pair in "repro_e4.out:$EXP/eisenstein_four_fourth_powers_20261027.out" \
            "repro_two_grids.out:$EXP/two_grids_class_group_20261027.out" \
            "repro_obstructions.out:$EXP/eisenstein_fourth_power_obstructions_20261027.out"; do
  mine=${pair%%:*}; theirs=${pair#*:}
  strip "$mine" > "$A/_a.tmp"; strip "$theirs" > "$A/_b.tmp"
  if diff -q "$A/_a.tmp" "$A/_b.tmp" > /dev/null; then echo "IDENTICAL (timing ignored): $mine vs $(basename "$theirs")";
  else echo "DIFFERENT: $mine vs $(basename "$theirs")"; diff "$A/_a.tmp" "$A/_b.tmp" | head -20; fi
done
rm -f "$A/_a.tmp" "$A/_b.tmp"
