#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
sha256sum -c SHA256SUMS >/dev/null
T=$(mktemp -d);trap 'rm -rf "$T"' EXIT
for mode in normal optimized;do
 if [[ $mode == normal ]];then py=(); opt=(-O2);else py=(-O);opt=(-O3 -DNDEBUG);fi
 g++ "${opt[@]}" -std=c++17 finite/authenticate_universe.cpp -o "$T/authenticate-$mode"
 python3 -B "${py[@]}" finite/universe_input.py | "$T/authenticate-$mode" > "$T/universe-$mode.json"
 cmp "$T/universe-$mode.json" "finite/universe-$mode.json"
 python3 -B "${py[@]}" finite/verify_witnesses.py > "$T/witness-$mode.json"
 cmp "$T/witness-$mode.json" "finite/witness-$mode.json"
 python3 -B "${py[@]}" finite/verify_counterexample.py > "$T/counterexample-$mode.json"
 cmp "$T/counterexample-$mode.json" "finite/counterexample-$mode.json"
 python3 -B "${py[@]}" family/check_family.py > "$T/family-$mode.json"
 cmp "$T/family-$mode.json" "family/$mode.json"
 python3 -B "${py[@]}" family/check_odd.py > "$T/odd-$mode.json"
 cmp "$T/odd-$mode.json" "family/odd-$mode.json"
 python3 -B "${py[@]}" survivors/check_survivors.py > "$T/survivor-$mode.json"
 cmp "$T/survivor-$mode.json" "survivors/$mode.json"
done
printf 'PASS: authenticated all ordinary classes n4..9; replayed every witness; all family and survivor controls agree.\n'
