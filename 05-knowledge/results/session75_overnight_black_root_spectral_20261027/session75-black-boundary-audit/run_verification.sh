#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
sha256sum -c SHA256SUMS >/dev/null
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
for mode in normal optimized;do
 if [[ $mode == normal ]];then py=();opt=(-O2);else py=(-O);opt=(-O3 -DNDEBUG);fi
 g++ -std=c++17 "${opt[@]}" audit_universe.cpp -o "$T/universe-$mode"
 python3 -B "${py[@]}" universe_input.py | "$T/universe-$mode" > "$T/universe-$mode.json"
 cmp "$T/universe-$mode.json" "universe-$mode.json"
 for topic in families survivors finite controls;do
  output=$topic; [[ $topic == families ]] && output=family
  python3 -B "${py[@]}" "audit_$topic.py" > "$T/$output-$mode.json"
  cmp "$T/$output-$mode.json" "$output-$mode.json"
 done
done
printf 'PASS: independent family, survivor, ordinary universe, every finite witness, and rejection controls.\n'
