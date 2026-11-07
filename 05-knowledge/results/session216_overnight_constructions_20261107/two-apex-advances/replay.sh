#!/bin/sh
set -eu
cd "$(dirname "$0")"
python3 -c 'import sys; assert sys.byteorder == "little", "Native C++ binary replay requires a little-endian host"'
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT HUP INT TERM
g++ -O3 -std=c++17 census_depths.cpp -o "$tmp/census"
"$tmp/census" 14 4000000 "$tmp/census.bin" > "$tmp/census.out"
cmp census14full.out "$tmp/census.out"
python3 reconstruct_certificate.py "$tmp/coverage2327673.bin.gz"
gzip -cd "$tmp/coverage2327673.bin.gz" > "$tmp/coverage.bin"
python3 - "$tmp" <<'PY'
from pathlib import Path
import sys
p=Path(sys.argv[1]); small=(p/'coverage.bin').read_bytes()
assert len(small)==9310696
assert small==(p/'census.bin').read_bytes()[:len(small)]
print('Full census receipt and complete packed certificate match')
PY
g++ -O3 -std=c++17 audit/check.cpp -o "$tmp/check"
"$tmp/check" "$tmp/coverage.bin" > "$tmp/check.out"
cmp audit/check.out "$tmp/check.out"
cat "$tmp/check.out"
g++ -O3 -std=c++17 audit/modes.cpp -o "$tmp/modes"
"$tmp/modes" > "$tmp/modes.out"
cmp audit/modes.out "$tmp/modes.out"
cat "$tmp/modes.out"
g++ -O3 -std=c++17 audit.cpp -o "$tmp/subset"
printf '2327675 120021002033010\n' | "$tmp/subset"
for name in global_map_wrappers two_letter_wrappers ternary_global_wrappers ternary_modes verify_modes; do
  PYTHONDONTWRITEBYTECODE=1 python3 "$name.py" > "$tmp/$name.out"
  python3 - "$name.out" "$tmp/$name.out" <<'PY'
import json,sys
old,new=[json.load(open(p)) for p in sys.argv[1:]]
old.pop('seconds',None);new.pop('seconds',None)
assert old==new,sys.argv[1]
print(sys.argv[1], 'semantic receipt matches')
PY
done
