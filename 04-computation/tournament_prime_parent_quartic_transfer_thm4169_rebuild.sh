#!/usr/bin/env bash
set -euo pipefail

THM4169_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
THM4169_SOURCE="$THM4169_ROOT/04-computation/tournament_prime_parent_quartic_transfer_thm4169_audit.cpp"
THM4169_PYTHON="$THM4169_ROOT/04-computation/tournament_prime_parent_quartic_transfer_thm4169_independent_audit.py"
THM4169_RESULTS="$THM4169_ROOT/05-knowledge/results"
THM4169_PARENT=111111101111111111111101111110110111110111111
THM4169_TMP="$(mktemp -d "${TMPDIR:-/tmp}/thm4169.XXXXXX")"
trap 'rm -rf -- "$THM4169_TMP"' EXIT

for THM4169_TOOL in clang++ gentourng python3 shortg sha256sum; do
  command -v "$THM4169_TOOL" >/dev/null || {
    echo "missing required tool: $THM4169_TOOL" >&2
    exit 2
  }
done

clang++ -std=c++20 -Wall -Wextra -Werror -O0 \
  "$THM4169_SOURCE" -o "$THM4169_TMP/audit_O0"
clang++ -std=c++20 -Wall -Wextra -Werror -O3 \
  "$THM4169_SOURCE" -o "$THM4169_TMP/audit"
clang++ -std=c++20 -Wall -Wextra -Werror -O1 -g \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  "$THM4169_SOURCE" -o "$THM4169_TMP/audit_san"

"$THM4169_TMP/audit" --selftest > "$THM4169_TMP/selftest.out"
"$THM4169_TMP/audit_O0" --selftest > "$THM4169_TMP/selftest_O0.out"
ASAN_OPTIONS=detect_leaks=0 UBSAN_OPTIONS=print_stacktrace=1 \
  "$THM4169_TMP/audit_san" --selftest > "$THM4169_TMP/selftest_san.out" \
  2> "$THM4169_TMP/selftest_san.err"
cmp "$THM4169_TMP/selftest.out" "$THM4169_TMP/selftest_O0.out"
cmp "$THM4169_TMP/selftest.out" "$THM4169_TMP/selftest_san.out"
[[ ! -s "$THM4169_TMP/selftest_san.err" ]]
cmp "$THM4169_TMP/selftest.out" \
  "$THM4169_RESULTS/tournament_prime_parent_quartic_transfer_thm4169_selftest.out"

"$THM4169_TMP/audit" --full-q10-poly "$THM4169_PARENT" \
  > "$THM4169_TMP/full_q10_poly.out"
cmp "$THM4169_TMP/full_q10_poly.out" \
  "$THM4169_RESULTS/tournament_prime_parent_quartic_transfer_thm4169_full_q10_poly.out"

python3 -B "$THM4169_PYTHON" > "$THM4169_TMP/independent_python.out"
python3 -B -O "$THM4169_PYTHON" > "$THM4169_TMP/independent_python_O.out"
PYTHONHASHSEED=4169 python3 -B "$THM4169_PYTHON" \
  > "$THM4169_TMP/independent_python_hashseed.out"
cmp "$THM4169_TMP/independent_python.out" \
  "$THM4169_TMP/independent_python_O.out"
cmp "$THM4169_TMP/independent_python.out" \
  "$THM4169_TMP/independent_python_hashseed.out"
cmp "$THM4169_TMP/independent_python.out" \
  "$THM4169_RESULTS/tournament_prime_parent_quartic_transfer_thm4169_independent_audit.out"

gentourng -q 10 | "$THM4169_TMP/audit" --stream-prime 10 \
  > "$THM4169_TMP/order10_prime_census.out"
cmp "$THM4169_TMP/order10_prime_census.out" \
  "$THM4169_RESULTS/tournament_prime_parent_quartic_transfer_thm4169_order10_prime_census.out"

"$THM4169_TMP/audit" --emit-prime-extensions "$THM4169_PARENT" \
  > "$THM4169_TMP/extensions.d6"
"$THM4169_TMP/audit" --list-prime-patterns "$THM4169_PARENT" \
  > "$THM4169_TMP/patterns.txt"
shortg -u "$THM4169_TMP/extensions.d6" >/dev/null \
  2> "$THM4169_TMP/shortg_count.err"
shortg -v -q "$THM4169_TMP/extensions.d6" \
  "$THM4169_TMP/extensions.canon.d6" >/dev/null \
  2> "$THM4169_TMP/shortg_mapping.err"

THM4169_READ="$(awk '/graphs read/{print $2}' "$THM4169_TMP/shortg_count.err")"
THM4169_PRODUCED="$(awk '/graphs produced/{print $2}' "$THM4169_TMP/shortg_count.err")"
awk -F: -v pairs="$THM4169_TMP/doubleton_indices.txt" '
  /^[[:space:]]*[0-9]+[[:space:]]*:/ {
    count=0; first=""; second="";
    fields=split($2, a, /[[:space:]]+/);
    for (i=1; i<=fields; ++i) if (a[i] ~ /^[0-9]+$/) {
      ++count; if (count==1) first=a[i]; if (count==2) second=a[i];
    }
    ++rows; inputs+=count;
    if (count==2) { ++doubletons; print first, second > pairs; }
    else if (count!=1) ++bad;
  }
  END { print rows+0, inputs+0, doubletons+0, bad+0; }
' "$THM4169_TMP/shortg_mapping.err" > "$THM4169_TMP/mapping_stats.txt"
read -r THM4169_ROWS THM4169_INPUTS THM4169_DOUBLETONS THM4169_BAD \
  < "$THM4169_TMP/mapping_stats.txt"

: > "$THM4169_TMP/doubleton_patterns.txt"
while read -r THM4169_LEFT THM4169_RIGHT; do
  THM4169_LEFT_PATTERN="$(awk -v key="$THM4169_LEFT" '$1==key{print $2}' \
    "$THM4169_TMP/patterns.txt")"
  THM4169_RIGHT_PATTERN="$(awk -v key="$THM4169_RIGHT" '$1==key{print $2}' \
    "$THM4169_TMP/patterns.txt")"
  echo "$THM4169_LEFT_PATTERN~$THM4169_RIGHT_PATTERN" \
    >> "$THM4169_TMP/doubleton_patterns.txt"
done < "$THM4169_TMP/doubleton_indices.txt"
THM4169_PAIRS="$(paste -sd, "$THM4169_TMP/doubleton_patterns.txt")"

[[ "$THM4169_READ" == 1002 && "$THM4169_PRODUCED" == 1000 ]]
[[ "$THM4169_ROWS" == 1000 && "$THM4169_INPUTS" == 1002 ]]
[[ "$THM4169_DOUBLETONS" == 2 && "$THM4169_BAD" == 0 ]]
[[ "$THM4169_PAIRS" == "336~432,368~400" ]]

THM4169_INPUT_HASH="$(sha256sum "$THM4169_TMP/extensions.d6" | awk '{print $1}')"
THM4169_CANON_HASH="$(sha256sum "$THM4169_TMP/extensions.canon.d6" | awk '{print $1}')"
THM4169_PATTERN_HASH="$(sha256sum "$THM4169_TMP/patterns.txt" | awk '{print $1}')"
[[ "$THM4169_INPUT_HASH" == 32248819c74a46d963685af521d36caf5609b519b8b3421d5b46ad91d32f8128 ]]
[[ "$THM4169_CANON_HASH" == 0974bf20679943a5c41ff10cf22d20c46703e3472fec9ea62ac37136a07df8df ]]
[[ "$THM4169_PATTERN_HASH" == 0a6f1ee2ddcc77c19c92a55f1ae78e316c21da612100c2ebe58e0edcb6238b6a ]]

{
  echo "parent=$THM4169_PARENT prime_patterns=1002 root_preserving_orbits=1002 unrooted_classes=1000"
  echo "doubletons=$THM4169_PAIRS"
  echo "extensions_d6_sha256=$THM4169_INPUT_HASH"
  echo "canonical_d6_sha256=$THM4169_CANON_HASH"
  echo "patterns_sha256=$THM4169_PATTERN_HASH"
  echo "status=PASS"
} > "$THM4169_TMP/unrooted.out"
cmp "$THM4169_TMP/unrooted.out" \
  "$THM4169_RESULTS/tournament_prime_parent_quartic_transfer_thm4169_unrooted.out"

sha256sum \
  "$THM4169_TMP/selftest.out" \
  "$THM4169_TMP/full_q10_poly.out" \
  "$THM4169_TMP/independent_python.out" \
  "$THM4169_TMP/order10_prime_census.out" \
  "$THM4169_TMP/unrooted.out"
echo "THM-4169 rebuild status=PASS"
