#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
for mode in normal optimized; do
  options=()
  if [[ $mode == optimized ]]; then options=(-O); fi
  for item in 'independent_counts.py:independent-counts.json' 'audit_comb.py:comb-normal.json' 'audit_covers.py:covers-normal.json' 'audit_deletions.py:deletion-normal.json'; do
    program=${item%%:*}
    receipt=${item#*:}
    python3 "${options[@]}" "$program" > "$tmp/result.json"
    cmp "$receipt" "$tmp/result.json"
  done
  python3 "${options[@]}" audit_support.py "${COMB_RESULT:-../cycle-to-path/RESULT.json}" > "$tmp/support.json"
  cmp support-independent.json "$tmp/support.json"
done
for optimization in O0 O2; do
  g++ -std=c++17 "-$optimization" -Wall -Wextra -pedantic audit_endpoint_equality.cpp -o "$tmp/endpoint"
  "$tmp/endpoint" "${CLASS7:-../../session212/cycle-levels/classes7.rows}" "${CLASS8:-../../session212/cycle-levels/classes8.rows}" > "$tmp/endpoint.jsonl"
  cmp endpoint-equality-audit.jsonl "$tmp/endpoint.jsonl"
done
echo 'PASS: all independent Session213 audit checks'
