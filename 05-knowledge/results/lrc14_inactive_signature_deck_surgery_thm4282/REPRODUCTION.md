# Reproduction

Run with Bash from the repository root. C++20, Python 3, Ruby, and
`sha256sum` are used by the full component suite.

```bash
set -euo pipefail
ROOT=$(git rev-parse --show-toplevel)
PACKET="$ROOT/05-knowledge/results/lrc14_inactive_signature_deck_surgery_thm4282"
SRC="$ROOT/04-computation/lrc14_inactive_signature_deck_surgery_thm4282"
TMP=$(mktemp -d)
CXX=${CXX:-clang++}
JOINT="$PACKET/components/index367/results/original_joint421_deck.txt"
BASE="$PACKET/components/surgery520/inputs/final8951_carrier.txt"
```

## 1. Root and component byte ledgers

```bash
cd "$ROOT"
sha256sum -c "$PACKET/SHA256SUMS"
for COMPONENT in index367 surgery520 carrier663 carrier645 surgery256; do
  (cd "$PACKET/components/$COMPONENT" && sha256sum -c SHA256SUMS)
done
```

The component documents route the full independent O0/O3/UBSan and literal
replays. Their old absolute `PACKET` or `REPO_ROOT` assignments should be
replaced by the corresponding directory under `$PACKET/components` and by
`$ROOT`.

## 2. Source-pinned primary carrier-band replay

```bash
"$CXX" -std=c++20 -O3 -DNDEBUG -I"$ROOT" \
  "$SRC/final8951_joint_exposure_scan.cpp" \
  -o "$TMP/final8996_primary"

"$TMP/final8996_primary" "$JOINT" "$BASE" \
  "$SRC/carrier_closed_663_645_pairs.csv" \
  "$SRC/carrier_through648_additions.txt" \
  "$TMP/failures.csv" > "$TMP/final8996.out" \
  2> "$TMP/final8996.err"

cmp "$TMP/final8996.out" "$PACKET/results/final8996_primary_audit.out"
cmp "$TMP/final8996.err" "$PACKET/results/final8996_primary_audit.err"
cmp "$TMP/failures.csv" "$PACKET/results/final8996_primary_failures.csv"
test ! -s "$TMP/final8996.err"
test "$(grep -c '^PAIR ' "$TMP/final8996.out")" = 90
grep -q '^PAIR_LEDGER_FNV 4ad413959cd58a9d TOTAL_FAILURES 0$' \
  "$TMP/final8996.out"
grep -q '^VERDICT PASS ALL_PAIRS_CLOSED$' "$TMP/final8996.out"
```

The pair-ledger FNV in the transcript includes each pair's complete audit
record. The input pair-only FNV pinned by the source is
`942995bee7469430`.

## 3. Hostile completeness controls

```bash
head -n 30 "$SRC/carrier_closed_663_645_pairs.csv" > "$TMP/short_pairs.csv"
head -n 42 "$SRC/carrier_through648_additions.txt" > "$TMP/short_additions.txt"

set +e
"$TMP/final8996_primary" "$JOINT" "$BASE" \
  "$TMP/short_pairs.csv" "$SRC/carrier_through648_additions.txt" \
  > "$TMP/reject_pairs.out" 2> "$TMP/reject_pairs.err"
PAIR_STATUS=$?
"$TMP/final8996_primary" "$JOINT" "$BASE" \
  "$SRC/carrier_closed_663_645_pairs.csv" "$TMP/short_additions.txt" \
  > "$TMP/reject_additions.out" 2> "$TMP/reject_additions.err"
ADDITION_STATUS=$?
set -e

test "$PAIR_STATUS" = 1
test "$ADDITION_STATUS" = 1
cmp "$TMP/reject_pairs.out" "$PACKET/controls/reject_short_pair_ledger.out"
cmp "$TMP/reject_pairs.err" "$PACKET/controls/reject_short_pair_ledger.err"
cmp "$TMP/reject_additions.out" \
  "$PACKET/controls/reject_short_addition_ledger.out"
cmp "$TMP/reject_additions.err" \
  "$PACKET/controls/reject_short_addition_ledger.err"
```

## 4. Endpoint-650 complete response and explicit greedy cover

```bash
"$CXX" -std=c++20 -O3 -DNDEBUG -I"$ROOT" \
  "$SRC/endpoint650_large_response_greedy.cpp" \
  -o "$TMP/endpoint650_greedy"

"$TMP/endpoint650_greedy" \
  "$PACKET/results/endpoint650_failures358.csv" "$BASE" \
  "$SRC/carrier_continuation_masks.txt" "$TMP/greedy31.txt" \
  > "$TMP/greedy31.out" 2> "$TMP/greedy31.err"

cmp "$TMP/greedy31.out" "$PACKET/results/endpoint650_greedy31.out"
cmp "$TMP/greedy31.txt" "$PACKET/results/endpoint650_greedy31_masks.txt"
test ! -s "$TMP/greedy31.err"
tail -n 31 "$SRC/carrier_through650_additions.txt" > "$TMP/stated31.txt"
cmp "$TMP/greedy31.txt" "$TMP/stated31.txt"
```

This proves a packing lower bound of 20 and an explicit upper bound of 31. It
does not prove the minimum is 31.

## 5. Exact proof-graph replay

```bash
python3 "$SRC/proof_graph_consequence.py" \
  --post4281 "$PACKET/components/index367/results/post_thm4281_residual.csv" \
  --combined586 \
    "$PACKET/components/surgery520/results/combined_gain_520_367.csv" \
  --postcombined \
    "$PACKET/components/surgery520/results/combined_residual_520_367.csv" \
  --surgery256 "$PACKET/components/surgery256/results/surgery_common.csv" \
  --emit-carrier "$TMP/carrier90.csv" \
  --emit-surgery256-net "$TMP/surgery256_net174.csv" \
  --emit-union "$TMP/union850.csv" \
  --emit-final "$TMP/final23373.csv" \
  > "$TMP/proof_graph.out"

cmp "$TMP/proof_graph.out" "$PACKET/results/proof_graph_consequence.out"
cmp "$TMP/carrier90.csv" "$PACKET/results/carrier90.csv"
cmp "$TMP/union850.csv" "$PACKET/results/deletion_union850.csv"
cmp "$TMP/final23373.csv" "$PACKET/results/final_residual23373.csv"
```

For a second implementation, run the Ruby overlap audit in
`components/surgery256/REPRODUCTION.md`.
