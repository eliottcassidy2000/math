# Reproduction

Run from the repository root with Bash, a C++20 compiler, and Python 3.

```bash
set -euo pipefail
ROOT=$(git rev-parse --show-toplevel)
SRC="$ROOT/04-computation/lrc14_endpoint638_carrier_response_thm4283"
OUT="$ROOT/05-knowledge/results/lrc14_endpoint638_carrier_response_thm4283"
P4282="$ROOT/05-knowledge/results/lrc14_inactive_signature_deck_surgery_thm4282"
S4282="$ROOT/04-computation/lrc14_inactive_signature_deck_surgery_thm4282"
JOINT="$P4282/components/index367/results/original_joint421_deck.txt"
BASE="$P4282/components/surgery520/inputs/final8951_carrier.txt"
ADDITIONS="$S4282/carrier_through648_additions.txt"
POST="$P4282/results/final_residual23373.csv"
TMP=$(mktemp -d)
CXX=${CXX:-g++}
```

Verify the byte ledger:

```bash
cd "$ROOT"
sha256sum -c "$OUT/SHA256SUMS"
```

Compile and run the primary response quotient:

```bash
"$CXX" -std=c++20 -O3 -DNDEBUG -I"$ROOT" \
  "$SRC/endpoint638_carrier_response_primary.cpp" \
  -o "$TMP/primary"
"$TMP/primary" "$JOINT" "$BASE" "$ADDITIONS" \
  > "$TMP/primary.out"
cmp "$TMP/primary.out" "$OUT/primary_response_audit.out"
```

Compile and run the exact activity-margin audit:

```bash
"$CXX" -std=c++20 -O3 -DNDEBUG -I"$ROOT" \
  "$SRC/endpoint638_activity_margins.cpp" \
  -o "$TMP/margins"
"$TMP/margins" "$JOINT" "$BASE" "$ADDITIONS" \
  > "$TMP/margins.out"
cmp "$TMP/margins.out" "$OUT/activity_margins.out"
```

Compile and run the detached literal audit:

```bash
"$CXX" -std=c++20 -O3 -DNDEBUG -I"$ROOT" \
  -I"$ROOT/04-computation" \
  "$SRC/endpoint638_detached_literal_audit.cpp" \
  -o "$TMP/literal"
"$TMP/literal" "$BASE" "$ADDITIONS" > "$TMP/literal.out"
cmp "$TMP/literal.out" "$OUT/detached_literal_audit.out"
grep -q '^LEDGER cbe0c99a6d552e23 TOTAL_FAILURES 0$' \
  "$TMP/literal.out"
```

Replay the exact proof-graph consequence:

```bash
python3 -B "$SRC/proof_graph_consequence.py" \
  --post4282 "$POST" --closed "$OUT/closed63.csv" \
  --base "$BASE" --additions "$ADDITIONS" \
  > "$TMP/consequence.out"
cmp "$TMP/consequence.out" "$OUT/proof_graph_consequence.out"

python3 -B -O "$SRC/proof_graph_consequence.py" \
  --post4282 "$POST" --closed "$OUT/closed63.csv" \
  --base "$BASE" --additions "$ADDITIONS" \
  > "$TMP/consequence-opt.out"
cmp "$TMP/consequence-opt.out" "$OUT/proof_graph_consequence.out"
```

For an additional compiler-mode control, repeat the three C++ builds at
`-O0`; all consequence-bearing transcripts must compare byte-for-byte.
