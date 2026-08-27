# Reproduction

Run from a checkout containing the promoted THM-4281 dependencies.  The packet
sources use repository-relative includes.  Only C++20 and Python 3 are needed.

```bash
ROOT=$(git rev-parse --show-toplevel)
PACKET=/private/tmp/thm4281-signature-atlas/thm4282-index367-fibre-promotion-packet
SRC="$PACKET/src"
RES="$PACKET/results"
TMP=$(mktemp -d)
CXX=${CXX:-clang++}
CXXFLAGS=(-std=c++20 -pthread -I"$ROOT" -I"$ROOT/04-computation")
DECK="$RES/original_joint421_deck.txt"
```

## 1. Full order-free inactive signatures

```bash
"$CXX" "${CXXFLAGS[@]}" -O3 -DNDEBUG \
  "$SRC/full421_inactive_signatures_primary.cpp" -o "$TMP/signatures"
"$TMP/signatures" "$RES/post_thm4281_residual.csv" "$DECK" \
  "$TMP/signatures.csv" > "$TMP/signatures.out" \
  2> "$TMP/signatures.progress.err"
cmp "$TMP/signatures.csv" \
  "$RES/post_thm4281_full421_inactive_signatures.csv"
cmp "$TMP/signatures.out" \
  "$RES/full421_inactive_signatures_primary_O3.out"
```

The progress stderr is deterministic but deliberately noncanonical.

```bash
python3 "$SRC/analyze_full421_signatures.py" \
  "$RES/post_thm4281_full421_inactive_signatures.csv" "$DECK" \
  > "$TMP/analysis.out"
cmp "$TMP/analysis.out" "$RES/analyze_full421_signatures.out"
```

The exact 216 lower bound is solver-free: 215 masks have weight-one witnesses,
they leave only `(70,302)`, and that row has signature exactly `{78,368}`.

## 2. Exact private-body fibre

```bash
"$CXX" -std=c++20 -O3 -DNDEBUG \
  "$SRC/index367_private_bodies.cpp" -o "$TMP/private367"
"$TMP/private367" "$DECK" > "$TMP/private367.out"
cmp "$TMP/private367.out" "$RES/index367_private_bodies_O3.out"
```

## 3. Complete active response quotient at `(366,663)`

```bash
"$CXX" "${CXXFLAGS[@]}" -O3 -DNDEBUG \
  -I"$ROOT/04-computation/lrc14_joint421_global_common_carrier_thm4281" \
  "$SRC/endpoint663_index367_response_atlas.cpp" -o "$TMP/atlas"
"$TMP/atlas" "$DECK" > "$TMP/atlas.out"
cmp "$TMP/atlas.out" \
  "$RES/endpoint663_index367_response_atlas_O3_NDEBUG.out"
```

## 4. Detached literal fibre audit and sanitizer control

```bash
for MODE in O3 UBSAN; do
  case "$MODE" in
    O3) FLAGS=(-O3 -DNDEBUG) ;;
    UBSAN) FLAGS=(-O1 -g -fsanitize=undefined \
                  -fno-sanitize-recover=all) ;;
  esac
  "$CXX" "${CXXFLAGS[@]}" "${FLAGS[@]}" \
    "$SRC/index367_singleton_fibre_literal_audit.cpp" \
    -o "$TMP/literal_$MODE"
  "$TMP/literal_$MODE" "$DECK" > "$TMP/literal_$MODE.out" \
    2> "$TMP/literal_$MODE.err"
done
cmp "$TMP/literal_O3.out" "$TMP/literal_UBSAN.out"
cmp "$TMP/literal_O3.out" \
  "$RES/index367_singleton_fibre_literal_audit_O3_NDEBUG.out"
test ! -s "$TMP/literal_UBSAN.err"
cmp "$TMP/literal_UBSAN.err" \
  "$RES/index367_singleton_fibre_literal_audit_UBSan.err"
```

## 5. Proof-graph consequence

```bash
python3 "$SRC/proof_graph_consequence.py" \
  "$RES/post_thm4281_residual.csv" \
  "$RES/post_thm4281_full421_inactive_signatures.csv" \
  "$RES/index367_singleton_inactive_fibre.csv" \
  "$RES/original_joint421_deck.txt" \
  "$RES/index367_fibre_rebuilt_deck.txt" \
  "$TMP/post_fibre.csv" > "$TMP/proof_graph.out"
cmp "$TMP/post_fibre.csv" \
  "$RES/post_index367_singleton_fibre_residual.csv"
cmp "$TMP/proof_graph.out" "$RES/proof_graph_consequence.out"
```

## 6. Artifact ledger

```bash
cd "$PACKET"
sha256sum -c SHA256SUMS
```
