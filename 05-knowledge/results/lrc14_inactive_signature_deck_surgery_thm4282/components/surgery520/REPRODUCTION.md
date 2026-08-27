# Reproduction

Run with Bash from a checkout containing the promoted THM-4281 dependencies.
The packet sources use repository-relative includes.  C++20 and Python 3 are
required.  The full active-response TSV is generated only under the temporary
directory.

```bash
set -euo pipefail
ROOT=$(git rev-parse --show-toplevel)
PACKET=/private/tmp/thm4281-signature-surgery-520-663/thm4282-520-663-signature-surgery-promotion-packet
SRC="$PACKET/src"
IN="$PACKET/inputs"
RES="$PACKET/results"
TMP=$(mktemp -d)
CXX=${CXX:-clang++}
BASE=(-std=c++20 -pthread -I"$ROOT")
```

## 1. Independent lost-body scan

```bash
"$CXX" -std=c++20 -O3 -DNDEBUG "$SRC/lost_body_probe.cpp" \
  -o "$TMP/lost_probe"
"$TMP/lost_probe" "$IN/original_joint421_deck.txt" \
  "$TMP/lost_probe.csv" > "$TMP/lost_probe.out" \
  2> "$TMP/lost_probe.err"
cmp "$TMP/lost_probe.csv" "$RES/lost_bodies_probe.csv"
cmp "$TMP/lost_probe.out" "$RES/lost_body_probe_O3.out"
cmp "$TMP/lost_probe.err" "$RES/lost_body_probe_O3.err"
test ! -s "$TMP/lost_probe.err"
```

## 2. Primary surgery atlas under O3, O0, and UBSan

The O3 active TSV is retained transiently for later independent checks.  O0
and UBSan overwrite a separate transient check path.

```bash
for MODE in O3 O0 UBSAN; do
  case "$MODE" in
    O3)
      FLAGS=(-O3 -DNDEBUG)
      ACTIVE="$TMP/active.primary.tsv"
      EXPECT_OUT="$RES/signature_surgery_atlas.out"
      EXPECT_ERR="$RES/signature_surgery_atlas.err"
      ;;
    O0)
      FLAGS=(-O0)
      ACTIVE="$TMP/active.check.tsv"
      EXPECT_OUT="$RES/signature_surgery_atlas_O0.out"
      EXPECT_ERR="$RES/signature_surgery_atlas_O0.err"
      ;;
    UBSAN)
      FLAGS=(-O1 -g -fsanitize=undefined -fno-sanitize-recover=undefined)
      ACTIVE="$TMP/active.check.tsv"
      EXPECT_OUT="$RES/signature_surgery_atlas_ubsan.out"
      EXPECT_ERR="$RES/signature_surgery_atlas_ubsan.err"
      ;;
  esac
  "$CXX" "${BASE[@]}" "${FLAGS[@]}" \
    "$SRC/signature_surgery_atlas.cpp" -o "$TMP/atlas_$MODE"
  "$TMP/atlas_$MODE" "$IN/full_signatures_primary.csv" \
    "$IN/original_joint421_deck.txt" "$TMP/lost.check.csv" \
    "$TMP/classes.check.tsv" "$ACTIVE" "$TMP/deck.check.txt" \
    > "$TMP/atlas_$MODE.out" 2> "$TMP/atlas_$MODE.err"
  cmp "$TMP/atlas_$MODE.out" "$EXPECT_OUT"
  cmp "$TMP/atlas_$MODE.err" "$EXPECT_ERR"
  test ! -s "$TMP/atlas_$MODE.err"
  cmp "$TMP/lost.check.csv" "$RES/lost_bodies_exact.csv"
  cmp "$TMP/classes.check.tsv" "$RES/response_classes.tsv"
  cmp "$TMP/deck.check.txt" "$RES/rebuilt_deck.txt"
  if [[ "$MODE" != O3 ]]; then
    cmp "$ACTIVE" "$TMP/active.primary.tsv"
  fi
done
test "$(sha256sum "$TMP/active.primary.tsv" | cut -d' ' -f1)" = \
  dbdd1dab8c1e0d4de03866d20155d2a985929bd9e279b3fb15eb643c351a851e
```

## 3. Detached literal response-universe completeness

This path does not call the primary primitive-pair/cocycle activity builder.

```bash
for MODE in O3 O0 UBSAN; do
  case "$MODE" in
    O3)
      FLAGS=(-O3 -DNDEBUG)
      EXPECT_OUT="$RES/literal_response_universe_audit.out"
      EXPECT_ERR="$RES/literal_response_universe_audit.err"
      ;;
    O0)
      FLAGS=(-O0)
      EXPECT_OUT="$RES/literal_response_universe_audit_O0.out"
      EXPECT_ERR="$RES/literal_response_universe_audit_O0.err"
      ;;
    UBSAN)
      FLAGS=(-O1 -g -fsanitize=undefined -fno-sanitize-recover=undefined)
      EXPECT_OUT="$RES/literal_response_universe_audit_ubsan.out"
      EXPECT_ERR="$RES/literal_response_universe_audit_ubsan.err"
      ;;
  esac
  "$CXX" "${BASE[@]}" "${FLAGS[@]}" \
    "$SRC/literal_response_universe_audit.cpp" -o "$TMP/literal_$MODE"
  "$TMP/literal_$MODE" "$RES/lost_bodies_exact.csv" \
    "$TMP/classes.literal.tsv" "$TMP/active.literal.tsv" \
    > "$TMP/literal_$MODE.out" 2> "$TMP/literal_$MODE.err"
  cmp "$TMP/literal_$MODE.out" "$EXPECT_OUT"
  cmp "$TMP/literal_$MODE.err" "$EXPECT_ERR"
  test ! -s "$TMP/literal_$MODE.err"
  cmp "$TMP/classes.literal.tsv" "$RES/response_classes.tsv"
  cmp "$TMP/active.literal.tsv" "$TMP/active.primary.tsv"
done
```

## 4. Rebuilt-deck target and neighbour checks

```bash
"$CXX" "${BASE[@]}" -O3 -DNDEBUG \
  "$SRC/rebuilt_pair_activity.cpp" -o "$TMP/pair"
for R in 663 589; do
  "$TMP/pair" "$RES/rebuilt_deck.txt" 520 "$R" \
    > "$TMP/pair_520_$R.out" 2> "$TMP/pair_520_$R.err"
  cmp "$TMP/pair_520_$R.out" "$RES/rebuilt_pair_520_$R.out"
  cmp "$TMP/pair_520_$R.err" "$RES/rebuilt_pair_520_$R.err"
  test ! -s "$TMP/pair_520_$R.err"
done
```

## 5. Global primary/literal scans

The primary program emits deterministic semantic files plus progress on
stderr.  Progress stderr is noncanonical and excluded.

```bash
"$CXX" "${BASE[@]}" -O3 -DNDEBUG \
  "$SRC/rebuilt_global_primary.cpp" -o "$TMP/global_primary"
"$TMP/global_primary" "$IN/post_thm4271_residual.csv" \
  "$RES/rebuilt_deck.txt" "$TMP/common.primary.csv" \
  "$TMP/layers.primary.csv" "$TMP/hostile.primary.csv" \
  > "$TMP/global_primary.out" 2> "$TMP/global_primary.progress.err"
cmp "$TMP/global_primary.out" "$RES/rebuilt_global_primary.out"
cmp "$TMP/common.primary.csv" "$RES/rebuilt_common_primary.csv"
cmp "$TMP/layers.primary.csv" "$RES/rebuilt_common_primary_layers.csv"
cmp "$TMP/hostile.primary.csv" "$RES/rebuilt_primary_hostile.csv"

"$CXX" "${BASE[@]}" -O3 -DNDEBUG \
  "$SRC/rebuilt_global_literal.cpp" -o "$TMP/global_literal"
"$TMP/global_literal" "$IN/post_thm4271_residual.csv" \
  "$RES/rebuilt_deck.txt" "$TMP/common.literal.csv" \
  > "$TMP/global_literal.out" 2> "$TMP/global_literal.err"
cmp "$TMP/global_literal.out" "$RES/rebuilt_global_literal.out"
cmp "$TMP/global_literal.err" "$RES/rebuilt_global_literal.err"
test ! -s "$TMP/global_literal.err"
cmp "$TMP/common.literal.csv" "$TMP/common.primary.csv"
```

## 6. Detached ledger/dual certificate

```bash
for PY in python3 "python3 -O"; do
  $PY "$SRC/verify_surgery_certificate.py" \
    --original-deck "$IN/original_joint421_deck.txt" \
    --rebuilt-deck "$RES/rebuilt_deck.txt" \
    --lost "$RES/lost_bodies_exact.csv" \
    --classes "$RES/response_classes.tsv" \
    --active-masks "$TMP/active.primary.tsv" \
    --primary-common "$TMP/common.primary.csv" \
    --literal-common "$TMP/common.literal.csv" \
    > "$TMP/verify.$(echo "$PY" | tr ' ' '_').out"
done
cmp "$TMP/verify.python3.out" "$RES/verify_surgery_certificate.out"
cmp "$TMP/verify.python3_-O.out" "$RES/verify_surgery_certificate_opt.out"
cmp "$TMP/verify.python3.out" "$TMP/verify.python3_-O.out"
```

## 7. Proof-graph consequence and conditional index-367 union

```bash
for OPT in normal optimized; do
  if [[ "$OPT" == normal ]]; then PY=(python3); else PY=(python3 -O); fi
  "${PY[@]}" "$SRC/combined_consequence.py" \
    --signatures "$IN/full_signatures_primary.csv" \
    --old-common "$IN/canonical_joint421_common.csv" \
    --new-common "$TMP/common.primary.csv" \
    --post4281 "$IN/post_thm4281_residual.csv" \
    --emit-520-gain "$TMP/gain520.$OPT.csv" \
    --emit-post520 "$TMP/post520.$OPT.csv" \
    --emit-367-fibre "$TMP/fibre367.$OPT.csv" \
    --emit-combined-gain "$TMP/combined.$OPT.csv" \
    --emit-final "$TMP/final.$OPT.csv" \
    > "$TMP/consequence.$OPT.out"
done
cmp "$TMP/consequence.normal.out" "$RES/combined_consequence.out"
cmp "$TMP/consequence.optimized.out" "$RES/combined_consequence_opt.out"
cmp "$TMP/consequence.normal.out" "$TMP/consequence.optimized.out"
cmp "$TMP/gain520.normal.csv" "$RES/proof_graph_new.csv"
cmp "$TMP/post520.normal.csv" "$RES/post_surgery_residual.csv"
cmp "$TMP/fibre367.normal.csv" "$RES/index367_fibre.csv"
cmp "$TMP/combined.normal.csv" "$RES/combined_gain_520_367.csv"
cmp "$TMP/final.normal.csv" "$RES/combined_residual_520_367.csv"
for STEM in gain520 post520 fibre367 combined final; do
  cmp "$TMP/$STEM.normal.csv" "$TMP/$STEM.optimized.csv"
done
```

## 8. Dormant-carrier reconstruction cross-filter

```bash
for OPT in normal optimized; do
  if [[ "$OPT" == normal ]]; then PY=(python3); else PY=(python3 -O); fi
  "${PY[@]}" "$SRC/cross_filter_dormant_carrier.py" \
    "$IN/final8951_carrier.txt" "$TMP/active.primary.tsv" \
    "$RES/lost_bodies_exact.csv" "$TMP/carrier.$OPT.tsv" \
    "$TMP/hits.$OPT.tsv" > "$TMP/cross.$OPT.out" \
    2> "$TMP/cross.$OPT.err"
done
cmp "$TMP/cross.normal.out" "$RES/cross_filter_dormant_carrier.out"
cmp "$TMP/cross.optimized.out" "$RES/cross_filter_dormant_carrier_opt.out"
cmp "$TMP/cross.normal.err" "$RES/cross_filter_dormant_carrier.err"
cmp "$TMP/cross.optimized.err" "$RES/cross_filter_dormant_carrier_opt.err"
test ! -s "$TMP/cross.normal.err"
test ! -s "$TMP/cross.optimized.err"
cmp "$TMP/carrier.normal.tsv" "$RES/dormant_carrier_active_responses.tsv"
cmp "$TMP/hits.normal.tsv" "$RES/dormant_carrier_obligation_hits.tsv"
cmp "$TMP/carrier.normal.tsv" "$TMP/carrier.optimized.tsv"
cmp "$TMP/hits.normal.tsv" "$TMP/hits.optimized.tsv"
```

## 9. Artifact ledger

```bash
cd "$PACKET"
sha256sum -c SHA256SUMS
```
