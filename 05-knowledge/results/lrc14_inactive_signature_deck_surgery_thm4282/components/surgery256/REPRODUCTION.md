# Reproduction

Run from this packet root with Bash, a C++20 compiler, Python 3, Ruby, and
`sha256sum` available.

```bash
ROOT="$PWD"
SRC="$ROOT/src"
IN="$ROOT/inputs"
RES="$ROOT/results"
TMP=$(mktemp -d)
CXX=${CXX:-clang++}
BASE_FLAGS=(-std=c++20 -pthread -I"$SRC" -I"$ROOT/deps")
```

## 1. Primary exact deck/body/common replay

```bash
for MODE in O0 O3 UBSAN; do
  case "$MODE" in
    O0) FLAGS=(-O0) ;;
    O3) FLAGS=(-O3 -DNDEBUG) ;;
    UBSAN) FLAGS=(-O1 -g -fsanitize=undefined \
                  -fno-sanitize-recover=all) ;;
  esac
  "$CXX" "${BASE_FLAGS[@]}" "${FLAGS[@]}" \
    "$SRC/verify_surgery_deck.cpp" -o "$TMP/primary_$MODE"
  "$TMP/primary_$MODE" "$IN/joint421_masks.txt" \
    "$IN/vulnerable_bodies.csv" "$IN/full_signatures_primary.csv" \
    "$TMP/deck_$MODE.txt" "$TMP/common_$MODE.csv" \
    > "$TMP/primary_$MODE.out" 2> "$TMP/primary_$MODE.err"
done

cmp "$TMP/primary_O0.out" "$TMP/primary_O3.out"
cmp "$TMP/primary_O3.out" "$TMP/primary_UBSAN.out"
cmp "$TMP/primary_O3.out" "$RES/primary_surgery_replay.out"
cmp "$TMP/deck_O0.txt" "$TMP/deck_O3.txt"
cmp "$TMP/deck_O3.txt" "$TMP/deck_UBSAN.txt"
cmp "$TMP/deck_O3.txt" "$RES/surgery422_masks.txt"
cmp "$TMP/common_O0.csv" "$TMP/common_O3.csv"
cmp "$TMP/common_O3.csv" "$TMP/common_UBSAN.csv"
cmp "$TMP/common_O3.csv" "$RES/surgery_common.csv"
test ! -s "$TMP/primary_O0.err"
test ! -s "$TMP/primary_O3.err"
cmp "$TMP/primary_UBSAN.err" \
  "$RES/primary_surgery_replay_UBSan.err"
```

This replay exhausts all 14,307,150 nine-bodies. Its common-set computation
uses the exact primitive-cocycle/atom implementation routed through the three
packet-local source dependencies.

## 2. Detached literal-wall replay

```bash
for MODE in O0 O3 UBSAN; do
  case "$MODE" in
    O0) FLAGS=(-O0) ;;
    O3) FLAGS=(-O3 -DNDEBUG) ;;
    UBSAN) FLAGS=(-O1 -g -fsanitize=undefined \
                  -fno-sanitize-recover=all) ;;
  esac
  mkdir "$TMP/literal_$MODE"
  "$CXX" -std=c++20 "${FLAGS[@]}" \
    "$SRC/surgery188_detached_literal_audit.cpp" \
    -o "$TMP/literal_bin_$MODE"
  "$TMP/literal_bin_$MODE" "$RES/surgery422_masks.txt" \
    "$RES/surgery_common.csv" "$IN/post_thm4281_residual.csv" \
    "$IN/vulnerable_obligations_full_carrier.csv" \
    "$TMP/literal_$MODE" > "$TMP/literal_$MODE.out" \
    2> "$TMP/literal_$MODE.err"
done

cmp "$TMP/literal_O0.out" "$TMP/literal_O3.out"
cmp "$TMP/literal_O3.out" "$TMP/literal_UBSAN.out"
cmp "$TMP/literal_O3.out" "$RES/detached_literal_surgery188.out"
for FILE in common_literal_rows.csv exposed_body_replacements.csv \
            residual_complement.csv verified_common.csv; do
  cmp "$TMP/literal_O0/$FILE" "$TMP/literal_O3/$FILE"
  cmp "$TMP/literal_O3/$FILE" "$TMP/literal_UBSAN/$FILE"
done
cmp "$TMP/literal_O3/common_literal_rows.csv" \
  "$RES/detached_literal_common_rows.csv"
cmp "$TMP/literal_O3/exposed_body_replacements.csv" \
  "$RES/detached_literal_exposed_body_replacements.csv"
cmp "$TMP/literal_O3/residual_complement.csv" \
  "$RES/detached_literal_residual_complement24035.csv"
cmp "$TMP/literal_O3/verified_common.csv" "$RES/surgery_common.csv"
test ! -s "$TMP/literal_O0.err"
test ! -s "$TMP/literal_O3.err"
cmp "$TMP/literal_UBSAN.err" \
  "$RES/detached_literal_surgery188_UBSan.err"
```

This source imports no project code. It directly checks all `188*422=79,336`
claimed common-deck activity cells, with zero equalities, and independently
reconstructs the exhaustive body decomposition.

## 3. Proof-graph union and independent set audit

```bash
for MODE in normal optimized; do
  mkdir "$TMP/graph_$MODE"
  if test "$MODE" = normal; then PY=(python3); else PY=(python3 -O); fi
  "${PY[@]}" "$SRC/proof_graph_consequence.py" \
    --surgery "$RES/surgery_common.csv" \
    --gain "$IN/prior_combined_gain586.csv" \
    --prior-residual "$IN/prior_combined_residual23637.csv" \
    --out "$TMP/graph_$MODE" > "$TMP/graph_$MODE.stdout"
done

cmp "$TMP/graph_normal.stdout" "$TMP/graph_optimized.stdout"
for FILE in audit.out carrier90.csv combined_union850.csv \
            exclusive_carrier85.csv exclusive_gain577.csv \
            exclusive_surgery174.csv final_residual23373.csv \
            overlap_gain_carrier0.csv overlap_surgery_carrier5.csv \
            overlap_surgery_gain9.csv overlap_triple0.csv; do
  cmp "$TMP/graph_normal/$FILE" "$TMP/graph_optimized/$FILE"
done

cmp "$TMP/graph_normal/audit.out" "$RES/proof_graph_consequence.out"
cmp "$TMP/graph_normal/carrier90.csv" "$RES/prior_carrier90.csv"
cmp "$TMP/graph_normal/combined_union850.csv" \
  "$RES/deletion_union850.csv"
cmp "$TMP/graph_normal/exclusive_surgery174.csv" \
  "$RES/surgery_novel174.csv"
cmp "$TMP/graph_normal/final_residual23373.csv" \
  "$RES/final_residual23373.csv"
cmp "$TMP/graph_normal/overlap_surgery_gain9.csv" \
  "$RES/surgery_overlap_prior_gain9.csv"
cmp "$TMP/graph_normal/overlap_surgery_carrier5.csv" \
  "$RES/surgery_overlap_carrier5.csv"

ruby "$SRC/independent_overlap_audit.rb" \
  "$RES/surgery_common.csv" "$IN/prior_combined_gain586.csv" \
  "$IN/prior_combined_residual23637.csv" "$TMP/graph_normal" \
  > "$TMP/independent_ruby.out"
cmp "$TMP/independent_ruby.out" "$RES/independent_overlap_audit.out"
```

## 4. Frozen byte ledger

```bash
sha256sum -c SHA256SUMS
```
