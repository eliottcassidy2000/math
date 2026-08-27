# THM-4281 reproduction

Run from a checkout containing promoted THM-4266, THM-4276, and THM-4277.
The commands below use only the C++ standard library, Python 3, and the frozen
inputs in this artifact directory.

```bash
ROOT=$(git rev-parse --show-toplevel)
SRC="$ROOT/04-computation/lrc14_joint421_global_common_carrier_thm4281"
RES="$ROOT/05-knowledge/results/lrc14_joint421_global_common_carrier_thm4281"
TMP=$(mktemp -d)
CXX=${CXX:-clang++}
CXXFLAGS=(-std=c++20 -pthread -I"$ROOT" -I"$ROOT/04-computation")
DECK="$RES/joint_iter39_shrunk_forward_masks.txt"
POST4271="$RES/post_thm4271_residual.csv"
```

## 1. Primary global common-activity scan

```bash
"$CXX" "${CXXFLAGS[@]}" -O3 -DNDEBUG \
  "$SRC/scout_joint421_post4277_primary.cpp" -o "$TMP/global_primary"
"$TMP/global_primary" "$POST4271" "$DECK" \
  "$TMP/common.csv" "$TMP/common_layers.csv" "$TMP/primary_hostile.csv" \
  > "$TMP/global_primary.out" 2> "$TMP/global_primary.err"
cmp "$TMP/common.csv" "$RES/post_thm4277_joint421_common_primary.csv"
cmp "$TMP/common_layers.csv" \
  "$RES/post_thm4277_joint421_common_primary_layers.csv"
cmp "$TMP/primary_hostile.csv" \
  "$RES/post_thm4277_joint421_primary_hostile_controls.csv"
cmp "$TMP/global_primary.out" "$RES/scout_joint421_post4277_primary_O3.out"
```

The stderr file contains deterministic progress only and is not canonical.

## 2. Detached global literal partition and body audit

Compile and run the same source four ways:

```bash
for MODE in O0 O3 O3N UBSAN; do
  case "$MODE" in
    O0) FLAGS=(-O0) ;;
    O3) FLAGS=(-O3) ;;
    O3N) FLAGS=(-O3 -DNDEBUG) ;;
    UBSAN) FLAGS=(-O1 -g -fsanitize=undefined \
                  -fno-sanitize-recover=all) ;;
  esac
  "$CXX" "${CXXFLAGS[@]}" "${FLAGS[@]}" \
    "$SRC/verify_joint421_literal_post4277.cpp" \
    -o "$TMP/global_literal_$MODE"
  "$TMP/global_literal_$MODE" "$DECK" \
    "$RES/post_thm4277_joint421_common_primary.csv" \
    "$RES/post_thm4277_joint421_primary_hostile_controls.csv" \
    "$TMP/literal_hostile_$MODE.csv" \
    > "$TMP/global_literal_$MODE.out" \
    2> "$TMP/global_literal_$MODE.err"
done
cmp "$TMP/global_literal_O0.out" "$TMP/global_literal_O3.out"
cmp "$TMP/global_literal_O3.out" "$TMP/global_literal_O3N.out"
cmp "$TMP/global_literal_O3.out" "$TMP/global_literal_UBSAN.out"
cmp "$TMP/global_literal_O3N.out" \
  "$RES/verify_joint421_literal_post4277_O3_NDEBUG.out"
cmp "$TMP/literal_hostile_O0.csv" "$TMP/literal_hostile_O3.csv"
cmp "$TMP/literal_hostile_O3.csv" "$TMP/literal_hostile_O3N.csv"
cmp "$TMP/literal_hostile_O3.csv" "$TMP/literal_hostile_UBSAN.csv"
cmp "$TMP/literal_hostile_O3N.csv" \
  "$RES/post_thm4277_joint421_literal_hostile_controls_O3_NDEBUG.csv"
test ! -s "$TMP/global_literal_UBSAN.err"
cmp "$TMP/global_literal_UBSAN.err" \
  "$RES/verify_joint421_literal_post4277_ubsan.err"
```

### 2a. The body-redundant 420 core has the same global common set

```bash
CORE420="$RES/joint420_core_without_003c900c_masks.txt"
"$CXX" "${CXXFLAGS[@]}" -O3 -DNDEBUG \
  "$SRC/joint420_global_literal_activity_scout.cpp" -o "$TMP/global420"
"$TMP/global420" "$POST4271" "$CORE420" "$TMP/common420.csv" \
  > "$TMP/global420.out" 2> "$TMP/global420.err"
cmp "$TMP/global420.out" \
  "$RES/joint420_global_literal_activity_scout.out"
cmp "$TMP/global420.err" \
  "$RES/joint420_global_literal_activity_scout.err"
cmp "$TMP/common420.csv" \
  "$RES/post_thm4277_joint420_common_literal.csv"
cmp "$TMP/common420.csv" "$RES/post_thm4277_joint421_common_primary.csv"
```

This detached literal replay checks all 420 masks at every common row and
short-circuits only after the first literal failure on a noncommon row.  The
last comparison is the exact set-equality certificate; no order-sensitive
auxiliary row ledger is a theorem dependency.

## 3. Exact response atlases

```bash
for NAME in endpoint669_full_response_atlas \
            endpoint667_full_response_atlas \
            endpoint665_single_response \
            endpoint664_full_response_atlas; do
  "$CXX" "${CXXFLAGS[@]}" -O3 -DNDEBUG -I"$SRC" \
    "$SRC/$NAME.cpp" -o "$TMP/$NAME"
  "$TMP/$NAME" > "$TMP/$NAME.out" 2> "$TMP/$NAME.err"
done
cmp "$TMP/endpoint669_full_response_atlas.out" \
  "$RES/endpoint669_full_response_atlas_joint421_O3.out"
cmp "$TMP/endpoint667_full_response_atlas.out" \
  "$RES/endpoint667_full_response_atlas_O3.out"
cmp "$TMP/endpoint665_single_response.out" \
  "$RES/endpoint665_single_response_O3.out"
cmp "$TMP/endpoint664_full_response_atlas.out" \
  "$RES/endpoint664_full_response_atlas_O3.out"
```

Each atlas enumerates the complete `C(30,8)=5852925` repair universe, or the
complete `C(21,8)=203490` disjoint fibre in the one-obligation endpoint-665
case.  No solver transcript is a theorem dependency.

## 4. Primary endpoint layers

The common argument vector is:

```bash
SCAN_ARGS=(
  "$ROOT/05-knowledge/results/lrc14_endpoint_cascade_thm4254"
  "$ROOT/05-knowledge/results/lrc14_three_round_learned_carrier_thm4266/full_discovery_416_704_O3.semantic.out"
  "$ROOT/05-knowledge/results/lrc14_three_round_learned_carrier_thm4266/full_discovery_416_700_O3.semantic.out"
  "$ROOT/05-knowledge/results/lrc14_three_round_learned_carrier_thm4266/full_discovery_520_700_O3.semantic.out"
  "$ROOT/05-knowledge/results/lrc14_three_round_learned_carrier_thm4266/full_discovery_384_694_O3.semantic.out"
  "$ROOT/05-knowledge/results/lrc14_fourth_round_learned_carrier_thm4271/full_discovery_520_688_O3.semantic.out"
  "$POST4271" "$DECK"
)
"$CXX" "${CXXFLAGS[@]}" -O3 -DNDEBUG \
  "$SRC/scan_augmented_carrier_endpoint669.cpp" -o "$TMP/scan669"
"$TMP/scan669" "${SCAN_ARGS[@]}" > "$TMP/scan669.out" 2> "$TMP/scan669.err"
cmp "$TMP/scan669.out" \
  "$RES/scan_augmented_carrier_endpoint669_joint421_O3.out"
```

Compile the generalized layer source at each frozen carrier stage:

```bash
build_scan() {
  local NAME=$1
  shift
  "$CXX" "${CXXFLAGS[@]}" -O3 -DNDEBUG "$@" \
    "$SRC/scan_augmented_carrier_endpoint668.cpp" -o "$TMP/$NAME"
}
build_scan scan668 -DTARGET_ENDPOINT=668
build_scan scan667 -DTARGET_ENDPOINT=667
build_scan scan666 -DTARGET_ENDPOINT=666 \
  -DINCLUDE_ENDPOINT667_MINIMUM=1
build_scan scan665 -DTARGET_ENDPOINT=665 \
  -DINCLUDE_ENDPOINT667_MINIMUM=1
build_scan scan664_hostile -DTARGET_ENDPOINT=664 \
  -DINCLUDE_ENDPOINT667_MINIMUM=1 -DINCLUDE_ENDPOINT665_MINIMUM=1
build_scan scan664_final -DTARGET_ENDPOINT=664 \
  -DINCLUDE_ENDPOINT667_MINIMUM=1 -DINCLUDE_ENDPOINT665_MINIMUM=1 \
  -DINCLUDE_ENDPOINT664_MINIMUM=1
for NAME in scan668 scan667 scan666 scan665 scan664_hostile scan664_final; do
  "$TMP/$NAME" "${SCAN_ARGS[@]}" \
    > "$TMP/$NAME.out" 2> "$TMP/$NAME.err"
done
cmp "$TMP/scan668.out" \
  "$RES/scan_augmented_carrier_endpoint668_joint421_O3_NDEBUG.out"
cmp "$TMP/scan667.out" \
  "$RES/scan_augmented_carrier_endpoint667_joint421_O3.out"
cmp "$TMP/scan666.out" \
  "$RES/scan_augmented_carrier_endpoint666_joint421_O3.out"
cmp "$TMP/scan665.out" \
  "$RES/scan_augmented_carrier_endpoint665_joint421_O3.out"
cmp "$TMP/scan664_hostile.out" \
  "$RES/scan_augmented_carrier_endpoint664_joint421_O3.out"
cmp "$TMP/scan664_final.out" \
  "$RES/scan_augmented_carrier_endpoint664_final_O3.out"
```

For the final endpoint-664 boundary, also compile `scan664_final` with
O3+NDEBUG and O1+UBSan.  Add `-DQUIET_ENDPOINT_PROGRESS=1` to the UBSan build
so stderr is reserved for sanitizer diagnostics.  Both semantic stdout files
must match the stored O3 representative; UBSan stderr must match
`scan_augmented_carrier_endpoint664_final_UBSan.err` and be empty.

## 5. Detached literal endpoint audit

```bash
LITERAL_ARGS=(
  "$ROOT/05-knowledge/results/lrc14_endpoint_cascade_thm4254/replay_band"
  "$ROOT/05-knowledge/results/lrc14_three_round_learned_carrier_thm4266/full_discovery_416_704_O3.semantic.out"
  "$ROOT/05-knowledge/results/lrc14_three_round_learned_carrier_thm4266/full_discovery_416_700_O3.semantic.out"
  "$ROOT/05-knowledge/results/lrc14_three_round_learned_carrier_thm4266/full_discovery_520_700_O3.semantic.out"
  "$ROOT/05-knowledge/results/lrc14_three_round_learned_carrier_thm4266/full_discovery_384_694_O3.semantic.out"
  "$ROOT/05-knowledge/results/lrc14_fourth_round_learned_carrier_thm4271/full_discovery_520_688_O3.semantic.out"
  "$DECK"
)
for MODE in O0 O3 O3N UBSAN; do
  case "$MODE" in
    O0) FLAGS=(-O0) ;;
    O3) FLAGS=(-O3) ;;
    O3N) FLAGS=(-O3 -DNDEBUG) ;;
    UBSAN) FLAGS=(-O1 -g -fsanitize=undefined \
                  -fno-sanitize-recover=all) ;;
  esac
  "$CXX" "${CXXFLAGS[@]}" "${FLAGS[@]}" \
    "$SRC/verify_joint421_literal_r670.cpp" \
    -o "$TMP/endpoint_literal_$MODE"
  "$TMP/endpoint_literal_$MODE" "${LITERAL_ARGS[@]}" \
    > "$TMP/endpoint_literal_$MODE.out" \
    2> "$TMP/endpoint_literal_$MODE.err"
done
cmp "$TMP/endpoint_literal_O0.out" "$TMP/endpoint_literal_O3.out"
cmp "$TMP/endpoint_literal_O3.out" "$TMP/endpoint_literal_O3N.out"
cmp "$TMP/endpoint_literal_O3.out" "$TMP/endpoint_literal_UBSAN.out"
cmp "$TMP/endpoint_literal_O3N.out" \
  "$RES/verify_joint421_literal_r670_O3_NDEBUG.out"
test ! -s "$TMP/endpoint_literal_UBSAN.err"
cmp "$TMP/endpoint_literal_UBSAN.err" \
  "$RES/verify_joint421_literal_r670_UBSan.err"
```

### 5a. Complete body-response census and endpoint role of `003c900c`

The response-fibre census is standalone:

```bash
for MODE in O0 O3 UBSAN; do
  case "$MODE" in
    O0) FLAGS=(-O0) ;;
    O3) FLAGS=(-O3 -DNDEBUG) ;;
    UBSAN) FLAGS=(-O1 -g -fsanitize=undefined \
                  -fno-sanitize-recover=all) ;;
  esac
  "$CXX" -std=c++20 -pthread "${FLAGS[@]}" \
    "$SRC/joint421_response_fibre_census.cpp" \
    -o "$TMP/fibre_$MODE"
  "$TMP/fibre_$MODE" "$DECK" > "$TMP/fibre_$MODE.out" \
    2> "$TMP/fibre_$MODE.err"
done
cmp "$TMP/fibre_O0.out" "$TMP/fibre_O3.out"
cmp "$TMP/fibre_O3.out" "$TMP/fibre_UBSAN.out"
cmp "$TMP/fibre_O3.out" "$RES/joint421_response_fibre_census_O3.out"
test ! -s "$TMP/fibre_UBSAN.err"
cmp "$TMP/fibre_UBSAN.err" \
  "$RES/joint421_response_fibre_census_UBSan.err"
```

The portable endpoint-role checker imports the endpoint literal source under
its library-only guard:

```bash
for MODE in O0 O3N UBSAN; do
  case "$MODE" in
    O0) FLAGS=(-O0) ;;
    O3N) FLAGS=(-O3 -DNDEBUG) ;;
    UBSAN) FLAGS=(-O1 -g -fsanitize=undefined \
                  -fno-sanitize-recover=all) ;;
  esac
  "$CXX" "${CXXFLAGS[@]}" "${FLAGS[@]}" -I"$SRC" \
    "$SRC/joint421_redundant_role.cpp" -o "$TMP/role_$MODE"
  "$TMP/role_$MODE" "${LITERAL_ARGS[@]}" \
    > "$TMP/role_$MODE.out" 2> "$TMP/role_$MODE.err"
done
cmp "$TMP/role_O0.out" "$TMP/role_O3N.out"
cmp "$TMP/role_O3N.out" "$TMP/role_UBSAN.out"
cmp "$TMP/role_O3N.out" "$RES/joint421_redundant_role_O3_NDEBUG.out"
test ! -s "$TMP/role_UBSAN.err"
cmp "$TMP/role_UBSAN.err" "$RES/joint421_redundant_role_UBSan.err"
```

## 6. Proof graph and independent set replay

```bash
python3 -B "$SRC/proof_graph_consequence.py" \
  --post4271 "$POST4271" \
  --common "$RES/post_thm4277_joint421_common_primary.csv" \
  --emit-final "$TMP/post_thm4281_residual.csv" \
  > "$TMP/proof_graph_consequence.out"
cmp "$TMP/proof_graph_consequence.out" "$RES/proof_graph_consequence.out"
cmp "$TMP/post_thm4281_residual.csv" "$RES/post_thm4281_residual.csv"

OVERLAP_ARGS=(
  --post4271 "$POST4271"
  --common "$RES/post_thm4277_joint421_common_primary.csv"
  --proof-graph "$RES/proof_graph_consequence.out"
  --literal670 "$RES/verify_joint421_literal_r670_O3_NDEBUG.out"
  --scan669 "$RES/scan_augmented_carrier_endpoint669_joint421_O3.out"
  --scan668 "$RES/scan_augmented_carrier_endpoint668_joint421_O3_NDEBUG.out"
  --scan667 "$RES/scan_augmented_carrier_endpoint667_joint421_O3.out"
  --scan666 "$RES/scan_augmented_carrier_endpoint666_joint421_O3.out"
  --scan665 "$RES/scan_augmented_carrier_endpoint665_joint421_O3.out"
  --scan664 "$RES/scan_augmented_carrier_endpoint664_joint421_O3.out"
  --scan664-final "$RES/scan_augmented_carrier_endpoint664_final_O3.out"
  --emit-complement "$TMP/complement.csv"
  --emit-complement-layers "$TMP/complement_layers.csv"
  --emit-overlap-layers "$TMP/overlap_layers.csv"
  --emit-carrier-not-common "$TMP/carrier_not_common.csv"
  --emit-final "$TMP/final_overlap_replay.csv"
)
python3 -B "$SRC/verify_common_complement_carrier_overlap.py" \
  "${OVERLAP_ARGS[@]}" > "$TMP/overlap.out"
python3 -O -B "$SRC/verify_common_complement_carrier_overlap.py" \
  "${OVERLAP_ARGS[@]}" > "$TMP/overlap_opt.out"
cmp "$TMP/overlap.out" "$TMP/overlap_opt.out"
cmp "$TMP/overlap.out" "$RES/verify_common_complement_carrier_overlap.out"
cmp "$TMP/complement.csv" \
  "$RES/post_thm4277_joint421_not_common_primary.csv"
cmp "$TMP/complement_layers.csv" \
  "$RES/post_thm4277_joint421_not_common_primary_layers.csv"
cmp "$TMP/overlap_layers.csv" \
  "$RES/joint421_common_carrier_overlap_layers.csv"
cmp "$TMP/carrier_not_common.csv" "$RES/thm4281_carrier_not_common.csv"
cmp "$TMP/final_overlap_replay.csv" "$RES/post_thm4281_residual.csv"
```

Finally, from the repository root:

```bash
cd "$ROOT"
sha256sum -c "$RES/SHA256SUMS"
```
