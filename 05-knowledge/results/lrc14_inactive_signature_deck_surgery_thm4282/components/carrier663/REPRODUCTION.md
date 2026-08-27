# Reproduction

Run with Bash from a checkout containing the promoted THM-4254, THM-4266,
THM-4271, THM-4276, and THM-4281 artifacts.  The frozen audit used Apple clang
17.0.0; any conforming C++20 compiler should work.

```bash
set -euo pipefail
PACKET_ROOT=$(pwd)
REPO_ROOT=/private/tmp/math-wt-cross-frontier-bridge
TMP_ROOT=$(mktemp -d)
CXX=${CXX:-clang++}

cd "$REPO_ROOT"
sha256sum -c "$PACKET_ROOT/controls/canonical_dependencies.sha256"
sha256sum -c "$PACKET_ROOT/controls/replay_band_files.sha256"

ARGS=(
  "$REPO_ROOT/05-knowledge/results/lrc14_endpoint_cascade_thm4254/replay_band"
  "$REPO_ROOT/05-knowledge/results/lrc14_three_round_learned_carrier_thm4266/full_discovery_416_704_O3.semantic.out"
  "$REPO_ROOT/05-knowledge/results/lrc14_three_round_learned_carrier_thm4266/full_discovery_416_700_O3.semantic.out"
  "$REPO_ROOT/05-knowledge/results/lrc14_three_round_learned_carrier_thm4266/full_discovery_520_700_O3.semantic.out"
  "$REPO_ROOT/05-knowledge/results/lrc14_three_round_learned_carrier_thm4266/full_discovery_384_694_O3.semantic.out"
  "$REPO_ROOT/05-knowledge/results/lrc14_fourth_round_learned_carrier_thm4271/full_discovery_520_688_O3.semantic.out"
  "$REPO_ROOT/05-knowledge/results/lrc14_joint421_global_common_carrier_thm4281/joint_iter39_shrunk_forward_masks.txt"
)

COMMON=(
  -std=c++20 -pthread
  -I"$REPO_ROOT"
  -I"$REPO_ROOT/04-computation"
  "$PACKET_ROOT/src/reconstruct_and_literal_audit.cpp"
)

"$CXX" "${COMMON[@]}" -O0 -o "$TMP_ROOT/audit_O0"
"$CXX" "${COMMON[@]}" -O3 -DNDEBUG -o "$TMP_ROOT/audit_O3"
"$CXX" "${COMMON[@]}" -O1 -g -fsanitize=undefined \
  -fno-sanitize-recover=undefined -o "$TMP_ROOT/audit_UBSan"

for MODE in O0 O3 UBSan; do
  "$TMP_ROOT/audit_$MODE" "${ARGS[@]}" \
    "$TMP_ROOT/carrier_$MODE.txt" \
    "$TMP_ROOT/vulnerable_$MODE.csv" \
    > "$TMP_ROOT/audit_$MODE.out" \
    2> "$TMP_ROOT/audit_$MODE.err"
done

cmp "$TMP_ROOT/audit_O0.out" "$TMP_ROOT/audit_O3.out"
cmp "$TMP_ROOT/audit_O3.out" "$TMP_ROOT/audit_UBSan.out"
cmp "$TMP_ROOT/carrier_O0.txt" "$TMP_ROOT/carrier_O3.txt"
cmp "$TMP_ROOT/carrier_O3.txt" "$TMP_ROOT/carrier_UBSan.txt"
cmp "$TMP_ROOT/vulnerable_O0.csv" "$TMP_ROOT/vulnerable_O3.csv"
cmp "$TMP_ROOT/vulnerable_O3.csv" "$TMP_ROOT/vulnerable_UBSan.csv"
test ! -s "$TMP_ROOT/audit_O0.err"
test ! -s "$TMP_ROOT/audit_O3.err"
test ! -s "$TMP_ROOT/audit_UBSan.err"

for MODE in O0 O3 UBSan; do
  cmp "$TMP_ROOT/audit_$MODE.out" \
    "$PACKET_ROOT/results/audit_$MODE.out"
  cmp "$TMP_ROOT/carrier_$MODE.txt" \
    "$PACKET_ROOT/results/reconstructed_final8951_carrier_$MODE.txt"
  cmp "$TMP_ROOT/vulnerable_$MODE.csv" \
    "$PACKET_ROOT/results/reconstructed_vulnerable_bodies_$MODE.csv"
  cmp "$TMP_ROOT/audit_$MODE.err" \
    "$PACKET_ROOT/results/audit_$MODE.err"
done

cd "$PACKET_ROOT"
sha256sum -c SHA256SUMS
```

On systems without `sha256sum`, substitute `shasum -a 256` and its compatible
check mode.  The scratch checkout path is not semantically significant; set
`REPO_ROOT` to the repository being audited.

