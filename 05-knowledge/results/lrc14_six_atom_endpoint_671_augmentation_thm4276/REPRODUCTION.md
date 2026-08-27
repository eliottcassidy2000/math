# THM-4276 reproduction

Run from the repository root after placing the packet files at their displayed
canonical-relative paths. The commands depend on canonical THM-4266 and
THM-4271 artifacts but do not edit the repository.

```bash
SRC=04-computation/lrc14_six_atom_endpoint_671_augmentation_thm4276
OUT=05-knowledge/results/lrc14_six_atom_endpoint_671_augmentation_thm4276
OLD=05-knowledge/results/lrc14_three_round_learned_carrier_thm4266
R4=05-knowledge/results/lrc14_fourth_round_learned_carrier_thm4271

env LC_ALL=C LANG=C shasum -a 256 -c "$OUT/SHA256SUMS"
```

## 1. Recompute the two complete discovery controls

Each run enumerates all `binom(30,8)=5,852,925` repairs. Only wall-clock
timings are removed before byte comparison.

```bash
clang++ -std=c++20 -O3 -DNDEBUG \
  04-computation/lrc14_three_round_learned_carrier_thm4266/cascade_pair_exhaustive_primary.cpp \
  -o /tmp/round5-candidate-discovery

/tmp/round5-candidate-discovery 256 671 \
  > /tmp/round5-candidate-256-671.raw
sed -E 's/ BUILD_SECONDS [^ ]+ SCAN_SECONDS [^ ]+$//' \
  /tmp/round5-candidate-256-671.raw \
  > /tmp/round5-candidate-256-671.semantic
cmp /tmp/round5-candidate-256-671.semantic \
  "$OUT/full_discovery_256_671.semantic.out"

/tmp/round5-candidate-discovery 384 671 \
  > /tmp/round5-candidate-384-671.raw
sed -E 's/ BUILD_SECONDS [^ ]+ SCAN_SECONDS [^ ]+$//' \
  /tmp/round5-candidate-384-671.raw \
  > /tmp/round5-candidate-384-671.semantic
cmp /tmp/round5-candidate-384-671.semantic \
  "$OUT/full_discovery_384_671.semantic.out"
```

## 2. Recompute the full-universe response-pattern atlas

This pass independently enumerates all 5,852,925 repairs, re-evaluates
activation at both retained THM-4271 rows, and freezes the least repair for
each of the 330 distinct nonempty response patterns on the 27 failed-body
obligations.

```bash
clang++ -std=c++20 -O3 -DNDEBUG -I. \
  "$SRC/full_universe_augmentation_atlas.cpp" \
  -o /tmp/round5-candidate-atlas
/tmp/round5-candidate-atlas > /tmp/round5-candidate-atlas.out
cmp /tmp/round5-candidate-atlas.out \
  "$OUT/full_universe_augmentation_atlas.out"
```

## 3. Verify the exact minimum-six set-cover certificate

The Python checker does not invoke an optimizer. It verifies all 330 atlas
patterns, a six-obligation packing that every realizable pattern hits at most
once, and a six-pattern cover of all 27 obligations.

```bash
python3 -B "$SRC/minimum_set_cover_certificate.py" \
  "$OUT/full_universe_augmentation_atlas.out" \
  "$OUT/compact_carrier_descent_O3.out" \
  "$OUT/full_discovery_256_671.semantic.out" \
  "$OUT/full_discovery_384_671.semantic.out" \
  > /tmp/round5-candidate-minimum.out
cmp /tmp/round5-candidate-minimum.out \
  "$OUT/minimum_set_cover_certificate.out"
```

## 4. Rebuild the post-THM-4271 proof graph and compact descent

The semantic checker reconstructs every residual ledger through canonical
THM-4271 and emits the exact 174,904-edge post-THM-4271 residual. The C++
replay then reconstructs the 8,518-mask THM-4271 carrier, appends the six
certified masks in frozen order, checks both endpoint-671 seeds, and scans all
163 current endpoint-670 rows.

```bash
python3 -B "$SRC/proof_graph_consequence.py" --repo . \
  --compact-transcript "$OUT/compact_carrier_descent_O3.out" \
  --emit-post4271 /tmp/round5-candidate-post4271.csv \
  > /tmp/round5-candidate-proof-graph.out
cmp /tmp/round5-candidate-proof-graph.out \
  "$OUT/proof_graph_consequence.out"

clang++ -std=c++20 -O3 -DNDEBUG -I. \
  "$SRC/compact_carrier_descent.cpp" \
  -o /tmp/round5-candidate-descent
/tmp/round5-candidate-descent \
  05-knowledge/results/lrc14_endpoint_cascade_thm4254 \
  "$OLD/full_discovery_416_704_O3.semantic.out" \
  "$OLD/full_discovery_416_700_O3.semantic.out" \
  "$OLD/full_discovery_520_700_O3.semantic.out" \
  "$OLD/full_discovery_384_694_O3.semantic.out" \
  "$R4/full_discovery_520_688_O3.semantic.out" \
  "$OUT/full_discovery_256_671.semantic.out" \
  "$OUT/full_discovery_384_671.semantic.out" \
  /tmp/round5-candidate-post4271.csv 670 670 \
  > /tmp/round5-candidate-descent.out
cmp /tmp/round5-candidate-descent.out \
  "$OUT/compact_carrier_descent_O3.out"
```

## 5. Replay the detached literal-joint-wall controls

This checker rebuilds direct joint walls and never calls the endpoint-atom or
cocycle activation routines. It checks the inherited carrier construction,
the selected response patterns, both newly closed seeds, both complete
order-minimal discovery decks, and both endpoint-670 hostile controls.

```bash
COMMON_ARGS="05-knowledge/results/lrc14_endpoint_cascade_thm4254/replay_band \
$OLD/full_discovery_416_704_O3.semantic.out \
$OLD/full_discovery_416_700_O3.semantic.out \
$OLD/full_discovery_520_700_O3.semantic.out \
$OLD/full_discovery_384_694_O3.semantic.out \
$R4/full_discovery_520_688_O3.semantic.out \
$OUT/full_discovery_256_671.semantic.out \
$OUT/full_discovery_384_671.semantic.out"

clang++ -std=c++20 -O0 -I. -I04-computation \
  "$SRC/literal_independent_audit.cpp" \
  -o /tmp/round5-candidate-literal-O0
clang++ -std=c++20 -O3 -DNDEBUG -I. -I04-computation \
  "$SRC/literal_independent_audit.cpp" \
  -o /tmp/round5-candidate-literal-O3
clang++ -std=c++20 -O1 -g -fsanitize=address,undefined \
  -fno-omit-frame-pointer -I. -I04-computation \
  "$SRC/literal_independent_audit.cpp" \
  -o /tmp/round5-candidate-literal-san

/tmp/round5-candidate-literal-O0 $COMMON_ARGS \
  > /tmp/round5-candidate-literal-O0.out
/tmp/round5-candidate-literal-O3 $COMMON_ARGS \
  > /tmp/round5-candidate-literal-O3.out
/tmp/round5-candidate-literal-san $COMMON_ARGS \
  > /tmp/round5-candidate-literal-san.out \
  2> /tmp/round5-candidate-literal-san.err

cmp /tmp/round5-candidate-literal-O0.out "$OUT/literal_O0.out"
cmp /tmp/round5-candidate-literal-O3.out "$OUT/literal_O3.out"
cmp /tmp/round5-candidate-literal-san.out "$OUT/literal_san.out"
cmp /tmp/round5-candidate-literal-O0.out \
  /tmp/round5-candidate-literal-O3.out
cmp /tmp/round5-candidate-literal-O3.out \
  /tmp/round5-candidate-literal-san.out
test ! -s /tmp/round5-candidate-literal-san.err
cmp /tmp/round5-candidate-literal-san.err "$OUT/literal_san.err"
```
