# THM-4271 reproduction

Run from the repository root after placing the packet files at their
canonical paths.

```bash
SRC=04-computation/lrc14_fourth_round_learned_carrier_thm4271
OUT=05-knowledge/results/lrc14_fourth_round_learned_carrier_thm4271
OLD=05-knowledge/results/lrc14_three_round_learned_carrier_thm4266

env LC_ALL=C LANG=C shasum -a 256 -c "$OUT/SHA256SUMS"
```

Recompute the complete `binom(30,8)=5,852,925` discovery at `(520,688)`.
Only wall-clock timings are removed before comparison.

```bash
clang++ -std=c++20 -O3 -DNDEBUG \
  04-computation/lrc14_three_round_learned_carrier_thm4266/cascade_pair_exhaustive_primary.cpp \
  -o /tmp/thm4271-discovery
/tmp/thm4271-discovery 520 688 > /tmp/thm4271-discovery.raw
sed -E 's/ BUILD_SECONDS [^ ]+ SCAN_SECONDS [^ ]+$//' \
  /tmp/thm4271-discovery.raw > /tmp/thm4271-discovery.semantic
cmp /tmp/thm4271-discovery.semantic \
  "$OUT/full_discovery_520_688_O3.semantic.out"
```

Reconstruct the exact post-THM-4267 residual from the canonical THM-4266
residual, audit the frozen layer transcript, and replay the live proof graph
through THM-4270.

```bash
python3 -B "$SRC/proof_graph_consequence.py" --repo . \
  --emit-post4267 /tmp/thm4271-post4267.csv \
  > /tmp/thm4271-proof-graph.out
cmp /tmp/thm4271-proof-graph.out "$OUT/proof_graph_consequence.out"
```

Rebuild the 8,319-mask THM-4266 carrier, append only the 199 novel masks,
and rescan every current residual row from endpoint 688 down to the first
resistant layer 671.

```bash
clang++ -std=c++20 -O3 -DNDEBUG -I. \
  "$SRC/round4_carrier_descent.cpp" -o /tmp/thm4271-descent
/tmp/thm4271-descent \
  05-knowledge/results/lrc14_endpoint_cascade_thm4254 \
  "$OLD/full_discovery_416_704_O3.semantic.out" \
  "$OLD/full_discovery_416_700_O3.semantic.out" \
  "$OLD/full_discovery_520_700_O3.semantic.out" \
  "$OLD/full_discovery_384_694_O3.semantic.out" \
  "$OUT/full_discovery_520_688_O3.semantic.out" \
  /tmp/thm4271-post4267.csv 688 671 \
  > /tmp/thm4271-descent.out
cmp /tmp/thm4271-descent.out "$OUT/round4_descent_688_671_O3.out"
```

Finally replay the detached literal-joint-wall audit. It never calls the
endpoint-atom or cocycle activation routines. Ordinary, optimized, and
ASan+UBSan transcripts must byte-match, and sanitizer stderr must be empty.

```bash
COMMON_ARGS="05-knowledge/results/lrc14_endpoint_cascade_thm4254/replay_band \
$OLD/full_discovery_416_704_O3.semantic.out \
$OLD/full_discovery_416_700_O3.semantic.out \
$OLD/full_discovery_520_700_O3.semantic.out \
$OLD/full_discovery_384_694_O3.semantic.out \
$OUT/full_discovery_520_688_O3.semantic.out"

clang++ -std=c++20 -O0 -I. -I04-computation \
  "$SRC/literal_independent_audit.cpp" -o /tmp/thm4271-literal-O0
clang++ -std=c++20 -O3 -DNDEBUG -I. -I04-computation \
  "$SRC/literal_independent_audit.cpp" -o /tmp/thm4271-literal-O3
clang++ -std=c++20 -O1 -g -fsanitize=address,undefined \
  -fno-omit-frame-pointer -I. -I04-computation \
  "$SRC/literal_independent_audit.cpp" -o /tmp/thm4271-literal-san

/tmp/thm4271-literal-O0 $COMMON_ARGS > /tmp/thm4271-literal-O0.out
/tmp/thm4271-literal-O3 $COMMON_ARGS > /tmp/thm4271-literal-O3.out
/tmp/thm4271-literal-san $COMMON_ARGS \
  > /tmp/thm4271-literal-san.out 2> /tmp/thm4271-literal-san.err

cmp /tmp/thm4271-literal-O0.out "$OUT/literal_O0.out"
cmp /tmp/thm4271-literal-O3.out "$OUT/literal_O3.out"
cmp /tmp/thm4271-literal-san.out "$OUT/literal_san.out"
cmp /tmp/thm4271-literal-O0.out /tmp/thm4271-literal-O3.out
cmp /tmp/thm4271-literal-O3.out /tmp/thm4271-literal-san.out
test ! -s /tmp/thm4271-literal-san.err
cmp /tmp/thm4271-literal-san.err "$OUT/literal_san.err"
```
