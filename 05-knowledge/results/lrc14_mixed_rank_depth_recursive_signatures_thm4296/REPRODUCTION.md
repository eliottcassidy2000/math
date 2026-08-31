# THM-4296 reproduction

Run from the repository root. The promotion replay used a C++20 compiler,
Python 3, NumPy, and SciPy. All mathematical predicates are checked with exact
integers or exact rational fractions. SciPy is used only to suggest dual
weights; every rationalized dual load is then verified exactly.

```bash
packet=05-knowledge/results/lrc14_mixed_rank_depth_recursive_signatures_thm4296
src=04-computation/lrc14_mixed_rank_depth_recursive_signatures_thm4296
incoming=05-knowledge/results/lrc14_endpoint636_exchange_recursive_ideals_thm4295
scratch_dir=$(mktemp -d /tmp/thm4296-replay.XXXXXX)
mkdir -p "$scratch_dir"/{endpoint636,r632,r629,singletons,j19,proof}
cxx=${CXX:-g++}

# The packet's exact text-byte convention is LF. MinGW may translate native
# stdout/ofstream newlines to CRLF, so normalize scratch text before cmp.
normalize_lf() {
  python3 - "$@" <<'PY'
from pathlib import Path
import sys
for argument in sys.argv[1:]:
    path = Path(argument)
    payload = path.read_bytes()
    path.write_bytes(payload.replace(b"\r\n", b"\n").replace(b"\r", b"\n"))
PY
}
cmp_frozen() {
  normalize_lf "$1"
  cmp "$1" "$2"
}
```

`cmp_frozen` changes only generated files under `scratch_dir`; it never
rewrites a packet input or result. This makes the raw-byte comparisons
portable across POSIX and MinGW hosts without weakening any comparison.

## 1. Endpoint-636 response, lower bound, and exchange

Generate the complete response atlases:

```bash
$cxx -std=c++20 -O3 -DNDEBUG -pthread -I. \
  "$src/endpoint636_response_atlas.cpp" -o "$scratch_dir/endpoint636/atlas"
"$scratch_dir/endpoint636/atlas" \
  "$packet/inputs/endpoint636_failures101.csv" \
  "$packet/inputs/reconstructed_final8951.txt" \
  "$packet/inputs/additions45.txt" \
  "$packet/inputs/endpoint638_response_witness9.txt" \
  "$scratch_dir/endpoint636/atlas_files" 8 \
  > "$scratch_dir/endpoint636/atlas.out"
cmp_frozen "$scratch_dir/endpoint636/atlas.out" \
  "$packet/results/endpoint/endpoint636_response_atlas.out"
for name in LOCAL_100 LOCAL_256 COMMON_ACTIVE CARRIER_UNION; do
  cmp_frozen "$scratch_dir/endpoint636/atlas_files/$name.tsv" \
    "$packet/results/endpoint/endpoint636_atlas/$name.tsv"
done
```

Replay the exact cover and the independent integer dual-gap audit:

```bash
python3 "$src/exact_cover.py" \
  "$scratch_dir/endpoint636/atlas_files" \
  LOCAL_100 LOCAL_256 COMMON_ACTIVE CARRIER_UNION \
  > "$scratch_dir/endpoint636/exact_cover.out"
cmp_frozen "$scratch_dir/endpoint636/exact_cover.out" \
  "$packet/results/endpoint/endpoint636_exact_cover.out"

$cxx -std=c++20 -O3 -DNDEBUG \
  "$src/endpoint636_dual_gap_independent.cpp" \
  -o "$scratch_dir/endpoint636/dual"
"$scratch_dir/endpoint636/dual" \
  "$scratch_dir/endpoint636/atlas_files/CARRIER_UNION.tsv" \
  > "$scratch_dir/endpoint636/dual.out"
cmp_frozen "$scratch_dir/endpoint636/dual.out" \
  "$packet/results/endpoint/endpoint636_dual_gap_independent.out"
```

Replay the literal exchange on the complete nine-row boundary:

```bash
$cxx -std=c++20 -O3 -DNDEBUG -pthread -I. \
  "$src/detached_exchange_audit.cpp" \
  -o "$scratch_dir/endpoint636/exchange"
"$scratch_dir/endpoint636/exchange" \
  "$packet/inputs/reconstructed_final8951.txt" \
  "$packet/inputs/additions45.txt" \
  "$packet/inputs/endpoint638_response_witness9.txt" \
  "$packet/inputs/joint421_masks.txt" \
  "$packet/inputs/current_residual22682.csv" \
  "$packet/inputs/endpoint636_failures101.csv" 8 \
  > "$scratch_dir/endpoint636/exchange.out"
cmp_frozen "$scratch_dir/endpoint636/exchange.out" \
  "$packet/results/endpoint/detached_exchange_audit.out"
```

Repeat the dual program at `-O0`; its output must match the corresponding
file under `results/controls`. The detached exchange transcript is frozen as
the primary exhaustive literal-wall audit, without a maintained O0 twin.

## 2. Endpoint-632 hostile and mixed response

The detached hostile program exhausts the rank-eight complement of all 66
`(256,632)` failures, checks the unique hostile, its rank-nine extensions,
and the exact ten-mask cover of the other 65 bodies:

```bash
$cxx -std=c++20 -O3 -DNDEBUG -I. \
  "$src/r632_detached_hostile_survivor.cpp" \
  -o "$scratch_dir/r632/hostile"
"$scratch_dir/r632/hostile" \
  "$packet/inputs/post_exchange_r632_failures72.csv" \
  "$scratch_dir/r632/rank8_responses.csv" \
  > "$scratch_dir/r632/hostile.out"
cmp_frozen "$scratch_dir/r632/hostile.out" \
  "$packet/results/endpoint/r632_hostile_survivor.out"
cmp_frozen "$scratch_dir/r632/rank8_responses.csv" \
  "$packet/results/endpoint/r632_survivor_responses.csv"
```

Generate the complete rank-eight response atlases and certify their separate
local minima:

```bash
mkdir -p "$scratch_dir/r632/rank8_atlas" "$scratch_dir/r632/r632_atlas"
$cxx -std=c++20 -O3 -DNDEBUG -I. \
  "$src/r632_response_atlas.cpp" \
  -o "$scratch_dir/r632/rank8"
"$scratch_dir/r632/rank8" \
  "$packet/inputs/post_exchange_r632_failures72.csv" \
  "$packet/inputs/reconstructed_final8951.txt" \
  "$packet/inputs/additions45.txt" \
  "$packet/inputs/endpoint638_response_witness9.txt" \
  "$scratch_dir/r632/rank8_atlas" \
  > "$scratch_dir/r632/rank8.out"
cmp_frozen "$scratch_dir/r632/rank8.out" \
  "$packet/results/endpoint/r632_rank8_response_atlas.out"
for name in R632_LOCAL_100 R632_LOCAL_256; do
  cmp_frozen "$scratch_dir/r632/rank8_atlas/$name.tsv" \
    "$packet/results/endpoint/r632_rank8_atlas/$name.tsv"
done

python3 "$src/exact_cover.py" "$scratch_dir/r632/rank8_atlas" \
  R632_LOCAL_100 R632_LOCAL_256 \
  > "$scratch_dir/r632/rank8_exact_cover.out"
cmp_frozen "$scratch_dir/r632/rank8_exact_cover.out" \
  "$packet/results/endpoint/r632_rank8_exact_cover.out"
```

Now generate the complete rank-nine and mixed response quotient, using the
class atlas rather than the hostile program's per-incidence CSV, then solve
the mixed problem:

```bash
$cxx -std=c++20 -O3 -DNDEBUG -I"$src" \
  "$src/r632_rank9_mixed_response_atlas.cpp" \
  -o "$scratch_dir/r632/rank9"
"$scratch_dir/r632/rank9" \
  "$packet/inputs/post_exchange_r632_failures72.csv" \
  "$scratch_dir/r632/rank8_atlas/R632_LOCAL_256.tsv" \
  "$scratch_dir/r632/r632_atlas/rank9.tsv" \
  "$scratch_dir/r632/r632_atlas/mixed89.tsv" \
  > "$scratch_dir/r632/rank9.out"
cmp_frozen "$scratch_dir/r632/rank9.out" \
  "$packet/results/endpoint/r632_rank9_mixed_response_atlas.out"
cmp_frozen "$scratch_dir/r632/r632_atlas/rank9.tsv" \
  "$packet/results/endpoint/r632_atlas/rank9.tsv"
cmp_frozen "$scratch_dir/r632/r632_atlas/mixed89.tsv" \
  "$packet/results/endpoint/r632_atlas/mixed89.tsv"

python3 "$src/exact_cover.py" "$scratch_dir/r632" r632_mixed89_atlas \
  > "$scratch_dir/r632/exact_cover.out"
cmp_frozen "$scratch_dir/r632/exact_cover.out" \
  "$packet/results/endpoint/r632_mixed_exact_cover.out"

$cxx -std=c++20 -O3 -DNDEBUG -I. \
  "$src/r632_mixed_witness_detached.cpp" \
  -o "$scratch_dir/r632/witness"
"$scratch_dir/r632/witness" \
  "$packet/inputs/post_exchange_r632_failures72.csv" \
  "$scratch_dir/r632/witness_responses.csv" \
  > "$scratch_dir/r632/witness.out"
cmp_frozen "$scratch_dir/r632/witness.out" \
  "$packet/results/endpoint/r632_mixed_witness.out"
cmp_frozen "$scratch_dir/r632/witness_responses.csv" \
  "$packet/results/endpoint/r632_mixed_witness_responses.csv"
```

The endpoint-local and mixed exact-cover paths are also summarized in
`results/endpoint/EXACT-COVER-REPRODUCTION.md`.

## 3. Later exact response rounds and final carrier scan

The maintained generators for the `630`, `629`, and `628` boundaries are:

```bash
$cxx -std=c++20 -O3 -DNDEBUG -I. \
  "$src/r630_mixed_response_detached.cpp" \
  -o "$scratch_dir/r629/r630"
"$scratch_dir/r629/r630" > "$scratch_dir/r629/r630.out"
cmp_frozen "$scratch_dir/r629/r630.out" \
  "$packet/results/endpoint/r630_mixed_response.out"

$cxx -std=c++20 -O3 -DNDEBUG -I. \
  "$src/r629_mixed_response_atlas_detached.cpp" \
  -o "$scratch_dir/r629/atlas"
"$scratch_dir/r629/atlas" \
  "$packet/inputs/r629_failures28.csv" \
  "$scratch_dir/r629/mixed.tsv" \
  > "$scratch_dir/r629/atlas.out"
cmp_frozen "$scratch_dir/r629/atlas.out" \
  "$packet/results/endpoint/r629_mixed_response_atlas.out"
cmp_frozen "$scratch_dir/r629/mixed.tsv" \
  "$packet/results/endpoint/r629_mixed_atlas.tsv"
python3 "$src/r629_exact_cover.py" "$scratch_dir/r629/mixed.tsv" \
  > "$scratch_dir/r629/exact_cover.out"
cmp_frozen "$scratch_dir/r629/exact_cover.out" \
  "$packet/results/endpoint/r629_exact_cover.out"

$cxx -std=c++20 -O3 -DNDEBUG -I. \
  "$src/r629_mixed_optimum_detached_audit.cpp" \
  -o "$scratch_dir/r629/independent"
"$scratch_dir/r629/independent" \
  "$packet/inputs/r629_failures28.csv" \
  "$scratch_dir/r629/mixed.tsv" \
  > "$scratch_dir/r629/independent.out"
cmp_frozen "$scratch_dir/r629/independent.out" \
  "$packet/results/endpoint/r629_mixed_optimum_detached_audit.out"

$cxx -std=c++20 -O3 -DNDEBUG -I. \
  "$src/r628_mixed_response_detached.cpp" \
  -o "$scratch_dir/r629/r628"
"$scratch_dir/r629/r628" > "$scratch_dir/r629/r628.out"
cmp_frozen "$scratch_dir/r629/r628.out" \
  "$packet/results/endpoint/r628_mixed_response.out"
```

The final carrier scan takes the seven post-632 additions as explicit masks:

```bash
$cxx -std=c++20 -O3 -DNDEBUG -pthread -I. \
  "$src/mixed_carrier_descending_detached.cpp" \
  -o "$scratch_dir/endpoint636/descent"
"$scratch_dir/endpoint636/descent" \
  "$packet/inputs/joint421_masks.txt" \
  "$packet/inputs/reconstructed_final8951.txt" \
  "$packet/inputs/additions45.txt" \
  "$packet/inputs/endpoint638_response_witness9.txt" \
  "$packet/inputs/current_residual22647.csv" 620 \
  "$scratch_dir/endpoint636/final_failures.csv" \
  0010e125 002ac4c0 3882a082 0041c325 08c28e40 02008327 0006e281 \
  > "$scratch_dir/endpoint636/final_descent.out"
cmp_frozen "$scratch_dir/endpoint636/final_descent.out" \
  "$packet/results/endpoint/final_mixed9019_descent.out"
cmp_frozen "$scratch_dir/endpoint636/final_failures.csv" \
  "$packet/results/endpoint/final_mixed9019_failures.csv"
```

Repeat at `-O0`; both files must match `results/controls` byte-for-byte.

## 4. Singleton-signature census and detached audit

Generate the complete list of 192 live singleton groups, then run every group:

```bash
python3 "$src/singleton_targets.py" \
  --live "$packet/inputs/current_residual22647.csv" \
  --signatures "$packet/inputs/full_signatures_primary.csv" \
  --output "$scratch_dir/singletons/targets.out"

$cxx -std=c++20 -O3 -DNDEBUG -pthread -I. \
  "$src/singleton_live_surgery.cpp" \
  -o "$scratch_dir/singletons/surgery"
while read -r index count target; do
  "$scratch_dir/singletons/surgery" \
    "$packet/inputs/joint421_masks.txt" \
    "$packet/inputs/full_signatures_primary.csv" \
    "$packet/inputs/current_residual22647.csv" \
    "$index" "$target" \
    > "$scratch_dir/singletons/singleton_v2_${index}.out"
done < "$scratch_dir/singletons/targets.out"

python3 "$src/aggregate_singletons.py" \
  --targets "$scratch_dir/singletons/targets.out" \
  --runs-dir "$scratch_dir/singletons" \
  --live "$packet/inputs/current_residual22647.csv" \
  --output-dir "$scratch_dir/singletons/aggregate" \
  --summary-output "$scratch_dir/singletons/aggregate.out"
cmp_frozen "$scratch_dir/singletons/aggregate.out" \
  "$packet/results/singletons/aggregate.out"
cmp_frozen "$scratch_dir/singletons/aggregate/singleton_success_union.csv" \
  "$packet/results/singletons/success_union1219.csv"
cmp_frozen "$scratch_dir/singletons/aggregate/singleton_success_remaining.csv" \
  "$packet/results/singletons/post_singleton_residual21428.csv"
```

The aggregate also emits the primary discovery witness ledger. Feed that
intermediate to the detached audit; only after every activity and body check
passes does the program write the corrected literal ledger:

```bash
$cxx -std=c++20 -O3 -DNDEBUG -pthread -I. \
  "$src/singleton_success_literal_audit.cpp" \
  -o "$scratch_dir/singletons/literal"
"$scratch_dir/singletons/literal" \
  "$packet/inputs/joint421_masks.txt" \
  "$packet/inputs/full_signatures_primary.csv" \
  "$packet/inputs/current_residual22647.csv" \
  "$scratch_dir/singletons/aggregate/singleton_success_witnesses.csv" \
  --literal-witness-output "$scratch_dir/singletons/witnesses_literal.csv" \
  > "$scratch_dir/singletons/literal.out"
cmp_frozen "$scratch_dir/singletons/literal.out" \
  "$packet/results/singletons/detached_literal_audit.out"
cmp_frozen "$scratch_dir/singletons/witnesses_literal.csv" \
  "$packet/results/singletons/witnesses_literal.csv"
```

Repeat the literal audit at `-O0`; its output must match the maintained O0
control.

## 5. Index-19 common deck

```bash
$cxx -std=c++20 -O3 -DNDEBUG -pthread -I. \
  "$src/j19_full_common_response.cpp" -o "$scratch_dir/j19/primary"
"$scratch_dir/j19/primary" \
  "$packet/inputs/joint421_masks.txt" \
  "$packet/inputs/full_signatures_primary.csv" \
  "$packet/inputs/current_residual22647.csv" \
  > "$scratch_dir/j19/primary.out"
cmp_frozen "$scratch_dir/j19/primary.out" "$packet/results/j19/primary.out"

$cxx -std=c++20 -O3 -DNDEBUG -pthread -I. \
  "$src/detached_j19_audit.cpp" -o "$scratch_dir/j19/detached"
"$scratch_dir/j19/detached" \
  "$packet/inputs/joint421_masks.txt" \
  "$packet/inputs/full_signatures_primary.csv" \
  "$packet/inputs/current_residual22647.csv" 8 \
  > "$scratch_dir/j19/detached.out"
cmp_frozen "$scratch_dir/j19/detached.out" \
  "$packet/results/j19/detached_literal_audit.out"
```

The primary file's second weakest-row annotation is noncanonical under
MISTAKE-532. The detached file is authoritative for exact cross-row margin
fractions. All primary sign, response, and body-cover fields reproduce.

## 6. Typed proof graph

Extract the carrier node directly from the raw final scan and consume the
three typed node sets:

```bash
python3 "$src/endpoint_success_from_raw_scan.py" \
  "$packet/results/endpoint/final_mixed9019_descent.out" \
  "$scratch_dir/proof/endpoint70.csv" \
  --failure-csv "$scratch_dir/proof/failed2.csv" \
  > "$scratch_dir/proof/extractor.out"
cmp_frozen "$scratch_dir/proof/endpoint70.csv" \
  "$packet/results/proof_graph/endpoint_success70.csv"
cmp_frozen "$scratch_dir/proof/failed2.csv" \
  "$packet/results/proof_graph/endpoint_failed_pairs2.csv"
cmp_frozen "$scratch_dir/proof/extractor.out" \
  "$packet/results/proof_graph/endpoint_extractor.out"

python3 "$src/typed_proof_graph_consumer.py" \
  "$packet/inputs/current_residual22647.csv" \
  "$packet/results/singletons/success_union1219.csv" \
  "$packet/results/j19/ideal36.csv" \
  "$scratch_dir/proof/endpoint70.csv" \
  --union-out "$scratch_dir/proof/union1321.csv" \
  --residual-out "$scratch_dir/proof/final21326.csv" \
  > "$scratch_dir/proof/typed.out"
cmp_frozen "$scratch_dir/proof/union1321.csv" \
  "$packet/results/proof_graph/typed_union1321.csv"
cmp_frozen "$scratch_dir/proof/final21326.csv" \
  "$packet/results/proof_graph/final_residual21326.csv"
cmp_frozen "$scratch_dir/proof/typed.out" \
  "$packet/results/proof_graph/typed_consumer.out"
```

Finally compile the independent C++ consumer. It derives the index-19 ideal,
parses the raw scan, verifies the complete inherited endpoint prefix, checks
all overlaps, and requires byte equality with the maintained ledgers:

```bash
$cxx -std=c++20 -O3 -DNDEBUG \
  "$src/proof_graph_union_independent.cpp" \
  -o "$scratch_dir/proof/independent"
"$scratch_dir/proof/independent" \
  "$packet/inputs/current_residual22647.csv" \
  "$packet/results/singletons/success_union1219.csv" \
  "$packet/inputs/full_signatures_primary.csv" \
  "$packet/results/proof_graph/endpoint_success70.csv" \
  "$packet/results/j19/ideal36.csv" \
  "$packet/results/proof_graph/typed_union1321.csv" \
  "$packet/results/proof_graph/final_residual21326.csv" \
  "$packet/results/endpoint/final_mixed9019_descent.out" \
  "$scratch_dir/proof/endpoint_independent.csv" \
  "$scratch_dir/proof/union_independent.csv" \
  "$scratch_dir/proof/final_independent.csv" \
  > "$scratch_dir/proof/independent.out"
cmp_frozen "$scratch_dir/proof/independent.out" \
  "$packet/results/proof_graph/proof_graph_union_independent.out"
```

The preceding replay preserves this theorem's own 1,321-row `S/J/K` subgraph.
Now join only its typed consequences with THM-4295's independently certified
carrier and signature-ideal row sets:

```bash
python3 "$src/incoming_thm4295_typed_union.py" \
  --union-out "$scratch_dir/proof/union1324.csv" \
  --residual-out "$scratch_dir/proof/final21323.csv" \
  --new-out "$scratch_dir/proof/incoming_new3.csv" \
  > "$scratch_dir/proof/incoming_join.out"
cmp_frozen "$scratch_dir/proof/union1324.csv" \
  "$packet/results/proof_graph/typed_union1324.csv"
cmp_frozen "$scratch_dir/proof/final21323.csv" \
  "$packet/results/proof_graph/final_residual21323.csv"
cmp_frozen "$scratch_dir/proof/incoming_new3.csv" \
  "$packet/results/proof_graph/incoming_new3.csv"
cmp_frozen "$scratch_dir/proof/incoming_join.out" \
  "$packet/results/proof_graph/incoming_thm4295_typed_union.out"

$cxx -std=c++20 -O3 -DNDEBUG \
  "$src/incoming_thm4295_typed_union_independent.cpp" \
  -o "$scratch_dir/proof/incoming_independent"
"$scratch_dir/proof/incoming_independent" \
  "$packet/inputs/current_residual22647.csv" \
  "$packet/results/singletons/success_union1219.csv" \
  "$packet/results/j19/ideal36.csv" \
  "$packet/results/proof_graph/endpoint_success70.csv" \
  "$incoming/inputs/carrier_layers636_633.csv" \
  "$incoming/inputs/signature19_fibre36.csv" \
  "$incoming/inputs/signature294_fibre21.csv" \
  "$incoming/inputs/signature372_fibre54.csv" \
  "$scratch_dir/proof/independent1324.csv" \
  "$scratch_dir/proof/independent21323.csv" \
  "$scratch_dir/proof/independent_new3.csv" \
  > "$scratch_dir/proof/incoming_independent.out"
cmp_frozen "$scratch_dir/proof/independent1324.csv" \
  "$packet/results/proof_graph/typed_union1324.csv"
cmp_frozen "$scratch_dir/proof/independent21323.csv" \
  "$packet/results/proof_graph/final_residual21323.csv"
cmp_frozen "$scratch_dir/proof/independent_new3.csv" \
  "$packet/results/proof_graph/incoming_new3.csv"
cmp_frozen "$scratch_dir/proof/incoming_independent.out" \
  "$packet/results/proof_graph/incoming_thm4295_typed_union_independent.out"
```

The two consumers agree that the incoming ten-row carrier node is contained
in the 70-row node, incoming `H19` equals `J`, and only three rows are new.
These are row-set statements only; no carrier or common-deck masks are merged.

## 7. Optimization controls

For each row below, repeat the preceding `-O3 -DNDEBUG` compile with `-O0`,
keep the same runtime arguments, write to a distinct scratch file, and compare
that transcript with `cmp_frozen`. The frozen controls are byte-identical to
their optimized primaries after the declared LF normalization.

| source | frozen `-O0` transcript |
|---|---|
| `endpoint636_dual_gap_independent.cpp` | `results/controls/endpoint636_dual_gap_independent_O0.out` |
| `r632_detached_hostile_survivor.cpp` | `results/controls/r632_hostile_survivor_O0.out` |
| `r630_mixed_response_detached.cpp` | `results/controls/r630_mixed_response_O0.out` |
| `r628_mixed_response_detached.cpp` | `results/controls/r628_mixed_response_O0.out` |
| `mixed_carrier_descending_detached.cpp` | `results/controls/final_mixed9019_descent_O0.out` (and `final_mixed9019_failures_O0.csv`) |
| `singleton_success_literal_audit.cpp` | `results/controls/singleton_detached_literal_O0.out` |
| `detached_j19_audit.cpp` | `results/controls/detached_j19_O0.out` |
| `proof_graph_union_independent.cpp` | `results/controls/proof_graph_union_independent_O0.out` |

The detached exchange audit has no maintained O0 twin; its primary exhaustive
transcript is the frozen literal-wall certificate.

The auxiliary r629 witness, carrier-identity, and post-9,017 consumer reports
are independent provenance checks documented in
`results/endpoint/R629-DETACHED-AUDIT.md`. They are not inputs to the typed
proof-graph consumer; the maintained 9,019 scan supersedes their intermediate
raw boundary. Likewise, `results/endpoint/pre_r630_failure_boundary.csv` is
frozen discovery provenance for the two r630 bodies; the load-bearing replay
uses those literal bodies in `r630_mixed_response_detached.cpp` and then
checks them again in the final carrier scan.

## 8. Packet identity

After the replays, verify every frozen packet byte:

```bash
(cd "$incoming" && sha256sum -c SHA256SUMS)
(cd "$packet" && sha256sum -c SHA256SUMS)
```

No command in this reproduction establishes physical entry or LRC(14).
