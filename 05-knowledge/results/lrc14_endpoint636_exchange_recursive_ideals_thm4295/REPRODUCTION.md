# THM-4295 reproduction

Run from the repository root. The promoted C++ outputs use exact integer
arithmetic. The two solver checks consume frozen exact binary response
atlases; their objective/status lines are load-bearing, while an optimal
witness and timing statistics may vary by solver version or thread schedule.

```bash
scratch_dir=$(mktemp -d /tmp/thm4295-replay.XXXXXX)
packet=05-knowledge/results/lrc14_endpoint636_exchange_recursive_ideals_thm4295
prior=05-knowledge/results/lrc14_repaired_carrier_endpoint637_descent_thm4287
src=04-computation/lrc14_endpoint636_exchange_recursive_ideals_thm4295
mkdir -p "$scratch_dir/carrier" "$scratch_dir/signatures" \
  "$scratch_dir/proof/primary" "$scratch_dir/proof/independent"
```

## 1. Inherited carrier boundary and 101-response minimum

```bash
clang++ -std=c++20 -O3 -DNDEBUG -pthread -I. \
  "$src/base9006_boundary_scan.cpp" -o "$scratch_dir/base9006_scan"
"$scratch_dir/base9006_scan" \
  "$prior/inputs/joint421_masks.txt" \
  "$prior/inputs/reconstructed_final8951.txt" \
  "$prior/inputs/additions45.txt" \
  "$prior/inputs/endpoint638_response_witness9.txt" \
  "$prior/results/proof_graph/final_residual.csv" 632 \
  "$scratch_dir/carrier/base9006_failures.csv" \
  > "$scratch_dir/carrier/base9006_boundary_scan.out"
cmp "$scratch_dir/carrier/base9006_boundary_scan.out" \
  "$packet/results/carrier/base9006_boundary_scan.out"
cmp "$scratch_dir/carrier/base9006_failures.csv" \
  "$packet/inputs/base9006_failures636_632.csv"

clang++ -std=c++20 -O3 -DNDEBUG -pthread -I. \
  "$src/failure_response_exact_cover.cpp" -o "$scratch_dir/failure_response"
"$scratch_dir/failure_response" \
  "$prior/results/primary/failures.csv" \
  "$scratch_dir/carrier/maximal_response.tsv" \
  > "$scratch_dir/carrier/failure_response.out"
cmp "$scratch_dir/carrier/failure_response.out" \
  "$packet/results/carrier/failure_response.out"
cmp "$scratch_dir/carrier/maximal_response.tsv" \
  "$packet/results/carrier/maximal_response.tsv"

python3 "$src/setcover_highs_independent.py" \
  "$packet/results/carrier/maximal_response.tsv" \
  > "$scratch_dir/carrier/setcover_highs.out"
rg -q '^STATUS 0 SUCCESS True OBJECTIVE 14$' \
  "$scratch_dir/carrier/setcover_highs.out"
rg -q '^MIP_GAP 0([.]0)? ' "$scratch_dir/carrier/setcover_highs.out"
```

The C++ depth-first search prints `EXACT_MINIMUM 14`; the independent HiGHS
transcript prints objective `14` and zero MIP gap.

## 2. Selected 14-mask carrier and descending boundary

```bash
clang++ -std=c++20 -O3 -DNDEBUG -pthread -I. \
  "$src/append14_descending_scan.cpp" -o "$scratch_dir/append14_scan"
"$scratch_dir/append14_scan" \
  "$prior/inputs/joint421_masks.txt" \
  "$prior/inputs/reconstructed_final8951.txt" \
  "$prior/inputs/additions45.txt" \
  "$prior/inputs/endpoint638_response_witness9.txt" \
  "$prior/results/proof_graph/final_residual.csv" \
  "$packet/inputs/append14_masks.txt" 632 \
  "$scratch_dir/carrier/endpoint632_failures.csv" \
  > "$scratch_dir/carrier/descending_scan.out"
cmp "$scratch_dir/carrier/descending_scan.out" \
  "$packet/results/carrier/descending_scan.out"
cmp "$scratch_dir/carrier/endpoint632_failures.csv" \
  "$packet/results/carrier/endpoint632_failures.csv"

clang++ -std=c++20 -O3 -DNDEBUG -pthread -I. \
  "$src/append14_detached_literal_audit.cpp" \
  -o "$scratch_dir/append14_literal"
"$scratch_dir/append14_literal" \
  "$prior/inputs/reconstructed_final8951.txt" \
  "$prior/inputs/additions45.txt" \
  "$prior/inputs/endpoint638_response_witness9.txt" \
  "$packet/inputs/append14_masks.txt" \
  "$prior/results/primary/failures.csv" \
  > "$scratch_dir/carrier/append14_detached_literal.out"
cmp "$scratch_dir/carrier/append14_detached_literal.out" \
  "$packet/results/carrier/append14_detached_literal.out"
```

Repeat the detached source with `-O0` and with
`-O1 -g -fsanitize=undefined -fno-sanitize-recover=undefined`. The stdout
files must match the two `controls/append14_*` files and stderr must be empty.

## 3. Four-row response quotient and downstream optimum

```bash
clang++ -std=c++20 -O3 -DNDEBUG -pthread -I. \
  "$src/rowaware176_response.cpp" -o "$scratch_dir/response176"
"$scratch_dir/response176" \
  "$packet/inputs/base9006_failures636_632.csv" \
  "$scratch_dir/carrier/rowaware176_classes.tsv" \
  > "$scratch_dir/carrier/rowaware176_response.out"
cmp "$scratch_dir/carrier/rowaware176_response.out" \
  "$packet/results/carrier/rowaware176_response.out"
cmp "$scratch_dir/carrier/rowaware176_classes.tsv" \
  "$packet/results/carrier/rowaware176_classes.tsv"

python3 "$src/rowaware176_optimize.py" \
  "$packet/results/carrier/rowaware176_classes.tsv" \
  > "$scratch_dir/carrier/rowaware176_optimize.out"
python3 "$src/rowaware176_highs_independent.py" \
  "$packet/results/carrier/rowaware176_classes.tsv" \
  > "$scratch_dir/carrier/rowaware176_highs.out"
```

Both solver transcripts must report an optimal objective/miss count of `37`,
cardinality `14`, full coverage of the first 101 bits, and zero gap. Do not
byte-compare the CP-SAT transcript: its optimal witness, timing, and search
statistics may vary. The stored HiGHS transcript records the promoted
environment and an independently produced optimum.

## 4. Absolute endpoint-632 obstruction

```bash
clang++ -std=c++20 -O3 -DNDEBUG -pthread -I. \
  "$src/endpoint632_zero_response.cpp" -o "$scratch_dir/zero_response"
"$scratch_dir/zero_response" \
  > "$scratch_dir/carrier/endpoint632_zero_response.out"
cmp "$scratch_dir/carrier/endpoint632_zero_response.out" \
  "$packet/results/carrier/endpoint632_zero_response.out"

clang++ -std=c++20 -O3 -DNDEBUG -pthread -I. \
  "$src/endpoint632_unresponsive_literal_audit.cpp" \
  -o "$scratch_dir/unresponsive_literal"
"$scratch_dir/unresponsive_literal" \
  > "$scratch_dir/carrier/endpoint632_unresponsive_literal.out"
cmp "$scratch_dir/carrier/endpoint632_unresponsive_literal.out" \
  "$packet/results/carrier/endpoint632_unresponsive_literal.out"
```

Repeat the detached source at O0 and under UBSan. The stdout files must match
`controls/unresponsive_O0.out` and `controls/unresponsive_ubsan.out`; stderr
must be empty.

## 5. Three recursive signature-ideal surgeries

The deletion selectors are the target signatures of the three top rows.

```bash
for target in 338 294 372; do
  clang++ -std=c++20 -O3 -DNDEBUG -pthread -I. \
    "$src/signature_ideal_surgery.cpp" -o "$scratch_dir/signature_surgery"
  "$scratch_dir/signature_surgery" \
    "$prior/inputs/joint421_masks.txt" \
    "$prior/inputs/full_signatures_primary.csv" \
    "$prior/results/proof_graph/final_residual.csv" \
    "@$target,636" \
    "$scratch_dir/signatures/signature${target}_fibre.csv" \
    "$scratch_dir/signatures/signature${target}_witness.txt" \
    > "$scratch_dir/signatures/signature${target}_surgery.out"
done

cmp "$scratch_dir/signatures/signature338_fibre.csv" \
  "$packet/inputs/signature19_fibre36.csv"
cmp "$scratch_dir/signatures/signature338_witness.txt" \
  "$packet/inputs/signature19_witness2.txt"
cmp "$scratch_dir/signatures/signature338_surgery.out" \
  "$packet/results/signature19/surgery.out"
cmp "$scratch_dir/signatures/signature294_fibre.csv" \
  "$packet/inputs/signature294_fibre21.csv"
cmp "$scratch_dir/signatures/signature294_witness.txt" \
  "$packet/inputs/signature294_witness4.txt"
cmp "$scratch_dir/signatures/signature294_surgery.out" \
  "$packet/results/signature294/surgery.out"
cmp "$scratch_dir/signatures/signature372_fibre.csv" \
  "$packet/inputs/signature372_fibre54.csv"
cmp "$scratch_dir/signatures/signature372_witness.txt" \
  "$packet/inputs/signature372_witness4.txt"
cmp "$scratch_dir/signatures/signature372_surgery.out" \
  "$packet/results/signature372/surgery.out"
```

Run the three detached literal audits:

```bash
for tag in 19 294 372; do
  clang++ -std=c++20 -O3 -DNDEBUG -pthread \
    "$src/signature${tag}_detached_literal_audit.cpp" \
    -o "$scratch_dir/signature${tag}_literal"
  "$scratch_dir/signature${tag}_literal" \
    "$prior/inputs/joint421_masks.txt" \
    "$packet/inputs/signature${tag}_fibre"*.csv \
    "$packet/inputs/signature${tag}_witness"*.txt \
    > "$scratch_dir/signatures/signature${tag}_detached.out"
  cmp "$scratch_dir/signatures/signature${tag}_detached.out" \
    "$packet/results/signature${tag}/detached_literal.out"
done
```

Repeat each detached source at O0 and under UBSan. The six stdout files must
match `controls/sig*_O0.out` and `controls/sig*_ubsan.out`; every `.err` file
must be empty.

## 6. Independent typed proof-graph consumers

```bash
python3 "$src/proof_graph_consequence.py" \
  "$prior/results/proof_graph/final_residual.csv" \
  "$packet/inputs/carrier_layers636_633.csv" \
  "$packet/inputs/signature19_fibre36.csv" \
  "$packet/inputs/signature294_fibre21.csv" \
  "$packet/inputs/signature372_fibre54.csv" \
  "$scratch_dir/proof/primary" \
  > "$scratch_dir/proof/primary.out"
env LC_ALL=C LANG=C ruby "$src/proof_graph_consequence_independent.rb" \
  "$prior/results/proof_graph/final_residual.csv" \
  "$packet/inputs/carrier_layers636_633.csv" \
  "$packet/inputs/signature19_fibre36.csv" \
  "$packet/inputs/signature294_fibre21.csv" \
  "$packet/inputs/signature372_fibre54.csv" \
  "$scratch_dir/proof/independent" \
  > "$scratch_dir/proof/independent.out"
cmp "$scratch_dir/proof/primary.out" "$scratch_dir/proof/independent.out"
cmp "$scratch_dir/proof/primary.out" \
  "$packet/results/proof_graph/primary.out"
cmp "$scratch_dir/proof/primary/carrier10_three_ideal_union.csv" \
  "$packet/results/proof_graph/primary/carrier10_three_ideal_union.csv"
cmp "$scratch_dir/proof/primary/post_carrier10_three_ideal_residual.csv" \
  "$packet/results/proof_graph/primary/post_carrier10_three_ideal_residual.csv"
cmp "$scratch_dir/proof/primary/carrier10_three_ideal_union.csv" \
  "$scratch_dir/proof/independent/carrier10_three_ideal_union.csv"
cmp "$scratch_dir/proof/primary/post_carrier10_three_ideal_residual.csv" \
  "$scratch_dir/proof/independent/post_carrier10_three_ideal_residual.csv"
```

Finally verify every frozen packet file:

```bash
(cd "$packet" && env LC_ALL=C LANG=C sha256sum -c SHA256SUMS)
```
