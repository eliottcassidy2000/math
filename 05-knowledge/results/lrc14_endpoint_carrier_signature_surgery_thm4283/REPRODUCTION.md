# THM-4283 reproduction

Run from the repository root. The recorded promotion used Apple clang 17.0.0
on arm64 macOS, Python 3.10.0, and Ruby 2.6.10. All programs are exact integer
computations; no solver, floating point, sampling, or randomized tie-break is
used.

Set a scratch directory outside the repository:

```bash
scratch=$(mktemp -d /private/tmp/thm4283-replay.XXXXXX)
packet=05-knowledge/results/lrc14_endpoint_carrier_signature_surgery_thm4283
src=04-computation/lrc14_endpoint_carrier_signature_surgery_thm4283
mkdir -p "$scratch/all127" "$scratch/lifts" "$scratch/common" \
  "$scratch/selected12" "$scratch/carrier_detached" "$scratch/proof"
```

## 1. Descending carrier scan

```bash
clang++ -std=c++20 -O3 -DNDEBUG -pthread -I. \
  "$src/endpoint_top_band_scan.cpp" -o "$scratch/endpoint_scan"
"$scratch/endpoint_scan" \
  "$packet/inputs/joint421_masks.txt" \
  "$packet/inputs/reconstructed_final8951.txt" \
  "$packet/inputs/additions45.txt" \
  "$packet/inputs/final_residual23373.csv" 630 \
  "$scratch/base_failures.csv" "$scratch/nested_failures.csv" \
  > "$scratch/top_band_scan.out"
cmp "$scratch/top_band_scan.out" \
  "$packet/results/carrier/top_band_scan.out"
cmp "$scratch/base_failures.csv" \
  "$packet/results/carrier/base_failures.csv"
cmp "$scratch/nested_failures.csv" \
  "$packet/results/carrier/nested_failures.csv"
```

The primary scan is the longest step. It deliberately completes the whole
endpoint-638 layer before stopping at the first nested-carrier failure.

Serialize and replay the exact boundary witness:

```bash
clang++ -std=c++20 -O3 -DNDEBUG -pthread -I. \
  "$src/endpoint638_exact_response_witness.cpp" \
  -o "$scratch/endpoint638_witness"
"$scratch/endpoint638_witness" \
  "$packet/inputs/joint421_masks.txt" \
  "$packet/inputs/reconstructed_final8951.txt" \
  "$packet/inputs/additions45.txt" \
  "$scratch/endpoint638_response_witness9.txt" \
  > "$scratch/endpoint638_response_witness.out"
cmp "$scratch/endpoint638_response_witness.out" \
  "$packet/results/carrier/endpoint638_response_witness.out"
cmp "$scratch/endpoint638_response_witness9.txt" \
  "$packet/results/carrier/endpoint638_response_witness9.txt"
```

## 2. All 127 Boolean-lattice covers

```bash
clang++ -std=c++20 -O3 -DNDEBUG -pthread -I. \
  "$src/all127_response_greedy.cpp" -o "$scratch/all127_greedy"
"$scratch/all127_greedy" \
  "$packet/inputs/joint421_masks.txt" \
  "$packet/inputs/full_signatures_primary.csv" \
  "$packet/inputs/final_residual23373.csv" \
  "$scratch/all127/witnesses" "$scratch/all127/summary.csv" \
  > "$scratch/all127/greedy_atlas.out"
cmp "$scratch/all127/greedy_atlas.out" \
  "$packet/results/all127/greedy_atlas.out"
cmp "$scratch/all127/summary.csv" \
  "$packet/results/all127/all127_summary.csv"
diff -ru "$scratch/all127/witnesses" \
  "$packet/results/all127/witnesses"
```

This step uses about 1.2 GB peak resident memory.

## 3. Four target-aware lifts

```bash
clang++ -std=c++20 -O3 -DNDEBUG -pthread -I. \
  "$src/targeted_response_lift.cpp" -o "$scratch/targeted_lift"
"$scratch/targeted_lift" \
  "$packet/inputs/joint421_masks.txt" \
  "$packet/inputs/full_signatures_primary.csv" \
  "$scratch/lifts/witnesses" \
  > "$scratch/lifts/targeted_response_lift.out"
cmp "$scratch/lifts/targeted_response_lift.out" \
  "$packet/results/target_lifts/targeted_response_lift.out"
diff -ru "$scratch/lifts/witnesses" \
  "$packet/results/target_lifts/witnesses"
```

This step also uses about 1.2 GB peak resident memory.

## 4. Exact common-family synthesis

```bash
clang++ -std=c++20 -O3 -DNDEBUG -pthread -I. \
  "$src/all127_common_family_audit.cpp" -o "$scratch/common_audit"
"$scratch/common_audit" \
  "$packet/inputs/joint421_masks.txt" \
  "$packet/inputs/full_signatures_primary.csv" \
  "$packet/inputs/final_residual23373.csv" \
  "$scratch/all127/witnesses" "$scratch/all127/summary.csv" \
  "$scratch/lifts/witnesses" "$scratch/common/results" \
  "$scratch/common/global_obligations.csv" \
  > "$scratch/common/common_family_audit.out"
cmp "$scratch/common/common_family_audit.out" \
  "$packet/results/common_families/common_family_audit.out"
cmp "$scratch/common/global_obligations.csv" \
  "$packet/results/common_families/global_obligations.csv"
for produced in "$scratch/common/results"/*; do
  cmp "$produced" \
    "$packet/results/common_families/$(basename "$produced")"
done
```

## 5. Detached twelve-deck audit

```bash
clang++ -std=c++20 -O3 -DNDEBUG -pthread \
  "$src/selected12_detached_literal_audit.cpp" \
  -o "$scratch/selected12_audit"
"$scratch/selected12_audit" \
  "$packet/inputs/joint421_masks.txt" \
  "$packet/inputs/full_signatures_primary.csv" \
  "$scratch/all127/witnesses" "$scratch/lifts/witnesses" \
  "$scratch/common/results" \
  "$scratch/common/results/selected_scenario_cover.txt" \
  "$scratch/selected12" > "$scratch/selected12.out"
cmp "$scratch/selected12.out" \
  "$packet/results/detached/detached_literal_audit.out"
cmp "$scratch/selected12/selected12_literal_rows.csv" \
  "$packet/results/detached/selected12_literal_rows.csv"
```

Recompile the same source with `-O0`, and with
`-O1 -g -fsanitize=undefined -fno-sanitize-recover=undefined`; both stdout
files must match `controls/detached_O0.out` and
`controls/detached_ubsan.out`, and both stderr files must be empty.

## 6. Detached carrier audit

```bash
clang++ -std=c++20 -O3 -DNDEBUG -pthread \
  "$src/carrier_band_detached_literal_audit.cpp" \
  -o "$scratch/carrier_detached_audit"
"$scratch/carrier_detached_audit" \
  "$packet/inputs/reconstructed_final8951.txt" \
  "$packet/inputs/additions45.txt" \
  "$packet/inputs/final_residual23373.csv" \
  "$scratch/base_failures.csv" "$scratch/nested_failures.csv" \
  "$scratch/endpoint638_response_witness9.txt" \
  "$scratch/carrier_detached/pair_audit.csv" 4 \
  > "$scratch/carrier_detached.out"
cmp "$scratch/carrier_detached.out" \
  "$packet/results/carrier_detached/detached_literal_audit.out"
cmp "$scratch/carrier_detached/pair_audit.csv" \
  "$packet/results/carrier_detached/pair_audit.csv"
```

The `-O0` and UBSan builds must match
`controls/carrier_detached_O0.out` and
`controls/carrier_detached_ubsan.out`; stderr is empty in both runs.

## 7. Typed proof graph

```bash
python3 "$src/proof_graph_consequence.py" \
  "$packet/inputs/final_residual23373.csv" \
  "$scratch/common/results/full_signature_fibre_common_union.csv" \
  "$scratch/top_band_scan.out" \
  "$scratch/endpoint638_response_witness.out" "$scratch/proof" \
  > "$scratch/proof_graph.out"
cmp "$scratch/proof_graph.out" "$packet/results/proof_graph/run.out"
for produced in "$scratch/proof"/*; do
  cmp "$produced" \
    "$packet/results/proof_graph/$(basename "$produced")"
done

ruby "$src/proof_graph_consequence_independent.rb" \
  "$packet/inputs/final_residual23373.csv" \
  "$scratch/common/results/full_signature_fibre_common_union.csv" \
  "$scratch/top_band_scan.out" \
  "$scratch/endpoint638_response_witness.out" \
  > "$scratch/proof_graph_ruby.out"
cmp "$scratch/proof_graph_ruby.out" "$packet/results/proof_graph/ruby.out"
```

Finally verify the frozen byte ledger:

```bash
(cd "$packet" && sha256sum -c SHA256SUMS)
```
