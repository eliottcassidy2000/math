# THM-4287 reproduction

Run from the repository root.  The recorded promotion used Apple clang 17
with exact integer arithmetic; no floating point, solver, sampling, or
randomized tie-break is used.

```bash
scratch_dir=$(mktemp -d /tmp/thm4287-replay.XXXXXX)
packet=05-knowledge/results/lrc14_repaired_carrier_endpoint637_descent_thm4287
src=04-computation/lrc14_repaired_carrier_endpoint637_descent_thm4287
mkdir -p "$scratch_dir/primary" "$scratch_dir/detached" \
  "$scratch_dir/signature" "$scratch_dir/proof"
```

## 1. Descending carrier scan

```bash
clang++ -std=c++20 -O3 -DNDEBUG -pthread -I. \
  "$src/repaired9006_descending_scan.cpp" \
  -o "$scratch_dir/repaired9006_scan"
"$scratch_dir/repaired9006_scan" \
  "$packet/inputs/joint421_masks.txt" \
  "$packet/inputs/reconstructed_final8951.txt" \
  "$packet/inputs/additions45.txt" \
  "$packet/inputs/endpoint638_response_witness9.txt" \
  "$packet/inputs/current_residual22682.csv" 630 \
  "$scratch_dir/primary/failures.csv" \
  > "$scratch_dir/primary/descending_scan.out"
cmp "$scratch_dir/primary/descending_scan.out" \
  "$packet/results/primary/descending_scan.out"
cmp "$scratch_dir/primary/failures.csv" \
  "$packet/results/primary/failures.csv"
```

The scan completes the entire endpoint-637 and endpoint-636 layers before
stopping.  The additional 8,996/8,997 prefix checks are performed only on the
three endpoint-637 rows and prevent attribution of their closure to the later
repair suffix.

## 2. Detached carrier audit

```bash
clang++ -std=c++20 -O3 -DNDEBUG -pthread \
  "$src/detached_endpoint637_literal_audit.cpp" \
  -o "$scratch_dir/detached/carrier_audit"
"$scratch_dir/detached/carrier_audit" \
  "$packet/inputs/reconstructed_final8951.txt" \
  "$packet/inputs/additions45.txt" \
  "$packet/inputs/endpoint638_response_witness9.txt" \
  "$packet/inputs/current_residual22682.csv" \
  "$scratch_dir/detached/failures.csv" 4 \
  > "$scratch_dir/detached/literal_audit.out"
cmp "$scratch_dir/detached/literal_audit.out" \
  "$packet/results/detached/literal_audit.out"
cmp "$scratch_dir/detached/failures.csv" \
  "$packet/results/detached/failures.csv"
cmp "$scratch_dir/detached/failures.csv" \
  "$packet/results/primary/failures.csv"
```

Repeat the detached source with `-O0`, and with
`-O1 -g -fsanitize=undefined -fno-sanitize-recover=undefined`.  The two
stdout files must match `controls/detached_O0.out` and
`controls/detached_ubsan.out`; stderr must be empty.

## 3. Exact signature-fibre surgery

```bash
clang++ -std=c++20 -O3 -DNDEBUG -pthread -I. \
  "$src/signature275_345_fibre_surgery.cpp" \
  -o "$scratch_dir/signature/surgery"
"$scratch_dir/signature/surgery" \
  "$packet/inputs/joint421_masks.txt" \
  "$packet/inputs/full_signatures_primary.csv" \
  "$packet/inputs/current_residual22682.csv" \
  "$scratch_dir/signature/fibre33.csv" \
  "$scratch_dir/signature/witness2.txt" \
  > "$scratch_dir/signature/surgery.out"
cmp "$scratch_dir/signature/surgery.out" \
  "$packet/results/signature_fibre/surgery.out"
cmp "$scratch_dir/signature/fibre33.csv" \
  "$packet/inputs/signature_fibre33.csv"
cmp "$scratch_dir/signature/witness2.txt" \
  "$packet/results/signature_fibre/witness2.txt"
```

This path enumerates the six deletion obligations, intersects all
`5,852,925` activity bits across all 33 rows, solves the complete 17-class
response quotient, rescans every labelled body for the rebuilt deck, and
checks all rebuilt-deck activity cells directly.

## 4. Detached signature-deck audit

```bash
clang++ -std=c++20 -O3 -DNDEBUG -pthread \
  "$src/detached_signature275_345_literal_audit.cpp" \
  -o "$scratch_dir/signature/detached_audit"
"$scratch_dir/signature/detached_audit" \
  "$packet/inputs/joint421_masks.txt" \
  "$packet/inputs/signature_fibre33.csv" \
  "$packet/results/signature_fibre/witness2.txt" \
  > "$scratch_dir/signature/detached_literal_audit.out"
cmp "$scratch_dir/signature/detached_literal_audit.out" \
  "$packet/results/signature_fibre/detached_literal_audit.out"
```

Repeat with `-O0`, and with
`-O1 -g -fsanitize=undefined -fno-sanitize-recover=undefined`.  The outputs
must match `controls/signature_detached_O0.out` and
`controls/signature_detached_ubsan.out`; stderr must be empty.

## 5. Independent proof-graph consumers

```bash
python3 "$src/proof_graph_consequence.py" \
  "$packet/inputs/current_residual22682.csv" \
  "$packet/inputs/carrier_endpoint637.csv" \
  "$packet/inputs/signature_fibre33.csv" \
  "$scratch_dir/proof/final_python.csv" \
  > "$scratch_dir/proof/primary.out"
env LC_ALL=C LANG=C ruby "$src/proof_graph_consequence_independent.rb" \
  "$packet/inputs/current_residual22682.csv" \
  "$packet/inputs/carrier_endpoint637.csv" \
  "$packet/inputs/signature_fibre33.csv" \
  "$scratch_dir/proof/final_ruby.csv" \
  > "$scratch_dir/proof/independent.out"
cmp "$scratch_dir/proof/primary.out" \
  "$packet/results/proof_graph/primary.out"
cmp "$scratch_dir/proof/independent.out" \
  "$packet/results/proof_graph/independent.out"
cmp "$scratch_dir/proof/primary.out" "$scratch_dir/proof/independent.out"
cmp "$scratch_dir/proof/final_python.csv" \
  "$packet/results/proof_graph/final_residual.csv"
cmp "$scratch_dir/proof/final_ruby.csv" \
  "$packet/results/proof_graph/independent_final_residual.csv"
cmp "$scratch_dir/proof/final_python.csv" "$scratch_dir/proof/final_ruby.csv"
```

Finally verify the artifact packet:

```bash
(cd "$packet" && env LC_ALL=C LANG=C sha256sum -c SHA256SUMS)
```
