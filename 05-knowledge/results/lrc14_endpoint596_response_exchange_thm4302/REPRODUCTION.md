# THM-4302 reproduction

Run from the repository root.  The commands below consume the frozen
THM-4296 and THM-4300 ledgers; they certify a fixed thirty-label carrier only.
They do not construct a physical entry and do not prove LRC(14).

```bash
src=04-computation/lrc14_endpoint596_response_exchange_thm4302
packet=05-knowledge/results/lrc14_endpoint596_response_exchange_thm4302
old=05-knowledge/results/lrc14_mixed_rank_depth_recursive_signatures_thm4296
prior=05-knowledge/results/lrc14_size_preserving_response_staircase_thm4300
scratch=$(mktemp -d)

common_args=(
  "$old/inputs/joint421_masks.txt"
  "$old/inputs/reconstructed_final8951.txt"
  "$old/inputs/additions45.txt"
  "$old/inputs/endpoint638_response_witness9.txt"
  "$old/inputs/current_residual22647.csv"
  "$prior/inputs/repairs76.txt"
)
```

## 1. Complete response quotient and prefix exchange

The primary program scans every rank-eight and rank-nine mask at `(210,596)`,
checks all 718 realized responses and both exact dual/cover pairs, and evaluates
all `9,088*363=3,298,944` augmented-carrier signs on the complete
endpoint-`>=596` prefix.  The O0 run is deliberately slower.

```bash
g++ -std=c++20 -O3 -DNDEBUG -I. \
  "$src/endpoint596_response_exchange_audit.cpp" \
  -o "$scratch/primary_O3"
g++ -std=c++20 -O0 -I. \
  "$src/endpoint596_response_exchange_audit.cpp" \
  -o "$scratch/primary_O0"

primary_tail=(
  "$prior/results/endpoint/augmented9088_failures.csv"
  "$packet/inputs/additions4.txt"
  "$packet/inputs/delete73.txt"
)
"$scratch/primary_O3" "${common_args[@]}" "${primary_tail[@]}" \
  > "$scratch/primary_O3.out"
"$scratch/primary_O0" "${common_args[@]}" "${primary_tail[@]}" \
  > "$scratch/primary_O0.out"

cmp "$scratch/primary_O3.out" "$packet/results/primary_O3.out"
cmp "$scratch/primary_O0.out" "$packet/results/controls/primary_O0.out"
cmp "$scratch/primary_O3.out" "$scratch/primary_O0.out"
```

## 2. Independent failure-universe and exchange replay

This program deliberately receives no failure CSV.  It reconstructs the
augmented 9,088-mask carrier, enumerates all `C(30,9)=14,307,150` labelled
nine-bodies at `(210,596)`, derives the 24 failures and their FNV, and repeats
the raw replay after the four additions and after the full exchange.  Its
response inventory is generated from complement subsets of those bodies,
rather than by the primary program's full-universe response scan.

```bash
g++ -std=c++20 -O3 -DNDEBUG -I. \
  "$src/endpoint596_independent_audit.cpp" \
  -o "$scratch/independent"

"$scratch/independent" "${common_args[@]}" \
  "$packet/inputs/additions4.txt" "$packet/inputs/delete73.txt" \
  > "$scratch/independent.out"

cmp "$scratch/independent.out" "$packet/results/independent_audit.out"
```

## 3. Typed row consequence

The set consumer derives the complete endpoint-`>=596` prefix from the
22,647-row universe and unions it with THM-4300's typed row set.  It checks the
old partition, all intersections, the new partition, and the complete top
layer.

```bash
python3 -B "$src/typed_union_consumer.py" \
  "$old/inputs/current_residual22647.csv" \
  "$prior/results/proof_graph/typed_union1624.csv" \
  "$prior/results/proof_graph/final_residual21023.csv" \
  > "$scratch/typed_union_consumer.out"

cmp "$scratch/typed_union_consumer.out" \
  "$packet/results/proof_graph/typed_union_consumer.out"
```

## 4. Packet identity

`SHA256SUMS` uses raw LF bytes and excludes itself.  On a system with GNU
Coreutils, verify it with:

```bash
(cd "$packet" && sha256sum -c SHA256SUMS)
```

On macOS the equivalent check is:

```bash
(cd "$packet" && shasum -a 256 -c SHA256SUMS)
```
