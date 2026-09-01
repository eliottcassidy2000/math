# THM-4305 reproduction

Run from the repository root. The exact-cover programs require Python 3 and
OR-Tools; the producing runs used OR-Tools `9.15.6755`. All C++ programs use
only the standard library and inherited audited sources.

```bash
packet=05-knowledge/results/lrc14_endpoint595_pair_tagged_response_exchange_thm4305
code=04-computation/lrc14_endpoint595_pair_tagged_response_exchange_thm4305
old=05-knowledge/results/lrc14_mixed_rank_depth_recursive_signatures_thm4296
t4300=05-knowledge/results/lrc14_size_preserving_response_staircase_thm4300
t4302=05-knowledge/results/lrc14_endpoint596_response_exchange_thm4302
t4303=05-knowledge/results/lrc14_endpoint595_twentyfive_closure_thm4303
scratch=$(mktemp -d)
```

## 1. Complete rank-eight/rank-nine atlas and exact covers

```bash
clang++ -std=c++20 -O3 -DNDEBUG -pthread -I. \
  "$code/pair_tagged_atlas.cpp" -o "$scratch/atlas_O3"
clang++ -std=c++20 -O0 -pthread -I. \
  "$code/pair_tagged_atlas.cpp" -o "$scratch/atlas_O0"

for opt in O3 O0; do
  "$scratch/atlas_$opt" \
    "$old/inputs/reconstructed_final8951.txt" \
    "$old/inputs/additions45.txt" \
    "$old/inputs/endpoint638_response_witness9.txt" \
    "$t4300/inputs/repairs76.txt" \
    "$t4302/inputs/delete73.txt" \
    "$t4302/inputs/additions4.txt" \
    "$t4303/results/endpoint595_failures.csv" \
    "$scratch/atlas_$opt.csv" > "$scratch/atlas_$opt.out"
  sed "s#$scratch/atlas_$opt.csv#/tmp/lrc595_pair_tagged_atlas.csv#" \
    "$scratch/atlas_$opt.out" > "$scratch/atlas_$opt.normalized.out"
done

cmp "$scratch/atlas_O3.csv" "$scratch/atlas_O0.csv"
cmp "$scratch/atlas_O3.csv" "$packet/results/pair_tagged_atlas.csv"
cmp "$scratch/atlas_O3.normalized.out" "$packet/results/pair_tagged_atlas_O3.out"
cmp "$scratch/atlas_O0.normalized.out" "$packet/results/pair_tagged_atlas_O0.out"

python3 -B "$code/exact_cover.py" "$packet/results/pair_tagged_atlas.csv" \
  > "$scratch/exact_cover.out"
cmp "$scratch/exact_cover.out" "$packet/results/exact_cover.out"
```

The atlas scans all `5,852,925 + 14,307,150` masks. The mixed lower bound is
the frozen integral dual in `mixed_dual_certificate.txt`; the rank-eight lower
bound is the deterministic infeasibility test in `exact_cover.py`.

## 2. Arbitrary-rank exact minimum five

The four-mask model contains every one of the `2^30` supports independently
in each slot. It retains all literal-wall failure classes and constrains both
activity and tagged hits as equivalences.

```bash
python3 -B "$code/rank_free_exact_cover.py" \
  "$t4303/results/endpoint595_failures.csv" \
  --mode general --masks 4 --seconds 180 --workers 4 \
  --first-solution --pure-feasibility --quiet \
  > "$scratch/no_four.out"
grep 'STATUS INFEASIBLE' "$scratch/no_four.out"

python3 -B "$code/rank_free_exact_cover.py" \
  "$t4303/results/endpoint595_failures.csv" \
  --mode general --masks 4 --seconds 180 --workers 1 \
  --first-solution --no-presolve --symmetry 0 --quiet \
  > "$scratch/no_four_hostile.out"
grep 'STATUS INFEASIBLE' "$scratch/no_four_hostile.out"

python3 -B "$code/rank_free_exact_cover.py" \
  "$t4303/results/endpoint595_failures.csv" \
  --mode general --masks 5 --seconds 180 --workers 4 \
  --first-solution --pure-feasibility --quiet \
  --hint 00612a76 00a183f6 024d8b32 0280a1ae 10110bf6 \
  > "$scratch/five.out"
grep 'JOINED 116,13,16 TOTAL 145 UNCOVERED 0' "$scratch/five.out"
grep 'WITNESS_FNV ebea2511eb7fa46f RESPONSE_FNV 8e13cc00c5bca917' \
  "$scratch/five.out"
```

Search telemetry can vary, so the decisive exact statuses and witness
identities are compared rather than wall time. The producing transcripts are
frozen in `rank_free_four_mask_infeasible.out`,
`rank_free_four_mask_hostile.out`, and `rank_free_five_mask_cover.out`.

The two positive controls are reproduced by:

```bash
python3 -B "$code/rank_free_exact_cover.py" \
  "$t4303/results/endpoint595_failures.csv" \
  --mode general --masks 4 --families 96 --seconds 180 --workers 4 \
  --first-solution --pure-feasibility --quiet \
  --hint 00310bf6 022183ae 0240a176 108122fe

python3 -B "$code/rank_free_exact_cover.py" \
  "$t4303/results/endpoint595_failures.csv" \
  --mode general --masks 4 --families 100 210 --seconds 180 --workers 4 \
  --first-solution --pure-feasibility --quiet \
  --hint 081a03d6 08e58a84 10410a3f 18270a72
```

The independent literal evaluator does not call OR-Tools:

```bash
clang++ -std=c++20 -O0 -I. "$code/rank_free_cover_independent.cpp" \
  -o "$scratch/five_O0"
clang++ -std=c++20 -O3 -DNDEBUG -I. \
  "$code/rank_free_cover_independent.cpp" -o "$scratch/five_O3"
"$scratch/five_O0" "$t4303/results/endpoint595_failures.csv" \
  > "$scratch/five_O0.out"
"$scratch/five_O3" "$t4303/results/endpoint595_failures.csv" \
  > "$scratch/five_O3.out"
cmp "$scratch/five_O0.out" "$scratch/five_O3.out"
cmp "$scratch/five_O3.out" "$packet/results/rank_free_five_mask_independent.out"
```

## 3. Independent arbitrary-rank two-mask controls

```bash
clang++ -std=c++20 -O3 -DNDEBUG -I. \
  "$code/rank_free_two_mask_intersection.cpp" -o "$scratch/intersection"
"$scratch/intersection" "$t4303/results/endpoint595_failures.csv" \
  > "$scratch/intersection.out"
cmp "$scratch/intersection.out" \
  "$packet/results/rank_free_two_mask_intersection_O3.out"

clang++ -std=c++20 -O3 -DNDEBUG -I. \
  "$code/rank_free_two_mask_nextclosure.cpp" -o "$scratch/nextclosure"
"$scratch/nextclosure" "$t4303/results/endpoint595_failures.csv" \
  > "$scratch/nextclosure.out"
cmp "$scratch/nextclosure.out" \
  "$packet/results/rank_free_two_mask_nextclosure.out"
```

The implementations enumerate the same 548 formal concepts by different
closure algorithms. Their FNVs differ because their traversal orders differ.

## 4. Protected 391-row common-inactive audit

```bash
common_args=(
  "$old/inputs/joint421_masks.txt"
  "$old/inputs/reconstructed_final8951.txt"
  "$old/inputs/additions45.txt"
  "$old/inputs/endpoint638_response_witness9.txt"
  "$old/inputs/current_residual22647.csv"
  "$t4300/results/proof_graph/typed_union1624.csv"
  "$t4300/inputs/repairs76.txt"
  "$t4302/inputs/additions4.txt"
  "$t4302/inputs/delete73.txt"
)

clang++ -std=c++20 -O3 -DNDEBUG -pthread -I. \
  "$code/protected391_activity.cpp" -o "$scratch/protected_O3"
clang++ -std=c++20 -O0 -pthread -I. \
  "$code/protected391_activity.cpp" -o "$scratch/protected_O0"

for opt in O3 O0; do
  "$scratch/protected_$opt" "${common_args[@]}" \
    "$scratch/protected_$opt.rows.csv" \
    "$scratch/protected_$opt.masks.csv" \
    "$scratch/protected_$opt.candidate.csv" \
    > "$scratch/protected_$opt.out"
done

cmp "$scratch/protected_O3.out" "$scratch/protected_O0.out"
cmp "$scratch/protected_O3.out" "$packet/results/protected391_O3.out"
cmp "$scratch/protected_O3.rows.csv" "$packet/results/protected391_rows.csv"
cmp "$scratch/protected_O3.masks.csv" "$packet/results/protected391_masks.csv"
cmp "$scratch/protected_O3.candidate.csv" \
  "$packet/results/protected391_candidate.csv"
```

The structurally independent replay consumes THM-4303's carrier and geometry
implementation:

```bash
clang++ -std=c++20 -O3 -DNDEBUG -pthread -I. \
  "$code/protected391_independent.cpp" -o "$scratch/protected_independent"
"$scratch/protected_independent" \
  "$old/inputs/joint421_masks.txt" \
  "$old/inputs/reconstructed_final8951.txt" \
  "$old/inputs/additions45.txt" \
  "$old/inputs/endpoint638_response_witness9.txt" \
  "$t4300/inputs/repairs76.txt" \
  "$t4302/inputs/delete73.txt" \
  "$t4302/inputs/additions4.txt" \
  "$old/inputs/current_residual22647.csv" \
  "$t4300/results/proof_graph/typed_union1624.csv" \
  > "$scratch/protected_independent.out"
cmp "$scratch/protected_independent.out" \
  "$packet/results/protected391_independent.out"
```

## 5. Typed row consequence

```bash
python3 -B "$code/typed_union_consumer.py" \
  "$old/inputs/current_residual22647.csv" \
  "$t4300/results/proof_graph/typed_union1624.csv" \
  "$t4300/results/proof_graph/final_residual21023.csv" \
  "$t4303/results/endpoint595_pair_audit.csv" \
  "$packet/results/rank_free_five_mask_independent.out" \
  "$scratch/typed_union1661.csv" \
  "$scratch/final_residual20986.csv" \
  "$scratch/final_residual_top594.csv" \
  > "$scratch/typed_union_consumer.out"

cmp "$scratch/typed_union_consumer.out" \
  "$packet/results/typed_union_consumer.out"
cmp "$scratch/typed_union1661.csv" "$packet/results/typed_union1661.csv"
cmp "$scratch/final_residual20986.csv" \
  "$packet/results/final_residual20986.csv"
cmp "$scratch/final_residual_top594.csv" \
  "$packet/results/final_residual_top594.csv"
```

Finally verify the frozen artifact manifest:

```bash
(cd "$packet" && LC_ALL=C shasum -a 256 -c SHA256SUMS)
```
