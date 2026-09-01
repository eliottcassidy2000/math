# THM-4303 reproduction

Run from the repository root.  Both raw programs reconstruct THM-4302's
9,019-mask carrier and derive the 28 audited rows from the frozen residual
universe and THM-4300 typed union.  These commands prove fixed-carrier row
consequences only; they do not construct physical entry and do not prove
LRC(14).

```bash
src=04-computation/lrc14_endpoint595_twentyfive_closure_thm4303
packet=05-knowledge/results/lrc14_endpoint595_twentyfive_closure_thm4303
old=05-knowledge/results/lrc14_mixed_rank_depth_recursive_signatures_thm4296
thm4300=05-knowledge/results/lrc14_size_preserving_response_staircase_thm4300
thm4302=05-knowledge/results/lrc14_endpoint596_response_exchange_thm4302
scratch=$(mktemp -d)
workers=8
```

## 1. Primary raw replay

The primary program reconstructs `C_596`, derives THM-4302's residual and its
complete maximum layer, and replays all `C(30,9)=14,307,150` labelled bodies
on each of the 28 rows.  It counts every nonjoint hit on every joint-exposed
body and freezes the complete pair and failure ledgers.

```bash
g++ -std=c++20 -O3 -DNDEBUG -pthread -I. \
  "$src/endpoint595_primary_raw_audit.cpp" \
  -o "$scratch/primary"

"$scratch/primary" \
  "$old/inputs/joint421_masks.txt" \
  "$old/inputs/reconstructed_final8951.txt" \
  "$old/inputs/additions45.txt" \
  "$old/inputs/endpoint638_response_witness9.txt" \
  "$old/inputs/current_residual22647.csv" \
  "$thm4300/results/proof_graph/typed_union1624.csv" \
  "$thm4300/inputs/repairs76.txt" \
  "$thm4302/inputs/additions4.txt" \
  "$thm4302/inputs/delete73.txt" \
  "$scratch/endpoint595_pair_audit.csv" \
  "$scratch/endpoint595_failures.csv" \
  "$workers" > "$scratch/primary_raw_audit.out"

cmp "$scratch/primary_raw_audit.out" \
  "$packet/results/primary_raw_audit.out"
cmp "$scratch/endpoint595_pair_audit.csv" \
  "$packet/results/endpoint595_pair_audit.csv"
cmp "$scratch/endpoint595_failures.csv" \
  "$packet/results/endpoint595_failures.csv"
```

## 2. Self-contained independent replay

The independent program does not include or call the primary raw routine.  It
implements the literal-wall geometry directly, reconstructs the carrier from
the component ledgers, generates nine-bodies recursively and sorts them into
numeric order, and stops at the first nonjoint witness for each exposed body.

```bash
g++ -std=c++20 -O3 -DNDEBUG -pthread -I. \
  "$src/endpoint595_independent_replay.cpp" \
  -o "$scratch/independent"

"$scratch/independent" \
  "$old/inputs/joint421_masks.txt" \
  "$old/inputs/reconstructed_final8951.txt" \
  "$old/inputs/additions45.txt" \
  "$old/inputs/endpoint638_response_witness9.txt" \
  "$thm4300/inputs/repairs76.txt" \
  "$thm4302/inputs/delete73.txt" \
  "$thm4302/inputs/additions4.txt" \
  "$old/inputs/current_residual22647.csv" \
  "$thm4300/results/proof_graph/typed_union1624.csv" \
  "$workers" > "$scratch/independent_replay.out"

cmp "$scratch/independent_replay.out" \
  "$packet/results/independent_replay.out"
```

The exact field consumer compares all eight common per-row values from the two
implementations on all 28 rows:

```bash
python3 -B "$src/primary_independent_crosscheck.py" \
  "$scratch/endpoint595_pair_audit.csv" \
  "$scratch/independent_replay.out" \
  > "$scratch/primary_independent_crosscheck.out"

cmp "$scratch/primary_independent_crosscheck.out" \
  "$packet/results/primary_independent_crosscheck.out"
```

## 3. Typed row-set consequence

This consumer re-derives THM-4302's partition, requires that the primary pair
ledger is exactly its complete 28-row top layer, and unions precisely the 25
zero-failure row consequences.

```bash
python3 -B "$src/typed_union_consumer.py" \
  "$old/inputs/current_residual22647.csv" \
  "$thm4300/results/proof_graph/typed_union1624.csv" \
  "$thm4300/results/proof_graph/final_residual21023.csv" \
  "$scratch/endpoint595_pair_audit.csv" \
  > "$scratch/typed_union_consumer.out"

cmp "$scratch/typed_union_consumer.out" \
  "$packet/results/proof_graph/typed_union_consumer.out"
```

## 4. Packet identity

`SHA256SUMS` uses raw LF bytes and excludes itself.  With GNU Coreutils:

```bash
(cd "$packet" && sha256sum -c SHA256SUMS)
```

On macOS:

```bash
(cd "$packet" && shasum -a 256 -c SHA256SUMS)
```
