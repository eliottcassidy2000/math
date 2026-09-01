# THM-4306 reproduction

Run from the repository root.  These commands certify a separate common deck
on the fixed singleton-signature ideal `H_265`; they do not merge that deck
with any carrier, construct physical entry, or prove LRC(14).

```bash
src=04-computation/lrc14_index265_recursive_ideal_thm4306
packet=05-knowledge/results/lrc14_index265_recursive_ideal_thm4306
old=05-knowledge/results/lrc14_mixed_rank_depth_recursive_signatures_thm4296
t4300=05-knowledge/results/lrc14_size_preserving_response_staircase_thm4300
scratch=$(mktemp -d)
workers=8
```

## 1. Complete response quotient

The primary computation traverses the complete rank-eight universe on all
367 rows.  It writes the ideal and complete response quotient.  This is the
long run (about nineteen minutes on the original host).

```bash
g++ -std=c++20 -O3 -DNDEBUG -pthread -I. \
  "$src/index265_response_quotient.cpp" -o "$scratch/primary"

"$scratch/primary" \
  "$old/inputs/joint421_masks.txt" \
  "$old/inputs/full_signatures_primary.csv" \
  "$old/inputs/current_residual22647.csv" \
  "$scratch/index265_ideal367.csv" \
  "$scratch/index265_response_quotient.csv" "$workers" \
  > "$scratch/primary.out"

cmp "$scratch/index265_ideal367.csv" \
  "$packet/results/index265_ideal367.csv"
cmp "$scratch/index265_response_quotient.csv" \
  "$packet/results/index265_response_quotient.csv"
cmp "$scratch/primary.out" "$packet/results/primary_response_quotient.out"
```

## 2. Independent direct-atom audit

This path does not consume the quotient.  It derives `H_265`, the private
bodies, the two replacement masks, and the lower-certificate candidate
universe directly from literal-wall atoms.  It exhausts all `C(17,8)=24,310`
possible one-mask lower candidates and scans all `C(30,9)=14,307,150` bodies
after rebuilding the deck.

```bash
g++ -std=c++20 -O3 -DNDEBUG -pthread -I. \
  "$src/index265_direct_audit.cpp" -o "$scratch/direct_O3"
g++ -std=c++20 -O0 -pthread -I. \
  "$src/index265_direct_audit.cpp" -o "$scratch/direct_O0"

for opt in O3 O0; do
  "$scratch/direct_$opt" \
    "$old/inputs/joint421_masks.txt" \
    "$old/inputs/full_signatures_primary.csv" \
    "$old/inputs/current_residual22647.csv" "$workers" \
    > "$scratch/direct_$opt.out"
done

cmp "$scratch/direct_O3.out" "$packet/results/direct_O3.out"
cmp "$scratch/direct_O0.out" "$packet/results/controls/direct_O0.out"
cmp "$scratch/direct_O3.out" "$scratch/direct_O0.out"
```

## 3. Typed row-set consequence

The consumer independently audits the response quotient, reconstructs the
THM-4303 and THM-4305 controls, unions only the proved `H_265` row
consequence, and writes canonical LF pair sets.

```bash
python3 -B "$src/typed_union_consumer.py" \
  --universe "$old/inputs/current_residual22647.csv" \
  --thm4300-union "$t4300/results/proof_graph/typed_union1624.csv" \
  --thm4300-residual "$t4300/results/proof_graph/final_residual21023.csv" \
  --ideal "$packet/results/index265_ideal367.csv" \
  --quotient "$packet/results/index265_response_quotient.csv" \
  --output-dir "$scratch/typed" \
  > "$scratch/typed.out"

cmp "$scratch/typed.out" "$packet/results/typed_union_consumer.out"
cmp "$scratch/typed/h265_overlap14.csv" \
  "$packet/results/proof_graph/typed_overlap14.csv"
cmp "$scratch/typed/h265_new353.csv" \
  "$packet/results/proof_graph/typed_new353.csv"
cmp "$scratch/typed/combined_union2014.csv" \
  "$packet/results/proof_graph/typed_union2014.csv"
cmp "$scratch/typed/combined_residual20633.csv" \
  "$packet/results/proof_graph/final_residual20633.csv"
cmp "$scratch/typed/combined_top594_22.csv" \
  "$packet/results/proof_graph/residual_top594.csv"
```

## 4. Packet identity

`SHA256SUMS` uses raw LF bytes and excludes itself:

```bash
(cd "$packet" && sha256sum -c SHA256SUMS)
```

