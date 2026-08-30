# Reproduction

Run with Bash from the repository root.  C++20, Python 3, Perl 5, and
`sha256sum` are required.

```bash
set -euo pipefail
ROOT=$(git rev-parse --show-toplevel)
SRC="$ROOT/04-computation/lrc14_signature_response_congruence_thm4286"
PACKET="$ROOT/05-knowledge/results/lrc14_signature_response_congruence_thm4286"
INHERITED="$ROOT/05-knowledge/results/lrc14_inactive_signature_deck_surgery_thm4282"
DECK="$INHERITED/components/index367/results/original_joint421_deck.txt"
SIG="$INHERITED/components/index367/results/post_thm4281_full421_inactive_signatures.csv"
PRIOR="$INHERITED/components/surgery256/results/final_residual23373.csv"
THM4283="$ROOT/05-knowledge/results/lrc14_endpoint_carrier_signature_surgery_thm4283/results/proof_graph"
THM4283_UNION="$THM4283/proof_union.csv"
POST4283="$THM4283/final_residual.csv"
TMP=$(mktemp -d)
CXX=${CXX:-g++}

cd "$ROOT"
sha256sum -c "$PACKET/SHA256SUMS"
```

## 1. Signature-response factorization census

```bash
python3 "$SRC/audit.py" > "$TMP/signature.out"
python3 -O "$SRC/audit.py" > "$TMP/signature.opt.out"
perl "$SRC/independent_audit.pl" > "$TMP/signature.independent.out"
cmp "$TMP/signature.out" "$TMP/signature.opt.out"
cmp "$TMP/signature.out" "$PACKET/signature_response_primary.out"
cmp "$TMP/signature.independent.out" \
  "$PACKET/signature_response_independent.out"
```

## 2. Compile the exact engines

```bash
for NAME in \
  singleton_fibre_cegar \
  singleton_fibre_literal_verify \
  singleton_fibre_literal_exhaustive \
  two_mask_fibre_cegar \
  two_mask_literal_verify \
  two_mask_full_gain \
  two_mask_net36_literal_verify
do
  "$CXX" -std=c++20 -O3 -DNDEBUG -I"$ROOT" \
    "$SRC/$NAME.cpp" -o "$TMP/$NAME"
done
```

## 3. Index-396 surgery

```bash
"$TMP/singleton_fibre_cegar" "$DECK" "$SIG" 396 366 644 \
  > "$TMP/singleton_cegar.out"
"$TMP/singleton_fibre_literal_verify" \
  "$DECK" "$SIG" 396 366 644 042022c9 \
  > "$TMP/singleton_literal.out"
"$TMP/singleton_fibre_literal_exhaustive" \
  "$DECK" "$SIG" 396 366 644 \
  > "$TMP/singleton_literal_exhaustive.out"
"$TMP/two_mask_full_gain" "$SIG" "$PRIOR" index396 \
  > "$TMP/full_gain_index396.out"

cmp "$TMP/singleton_cegar.out" "$PACKET/singleton_fibre_cegar.out"
cmp "$TMP/singleton_literal.out" "$PACKET/singleton_fibre_literal.out"
cmp "$TMP/singleton_literal_exhaustive.out" \
  "$PACKET/singleton_fibre_literal_exhaustive.out"
cmp "$TMP/full_gain_index396.out" "$PACKET/full_gain_index396.out"

awk 'NR!=397 {print} END {print "042022c9"}' "$DECK" \
  > "$TMP/index396_rebuilt421.txt"
cmp "$TMP/index396_rebuilt421.txt" "$PACKET/index396_rebuilt421.txt"
```

## 4. Two-mask minimum-three surgery and exact shrink

```bash
"$TMP/two_mask_fibre_cegar" \
  "$DECK" "$SIG" 107 374 512 644 fibre \
  > "$TMP/two_mask_cegar.out"
"$TMP/two_mask_literal_verify" \
  "$DECK" "$SIG" 107 374 512 644 \
  32043014 20807016 128c8012 \
  > "$TMP/two_mask_literal.out"
"$TMP/two_mask_full_gain" "$SIG" "$PRIOR" two-mask \
  > "$TMP/full_gain_two_mask.out"
"$TMP/two_mask_net36_literal_verify" \
  "$DECK" "$SIG" "$PRIOR" 32043014 20807016 128c8012 \
  > "$TMP/net36_literal.out"

cmp "$TMP/two_mask_cegar.out" "$PACKET/two_mask_fibre_cegar.out"
cmp "$TMP/two_mask_literal.out" "$PACKET/two_mask_fibre_literal.out"
cmp "$TMP/full_gain_two_mask.out" "$PACKET/full_gain_two_mask.out"
cmp "$TMP/net36_literal.out" "$PACKET/two_mask_net36_literal.out"

sed -n 's/^NET \([0-9][0-9]*,[0-9][0-9]*\).*/\1/p' \
  "$TMP/net36_literal.out" > "$TMP/two_mask_net36.csv"
cmp "$TMP/two_mask_net36.csv" "$PACKET/two_mask_net36.csv"

awk 'NR!=108 && NR!=319 && NR!=375 {print}
     END {print "32043014"; print "20807016"; print "128c8012"}' \
  "$DECK" > "$TMP/two_mask_rebuilt421.txt"
cmp "$TMP/two_mask_rebuilt421.txt" "$PACKET/two_mask_rebuilt421.txt"
```

## 5. Exact proof graph

```bash
python3 "$SRC/proof_graph_consequence.py" \
  --post4282 "$PRIOR" --thm4283-union "$THM4283_UNION" \
  --post4283 "$POST4283" --signatures "$SIG" \
  --two-mask-net36 "$PACKET/two_mask_net36.csv" \
  --emit-fibre36 "$TMP/index396_fibre36.csv" \
  --emit-redundant-union72 "$TMP/alternate_node_union72.csv" \
  > "$TMP/proof_graph.out"

python3 -O "$SRC/proof_graph_consequence.py" \
  --post4282 "$PRIOR" --thm4283-union "$THM4283_UNION" \
  --post4283 "$POST4283" --signatures "$SIG" \
  --two-mask-net36 "$PACKET/two_mask_net36.csv" \
  > "$TMP/proof_graph.opt.out"

perl "$SRC/proof_graph_subsumption_independent.pl" \
  > "$TMP/proof_graph.independent.out"

cmp "$TMP/proof_graph.out" "$TMP/proof_graph.opt.out"
cmp "$TMP/proof_graph.out" "$PACKET/proof_graph_consequence.out"
cmp "$TMP/proof_graph.independent.out" \
  "$PACKET/proof_graph_subsumption_independent.out"
cmp "$TMP/index396_fibre36.csv" "$PACKET/index396_fibre36.csv"
cmp "$TMP/alternate_node_union72.csv" "$PACKET/alternate_node_union72.csv"
```

For an undefined-behavior control, rebuild the seven C++ programs with
`-O1 -g -fsanitize=undefined -fno-sanitize-recover=all` and rerun the same
commands.  The maintained sources pin all consequence-bearing counts, FNVs,
strict margins, deck identities, and body-scan ledgers before printing PASS.
