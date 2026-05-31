---
id: HYP-1821
status: EXPLORATORY
source: codex-2026-05-31-S365b
related:
  - THM-361
---

# HYP-1821: Natural mode graphs meet along the product-sum defect interface

## Statement

The additive and multiplicative natural-number mode graphs should be studied
through the defect

```text
D(C) = product(C) - sum(C)
```

on cores `C` whose entries are at least `2`.

The additive mode graph is dense: from the seeds `{2,3}` it reaches every
positive integer except `{1,4,6}` under the old distinct-summand rule.  The
multiplicative mode graph is sparse: its atoms are primes and prime squares.
The family

```text
x_1 + ... + x_k = x_1 * ... * x_k
```

is the interface between these regimes.  By THM-361, stripping all `1`s turns
every such equation into the exact slack law

```text
D(C) = number of stripped 1s.
```

## Evidence

`natural_mode_graph_s365.py` rechecked the older summand/product graph picture.
Up to `80`, the additive graph has `78` non-atom nodes while the product graph
has only `53` non-atom nodes.  The product atoms up to `80` are exactly primes
and prime squares:

```text
2,3,4,5,7,9,11,13,17,19,23,25,29,31,37,41,43,47,49,...
```

The same run found:

```text
binary distinct x+y=x*y: none
binary diagonal allowed: (2,2)->4
ternary distinct x+y+z=xyz: (1,2,3)->6
```

The product-sum ledger through arity `10` begins:

```text
k=2:  (2,2)
k=3:  (1,2,3)
k=4:  (1,1,2,4)
k=5:  (1,1,1,2,5), (1,1,1,3,3), (1,1,2,2,2)
```

and every row obeys the defect normal form.

The S366 Lean session formalized the list-level normal form in
`TournamentH7.ProductSum`: deleting all `1`s preserves product and records the
sum loss as `ones`, adjoining exactly `d` ones repairs any core with
`prod=sum+d`, and the ordered positive two-entry case is forced to `(2,2)`.
The full build `lake build TournamentH7` succeeds with these audits.

## Interpretation

The old exceptional additive module `{1,4,6}` becomes legible as a three-part
resonance boundary:

```text
1 = source / multiplicative identity
4 = hidden binary diagonal resonance 2+2=2*2
6 = first visible distinct ternary resonance 1+2+3=1*2*3
```

This suggests a general tactic for arithmetic graph analogies in the repo:
look for a dense closure process, a sparse factorization process, and the
defect layer where identities or diagonal collisions compensate for the gap.

## Predictions

1. Other previously observed "special numbers" should often decompose into
   source, diagonal, and first-distinct-resonance roles for an appropriate mode
   graph.
2. Product-sum cores of fixed defect should form a small finite atlas useful
   for explaining why `6`, `42`, and related self-product constants recur.
3. Tournament forbidden-value phenomena should be compared against defect
   atlases, not raw sums or products alone.

## Sources

- THM-361.
- `04-computation/lean/TournamentH7/TournamentH7/ProductSum.lean`.
- `05-knowledge/results/lean_product_sum_s366.out`.
- `04-computation/natural_mode_graph_s365.py`.
- `05-knowledge/results/natural_mode_graph_s365.out`.
- `07-reflections/summand-graph-fermat-zeckendorf.md`.
- `07-reflections/product-graph-sc-spine-fractal-dimensions.md`.
