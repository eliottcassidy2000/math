---
date: 2026-05-31
session: codex-2026-05-31-S366
tags: [formalization, lean, product-sum, natural-modes]
related:
  - THM-361
  - HYP-1821
  - 04-computation/lean/TournamentH7/TournamentH7/ProductSum.lean
  - 05-knowledge/results/lean_product_sum_s366.out
---

# Product-Sum Defect in Lean

The natural-mode graph session produced a clean theorem:

```text
product(M) = sum(M)
```

is equivalent to a defect ledger after stripping the `1`s.  The formalization
landed the list-level heart of that statement in `TournamentH7.ProductSum`.

## Formal Core

The Lean module defines:

```lean
ones xs := (xs.filter (fun x => x = 1)).length
core xs := xs.filter (fun x => x ≠ 1)
IsProductSum xs := xs.prod = xs.sum
```

and proves:

```lean
xs.prod = (core xs).prod
xs.sum = (core xs).sum + ones xs
IsProductSum xs <-> (core xs).prod = (core xs).sum + ones xs
```

This is exactly THM-361 without the multiset quotient: lists keep order, but
the theorem itself is invariant under permutation because the product and sum
folds are.

## Converse

The repair direction is also formal:

```lean
c.prod = c.sum + d
```

implies:

```lean
IsProductSum ((List.replicate d 1) ++ c)
```

So the entry `1` is now a literal theorem-prover object: one unit of additive
slack that does not touch the product.

## First Resonance

The session also formalized the ordered positive two-entry clause:

```lean
0 < a -> 0 < b -> IsProductSum [a,b] -> a = 2 /\ b = 2
```

This pins down the hidden diagonal resonance `2+2=2*2`.  The distinct ternary
resonance `{1,2,3}` is not yet formalized, but the route is visible: prove the
finite sorted/distinct core reduction, then show every distinct no-one core of
size at least three has positive defect too large for the available number of
ones except the two-core `(2,3)` with one repaired unit.

## Why This Matters

The old additive graph on natural modes collapses to the transitive completion
when one forgets labels, while multiplication stays sparse as divisibility.
The product-sum defect theorem gives a formal bridge between them:

```text
additive equality = multiplicative excess + identity slack.
```

That bridge is now not just a story in the canon; it has a Lean surface that
future sessions can import, audit, and extend.
