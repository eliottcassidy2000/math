---
id: THM-410
title: Interval-reversal triangle count
status: PROVED
source: monad-formalizer-2026-06-04-S1
depends_on:
  - THM-399
related:
  - INV-190
  - THM-316
  - THM-317
  - THM-319
---

# THM-410 - Interval-Reversal Triangle Count

## Statement

Let `V` be a finite linearly ordered vertex set.  Start with the transitive
tournament oriented by the order:

```text
x -> y  iff  x < y.
```

Let `M` be a matching of ordered pairs `(a,b)` with `a < b`; that is, no vertex
appears in two pairs of `M`.  Form a new tournament `T_M` by reversing exactly
the arcs in `M`.

Then a 3-element set `{x,y,z}` with `x < y < z` is a directed 3-cycle in `T_M`
if and only if `(x,z) in M`.

Consequently, the number of directed 3-cycles in `T_M` is

```text
sum_{(a,b) in M} #{v in V : a < v < b}.
```

## Proof

Fix a triple `x < y < z`.  Before any reversals, its arcs are

```text
x -> y,   y -> z,   x -> z,
```

so the induced subtournament is transitive.

Because `M` is a matching, at most one of the three pair-arcs

```text
(x,y), (y,z), (x,z)
```

can be reversed.

There are four cases.

1. No edge of the triple is reversed.  The triple remains transitive.
2. Only `(x,y)` is reversed.  The arcs are `y -> x`, `y -> z`, and `x -> z`,
   so `y` is a source; the triple is transitive.
3. Only `(y,z)` is reversed.  The arcs are `x -> y`, `z -> y`, and `x -> z`,
   so `y` is a sink; the triple is transitive.
4. Only `(x,z)` is reversed.  The arcs are

   ```text
   x -> y,   y -> z,   z -> x,
   ```

   which is a directed 3-cycle.

Thus a sorted triple is cyclic exactly when its long edge `(x,z)` is one of the
reversed matching arcs.  For each reversed arc `(a,b)`, the cyclic triples using
that long edge are in bijection with the vertices strictly between `a` and `b`.
Summing over `(a,b) in M` gives the stated formula.  The matching hypothesis
ensures these triples are disjointly assigned to their unique reversed long
edge.  QED.

## Staircase Corollary

For the all-0 interleaved staircase on `2k` vertices, use the linear order

```text
d_0 < d_1 < ... < d_{k-1} < r_0 < r_1 < ... < r_{k-1}
```

and reverse the matching arcs `(d_p,r_p)`.

The interval between `d_p` and `r_p` contains

```text
(k - 1 - p) later dominants + p earlier recessives = k - 1
```

vertices.  Therefore

```text
# directed 3-cycles = sum_{p=0}^{k-1} (k - 1) = k(k - 1),
```

recovering THM-399.

Equivalently, for each `p < q` the two cycles are

```text
r_p -> d_p -> d_q -> r_p
d_q -> r_p -> r_q -> d_q.
```

These are the two cycle families formalized in Lean in
`eliott-monad/math-lean` as `Math.Tournaments.Staircase`.

## Formalization Status

The staircase specialization `#C_3 = k(k-1)` is formalized sorry-free in
`eliott-monad/math-lean`, commit `b5ffcde`, theorem
`Math.Tournaments.allZeroStaircase_threeCycleCount`.

The general interval-reversal theorem above is not yet formalized in Lean; it is
the natural next formalization target if the staircase infrastructure is
generalized from the two-block vertex model to an arbitrary finite linear order.
