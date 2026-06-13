---
source: codex-2026-06-03-S600
status: synthesis refinement of HYP-2089 after the cap-face and determinant branches
tags: [LRC, strong-tournament, round-tournament, SCC, n14, HYP-2089, HYP-2140, HYP-2144]
---

# The LRC-relevant subset of strong tournaments

The useful answer is not "strong tournaments matter."  That set is too large.
The useful subset is:

```text
LRC-accessible round strong runner blocks
  that survive the easy source/semicircle route
  and whose proof-obligation quotients have no peel order.
```

This is the strong-tournament version of the current n=14 proof split.

## Three nested sets

Start with all tournament isomorphism classes on `m=n-1` runners.  LRC does not
see most of them.  The half-turn phase comparator only produces round
tournaments: the out-neighborhood of each runner is a contiguous circular arc.
So Paley-like or generic strong tournaments may be beautiful, but they are not
the geometric image of an LRC time cell.

Inside round tournaments, THM-374 gives the first cut:

```text
points in an open semicircle  <=>  transitive tournament,
points not in an open semicircle  =>  strong encirclement in the generic source cases.
```

At a lonely time, S566 found the runner block in n=14 samples has only the two
extremes:

```text
transitive semicircle block, or one strong block;
no intermediate SCC decompositions in the 152-row probe.
```

Thus the LRC source event is either a clustering event or an encirclement event.
The transitive side is an easy source: all runners fit on one side.  The strong
side is the generic source: the runners Hamiltonian-encircle the observer, and
the observer escapes through a gap.

## Why only a subset of strong blocks is hard

Moon/Camion says a strong tournament has a Hamiltonian cycle.  In the LRC
half-turn setting that cycle is not merely graph-theoretic: it is a circular
encirclement of the observer.  But most such encirclements are uneven.  Uneven
means at least one of the following proof handles tends to appear:

```text
large observer gap,
semicircle-near transitive cut,
load imbalance in n-clock cells,
endpoint-owner private pivot,
small determinant wall among component languages.
```

The hard subset is therefore not "strong" but "regular strong" or
"certificate-resistant strong."  For `m=13` this is the rotational round
tournament `R_13`: every runner has out-degree `6`, the circle is balanced, the
class is self-converse, and the AP/V* boundary examples realize it.  This is why
the same object appears under many names:

```text
regular rotational encirclement,
round chi=2 tight orbit,
self-converse Burnside fixed side,
perfect antipodal transversal / prime-3 shell residual,
measure-zero floor-tight source event.
```

The subset matters because it filters away almost all strong tournament
complexity and leaves a single proof target:

```text
a regular round strong encirclement cannot perfectly cover the observer's
1/n-neighborhood.
```

That is the strong-language form of "the danger arcs never completely cover."

## Relation to the cap and determinant branches

HYP-2140 and HYP-2144 sharpen what "certificate-resistant" means.

The dual cap face chooses primary `n`-clock cells as tournament vertices.  If
one cell has too much safe mass, the cell-load tournament has a source and the
row is loose by capacity.  The origin-bisection refinement then tries upper and
lower half-cells.  These are mostly transitive proof tournaments: one
certificate dominates and gives a peel.

If all cap cells are under capacity, the row enters the large-owner determinant
layer.  There the natural vertices are component languages in the bounded
multiplier `w`.  A singleton or pair Helly wall is again a peel/source
certificate.  A genuinely hard residual would be a strong proof-obligation
tournament: every component language helps block every other, and no small
subfamily exposes emptiness.

So there are two strong notions:

```text
geometric strong: runner block encircles the observer;
proof strong: certificate obligations form a non-peelable SCC.
```

The LRC n=14 route wants to prove that the intersection is empty except for the
known floor-tight boundary, where the gap is exactly `1/n` and not a
counterexample.

## Assumption challenge

Candidate tournament vertices considered:

```text
runners,
gaps,
cap centers,
n-clock cells,
safe components,
endpoint owners,
component languages,
residue channels,
proof obligations.
```

The right vertex set depends on the question.  Runner vertices preserve the
geometric predicate "does the observer escape a strong encirclement?"  Clock
cells preserve the cap predicate "does capacity cover this safe mass?"
Component languages preserve the determinant predicate "does any bounded
multiplier cover all components?"

The challenged assumption is that a single raw runner tournament should solve
the proof.  Strong runner blocks identify the geometric worry set, while strong
proof-obligation blocks identify the algebraic residual after cap and Helly
certificates have failed.

## Practical proof use

For n=14, strong tournaments matter as a sieve, not as an endpoint:

1. discard non-round strong tournaments because LRC cannot realize them;
2. discard transitive/semicircle source events because they are easy;
3. discharge irregular round strong events by gap, cap, owner, or determinant
   certificates;
4. reduce the hard face to the regular rotational strong block and prove its
   danger cover has a gap.

That is why this subset matters: it is the place where all cheap asymmetry has
been removed.  If the proof works there, the surrounding irregular strong cases
should fall by perturbation, capacity, or private-pivot descent.
