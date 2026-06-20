---
id: HYP-2679
title: LRC(14) true-wide boundary-curvature ledger
status: OPEN; namespace claimed, exact scan in progress
source: codex-2026-06-20-S50
depends_on:
  - HYP-2678
  - HYP-2677
  - HYP-2676
  - HYP-2675
  - THM-547
  - THM-546
  - THM-531
  - HYP-2648
  - HYP-2639
related:
  - HYP-2671
  - HYP-2672
  - HYP-2657
  - HYP-2637
  - HYP-2606
  - OPEN-Q-108
---

# HYP-2679 - True-Wide Boundary Curvature

## Claim Being Tested

The remaining `LRC(14)` sector crux after THM-547 is the true-wide case
`second-largest(E)>14`, where at least two far speeds interact.  The
boundary-function analogy suggests that this is not a one-far boundary limit
problem any more: different far-speed approach arcs can have different first
derivatives even though the final direct quantity `p0(E)` is order-independent.

For a core `B` and two far speeds `u<v`, define the exact direct-p0 curvature

```text
I_B(u,v) =
  p0(B union {u,v}) - p0(B union {u}) - p0(B union {v}) + p0(B).
```

Equivalently, `I_B(u,v)` is the change in the `v`-increment after `u` has
already been added.  HYP-2679 tests whether the true-wide leaders are explained
by a finite low-rank list where this two-far curvature is large and positive,
while high-growth or higher-dimensional Freiman models have small, negative, or
packet-cancelling curvature.

This is meant to merge with HYP-2678, not replace it:

```text
d=1 GAP / dilated AP core
  -> scale-invariance (THM-531) + finite base + O(1) curvature ledger

d>=2 GAP / higher relation lattice rank
  -> signed dimension penalty + small two-far curvature
```

## Boundary-Function Analogy

External compass checked during the session:

- Kaczynski's "Boundary functions and sets of curvilinear convergence for
  continuous functions" records that curved approach sets and boundary
  functions are separate data.
- Fatou-type theorems for bounded harmonic functions distinguish controlled
  nontangential approach regions from wider/tangential approach behavior.
- Continuous self-maps of the Riemann sphere carry a topological degree; the
  repo's Cayley-monad/Riemann-sphere material is the local analogue where
  boundary winding is an invariant of rational transform data.

The LRC translation is intentionally modest: the core state word is the
boundary datum, and far-speed additions are approach arcs.  The invariant to
measure is not only the endpoint `p0(B union F)`, but the exact curvature of
the approach path in the finite wall arrangement.

## Planned Exact Test

The first computation should scan exact `k=9` two-far rows

```text
E = B union {u,v},  B subset {0,...,14}, 14 < u < v,
```

and store, in exact Fractions:

- `p0(B)`, one-far values, two-far value, LRC cap margin;
- `I_B(u,v)` and the two ordered increments;
- additive/Freiman fingerprints (`K2`, excess, energy, longest AP/run);
- squarefree divisor profile and state-word support/entropy;
- a Tournament Analysis whose vertices are rows or far pairs, not runners.

The pairwise observable for the first tournament is direct risk `p0(E)/cap_k`,
with tie path lexicographic by `(curvature, additive excess, row)`.  A second
far-pair tournament can orient `{u,v}` by the larger first increment
`p0(B union {u})-p0(B)` versus `p0(B union {v})-p0(B)`.

## Current Status

This file claims `HYP-2679/T918` for the true-wide boundary-curvature target.
No proof is claimed yet.  The immediate deliverable is an exact atlas that
decides whether the HYP-2675 true-wide leader
`(0,4,6,8,10,12,14,15,16)` is a finite curvature exception over the dilated
bounded core `(0,4,6,8,10,12,14)`, or whether the same curvature pattern
persists across enough higher-dimensional rows that a new analytic bound is
needed.

## Assumption Challenge

Alternate vertex sets considered: runners, far speeds, far-speed pairs, rows,
approach orders, missed-sector state words, packet signs, K4 edge pairs,
relation-lattice supports, and proof obligations.  The chosen quotient preserves
the direct LRC predicate `p0(E)<=cap_k` and the order-sensitive two-far
interaction.  It destroys detailed endpoint geometry inside each wall atom, so
any final theorem must reattach HYP-2648 state words and HYP-2639
relation-shell visibility.
