---
id: HYP-2679
title: LRC(14) true-wide boundary-curvature ledger
status: OPEN; exact two-far atlas supports d=1 finite-curvature route
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

## Exact Test

Script:

- `04-computation/lrc14_truewide_boundary_curvature_codex_s50.py`
- output: `05-knowledge/results/lrc14_truewide_boundary_curvature_codex_s50.out`

The first computation scans exact `k=9` two-far rows

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

## Exact Findings

The scan over `core=(0)+6-subsets of [1,14]` and far pairs `15..24` checked
`135135` rows, `135065` primitive.  Curvature signs:

```text
positive: 131003
negative:   3961
zero:        101
```

The direct-risk leader is still the HYP-2675 true-wide leader:

```text
E=(0,4,6,8,10,12,14,15,16)
p0=321/980
cap-p0=11681/70070
core=(0,4,6,8,10,12,14)=2*(0,2,3,4,5,6,7)
I_B(15,16)=-13/1470
```

This is the important correction: the top k=9 true-wide row is not a positive
two-far synergy.  It is a dilated bounded-core row with overlapping one-far
mass.  The large increment is the `16` peel:

```text
p0(B)=67/1470
p0(B union {15})=13/120       increment=123/1960
p0(B union {16})=1609/5880    increment=447/1960
p0(B union {15,16})=321/980
```

The largest positive-curvature row in the scan is

```text
E=(0,1,4,8,10,12,14,16,20)
p0=89/336
cap-p0=11021/48048
I_B(16,20)=307/1960
```

It is safely below the top-risk rows.  Thus positive curvature is a genuine
object, but it is not yet the cap-threatening one.

Bucket maxima by additive excess reinforce the split:

```text
small excess:  count=4188,   max p0=321/980 at the HYP-2675 leader
medium excess: count=130399, max p0=361/1260
large excess:  count=478,    max p0=159017/966966
```

Core-gcd maxima are also revealing:

```text
gcd(core)=1: count=134820, max p0=1729558/5360355
gcd(core)=2: count=245,    max p0=321/980
```

The unique direct-risk champion in this box lives in the tiny `gcd(core)=2`
dilated-core bank, exactly the HYP-2678 `d=1`/scale-invariance route.

The KPS third-pocket path

```text
(0,3,5,16,28,30,33,35)
```

behaves differently: over the tiny core `(0,3,5)`, every two-far pair curvature
among the far speeds is `0`; the sorted path has `p0=0` until the `33` and `35`
steps.  That row is therefore a multi-far threshold / packet-cancellation
object, not a two-far curvature exception.

## Tournament Analysis

Row-risk tournament:

```text
vertices: true-wide two-far rows / proof obligations
pairwise observable: larger exact p0/cap wins
tie path: curvature, then additive excess, then lexicographic row
score_hist={0:1,1:1,...,11:1}
directed_3cycles=0
Hamiltonian path leader=(0,4,6,8,10,12,14,15,16)
```

Far-speed first-increment tournament on vertices `15..24`:

```text
pairwise observable: aggregate_B [(p0(B+u)-p0(B))-(p0(B+v)-p0(B))]
Hamiltonian path=[15,18,16,20,17,24,21,22,19,23]
directed_3cycles=0
ties=0
```

This far-speed tournament is not monotone in numeric speed; it is a residue and
core-interaction order.  It is a useful finite gauge for the HYP-2678 peel
ledger.

## Current Status

No proof is claimed.  The exact atlas suggests the next sharp theorem is:

```text
d=1 / dilated-core true-wide rows:
  reduce by scale-invariance to a finite curvature ledger.

d>=2 / high-growth rows:
  positive curvature can occur, but observed direct risk is lower;
  prove a signed dimension-penalty or excess/curvature cap.
```

So HYP-2679 narrows HYP-2678 rather than broadening it.  The live k=9 leader is
a finite dilated-core overlap row, not a high-dimensional positive-curvature
row.

## Assumption Challenge

Alternate vertex sets considered: runners, far speeds, far-speed pairs, rows,
approach orders, missed-sector state words, packet signs, K4 edge pairs,
relation-lattice supports, and proof obligations.  The chosen quotient preserves
the direct LRC predicate `p0(E)<=cap_k` and the order-sensitive two-far
interaction.  It destroys detailed endpoint geometry inside each wall atom, so
any final theorem must reattach HYP-2648 state words and HYP-2639
relation-shell visibility.
