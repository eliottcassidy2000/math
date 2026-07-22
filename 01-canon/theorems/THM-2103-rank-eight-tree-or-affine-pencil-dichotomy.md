---
id: THM-2103
title: "The rank-eight mixed-torus tree-or-affine-pencil dichotomy"
status: >
  RESERVED / CONJECTURAL TARGET. For eight terminal characters transverse to
  one guard, the proposed dichotomy is: either the restricted-overlap maximum
  spanning tree has weight greater than 5/49, or some choice of terminal
  signs places the eight characters in a nondegenerate affine pencil. The
  first branch contradicts THM-2098 under cover; the second escapes by
  THM-2099. Exact bounded exhaustion and a structural proof are both pending.
  No rank-eight closure or LRC(14) conclusion is claimed.
source: codex-2026-07-22-LRC-rank-eight-tree-pencil
depends_on:
  - THM-2098
  - THM-2099
related:
  - THM-1221
  - THM-2096
---

# THM-2103 -- reserved tree-or-pencil dichotomy

## Exact target

Let `g,c_1,...,c_8 in Z^2`, with every `c_i` transverse to `g`, and suppose
some integer direction gives an odd positive guard and distinct positive
terminal specializations. With

```text
w_ij=measure({||c_i.X||<1/14,
              ||c_j.X||<1/14,
              ||g.X||>1/7}),
tau=max_(spanning trees T) sum_(ij in T) w_ij,
```

prove or refute

```text
tau>5/49,
or there are signs sigma_i and independent p,q with
sigma_i c_i-q in Q p for every i.                      (T)
```

THM-2098 makes the first alternative incompatible with a mixed cover.
THM-2099 makes the second alternative safely escapable.

## Current evidence and missing step

THM-2099's infinite dyadic affine pencils occupy the second branch even when
their tree weight is below or equal to `5/49`. Triadic, consecutive,
Fibonacci, alternating, all single-site dyadic perturbations, and the seeded
bounded controls tested there occupy the first branch. This is evidence only.

The missing theorem is an inverse statement for the exact primitive-relation
weights: a nonpositive tree surplus must force signed affine rank one. A
finite coefficient scan cannot replace its height tail.

## Assumption challenge and Tournament Analysis

The challenged assumption is that “non-pencil” should mean the unsigned
columns are not collinear. Character danger bands are sign-invariant, so all
`2^7` relative sign gauges must be checked. The faithful carrier retains the
edge weights, primitive relation labels, graphic-matroid maximum, and signed
affine rank. A tournament orientation by weight is only a search scheduler;
ties and label choices destroy the proposed inverse theorem.
