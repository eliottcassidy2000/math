---
id: THM-1129
title: The sharp bounded-offset polygon floor and the K at least 832 four-comb needle tail
status: PROVED / COMPUTER-ASSISTED — exact finite-candidate theorem on all 495 eight-speed cores and all 4,060 normalized offset shapes of span at most 30; closes their sharp four-comb target uniformly for K at least 832
source: codex-2026-07-18-S73 r5 covariogram/toothpick continuation
depends_on: [THM-1097, THM-1123, THM-1126, THM-1127, THM-1128]
related: [THM-1101, MISTAKE-164]
verification:
  - 04-computation/lrc14_r5_offset_polygon_floor_exact_codex_S73.py
  - 05-knowledge/results/lrc14_r5_offset_polygon_floor_exact_codex_S73.out
---

# THM-1129 — bounded-offset polygon floor and needle tail

Let `P` be an eight-element subset of `{1,...,12}` and let

```text
A={0<a<b<c},          c<=30.
```

For `t in R/Z`, place four offset centers at `-alpha*t`, `alpha in A`, and
write

```text
w_A(t)=max circular gap between the four centers - 1/7.       (1)
```

Thus `w_A(t)` is the longest vertical safe interval in the fixed torus
polygon

```text
H_A={(t,x): ||x+alpha*t||>=1/14 for every alpha in A}.         (2)
```

Let `G_P^+` be the closure of the union of the positive-length components of
the core-safe set

```text
G_P={t: ||pt||>=1/14 for every p in P}.                        (3)
```

The superscript excludes an isolated safe point, which cannot support a
one-dimensional component.

## The exact polygon theorem

For all

```text
C(12,8)*C(30,3)=495*4060=2,009,700                       (4)
```

core/shape pairs,

```text
M(P,A):=max_{t in G_P^+} w_A(t) >= 2/5.                  (5)
```

The bound is sharp in this atlas.  There are exactly `32` attaining pairs,
with offset-shape histogram

```text
(0,1,8,10): 1,   (0,2,3,7): 15,
(0,2,9,10): 1,  (0,4,5,7): 15.                         (6)
```

In particular, the hard core from THM-1123,

```text
P={1,2,4,5,7,9,11,12},
```

attains equality for both `(0,2,3,7)` and `(0,4,5,7)`.  Their reflected
maximizers are `t=27/70,43/70`.

## The finite-candidate lemma

The computation in (4) is exhaustive because it rests on a short exact
lemma, rather than a time mesh.  Split `[0,1]` at

```text
B = {(14j+-1)/(14p): 1<=p<=12}
    union {j/(beta-alpha): alpha<beta in A}.             (7)
```

The first set contains every possible core-tooth endpoint.  The second
contains every collision of two offset centers.  On an open cell of (7),
the core predicate is constant, the cyclic center order is fixed, and each
labelled circular gap is affine in `t`.  Hence `w_A`, the maximum of those
affine functions minus `1/7`, is convex on the cell.  Its maximum over the
closure of a core-safe cell occurs at an endpoint.  Testing the two endpoints
of every live cell therefore proves (5) exactly.

The referee uses `fractions.Fraction` for all `107,434,406` eligible
core/order-cell incidences.  No floating-point number enters a certificate
decision.

## A common rectangle and the all-K tail

The same atlas records more than the vertical maximum.  For every pair
`(P,A)`, it selects a one-sided core-safe order cell adjacent to a point `t0`
where

```text
w_A(t0)>=2/5.                                           (8)
```

Choose compatible lifts of the two walls bounding this gap.  Every wall has
slope in `{-alpha:alpha in A}`, hence in `[-c,0]`.  Over a one-sided interval
of time length `d`, the intersection of the lifted vertical gaps loses at
most `c*d`: on the right side only the upper extreme can retreat, and on the
left side only the lower extreme can advance.  Consequently the exact choice

```text
d=min(cell length,(w_A(t0)-1/5)/c)                      (9)
```

produces a rectangle `I times X subset G_P times H_A` with `|X|=1/5`.
Maximising (9) over the eligible endpoint/safe-side certificates gives the
uniform exact bound

```text
|I|>=1/728.                                             (10)
```

There are `73` pairs attaining the certificate floor.  A canonical one is

```text
P={1,3,4,5,8,9,11,12}, A=(0,3,8,13),
I=[3/13,13/56], w_A(3/13)=50/91,
|I|=1/728.                                              (11)
```

Now put `k_i=K+alpha_i`, so `k4=K+c`, and sample (2) along the Kakeya needle

```text
x=K*t mod 1.                                            (12)
```

The arc `X` contains a subarc `Y` of length `1/7`.  Its preimages under
`t -> Kt mod 1` are intervals of length `1/(7K)`, with left endpoints spaced
by `1/K`.  An interval of length `1/728` contains a complete such preimage
whenever

```text
1/728 - 1/(7K) >= 1/K,
```

which is equivalent to

```text
K>=832.                                                 (13)
```

That preimage lies in `G_P` and in the four-comb polygon, so the final safe
set has a component `J` with

```text
|J|=1/(7K)>1/(7(K+c))=1/(7k4).                         (14)
```

Thus the sharp four-comb target required before a fifth killer holds for
every eight-speed core, every normalized offset shape of span at most `30`,
and every `K>=832`.

## The exact finite complement

The uniform threshold (13) is a worst-case bound, not the finite workload.
For each of the `2,009,700` pairs, the referee also maximises the one-sided
rectangle supporting a `1/7` arc and records the resulting individual tail
start.  The exact ledger is

```text
pair rays whose tail begins by their legal first K:       1,802,872
raw legal K rows below the individual tails:               6,040,056
rows already covered by THM-1123's first-forty bank:        2,500,120
remaining bounded-offset finite rows:                       3,539,936. (15)
```

This is a finite target of the correct size, but it has not been run here.
Equation (15) is a workload reduction, not a silent promotion of those rows.

## Relation to the universal thirteen-grid cone

THM-1128 applies to arbitrary offset span.  In its centred coordinates it
proves the all-shape cone

```text
B>=53*max(A,24),
```

where `A` is the maximum centred offset.  Its midpoint corollary is

```text
K>=max(1272,26(c+1)).                                   (16)
```

THM-1129 does not replace this theorem: it restricts to `c<=30` and spends an
exact finite atlas.  In return, it improves the bounded-shape onset to the
constant `K>=832` and reduces the honest finite complement to (15).  The two
theorems cover different compactifications of the same torus object:
THM-1128 uses one universal centred thirteen-grid rectangle, while THM-1129
optimises all labelled bounded-shape cells.

## Tournament and carrier audit

Runner, comb, core-cell, section-boundary, wall-event, residue, and
proof-obligation vertices were all challenged.  Exact candidate-boundary
order gives a transitive tournament: after ties are coalesced there are no
directed cycles, every SCC is a singleton, and the sorted order is the unique
Hamiltonian path.  This quotient is not faithful.  It destroys the cyclic
gap lengths, adjacent owners, wall slopes, safe side, and rectangle width.

The proof-bearing carrier is instead

```text
(labelled cyclic gap, endpoint coordinate, safe-side cell, owner slopes). (17)
```

It preserves the rectangle and needle-crossing predicate.  The naked order
tournament does not.  This is another instance where the useful Kakeya
object is a labelled wall section rather than a tournament on runners.

## Verification and scope

Normal and optimized Python executions are byte-identical to the frozen
output.  The SHA-256 hashes are

```text
source d55b42f89a6cf06040954a25f61277d334e663fda2909630896380b0253c8518
output d2a84a3097b2da1e27eb709ade85edb9c4d6a5679a535efd304ff104af2816a9
```

THM-1129 closes only the tail (13) for `c<=30`; the `3,539,936` finite rows in
(15), all shapes with `c>30` outside THM-1128's cone, and uniform `r=5` remain
open.
