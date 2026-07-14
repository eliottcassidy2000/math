---
id: THM-773
title: Prime-seven sheet monodromy and the exact tournament-atlas fibre
status: CLAIMED (number reserved; elementary proof written in outline and exact atlas/monodromy audit in progress)
source: codex-2026-07-14-S6
depends_on:
  - THM-771   # exact seven-owner sheet defect and endpoint convention
  - HYP-6825  # canonical merged-metagraph addresses and inverse tiling fibres
related: [THM-754, HYP-3802, HYP-6835]
---

# THM-773 — Prime-seven sheet monodromy and tournament fibre

This stub reserves the theorem number for the exact `c=7`, seven-owner bridge.
Let `W={w_0,...,w_6}` with `7` dividing no `w_a`, and work away from the
endpoint set `w_a x in Z+1/2`.  Each owner has one strict bad sheet, namely

```text
k_a(x) = -w_a^(-1) round(w_a x)  (mod 7).
```

The sheet cover is exact precisely when `a -> k_a(x)` is a bijection of
`F_7`.  Equivalently,

```text
product_a (X-k_a)=X^7-X  in F_7[X],
```

or its six nontrivial finite-field power moments have the full-grid values.
At an event of owner `a`, its token changes by `-w_a^(-1)`.  Hence return
between exact-cover chambers requires the first-moment holonomy

```text
sum_(events a) w_a^(-1) = 0  (mod 7),
```

with the higher moment changes giving the complete prime-sheet obstruction.

For an exact owner assignment, two natural tournament gauges deliberately
forget different data:

1. orient owners by the marked linear sheet order; every assignment becomes
   the transitive class `n7-a000`;
2. orient `a -> b` when `k_b-k_a in {1,2,3}`; every assignment becomes the
   rotational heptagon class `n7-a267`.

The gauges differ on exactly the six long chords.  The exact HYP-6825 atlas
places `n7-a267` at local depth six with address word `SSRRSS`, blue/black
root word `BBBBBB`, and a 25-mask inverse staircase fibre.  The intended
finite-exact result is that lexicographically choosing a Hamiltonian path in
the labelled circular tournament maps all `7!=5040` owner assignments onto
exactly those 25 masks.  Thus the ordinary isomorphism-class node is constant;
the owner assignment, event steps, and endpoint schedule are a genuine stalk,
not optional decoration.

The companion computation will record the full map, fibre multiplicities,
the two example chamber movies from THM-771, and Tournament Analysis.  Until
that audit and the endpoint-side proof are committed, cite this file only as
a claimed theorem.
