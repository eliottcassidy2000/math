---
id: THM-968
title: T5 plane decomposition — exact affine-plane fibers, open uniform slab estimate
status: OPEN SCHEME + FINITE FLOATING-POINT PROBE — the slab fibers and coordinate-line loci are exact, but T5(H) = O(L^4/H) is unproved; the probe uses BOX=60 and capped K0=3000 with no tail estimate, and every advertised dissociated quintuple has a short exact relation
source: kind-pasteur-2026-07-17-S128 cont.42; rigor correction after upstream audit
depends_on:
  - THM-946 (corrected two-pole atom and the slab problem statement)
  - THM-952 (proved support-three near-pole treatment)
  - THM-965 (exact T4 line fibers, but not a proved T4 tail bound)
related:
  - THM-935/948 (the E_s frame and exact packet audits)
script: 04-computation/t5_slab_referee_kps_S128c42.py -> 05-knowledge/results/t5_slab_referee_kps_S128c42.out
---

# THM-968 — the T5 plane scheme

## Exact structural reduction

For an outer triple `(u,t,r)`, put

```text
k = w3*u + w4*t + w5*r.
```

For fixed `k`, the integer solutions are either empty or an affine coset of
the rank-two kernel lattice.  On such a plane the coordinate restrictions
`u`, `t`, and `r` are affine integer functionals.  Each zero locus is therefore
empty, the whole plane in a degenerate case, or an affine lattice line.  These
facts rigorously identify the plane and line carriers of the slab.

## What is not proved

The proposed estimate

```text
T5(H) <= C5 * L^4 / H
```

requires substantially more than this carrier decomposition.  In particular,
one needs uniform covolume and successive-minimum control for shifted planes,
uniform near-line collar estimates, treatment of coordinate-zero and gcd
degeneracies, an honest structured-relation branch, and an infinite-tail
bound.  THM-965 does not supply the required rank-one price: its T4 estimate
is itself open.  Consequently the recursive “price each plane by T4” step is
circular as a proof of the T5 bound.

The plane/collar/bulk picture is therefore an **open scheme**, not a completed
coarse-rigor pass.  At present the proved analytic ladder stops with the T3
case in THM-952.

## Exact scope of the finite probe

`t5_slab_referee_kps_S128c42.py` is a finite floating-point experiment:

- all five coefficients are restricted to `[-60,60]` (`BOX=60`), with no
  estimate for the omitted tail;
- `K0` is replaced by `min(K0,3000)` for all three displayed tuples, so the
  measured slab is not the full proposed resonance slab;
- the Fourier weights and the reported near-line percentages use ordinary
  floating-point arithmetic without certified error bounds;
- for `H=40`, only the thin shell `40 < max|h_i| <= 60` is measured.

All three tuples labelled “dissociated quintuples” have exact small relations:

```text
 307 - 2*425 - 541 + 4*671 - 2*800 = 0,
2*541 - 3*800 + 3*1087 - 1943       = 0,
2*671 + 4*2147 - 3*3310             = 0.
```

(The omitted coordinates in the last two equations have coefficient zero.)
Their coefficient heights are at most four, so none of the three tuples is
`H`-dissociated for either tested horizon `H=10,40`.  The observed bounded
truncated envelopes and `37–60%` near-line concentration are useful probes of
the finite box, but they cannot justify a dissociation floor or a uniform tail
theorem.

## Remaining proof obligations

- prove a uniform affine-plane lattice count with explicit covolume and shift
  dependence;
- control all three near-line collars and their intersections without assuming
  the desired T4 estimate;
- dispatch short relations to the structured branch before using `1/H`;
- remove both the `K0` and `BOX` caps and bound the infinite complement.
