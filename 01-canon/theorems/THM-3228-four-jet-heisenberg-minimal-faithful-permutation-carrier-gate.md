---
id: THM-3228
title: "Four-jet Heisenberg minimal faithful permutation carrier gate"
status: >
  PROVED + INDEPENDENTLY HOSTILE-AUDITED.  For every prime p, any
  permutation action of the mod-p tangent four-jet group on fewer than p^2
  points kills its central commutator.  The bound is sharp: a noncentral
  order-p coset action, equivalently the standard affine action on F_p^2,
  is faithful on p^2 points.  Thus a p-root permutation chart cannot retain
  the full four-jet curvature, while p=13 points toward a cardinality-minimal
  169-point carrier.  No lawful LRC carrier map is constructed.
source: root/2026-08-02
audit: >
  The proof composes THM-3220's exact four-jet/Heisenberg identification with
  THM-2779's sharp permutation-degree theorem and restates the orbit-size
  argument in center-faithful form.  An independent hostile audit rechecked
  every group/action sign, the orbit and core arguments, the D8 boundary, and
  the cardinality-only LRC scope.
depends_on:
  - THM-2779-bockstein-symplectic-decoder-frame-torsor-and-heisenberg-root-degree-gate
  - THM-3220-root-four-jet-schwarzian-heisenberg-transgression-and-oriented-discriminant-holonomy
related:
  - THM-2309-owner-aligned-pivot-packets-and-visible-height-separation
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2419-valuation-normalized-homogenization-of-affine-sideband-shells
---

# THM-3228 -- four-jet Heisenberg minimal faithful permutation carrier gate

**PROVED + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3220 identifies the selected-root tangent four-jet group over an odd
residue field with the finite Heisenberg group, and identifies its raw
characteristic-two boundary with `D_8=H_2`.  THM-2779 proves that the minimum
degree of a faithful permutation action of `H_p` is `p^2`.  The useful joint
conclusion is sharper than a cardinality slogan: every smaller permutation
model kills the central commutator that carries the four-jet area (and its
kernel may be larger).

## 1. The four-jet group

Let `J_4(p)` be the group of tangent-to-identity coefficient jets

```text
u+alpha u^2+beta u^3+gamma u^4              (mod u^5)    (1)
```

over `F_p`, under functional composition.  For odd `p`, THM-3220's
logarithmic coordinates and exact dictionary give

```text
J_4(p) isomorphic H_p,
H_p={(x,y,z):x,y,z in F_p},
(x,y,z)(x',y',z')=(x+x',y+y',z+z'-yx').                  (2)
```

For `p=2`, the raw coefficient law is `D_8`, with quotient square form
`q(x,y)=xy`; this is exactly THM-2779's `H_2` boundary.  Uniformly,

```text
Z(J_4(p))=[J_4(p),J_4(p)] isomorphic F_p.                (3)
```

Under THM-3220, this center is the fourth-jet half-symplectic fill and the
oriented discriminant/commutator line.  A carrier which makes that central
curvature observable must therefore act nontrivially on `(3)`.

## 2. Every smaller permutation action kills the center

Let `J_4(p)` act on a finite set `X` with

```text
|X|<p^2.                                                  (4)
```

Because `J_4(p)` is a `p`-group, every orbit has `p`-power size.  Under
`(4)` the only possible orbit sizes are `1` and `p`.  A fixed point plainly
kills the center on that orbit.  On an orbit of size `p`, a point stabilizer
has index `p`.  Every index-`p` subgroup of a finite `p`-group is normal, and
the quotient has order `p`, hence is abelian.  Therefore the stabilizer
contains the commutator subgroup `(3)`, so the center again fixes the entire
orbit.

Since this holds on every orbit,

```text
Z(J_4(p)) <= kernel(J_4(p) acting on X).                  (5)
```

Thus degree below `p^2` is impossible not only for a faithful action but for
any **center-faithful** permutation action.  In particular, a permutation of
only the `p` selected roots necessarily erases the fourth-jet curvature.

## 3. Sharpness at p squared

Choose a noncentral subgroup `K<=H_p` of order `p`.  Then `K intersect Z=1`.
Its normal core is trivial: any nontrivial normal subgroup of a finite
`p`-group meets the center nontrivially, whereas the core lies inside `K`.
Consequently the coset action on

```text
H_p/K,                         |H_p/K|=p^2,               (6)
```

is faithful.  Equivalently, THM-2779 gives the explicit action on `F_p^2`

```text
(x,y,z):(r,w) |-> (r+x,w+z-yr).                          (7)
```

The stabilizer calculation in THM-2779 proves `(7)` faithful, including the
`D_8` case `p=2`.  Combining `(5)--(7)` gives the exact equality

```text
minimum center-faithful permutation degree of J_4(p)
 =minimum faithful permutation degree of J_4(p)=p^2.      (8)
```

The correct-size hostile is already instructive.  The abelian prefix quotient

```text
J_4(p)/Z isomorphic F_p^2                                  (8a)
```

also has `p^2` elements, but pulling back its regular action kills `Z` by
construction.  The faithful minimal sets are instead the coset sets `J_4(p)/K`
with `K` a noncentral order-`p` subgroup.  Such `K` is necessarily nonnormal.
THM-2779 classifies `p+1` minimal transitive classes for odd `p` and the two
reflection-stabilizer classes for `D_8`.  Thus correct cardinality plus the
obvious `(A,B)` prefix labels is still insufficient: a minimal carrier needs
an oblique nonnormal stabilizer-line clutch.

## 4. The 13-root/169-target gate

At the LRC prime `p=13`, equation `(8)` reads

```text
13 roots:   too small; every permutation model kills the center,
169 points: cardinality-minimal; an abstract faithful model exists.       (9)
```

This makes the two-coordinate `F_13^2` target/address grid a structurally
correct **candidate size** for a four-jet carrier.  It also explains why a
root-only selector can preserve the affine two-jet reflection while losing
the higher central phase: the latter cannot fit faithfully in the same
permutation chart.

The connection contract is

```text
source:      the full mod-p selected-root four-jet group;
map:         a proposed permutation action on a finite physical carrier;
preserved:   the central commutator / oriented-discriminant curvature;
obstruction: every carrier of size <p^2 kills that center;
sharp model: abstract cosets H_p/K, equivalently F_p^2;
needed:      a lawful same-ancestry observable realizing that action.      (10)
```

Equation `(9)` does **not** identify THM-3220's root-jet center with an LRC
target coefficient, owner current, endpoint phase, or relation address.
Cardinality and abstract group isomorphism do not provide such a map.

## 5. Sharp boundaries

1. The result concerns permutation carriers.  THM-2779 also has a
   `p`-dimensional complex coefficient representation in which the center
   acts by phase; `(8)` does not rule out phase/Fourier carriers.
2. Faithfulness of the whole group is stronger than necessary.  The proof
   actually establishes the center-faithful lower bound `(5)`, which is the
   load-bearing statement for curvature.
3. A disjoint collection with at least `p^2` labels is not automatically a
   carrier.  Even the regular quotient action on `(8a)` kills the center.
   The nonnormal stabilizer clutch, common ancestry, and physical typing must
   still be supplied.
4. Multiple-root and root-at-infinity charts remain outside THM-3220 and
   therefore outside this corollary.

An independent hostile audit rechecked THM-2779's multiplication and action
signs, the `1/p` orbit dichotomy below `p^2`, normality and commutator
containment of every index-`p` stabilizer, the trivial-core coset action, and
the reflection-stabilizer `D_8` endpoint.  It also confirmed that `(8)` is a
center-faithful statement and that `(9)` asserts only a minimal cardinality,
not a physical LRC realization.

QED.
