---
id: THM-3695
title: "Danielewski embedding and seven-piece floor for the y=0 collision ring"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The y=0 collision
  ring of THM-3686 is a graded subalgebra of the squarefree exponent-two
  Danielewski ring with Sigma(b)=1-b^2, whose embedding into the source is
  Poisson.  Consequently the universal nonentry theorems
  THM-3569, THM-3583, and THM-3592 apply without relaxation: after scalar
  removal every Darboux pair in the collision ring has at least seven active
  homogeneous pieces.  At total width seven only the partitions 2x5 and 3x4,
  up to output exchange, are not excluded by this corollary.  Existence in
  those cells, a Darboux pair, and a counterexample to JC(2) remain OPEN.
source: root / jc-quartic-c3-construction, 2026-08-22
audit: >
  PASS after terminology repair.  The independent audit checked injectivity,
  all Poisson signs, grading-support preservation, the exact scopes of the
  three inherited theorems, and every support partition through total seven.
  It caught and removed the false preliminary phrase "Poisson subalgebra": R
  is not bracket-closed, and only bracket compatibility inside D is used.
depends_on:
  - THM-3686-y0-collision-normalization-and-bracket-anatomy
  - THM-3569-danielewski-two-by-three-weight-darboux-nonentry
  - THM-3583-universal-exponent-two-two-by-four-weight-darboux-nonentry
  - THM-3592-universal-exponent-two-three-by-three-weight-darboux-nonentry
related:
  - THM-3691-y0-collision-ring-two-weight-darboux-no-go
  - THM-3693-y0-collision-ring-two-by-three-weight-darboux-no-go
script: 04-computation/jacobian_y0_collision_danielewski_embedding_thm3695.py
output: 05-knowledge/results/jacobian_y0_collision_danielewski_embedding_thm3695.out
script_sha256: e3303baa68e110cffd9db828781ba3588afccc28fde12b68b1c9dfca67a863aa
output_sha256: 63551951598e0086f32748b640f28059bad8d55d91099c44d2d0286be0446173
hash_basis: LF-normalized bytes
---

# THM-3695 -- the collision ring already lives in the universal Danielewski model

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  This theorem is
an inheritance result.  It identifies the exact
old object containing the new collision ring and then applies three already
audited universal theorems.  It supersedes the narrower planned
target-degree-four arithmetic `3 x 3` computation that originally reserved
this theorem number.

All rings are over `C`.  On the source plane with coordinates `(x,z)`, put

```text
b=1-x^2z,             c=x,             e=z(2-x^2z).    (1)
```

Then

```text
c^2e=1-b^2.                                             (2)
```

Thus the assignments in `(1)` define a homomorphism

```text
D=C[b,c,e]/(c^2e-(1-b^2))  -->  C[x,z].                 (3)
```

It is injective.  Indeed, `c^2e+b^2-1` is primitive and linear in `e`, hence
irreducible, so `D` is a domain.  After inverting `c`, `(3)` is the
isomorphism

```text
D[c^-1]=C[b,c,c^-1]  ~=  C[x,x^-1,z],
x=c,                 z=(1-b)/c^2.                      (4)
```

The localization map from the domain `D` is injective, proving the claim.

## 1. The ambient embedding preserves the Poisson object

Use the source bracket `{F,G}=F_xG_z-F_zG_x`.  Direct differentiation gives

```text
{b,c}=c^2,             {c,e}=2b,             {b,e}=-2ce. (5)
```

For `Sigma(b)=1-b^2`, these are exactly

```text
{b,c}=c^2,       {c,e}=-Sigma'(b),       {b,e}=-2ce,   (6)
```

the bracket used in THM-3569/3583/3592.  Moreover `Sigma` is squarefree of
degree two, so every hypothesis of those universal theorems holds.  Give `D`
their grading

```text
wt(b,c,e)=(0,1,-2).                                    (7)
```

The three collision-ring generators of THM-3686 become

```text
A=3e,                 B=2ce,                 C=bc.     (8)
```

Consequently

```text
R=C[A,B,C]=C[3e,2ce,bc]  subset D                      (9)
```

is a graded subalgebra, with `wt(A,B,C)=(-2,-1,1)`.  It is important that
`R` is not asserted to be closed under the bracket; for example `{A,B}` need
not lie in `R`.  What `(5)--(6)` prove is that the injection `D -> C[x,z]` is
Poisson.  Hence for `P,Q in R subset D`, their ambient bracket maps exactly to
their source bracket.  In particular, a source equation `{P,Q}=1` is the same
equation in `D`, by injectivity.

Because the inclusion `R subset D` is injective and grading-preserving, the
nonzero homogeneous pieces
of an element of `R` are exactly the pieces seen when that element is regarded
in `D`; neither quotienting nor cancellation changes the support count.

## 2. The inherited seven-piece floor

Let `P,Q in R` and suppose `{P,Q}=1`.  Regard them as elements of `D`, and
subtract scalar weight-zero components before counting.

* THM-3569 excludes a one-piece output against an arbitrary mate and excludes
  every `2 x 3` or smaller pair.
* THM-3583 excludes `2 x 4`; hence a two-piece output needs at least five
  pieces in its mate.
* THM-3592 excludes `3 x 3`.

These cases exhaust every positive support partition of total size at most
six.  Therefore

```text
#supp(P)+#supp(Q) >= 7.                                (10)
```

At equality, the one-by-six boundary is already excluded.  The only
partitions not decided by the cited theorems are

```text
(2,5), (3,4), (4,3), (5,2).                           (11)
```

Thus, up to exchanging the outputs, the first live collision-ring cells are
`2 x 5` and `3 x 4`.  This is strictly stronger than separately proving the
finite target-degree-four arithmetic `3 x 3` gate and also reveals the right
next objects: the universal Danielewski arguments permit arbitrary
coefficients in `C[b]`, whereas membership in the proper subalgebra `(9)`
adds simultaneous gluing constraints at `b=-1,0,1`.  Those extra constraints
are not used in `(10)` and are the natural sidecar for attacking `(11)`.

Finally, `(8)` retains the explicit source collision

```text
(x,z)=(1,0),(-1,2)  -->  (A,B,C)=(0,0,1).              (12)
```

Hence any Darboux pair found in `R` would still give a noninjective planar
Keller map.  This theorem proves only the support floor; it does not exhibit
such a pair and does not settle `JC(2)`.

## Reproduction

```bash
python3 04-computation/jacobian_y0_collision_danielewski_embedding_thm3695.py
python3 -O 04-computation/jacobian_y0_collision_danielewski_embedding_thm3695.py
```

Both commands must agree byte for byte with the frozen output.
