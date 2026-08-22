---
id: THM-3545
title: "Catalan self-intersection Keller thickening and polynomial boundary"
status: >
  PROVED + VERIFIED-EXACT.  In the separated thickening ansatz
  P=v^2+A(w), Q=v^3-v+vB(w), A(0)=B(0)=0 over a characteristic-zero field,
  Jac(P,Q)=kappa in k* has exactly one formal solution:
  B=3A/2 and A-3A^2/4=kappa w.  Hence
  A=2(1-sqrt(1-3kappa w))/3 has nonzero Catalan coefficients in every
  positive degree and is never polynomial.  The formal/holomorphic map has
  the transverse double collision (+/-1,0)->(1,0) while retaining constant
  Jacobian.  This is an exact algebraic near-counterexample and a no-go only
  for the stated separated ansatz; it is not a polynomial endomorphism of
  A^2 and proves neither a counterexample nor JC(2).
source: boxeph-2026-08-18-planar-jacobian-dephasing
depends_on: []
related:
  - THM-2063-one-fiber-linear-planar-keller-pairs
  - THM-2102-power-free-weight-face-and-first-defect-descent
  - THM-3543-torus-quotient-ramification-square-no-go
script: 04-computation/catalan_self_intersection_keller_thickening_thm3545.py
output: 05-knowledge/results/catalan_self_intersection_keller_thickening_thm3545.out
script_sha256: 70bcccaf01bdc03de1833f21b908ca9af0b4d7832330ce6faba32d107ae26dd4
output_sha256: d7cde63e4c7d30759ae5713695ca3d692563be993728f3d895b9c7101e478880
hash_basis: LF-normalized bytes
---

# THM-3545 -- Catalan self-intersection Keller thickening

**PROVED + VERIFIED-EXACT.**

## 1. The separated self-intersection ansatz

Let `k` be a characteristic-zero field, let `kappa in k*`, and work in
`k[v][[w]]`.  For series `A,B in w k[[w]]`, put

```text
P(v,w)=v^2+A(w),
Q(v,w)=v^3-v+vB(w).                                   (1)
```

Then

```text
Jac_(v,w)(P,Q)=kappa                                  (2)
```

if and only if

```text
B=(3/2)A,
A-(3/4)A^2=kappa w.                                   (3)
```

Consequently `(2)` has exactly one solution in the normalized ansatz:

```text
A(w)=(2/3)(1-sqrt(1-3kappa w)),
B(w)=1-sqrt(1-3kappa w).                              (4)
```

Here the square root is the unique formal branch with constant term one.
For `k=C`, `(4)` also defines a single-valued holomorphic branch on a small
disk about `w=0`.

At `w=0`, the restriction of `(1)` is the self-intersecting polynomial curve

```text
gamma(v)=(v^2,v^3-v).                                 (5)
```

The two distinct points satisfy

```text
(P,Q)(1,0)=(P,Q)(-1,0)=(1,0),                         (6)
```

and their boundary tangent vectors `(2,2)` and `(-2,2)` have determinant
`8`.  Thus the collision is transverse while the full formal/holomorphic
Jacobian remains the nonzero constant `kappa`.

## 2. Rigidity proof

Differentiate `(1)` directly:

```text
P_v=2v,                 P_w=A',
Q_v=3v^2-1+B,           Q_w=vB'.                      (7)
```

Therefore

```text
Jac(P,Q)=(2B'-3A')v^2+A'(1-B).                        (8)
```

Identity `(2)` is equivalent to the two coefficient equations

```text
2B'=3A',                 A'(1-B)=kappa.               (9)
```

The first equation and `A(0)=B(0)=0` give `B=3A/2`.  Substitution into the
second gives

```text
A'(1-(3/2)A)=kappa.                                   (10)
```

Integration in the characteristic-zero formal power-series ring, followed
by `A(0)=0`, yields the second equation in `(3)`.  Conversely, differentiating
that equation gives `(10)`, and `(8)` then gives `(2)`.  Solving the quadratic
with the normalized square-root branch proves `(4)` and uniqueness.

Equation `(6)` follows immediately from `A(0)=B(0)=0`; the tangent
calculation proves transversality.  This completes the proof.

## 3. Catalan walk expansion and the polynomial obstruction

Rewrite `(3)` as

```text
A=kappa w+(3/4)A^2.                                   (11)
```

The standard quadratic generating-function recurrence gives

```text
A(w)=sum_(n>=1) C_(n-1)(3/4)^(n-1) kappa^n w^n,       (12)
C_j=(1/(j+1))*binom(2j,j).                            (13)
```

One can also obtain `(12)` by expanding `(4)`.  Every displayed coefficient
is nonzero in characteristic zero, so the solution does not terminate.
Equivalently, if a nonconstant polynomial `A` has degree `d`, then

```text
deg(A'(1-(3/2)A))=2d-1>0,                             (14)
```

whereas `(10)` asks for a unit.  Constant `A` is excluded by `A(0)=0` and
`kappa!=0`.  Hence there is no polynomial solution in the separated ansatz.

The Catalan numbers in `(12)` count rooted binary trees, or equivalently
Dyck excursions.  This is an exact positive walk expansion, not an assertion
that `(1)` is itself a finite electrical network.  Its relevance is that the
constant-Jacobian thickening accumulates a nonterminating positive walk tail;
a genuine polynomial counterexample would have to leave this ansatz and
cancel or reroute that tail.

## 4. What the theorem does not say

The map `(1)` with `(4)` is not a polynomial endomorphism of `A^2`.  It lives
on a formal branch, or analytically on a proper domain; its collision does not
contradict the planar Jacobian conjecture.  The theorem excludes only the
separated form `(1)`.  It does not exclude mixed corrections such as

```text
sum_(i,j>=1) c_(ij) v^i w^j                           (15)
```

in either component, algebraic elimination through extra variables, or a
different self-intersecting boundary curve.  Those are precisely the live
counterexample-construction directions exposed by this near miss.

The closest polynomial hostile is the torus quotient reserved in
`THM-3543-torus-quotient-ramification-square-no-go.md`: it retains
polynomiality and a large collision but loses constant Jacobian through a
square ramification factor.  THM-3545 retains the collision and constant
Jacobian but loses polynomiality.  This comparison is a typed pincer, not a
proof that the two defects exhaust all constructions.

## 5. Exact verification

Reproduce with

```bash
python3 04-computation/catalan_self_intersection_keller_thickening_thm3545.py
python3 -O 04-computation/catalan_self_intersection_keller_thickening_thm3545.py
```

The dependency-free algebra is checked symbolically along two paths: direct
Jacobian differentiation of `(4)`, and the generic coefficient identity
`(8)`.  The companion verifies the quadratic equation, normalization,
collision, transverse tangents, and the Catalan formula through degree `13`.
For every truncation degree `1<=N<=12`, it also locates the first nonconstant
Jacobian defect exactly at degree `N`.  Ordinary and optimized transcripts
agree byte-for-byte.

```text
script_sha256 = 70bcccaf01bdc03de1833f21b908ca9af0b4d7832330ce6faba32d107ae26dd4
output_sha256 = d7cde63e4c7d30759ae5713695ca3d692563be993728f3d895b9c7101e478880
hash_basis    = LF-normalized bytes
```

**QED.**
