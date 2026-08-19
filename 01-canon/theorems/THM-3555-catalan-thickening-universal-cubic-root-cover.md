---
id: THM-3555
title: "Catalan thickening as the universal cubic marked-root cover"
status: >
  PROVED + VERIFIED-EXACT.  After adjoining
  r=sqrt(1-3*kappa*w), THM-3545's Catalan map becomes polynomial and is
  affinely equivalent to the universal marked-root cover
  (t,p)->(p,-t^3-pt).  The Catalan constant Jacobian is exactly the
  cancellation of this cubic cover's ramification factor by dr/dw.  The
  cover is finite flat of degree three, etale off the cubic discriminant,
  and its fiber over the Catalan collision has three simple points.  Any
  polynomial correction that fixes the ramification line pointwise still
  has zero Jacobian at its cusp point.  This is a construction template and
  obstruction, not a planar Keller counterexample.
source: kps-s188
depends_on:
  - THM-3545-catalan-self-intersection-keller-thickening-boundary
  - THM-1300-jacobian-counterexample-dixmier-A3-explicit
  - THM-1310-jacobian-counterexample-fiber-geometry-S3-resolvent-jelonek
  - THM-3554-punctured-kummer-collision-surface-normal-form
companion: 04-computation/catalan_universal_cubic_root_cover_kps_s188.py
output: 05-knowledge/results/catalan_universal_cubic_root_cover_kps_s188.out
script_sha256: cd9f593fa8efe94b839e2175f63b7b19724d36815496042904daffd7a38e58ef
output_sha256: 4b72a07e5914c114defedea2b1dfc6571005bb9dd27ec0340a5a71f3331533d8
hash_basis: LF-normalized bytes
---

# THM-3555 -- Catalan thickening as the universal cubic marked-root cover

**PROVED + VERIFIED-EXACT.**  The nonterminating square root in THM-3545 is
not merely a generating-function accident.  It is the local coordinate that
divides out the ramification of the universal degree-three root cover.

The field `k` has characteristic zero, and `kappa` is a nonzero scalar.

## 1. Polynomializing the Catalan branch

THM-3545's unique normalized formal solution is

```text
A=(2/3)(1-sqrt(1-3kappa w)),
B=1-sqrt(1-3kappa w),                                  (1)
P=v^2+A,                 Q=v^3-v+vB.                   (2)
```

Adjoin the algebraic coordinate

```text
r^2=1-3kappa w                                         (3)
```

and use the branch `r=1+O(w)` at `w=0`.  Equations `(1)`--`(2)` become the
polynomial map

```text
H(v,r)=(v^2+(2/3)(1-r), v^3-vr).                       (4)
```

Direct differentiation gives

```text
det dH/d(v,r)=-(2/3)r.                                 (5)
```

On the formal branch of `(3)`,

```text
dr/dw=-3kappa/(2r).                                    (6)
```

The chain rule now explains THM-3545's constant determinant without any
coefficient recurrence:

```text
det dH/d(v,w)=(-(2/3)r)(-3kappa/(2r))=kappa.           (7)
```

Thus the formal square root contributes exactly the inverse of the
ramification factor in `(5)`.  At `r=0` the coordinate change `(3)` ramifies;
the cancellation cannot be made by a polynomial automorphism of `A^2`.

## 2. Affine equivalence to the universal cubic root cover

Make the polynomial source automorphism

```text
t=v,                    p=2r-3v^2,                    (8)
```

with inverse `r=(p+3t^2)/2`, and the affine target automorphism

```text
p=2-3P,                 q=2Q.                          (9)
```

Their Jacobians are `2` and `-6`, respectively.  Substitution in `(4)` gives

```text
q=-t^3-pt.                                               (10)
```

Therefore `H` is left-right affinely equivalent to

```text
G:A^2_(t,p) -> A^2_(p,q),
G(t,p)=(p,-t^3-pt).                                    (11)
```

The source point `t` is a marked root of

```text
X^3+pX+q=0.                                             (12)
```

Conversely every root of `(12)` gives one point of the fiber of `(11)`, so
`G` is the universal marked-root cover for depressed cubics.

Algebraically,

```text
k[t,p] ~= k[p,q,t]/(t^3+pt+q)                          (13)
```

is free of rank three over `k[p,q]`, with basis `{1,t,t^2}`.  Hence `G` is
finite flat of degree three.  Its Jacobian and the cubic discriminant are

```text
R=p+3t^2,
Delta=-4p^3-27q^2.                                     (14)
```

It is finite etale over `Delta!=0`.  The generic cubic is irreducible and
`Delta` is not a square, so its generic Galois closure has group `S_3`.
The pullback factorization is

```text
G^*Delta=-(p+3t^2)^2(4p+3t^2).                         (15)
```

The squared factor is the ramification line `R=0`; the second factor records
where the two *unmarked* roots collide.  In `(v,r)` coordinates, `(15)` is

```text
G^*Delta=-4r^2(8r-9v^2).                               (16)
```

The line `r=0` maps to the cuspidal discriminant by

```text
(P,Q)=(v^2+2/3,v^3),
(3P-2)^3=27Q^2.                                        (17)
```

## 3. The Catalan pair is two points of a three-point fiber

The two points retained in THM-3545 are `(v,w)=(+/-1,0)`.  On the chosen
branch `r=1`, they are

```text
(v,r)=(1,1),(-1,1).                                    (18)
```

Polynomialization supplies a third point,

```text
(v,r)=(0,-1/2).                                        (19)
```

All three map under `(4)` to `(P,Q)=(1,0)`.  Their Jacobians are,
respectively,

```text
-2/3, -2/3, 1/3,                                      (20)
```

so this is a completely split simple fiber, not a branch fiber.  Under
`(8)`--`(9)` it is the depressed cubic

```text
X^3-X=X(X-1)(X+1).                                     (21)
```

This makes the recurring `1+2` collision anatomy exact: THM-3545 displayed
the symmetric pair on one local square-root branch, while the universal
cubic cover reveals the third marked root.

## 4. First-jet gate for ramification surgery

One tempting counterexample program is to correct `(11)` while fixing its
ramification line.  The first jet already rules this out.  Let

```text
G_tilde=(p+R A(t,p), -t^3-pt+R B(t,p)),
R=p+3t^2.                                               (22)
```

Thus `G_tilde` agrees with `G` pointwise on `R=0`.  Restricting its Jacobian
to that line gives

```text
det JG_tilde|_(R=0)=-6t(t A+B)|_(R=0).                 (23)
```

Indeed, on `R=0` the two differential rows are

```text
(6tA,1+A),                 (6tB,-t+B),                 (24)
```

whose determinant is `(23)`.  At the cusp preimage `(t,p)=(0,0)`, equation
`(23)` vanishes for every polynomial `A,B`.  Therefore no correction of the
form `(22)` has a nonzero constant Jacobian.

This is stronger than saying that second-order corrections are too small:
any polynomial surgery that leaves the entire ramification line fixed is
impossible.  A viable deformation must move that line already at order zero,
while separately retaining a chosen collision fiber.

## 5. The exact `1 plus 2` / connected-cubic design gap

THM-3554 supplies a complementary model inside the fixed THM-1300 Keller
map.  Since `F_3=xC`, with `C=2-3xy-x^2z`, the preimage of the target plane
`F_3=0` is the disjoint scheme union

```text
V(xC)=V(x) disjoint-union V(C);                         (25)
```

the ideals are comaximal because `C mod x=2`.  On `V(x)`,

```text
(y,z)->(F_1,F_2)=(z+4y^2,y)                            (26)
```

is a polynomial-plane isomorphism.  On `V(C)`, THM-3554 gives the degree-two
Kummer cover `(s,b)->(b,4s^2)`.  Away from its target parabola, `(25)` is
therefore a disconnected etale cover of degree `1+2`.

The universal root cover `(11)` has the opposite defect: its three sheets
are connected by generic `S_3` monodromy, but it ramifies over the
discriminant cusp.  The full three-dimensional THM-1300 map combines the two
desired properties--connected generic `S_3` fibers and no finite
ramification--by exporting failure of properness to its Jelonek surface
(THM-1310).

This triangulates a concrete planar construction problem:

```text
disconnected 1+2 Kummer slice: etale, but punctured/disconnected;
connected universal cubic:     three sheets, but finitely ramified;
desired planar counterexample: connected three sheets, with all branching
                               replaced by nonproper escape at infinity. (27)
```

Equation `(23)` says that merely correcting the normal derivative of the
universal cover cannot perform this surgery.  The branch curve itself must
move, and the collision constraints must be imposed only on selected fibers,
not on the whole branch line.

## 6. Scope

This theorem does not construct a planar Keller map.  It does not exclude
mixed polynomial corrections that move `R=0`, higher cover degree, a
different Jelonek curve, or affine modifications that change the function
field.  It proves an affine normal form, its finite-cover and discriminant
structure, and the fixed-ramification-line no-go `(23)`.  The final line of
`(27)` remains the open counterexample architecture.

This identifier is reserved for an audited affine-equivalence theorem between
the algebraic square-root polynomialization of THM-3545's Catalan thickening
and the universal depressed-cubic marked-root cover.  No theorem, proof, or
dependency is asserted by this stub.
