---
id: THM-3745
title: "Monomial-plane branch conductor, triangular delta, and Pell square selector"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  For every monic
  polynomial F of degree m>=2 over a field, the normalization of
  A=k[F(b),bF(b)] is B=k[b], B/A is the direct sum of k[X]/(X^i) for
  1<=i<m, its delta invariant is T_(m-1), and its conductor is
  F(b)^(m-1)B.  If F is separable, the unique singular image is geometrically
  an ordinary m-fold point.  The degrees with square delta are exactly a Pell orbit.
  At m=9 the numbers 36,72,108 form a typed planar-JC near miss only; no
  Keller pair, LRC(14) result, or cross-conjecture reduction follows.
source: root + jc_arithmetic_archaeology + cross_domain_skeptic / 2026-08-23
audit: >
  PASS.  Independent review checked normalization, module decomposition,
  finite support, Gorenstein conductor duality, arbitrary-characteristic
  algebraic scope, the separability repair for the ordinary-fold label, Pell
  completeness, and JC typing.  The exact companion checks 10,186 gates,
  computes conductor quotients for degrees 2..10, includes repeated-root
  hostiles, and passes normal/-O byte identity.
depends_on: []
related:
  - THM-3335-square-triangular-pell-markov-pythagorean-selector
  - THM-3586-nodal-cylinder-cap38-width-period-and-second-conductor-keller-gates
  - THM-3696-y0-collision-ring-three-branch-conductor-and-graded-modules
  - THM-3734-automorphic-cohn-diagonal-binomial-divided-power-towers
  - THM-3770-vertical-principal-part-equalizer-and-log-canonical-dressing-gate
  - THM-3742-square-triangular-pell-mod13-central-sign-projective-cycle
  - THM-3744-pell-prefix-loneliness-constant-carry-exact-formula
script: 04-computation/monomial_plane_conductor_triangular_pell_thm3745.py
output: 05-knowledge/results/monomial_plane_conductor_triangular_pell_thm3745.out
script_sha256: 7bc16ec6ae49e5fc5050d0559e21ab1668512208f3b216ce76c86f140f87fe14
output_sha256: fccbe03773359e89002167953e8b4043e3650db63f83156161b56eef54fad581
hash_basis: raw LF bytes
---

# THM-3745 -- triangular normalization defect and its Pell selector

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

Because the legacy numeric namespace also contains `HYP-3745`, cite this
result as **THM-3745 plus its full theorem slug/file path**, never by the bare
number.

No literature-priority claim is made.  The `m=3` collision ring was already
proved in
[THM-3696](THM-3696-y0-collision-ring-three-branch-conductor-and-graded-modules.md).
The theorem below isolates and proves the all-degree mechanism.

## 1. Statement and characteristic scope

Let `k` be any field, let

```text
F(T)=T^m+c_(m-1)T^(m-1)+...+c_1T+c_0,      m>=2,      (1)
```

be monic, and in `B=k[b]` put

```text
X=F(b),             Y=bF(b),
A=k[X,Y] subset B.                                      (2)
```

### Theorem 1.1 (normalization, quotient, and conductor)

The fraction fields of `A` and `B` agree, `B` is the normalization of `A`,
and, as `k[X]`-modules,

```text
B/A = direct_sum_(i=1)^(m-1) k[X]/(X^i).               (3)
```

Consequently

```text
delta(A)=dim_k(B/A)=sum_(i=1)^(m-1)i=m(m-1)/2=T_(m-1). (4)
```

The conductor of `A` in `B` is

```text
c_(B/A)={z in B:zB subset A}=X^(m-1)B=F(b)^(m-1)B.     (5)
```

Equations `(3)--(5)` require only that `F` be monic.  If `F` is separable,
equivalently `gcd(F,F')=1`, the sole nonnormal image `(X,Y)=(0,0)` becomes an
ordinary `m`-fold point after a splitting-field extension.  Separability is
load-bearing only for this geometric branch interpretation, not for the
algebraic formulas.

### Corollary 1.2 (Pell selector for square delta)

For integers `m>=1` and `q in Z_(>=0)`,

```text
binom(m,2)=q^2  <=>  (2m-1)^2-8q^2=1.                  (6)
```

For the ring theorem retain `m>=2`; `m=1` is the degenerate Pell base point.
The complete nonnegative sequence is

```text
m=1,2,9,50,289,1682,9801,...,
q=0,1,6,35,204,1189,6930,... .                         (7)
```

## 2. Plane equation and normalization

Substitute `b=Y/X` into `F(b)=X` and clear `X^m`.  The image curve obeys

```text
G(X,Y)=Y^m+c_(m-1)XY^(m-1)+...+c_1X^(m-1)Y+c_0X^m
       -X^(m+1)=0.                                     (8)
```

The polynomial `F(T)-X` is irreducible over `k(X)`: before localization it is
primitive and linear in the indeterminate `X`.  Hence it is the minimal
polynomial of `b` over `k(X)`, and `(8)` is the degree-`m` monic relation for
`Y=Xb`.  Therefore

```text
A is isomorphic to k[X,Y]/(G),
Frac(A)=k(X,Y)=k(b)=Frac(B).                            (9)
```

The element `b` is integral over `A` because it satisfies the monic equation
`F(b)-X=0`.  Thus `B` is integral over `A`.  Conversely, any element of
`k(b)` integral over `A` is integral over the integrally closed PID `B`, hence
lies in `B`.  This proves the normalization assertion.

Away from `X=0`, equation `b=Y/X` makes the normalization map an isomorphism.
Every defect is therefore supported at the single image `(X,Y)=(0,0)`.

## 3. Exact module decomposition

The monic equation `F(b)=X` gives

```text
B=direct_sum_(i=0)^(m-1) k[X]b^i.                     (10)
```

Likewise, because `(8)` is monic of degree `m` in `Y`,

```text
A=direct_sum_(i=0)^(m-1) k[X]Y^i
 =direct_sum_(i=0)^(m-1) X^i k[X]b^i.                 (11)
```

Taking the coordinatewise quotient of `(10)` by `(11)` proves `(3)`, and
summing its finite `k`-dimensions proves `(4)`.

This is the source of the triangular number.  It is not inferred from a
finite list: the successive module defects have exact lengths
`1,2,...,m-1`.

## 4. The conductor is exactly F^(m-1)

Put

```text
I=X^(m-1)B.                                             (12)
```

For `0<=i<=m-1`,

```text
X^(m-1)b^i=X^(m-1-i)Y^i in A.                          (13)
```

Thus `I` is a `B`-ideal contained in `A`, so `I subset c_(B/A)`.

It remains to exclude a larger conductor.  Localize at the only defect
`(X,Y)`.  By `(8)`, the local ring `R` is a one-dimensional hypersurface,
hence Gorenstein.  Let `S` be its finite normalization and `c` its conductor.
For a one-dimensional Gorenstein local domain, conductor duality gives

```text
length_R(R/c)=length_R(S/R).                           (14)
```

For completeness, identify `Hom_R(S,R)` with `c`: every `R`-linear map
`S -> R` (where `S` is rank-one birational) is multiplication in the common
fraction field by its value at `1`, and the multiplier is exactly an element
of `S` carrying `S` into `R`.  Apply `Hom_R(-,R)` to
`0 -> R -> S -> S/R -> 0`.
The torsionfree module `S` is maximal Cohen--Macaulay, so its positive Ext
vanishes over the one-dimensional Gorenstein ring; local duality identifies
the remaining `Ext^1(S/R,R)` with the Matlis dual of `S/R`.  Taking lengths
gives `(14)`.

Equations `(3)` and `(14)` now give

```text
length_k(B/c)=length_k(B/A)+length_k(A/c)
             =2 delta(A)=m(m-1).                      (15)
```

On the other hand,

```text
length_k(B/I)=deg(F^(m-1))=m(m-1).                    (16)
```

Since `I subset c` and the two quotients have the same finite length,
`I=c`.  This proves `(5)` in every characteristic.

## 5. Separable branches and the repeated-root hostile

The degree-`m` tangent cone of `(8)` is

```text
X^m F(Y/X).                                            (17)
```

If `F` is separable, then over a splitting field

```text
X^mF(Y/X)=product_(F(alpha)=0)(Y-alpha X),             (18)
```

a product of `m` distinct lines.  Hence the origin is an ordinary `m`-fold
point: its normalization has `m` smooth branches and each pair is transverse.
The conductor has order `m-1` on every branch, exactly as `(5)` says.

Without separability, `(3)--(5)` survive but `(18)` can fail.  The minimal
hostile is

```text
F=b^2,
A=k[b^2,b^3],
c=b^2k[b],
delta=1.                                               (19)
```

This is a cusp with repeated tangent line `Y^2`, not an ordinary double point.
Any statement that calls all monic rows ordinary is therefore false.

## 6. The square-triangular Pell orbit

For integers `m>=1` and `q>=0`, equation `(4)` gives

```text
delta(A)=q^2
 <=> m(m-1)/2=q^2
 <=> (2m-1)^2-8q^2=1,                                 (20)
```

which proves `(6)`.  Starting from `(m,q)=(1,0)`, all solutions with
`m>=1` and `q>=0` are generated by

```text
(2m'-1)+q'sqrt(8)=(3+sqrt(8))((2m-1)+qsqrt(8)).        (21)
```

Indeed, put `x=2m-1`.  If `q>0`, inverse multiplication by `3-sqrt(8)` gives

```text
(x,q) |-> (3x-8q, 3q-x).                              (21a)
```

The Pell identity gives `sqrt(8)q<x<=3q` (with equality only at
`(x,q)=(3,1)`), so both entries in `(21a)` are nonnegative; also
`3x-8q<x` because `x<4q`.  Thus `(21a)` strictly descends in `x` and
terminates at the unique `q=0` solution `(x,q)=(1,0)`.  Reversing the descent
proves completeness of `(21)`.  This yields exactly `(7)` and is the
square-triangular orbit of
[THM-3335](THM-3335-square-triangular-pell-markov-pythagorean-selector.md),
with triangular index `m-1`.

The map is precisely typed:

```text
source:    normalization defect delta=T_(m-1)
target:    Pell conic x^2-8q^2=1
map:       (m,q) |-> (x,q)=(2m-1,q)
preserved: the scalar equality delta=q^2
destroyed: F, branch labels/incidence, conductor ideal, ring multiplication
sidecar:   the complete module packet (3) and conductor (5).              (22)
```

## 7. The m=9 shadow of (72,108), and why it is not JC

At the first nontrivial square row after `m=2`,

```text
m=9,          delta=36,
length(B/c)=2delta=72.                                 (23)
```

The cited sub-`125` classification leaves `(72,108)` as the only hypothetical
reduced counterexample degree pair; no such pair is known to exist
(THM-3586 and MISTAKE-427).  There is an exact syntactic leading-form shadow.
In characteristic zero, homogenize

```text
H(X,Y)=X^9 F(Y/X).
```

Then `c=F^8B`; if `K=H^4`, the two appended homogeneous forms satisfy

```text
K^2=H^8 has degree 72,
K^3=H^12 has degree 108.                               (24)
```

But `(24)` appends `H^12`; the conductor theorem selects only `F^8` in the
normalization ring.
Normalization length, conductor colength, and polynomial map degree are
different types.  Worse, the apparent pair is the hostile

```text
V^2=U^3             for (U,V)=(H^8,H^12),             (25)
```

so the displayed pair is algebraically dependent and has Jacobian zero.  A
zero top Jacobian is compatible with the leading forms of a Keller pair, so
this excludes the displayed common-power pair, not every possible lower-term
completion.

The squarefree tangent cone itself also has a strict log-canonical nonentry.
For separable `F`, write `H=prod_i H_i` over a splitting field.  If
`J(H,W)=lambda H` with `lambda!=0`, reduction modulo each line `H_i` makes
`W|_(H_i=0)` constant.  All lines meet at the origin, so all constants equal
`W(0)` and squarefreeness gives `W-W(0)=HV`.  Cancellation then gives
`J(H,V)=lambda`, impossible at the origin because `grad H(0)=0` for `m>=2`.
This is the direct component-equalization boundary abstracted in THM-3770;
it still does not exclude arbitrary inhomogeneous Keller leading forms.
Thus `(72,108)` is an **exact near-miss leading-form shadow**, never a JC
consequence.

There is a second component-count near miss with
[THM-3734](THM-3734-automorphic-cohn-diagonal-binomial-divided-power-towers.md)
at `r=8`: both displays have nine labelled components when `F=b^9-1`.
Here all nine branches meet at one ordinary point; THM-3734 has an axis plus
eight pairwise-disjoint hyperbolas and proves no polynomial mate.  Incidence
and conductor are destroyed, already at the three-component row.

## 8. LRC boundary and exact reproduction

The Pell selector `(20)` is only a scalar invariant.  No map sends
`A=k[F,bF]` to a thirteen-speed set while preserving maximum loneliness.
Routing `(20)` through the mod-13 clock of THM-3742 loses branch/conductor data
and still lacks a physical runner owner.  Nothing here proves LRC(14).

Reproduce the exact audit with

```bash
python3 -B 04-computation/monomial_plane_conductor_triangular_pell_thm3745.py
python3 -B -O 04-computation/monomial_plane_conductor_triangular_pell_thm3745.py
```

The companion performs `10,186` exact gates.  For split separable polynomials
of degrees `2..10` over `GF(1000003)`, it builds `B/I`, `A/I`, and the complete
linear system defining `c/I`, obtaining zero extra conductor.  It repeats the
audit for `F=b^m`, checks every square triangular row through `m=10000`, and
records the Jacobian-zero `(72,108)` hostile.  Normal and optimized outputs
are byte-identical.

The theorem proves an all-degree conductor mechanism and a lawful scalar Pell
selector.  It proves neither a bivariate Keller pair, JC(2), nor LRC(14).
