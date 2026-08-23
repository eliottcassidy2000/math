---
id: THM-3854
title: "The integrated three-cusp front is naturally an S5 quintic discriminant, not an S3 cubic"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The integrated
  squarefree-cubic front is an
  irreducible rational quintic with three A2 cusps, three transverse A1
  nodes, and one smooth normalization place at infinity.  It is exactly the
  ramification image of a finite polynomial-plane quintic whose generic
  geometric monodromy is S5.  Its discriminant quadratic layer has group A5
  and therefore no cyclic cubic quotient.  Two explicit seminormal
  square/cube descents have nonsquare residual discriminants.  The cuspidal
  normalization also cannot be the restriction of a Keller map to any
  interior smooth source arm, in arbitrary normal degree.  These statements
  still leave the three-torsion of the quadratic-resolvent surface,
  alternative cubic orders, boundary/Jelonek realizations, Keller atlases,
  and JC(2) OPEN.
source: jc_zero_debt_lift / integrated cusp-front completion lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (root, 2026-08-23).  The audit separately
  checked the finite birational normalization, singularity exhaustion by
  the collision resultant and genus count, the smooth unique infinity
  place, the finite polynomial-plane quintic, and the discriminant identity.
  It also checked that the two good specialization cycle types generate S5,
  that nonsquareness over the algebraic closure forces geometric S5, that
  A5 has no C3 quotient, and that both residual nonsquare arguments are
  valid in the rational function field.  The companion verifies the
  normalization and rational inverse, complete collision packet, all cusp
  and node types, smooth fivefold infinity contact, polynomial quintic map,
  exact discriminant constant, the all-degree constant-bracket arm factor,
  two good-prime cycle types, and both seminormal residual quotients.  Normal
  and optimized runs byte-match the frozen 55-gate transcript and both
  recorded hashes.  A second 74-gate audit recovers the complete singular
  scheme directly, computes the specialized Galois group by an exact
  number-field algorithm, checks A5 perfection, proves the seminormal
  elements proper, and reconstructs the derivative-ideal arm obstruction.
related:
  - THM-3811-ramification-class-unit-criterion-and-nonlinear-cubic-packet
  - THM-3841-deleted-ramification-three-puncture-jelonek-nonentry
  - THM-3844-two-cusp-polynomial-branch-quadratic-resolvent-design-gate
  - THM-3851-tricuspidal-quartic-rank-two-two-place-tradeoff
  - THM-3856-quadratic-normal-strip-keller-pairs-are-automorphisms
script: 04-computation/jc2_integrated_three_cusp_quintic_s5_thm3854.py
output: 05-knowledge/results/jc2_integrated_three_cusp_quintic_s5_thm3854.out
script_sha256: ffe82ca9bd4c147b685d75b16bfa94e18269b143fc9cbe629f5edf227538ae5e
output_sha256: 0a4598e2f042a08afd7b1d046c3cea2bc345e48fc6282980e1d2f265a619847e
semantic_sha256: 1f9a5b7d00f3883d38bd65c8a2f0546b5627e7a341cd29560e4f80993fc92f24
independent_audit_script: 04-computation/jc2_integrated_three_cusp_quintic_s5_independent_audit_thm3854.py
independent_audit_output: 05-knowledge/results/jc2_integrated_three_cusp_quintic_s5_independent_audit_thm3854.out
independent_audit_script_sha256: 218b7f99ce8ab062bcc4fa1d20db0f9262a9991d365ab55e8284acfd65ab45b2
independent_audit_output_sha256: d720392427c5837100ce9bcdff939147aa8a304c7188e5528187179c74676f5c
independent_audit_semantic_sha256: e6b1f10501fb4bf74f8b7e8c92037e68016e1feadb45501d5ed1b6fad2f6e8ae
hash_basis: raw LF bytes
---

# THM-3854 -- the natural completion is quintic `S5`, not cubic `S3`

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over an
algebraically closed field `k` of characteristic zero.  Put

```text
g(t)=t(t^2-1),
x(t)=t^4-2t^2,                       x'(t)=4g(t),
y(t)=3t^5-5t^3,                     y'(t)=15t g(t).            (1)
```

The image is the irreducible quintic

```text
Delta = 81x^5+90x^4+25x^3+30x^2y^2+30xy^2-y^4+8y^2.          (2)
```

It has precisely three `A2` cusps and three transverse `A1` nodes in the
affine plane, and its projective closure has a unique point at infinity.
That point is smooth and is the unique normalization place at infinity.

More strongly, `(2)` is the discriminant front of the finite polynomial map

```text
Phi:A2_(x,T) -> A2_(x,y),
(x,T) |-> (x,(-3T^5+10T^3+15xT)/4).                           (3)
```

Its monic quintic equation, up to the unit `3`, is

```text
F(T)=3T^5-10T^3-15xT+4y,                                    (4)
disc_T(F)=-64,800,000 Delta.                                 (5)
```

The generic geometric Galois group of `(4)` is `S5`.  After adjoining
`sqrt(Delta)`, the natural splitting-field group is `A5`, which has no
nontrivial cyclic cubic quotient.  Thus the most immediate completion of
this one-place three-cusp front supplies the wrong cover grammar for the
desired `S3/C3` counterexample lane.

There is an independent all-degree arm obstruction.  No polynomial pair on
a smooth source graph `r=0` can both restrict to `(x(t),y(t))` and have
nonzero constant Jacobian.  Thus this cuspidal normalization can occur only
as boundary/Jelonek/discriminant data in a Keller design, not as an interior
smooth Keller arm.

This is an obstruction to **the natural quintic splitting field only**.  It
does not compute `Cl(k[x,y,W]/(W^2-Delta))[3]`, exclude another binary-cubic
order with discriminant `(2)`, construct a Keller atlas, or settle `JC(2)`.

## 1. Normalization and irreducibility

Direct substitution of `(1)` into `(2)` gives zero.  On a dense open set the
inverse normalization coordinate is

```text
t = y(9x^2+13x+4)/(27x^3+33x^2+10x+y^2).                     (6)
```

Indeed the denominator and numerator pull back respectively to

```text
E(t)=t^2(t-1)^2(t+1)^2(3t^2-5)(9t^4-18t^2+4),
tE(t).                                                        (7)
```

Thus `(1)` is birational onto its image.  Its projective form is the
basepoint-free degree-five map

```text
[t:s] |->
[t^4s-2t^2s^3 : 3t^5-5t^3s^2 : s^5].                        (8)
```

If `s!=0`, the last coordinate is nonzero; if `s=0`, the second
coordinate is `3t^5!=0`.  Hence `(8)` has no base point.
The image consequently has degree five.  The nonzero quintic `(2)` vanishes
on that irreducible image, so it is the image equation and is irreducible.
Equation `(8)` is therefore the projective normalization.

## 2. Complete singularity ledger

For a second normalization address `u`, divide both coordinate differences
by `t-u`.  Their resultant in `u` is exactly

```text
Res_u((x(t)-x(u))/(t-u),(y(t)-y(u))/(t-u))
  =t^2(t-1)^2(t+1)^2(3t^2-5)(9t^4-18t^2+4).                  (9)
```

This is both the ramification and off-diagonal collision packet.  At
`t=0,1,-1`, both derivatives in `(1)` vanish.  Put `t=a+v`.  At the three
addresses one may take the following local coordinate pairs:

```text
a=0:   x,                         y;
a=1:   x+1,                       y+2-(15/4)(x+1);
a=-1:  x+1,                       y-2+(15/4)(x+1).             (10)
```

Their orders in `v` are `(2,3)` in every row, with nonzero leading
coefficients.  Hence the corresponding image points

```text
(0,0),                         (-1,-2),                    (-1,2) (11)
```

are `A2` cusps.

The divided first-coordinate difference factors as

```text
(t+u)(t^2+u^2-2).                                           (12)
```

On `u=-t`, the second divided difference becomes
`t^2(3t^2-5)`.  The two roots of `3t^2-5` therefore collide at

```text
(-5/9,0).                                                     (13)
```

On the other factor of `(12)`, reduction of the second difference is

```text
3t^4-6t^2+tu+2.                                              (14)
```

The remaining four addresses obey

```text
9t^4-18t^2+4=0,               u=-2/(3t),                    (15)
```

and form two collision pairs at

```text
x=-4/9,                         y^2=8/27.                     (16)
```

Away from the cusp addresses the tangent slope is

```text
(dy/dt)/(dx/dt)=15t/4.                                      (17)
```

Every pair in `(13)` and `(15)` has distinct addresses, hence distinct
slopes.  These three points are transverse `A1` nodes.

There are no omitted affine singularities: a finite birational
normalization fails to be a local isomorphism precisely at a ramified or
multiply-addressed point, and `(9)` exhausts both.  Equivalently, the three
cusps and three nodes contribute total delta invariant six, equal to the
arithmetic genus `(5-1)(5-2)/2` of a quintic whose normalization is `P1`.

Homogenizing `(2)` to `Delta_h(X,Y,Z)`, one gets

```text
Delta_h(X,Y,0)=81X^5.                                       (18)
```

Thus `[0:1:0]` is the only projective point at infinity.  Moreover

```text
(partial Delta_h/partial Z)(0,1,0)=-1,                       (19)
```

so it is smooth.  Its pullback under `(8)` is the single address `s=0`,
with contact five.  The front therefore passes the one-place branch test
that the tricuspidal quartic of THM-3851 necessarily fails.

## 3. The tautological polynomial-plane completion

Equation `(4)` is linear in `y` with unit coefficient `4`.  Consequently

```text
k[x,y,T]/(F) isomorphic to k[x,T].                            (20)
```

After division by the unit `3`, `(4)` is monic of degree five in `T`.
Hence `(3)` is a finite degree-five polynomial map from an actual affine
plane.  Its Jacobian is

```text
J(Phi)=partial y/partial T
      =-(15/4)(T^4-2T^2-x).                                  (21)
```

The ramification curve is the polynomial line

```text
x=T^4-2T^2.                                                  (22)
```

Substitution of `(22)` into the second coordinate of `(3)` gives exactly
`y=3T^5-5T^3`; thus its ramification image is `(1)-(2)`.  Finally

```text
F_T=15(T^4-2T^2-x),
disc_T(F)=-64,800,000 Delta,                                 (23)
```

which proves `(5)` and shows that the entire front, including its
self-collisions, is the natural quintic discriminant.

This attractive completion is not Keller: `(21)` vanishes along `(22)`.
Deleting the ramification curve also makes its defining equation a
nonconstant unit, the same basic deleted-different debt isolated throughout
THM-3811 and THM-3841.

### 3.1. All normal degrees fail for an interior arm

The preceding failure is not repairable by simply adding higher powers of a
normal coordinate while fixing the arm.  Let `A(t,r),C(t,r) in k[t,r]` be
arbitrary and suppose

```text
A(t,0)=x(t),                         C(t,0)=y(t).               (23a)
```

Write

```text
alpha(t)=A_r(t,0),                   beta(t)=C_r(t,0).          (23b)
```

The constant `r`-coefficient of the bracket is necessarily

```text
J_(r,t)(A,C)|_(r=0)
 =alpha(t)y'(t)-x'(t)beta(t)
 =g(t)(15t alpha(t)-4 beta(t)).                                (23c)
```

Because `g=t(t^2-1)` is a nonunit, `(23c)` cannot equal an element of
`k*`.  Terms of order `r^2` and higher do not enter this constant bucket, so
the obstruction is independent of the normal degree and survives formal
power-series corrections as well.

Geometrically, the differential of an etale map cannot kill the nonzero
tangent vector of a smooth source arm.  But `(x'(t),y'(t))=(0,0)` at all
three cusp addresses.  This explains why the immersed, derivative-unimodular
Russell arm of THM-3843 is a different object.  THM-3856's quadratic-strip
classification becomes relevant only after replacing the cuspidal carrier;
it is not a degree-three escape for `(1)`.

## 4. Exact `S5` monodromy

The polynomial `F` is irreducible in `k[x,y,T]`: it is degree one in `y`
with unit leading coefficient, and its quotient is the domain `(20)`.
It is primitive as a polynomial in `T`, so Gauss's lemma gives generic
irreducibility over `k(x,y)`.

To determine the group rather than merely its transitivity, first work over
`Q(x,y)` and specialize at `(x,y)=(-3,1)`:

```text
f(T)=3T^5-10T^3+45T+4,
disc(f)=834,688,800,000.                                     (24)
```

The primes `29` and `67` do not divide `(24)`.  Modulo `29`, after making
the polynomial monic,

```text
f/3 = T^5-13T^3-14T+11                                      (25)
```

is irreducible.  The companion checks the degree-five Frobenius criterion
`T^(29^5)=T mod f` together with
`gcd(f,T^29-T)=1`.  Thus the specialized Galois group contains a five-cycle.

Modulo `67`, one has

```text
f = 3(T+7)(T+21)(T-23)(T^2-5T+5),                           (26)
```

and `5` is a nonsquare.  This is the cycle type `(1,1,1,2)`, so the same
group contains a transposition.

A five-cycle `c` and any transposition `(a b)` generate `S5`: conjugating
the latter by the powers of `c` gives all translated edges of a nonzero
difference on `Z/5`; because five is prime, that graph is connected, and
edge transpositions of a connected graph generate the full symmetric
group.  Therefore the specialized group is `S5`.  Good specialization
embeds it in the generic arithmetic group, which must also be `S5`.

The geometric generic group is a normal subgroup of this `S5`.  The
irreducible quintic `Delta` occurs to odd order in `(5)`, so the discriminant
is not a square even over the algebraic closure of the constant field.
The geometric group is therefore not contained in `A5`; the normal-subgroup
classification of `S5` forces it to be `S5`.  This remains true after base
change to any algebraically closed characteristic-zero `k`.

The sign kernel after adjoining `sqrt(Delta)` is consequently `A5`.  Since
`A5` is simple nonabelian (in particular perfect), it has no nontrivial
`C3` quotient.  Equivalently, the `S5` splitting field has neither the
cyclic cubic quotient over its discriminant field nor the `S3` quotient over
the target that the natural binary-cubic/Kummer construction would require.

## 5. Two seminormal cubic controls

The failure is visible before the full quintic Galois computation.  On the
normalization define

```text
h_e=t^2(t^2-1)(9t^2-14).                                    (27)
```

Its square and cube descend to the branch coordinate ring as

```text
P_e=81x^3+49x^2+8y^2,
Q_e=-243x^4-143x^3+81x^2y^2+120xy^2+64y^2,                  (28)
```

and exact division gives

```text
P_e^3-Q_e^2
 =Delta(6561x^4+3888x^3-512y^2).                            (29)
```

An independent odd-address control is

```text
h_o=t(t^2-1)(9t^2-4)(3t^2-5),                               (30)
P_o=-648x^3-720x^2+81xy^2-200x+49y^2,
Q_o=6561x^4y+17010x^3y+18009x^2y
    +243xy^3+8760xy+143y^3+1600y.                           (31)
```

Again `P_o=h_o^2` and `Q_o=h_o^3` on the normalization, but

```text
P_o^3-Q_o^2
 =-Delta(9x+5)^2
   (41472x^2+6561xy^2+46080x+3888y^2+12800).                (32)
```

For the depressed cubic `S^3-3PS+2Q`, the discriminant is
`108(P^3-Q^2)`.  Neither residual factor in `(29)` and `(32)` is a square in
`k[x,y]`.  Indeed each has `y`-degree two, a nonzero `y^2` coefficient, no
linear `y` term, and a nonzero `y`-independent part.  If it were
`(a(x)y+b(x))^2`, the missing cross term would force `a(x)b(x)=0`,
contradicting the other two facts.  In `(32)` the displayed `(9x+5)^2` is
already separated off.  Since `k[x,y]` is a UFD and every nonzero constant
of `k` is a square, the same valuation argument rules out a square in
`k(x,y)`: a polynomial square in the fraction field would already be a
constant times a polynomial square.

Thus two explicit low-degree elements of the cusp seminormalization do have
the required square/cube descent, but neither makes the discriminant a
square multiple of `Delta`.  These are corroborating controls, not an
exhaustion of all binary cubic orders.

## 6. The exact remaining design problem

The affine quadratic-resolvent surface

```text
Q: W^2=Delta.                                                 (33)
```

is not assigned a class-group conclusion here.  The natural double-plane
compactification branches over the sextic `Z Delta_h=0`.  At infinity the
line and quintic are smooth and meet with multiplicity five by `(18)-(19)`,
so the double cover has an `A9` packet there; the finite packets are
`3A2+3A1`.  Its minimal resolution is therefore a K3 surface, not the weak
degree-two del Pezzo used for the finite lattice computations in THM-3844
and THM-3851.  In particular their root-lattice quotient argument does not
compute `Cl(Q)[3]` here.

The live positive question is consequently precise:

> Does `(33)` carry anti-invariant three-torsion not generated by the
> natural `A5` splitting field, and can such a class be realized by a
> nonmonogenic binary-cubic order while retaining the one-place branch?

THM-3854 proves only that the tautological degree-five completion and the two
displayed seminormal cubic descents do not answer that question.

## 7. Reproduction and scope

Run

```bash
python3 04-computation/jc2_integrated_three_cusp_quintic_s5_thm3854.py
python3 -O 04-computation/jc2_integrated_three_cusp_quintic_s5_thm3854.py
python3 04-computation/jc2_integrated_three_cusp_quintic_s5_independent_audit_thm3854.py
python3 -O 04-computation/jc2_integrated_three_cusp_quintic_s5_independent_audit_thm3854.py
```

and compare both outputs byte-for-byte with
`05-knowledge/results/jc2_integrated_three_cusp_quintic_s5_thm3854.out`.
The companion uses exact integer/rational polynomial arithmetic.  The two
finite-field reductions are good-prime Frobenius witnesses inside an exact
characteristic-zero specialization argument; they are not used as rank
proxies or evidence for `Cl(Q)`.

The independent companion performs 74 exact gates by direct resultants,
singular-ideal elimination and an exact specialization Galois computation.
Its normal and optimized streams byte-match the frozen output.

No planar Keller counterexample, cubic completion, class-group computation,
or general obstruction to one-place branch curves is claimed.
