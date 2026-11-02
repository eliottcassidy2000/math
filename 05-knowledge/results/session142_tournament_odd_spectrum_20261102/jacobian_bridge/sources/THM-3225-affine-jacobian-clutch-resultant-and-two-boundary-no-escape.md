---
id: THM-3225
title: "Affine Jacobian clutch resultant and two-boundary no-escape"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  On both THM-3212 accessory algebras and all four heptic covers, every
  affine clutch B=b0+b1*x leaves a critical point in
  P_B=(Vz^2+Bz)^2+A_src*z+x.  A clutch meeting ST gives an explicit finite
  critical point.  Otherwise the saturated resultant has degree 52, no
  multiplicity escapes through T, and its sole S wall loses at most two,
  leaving at least 50.  Generic affine clutches have 52 reduced Morse points.
  This supplies no Keller mate, owner, or proof of JC(2).
source: root/frontier-sidecars-cont-2026-08-02
audit: >
  A self-contained exact companion derives the new covariant
  J_B=2VB'-BV', its universal y-resultant, all six affine-parameter boundary
  quotients, the T-leading coefficient, and the S-wall jets.  Exact PRS
  calculations in both cubic accessory fields give degrees
  (5,3),(3,2),(2,1),(1,0) and gcd one.  Good reductions prove the B=1+x
  control squarefree and boundary-disjoint and audit finite-collision and
  exact S-wall hostiles.  Normal, optimized, and stored transcripts agree
  byte-for-byte.  An independent native number-field implementation rederived
  the gradient reduction, resultant, saturation, local valuations, PRS,
  squarefree control, Morse interpretation, and every scope boundary.
depends_on:
  - THM-3212-centered-heptic-source-morse-obstruction-and-offcenter-clutch
related:
  - THM-3123-heptic-e3-remaining-accessory-classification-and-s7-monodromy
  - THM-3172-shear-invariant-differential-owner-filtration-and-transverse-recurrence
  - THM-3133-common-simple-zero-faber-exclusion-and-odd-bipole-boundary
script: 04-computation/jc_heptic_affine_B_source_obstruction_thm3225.py
output: 05-knowledge/results/jc_heptic_affine_B_source_obstruction_thm3225.out
script_sha256: 0ada8f35c1523ee802c52ca5e090f990909b829c7ebb767ac1b1f744607ad631
output_sha256: 7237121880bb8e38841f969dce420595a05980437eaaae3682395b7a6b854652
hash_basis: LF-normalized bytes
---

# THM-3225 -- affine Jacobian clutch resultant and two-boundary no-escape

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Statement

Let `K_i` be either of the two cubic accessory fields in THM-3212, let
`K_i -> K_0` be any embedding into an algebraically closed characteristic-zero
field, and retain the response polynomials

```text
V=4SDT^2/C^2,                 A=A_src=2SET/C,
g=ST,                         2VA'-AV'=2V.
```

The embeddings cover the one `(4,1,1,1)` and three `(3,2,1,1)` unmarked
`S_7` covers.  For every affine polynomial

```text
B(x)=b_0+b_1 x in K_0[x]
```

the polynomial

```text
P_B(x,z)=(V(x)z^2+B(x)z)^2+A(x)z+x                         (1)
```

has a critical point.  More precisely, either:

1. `B` meets `g`, in which case `(1)` has an explicit critical point over
   that common root; or
2. `B` is a unit modulo `g`, in which case at least `50` units of critical
   resultant multiplicity remain on `V!=0`.

On a nonempty Zariski-open set of affine `B`, there are exactly `52` distinct
reduced Morse critical points.  The rational control `B=1+x` lies in that open
for both accessory fields.  Therefore no `P_B` is one coordinate of a
polynomial Keller pair.

The count `>=50` in the exceptional lane is an intersection-multiplicity
count, not a distinct-point count; only existence is used there.

## 2. The genuinely new covariant

On `V!=0` put

```text
y=Vz,                  L=y^2+B y,
J_B=2VB'-BV'.                                               (2)
```

Multiplying `P_z` by `V` gives

```text
R_1=2L(2y+B)+VA.                                            (3)
```

Multiplying `P_x` by `V^3` gives the raw second gradient equation.  The
response identity and `(3)` reduce it exactly to

```text
R_2=V^3+V^2y+J_B yL.                                       (4)
```

Indeed, raw minus `(4)` is `(V'y/2)R_1`.  Thus `(3),(4)` generate the gradient
ideal after localizing at `V`.  Constant clutches have `J_B=-bV'`; the affine
deformation changes the universal resultant precisely through the covariant
`J_B`.

Exact symbolic elimination gives

```text
Res_y(R_1,R_2)=V^3 K_B,                                    (5)
```

where

```text
K_B=
 -A^3J^3+12A^2J^2V^2-4AB^3J^2V+4AB^2JV^2
 +24ABJV^3-48AJV^4-16AV^4-8B^4JV^2+8B^3JV^3
 +32B^2V^4-96BV^5+64V^6.                                  (6)
```

Here and below `J=J_B`.

## 3. Finite `g`-fibre clutch

At a root `alpha` of `g`, THM-3212 gives

```text
V(alpha)=A(alpha)=0,                 A'(alpha)!=0.          (7)
```

If `B(alpha)!=0`, then `P_z(alpha,z)=2B(alpha)^2z`, whose
unique zero `z=0` has `P_x(alpha,0)=1`; this fibre has no critical point.

If `B(alpha)=0`, the squared term itself vanishes identically on the fibre and

```text
P_z(alpha,z)=0,                     P_x(alpha,z)=A'(alpha)z+1.
```

Hence

```text
(alpha,-1/A'(alpha))                                      (8)
```

is a critical point.  Thus the finite clutch passes if and only if
`gcd(B,g)=1`.  The hostile `B=1-x` meets the `x=1` root and realizes `(8)`.

This identifies the first failed repair mechanism: making `B` nonconstant
does not permit it to vanish harmlessly on the old owner divisor.  Every such
collision reinstates a finite critical point immediately.

## 4. Passport saturation and uniform degree

Let

```text
boundary_(4111)=S^3T^8x^9,
boundary_(3211)=S^3T^8x^6(x-1)^3.                         (9)
```

These both have degree `44`.  At a root of `T`, write

```text
V=v t^m+...,                 A=a t+...,
a=2/(2-m),                   m in {3,4,5,6}.              (10)
```

For `B=b+e t`, direct substitution in `(6)` gives order at least `3m-1`, and

```text
[t^(3m-1)]K_B
 =16m(m-1)/(m-2) b^5v^3.                                  (11)
```

The exponents `3m-1` are exactly those encoded in `(9)`.  At the simple root
of `S`, the response identity cancels the apparent order-two term and gives
order at least three.  Therefore `(9)` divides `(6)` for every affine `B`:

```text
K_B=boundary_i H_(i;B).                                   (12)
```

Since `deg(V)=16`, `deg(A)=8`, `deg(B)<=1`, and `deg(J_B)<=16`, the unique
degree-`96` term in `(6)` is `64V^6`.  Hence, uniformly in both affine
parameters,

```text
deg_x H_(i;B)=96-44=52.                                   (13)
```

The exact companion independently expands the subfamily `B=1+t x` in `t`.
All six parameter coefficients divide `(9)` in both accessory fields; their
residual `x`-degrees are

```text
(52,44,36,28,23,15).                                     (14)
```

## 5. No `T` escape after the finite clutch

Assume now `gcd(B,g)=1`.  Then `b=B(alpha)!=0` at every root of `T`, so `(11)`
is nonzero.  Consequently

```text
H_(i;B)(alpha)!=0                                         (15)
```

at all four roots of `T`.  Thus every post-clutch boundary root must lie at
the single simple `S`-root.

## 6. The `S` wall and the two-jet escape cap

Let `s` be the root of `S`, use `t=x-s`, and write

```text
V=v_1t+v_2t^2+v_3t^3+v_4t^4+...,
B=b+e t,                   b=B(s)!=0.                     (16)
```

The response identity fixes the corresponding jet of `A`.  Substitution in
`(6)` gives

```text
[t^3]K_B=-(8/3)b^4v_1^2 W,
W=4bv_2-6ev_1+3v_1^2.                                    (17)
```

If `W!=0`, then `H_(i;B)` does not meet `S`.  On the unique escape wall
`W=0`, solve

```text
e=(4bv_2+3v_1^2)/(6v_1).                                 (18)
```

The next two coefficients are

```text
[t^4]K_B=(16/45)b^2v_1 Q_4(b),
[t^5]K_B=-(8/945)Q_5(b),                                 (19)
```

with

```text
Q_4=b^3(-54v_1v_3+8v_2^2)-30b^2v_1^2v_2+45v_1^3,

Q_5=b^5(3240v_1^2v_4+11016v_1v_2v_3-1552v_2^3)
   +b^4(3969v_1^3v_3+6132v_1^2v_2^2)
   +b^3(1260v_1^4v_2-4536v_1^2v_3+672v_1v_2^2)
   -12600b^2v_1^3v_2+945bv_1^5+3780v_1^4.                (20)
```

Exact Euclidean division in each cubic accessory field gives the identical
PRS degree profile

```text
(5,3) -> (3,2) -> (2,1) -> (1,0),
gcd(Q_4,Q_5)=1.                                           (21)
```

Thus `(19)` cannot vanish simultaneously.  Since the passport boundary
contains exactly `S^3`,

```text
ord_S H_(i;B)<=2.                                         (22)
```

The structured subfamily `B=1+t x` has one exact `S`-wall value in each
accessory field.  It remains a unit modulo `g`, and direct exact division gives
`ord_S H=1`, leaving `51` units away.  This is the hostile control showing that
the wall is genuine rather than an algebraic artifact.

## 7. Surviving critical points

If `B` meets `g`, `(8)` is an explicit critical point.  Otherwise `(13)`,
`(15)`, and `(22)` leave at least

```text
52-2=50                                                    (23)
```

units of resultant multiplicity at roots `x_0` with `V(x_0)!=0`.
Because the leading `y`-coefficient of `R_1` is the constant `4`, each such
resultant root supports an affine common zero of `(3),(4)`.  By the exact
gradient-ideal reduction, this is a critical point of `P_B`.

At any critical point of `P_B`, `Jac(P_B,Q)` vanishes for every polynomial
`Q`; hence no constant-Jacobian mate exists.

The good reductions

```text
(passport,p,u)=(4111,113,85),             (3211,101,64)
```

show that the rational genuinely nonconstant control `B=1+x` has `H`
squarefree and disjoint from `g` in both accessory fields.  It therefore has
exactly `52` reduced Morse critical points on all four covers.  The same
nonvanishing conditions define a nonempty Zariski-open affine-`B` locus.

## 8. Failure anatomy and strongest survivor

The bold reframe was that the new covariant `J_B=2VB'-BV'` might change the
infinity balance enough to remove the constant-clutch resultant.  It does
change the local `S` wall, but the first hoped-for implication fails:

```text
nonconstant clutch -> degree loss of the saturated resultant              FALSE.
```

The unique `64V^6` term freezes the residual degree at `52`; affine motion can
move at most two units through `S`.  The strongest surviving repair principle
is therefore:

> A source deformation must change more than the affine `B` covariant.  It
> must alter the degree-`96` bulk term or create a new boundary channel capable
> of absorbing the residual degree.

The next minimal lawful probes are a nonconstant `C_0`, a nonlinear `E_0`, or
quadratic `B`; each changes `(4)` or the infinity degree ledger.  No point in
the present family supplies the complete marked inverse pair required by the
THM-3172 owner filtration.

## 9. Exact reproduction and scope

Run

```text
python3 04-computation/jc_heptic_affine_B_source_obstruction_thm3225.py
python3 -O 04-computation/jc_heptic_affine_B_source_obstruction_thm3225.py
```

and compare LF-normalized bytes with the declared output.  The companion is
self-contained and uses exact rational, polynomial, resultant, gcd, and
accessory-algebra arithmetic.  Finite-field reductions are used only as
independent nonvanishing and squarefreeness certificates.  Every logical gate
uses an explicit exception, so optimized execution retains it; normal,
optimized, and stored transcripts agree byte-for-byte.

The exceptional lower bound `>=50` is resultant/intersection multiplicity,
not a distinct-Morse-point count.  Distinctness and Morse reduction are
asserted only on the audited nonempty open locus containing `B=1+x`.
No member supplies a marked inverse pair or Keller mate, and no claim proves
`JC(2)` or `DC(2)`.

**QED.**
