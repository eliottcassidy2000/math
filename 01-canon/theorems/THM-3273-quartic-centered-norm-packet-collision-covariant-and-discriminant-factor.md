---
id: THM-3273
title: "Quartic centered-norm packet collision covariant and discriminant factor"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For a depressed separable quartic
  f=X^4+pX^2+qX+r and its THM-3271 centered-norm packet
  N(X)=-(20X^3+18pX+27q)/27, the packet characteristic polynomial satisfies
  disc(P_N)=(2^12/3^36) C(p,q,r)^2 disc(f), where
  C=9p^6-1900p^4r+1000p^3q^2+110000p^2r^2-1000000r^3.
  Directly, each packet-value difference is the corresponding root
  difference times a quadratic pair quotient, and the product of all six
  quotients is -(2^6/3^18)C.  Hence for separable f, C=0 is exactly the new
  sheet-collision divisor introduced by scalarizing the algebra-valued
  packet.  After depression C is an affine relative covariant of weight
  twelve, not an absolute invariant in arbitrary coefficients.  On the tame
  special fibre with roots (a,m,m,m), C=27(a-m)^12/8, a unit when 2,3 and
  a-m are units; this explains THM-3272's fixed/moving separation.  No
  Keller cofactor, forbidden collision, C3/S4 exclusion, or JC(2) theorem
  follows.
source: jc-centered-norm-collision-covariant-2026-08-02
audit: >
  An independent audit checked the pair quotient, resultant discriminant
  convention, affine weights, tame special fibre, and separable collision
  equivalence.  Fresh normal and optimized symbolic replays byte-match the
  archived transcript.
depends_on:
  - THM-3271-universal-quartic-centered-norm-packet-and-local-singleton-projector
  - THM-3272-tame-pure-c3-centered-norm-packet-integral-fixed-projector
related:
  - THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary
  - THM-2993-quartic-signed-edge-star-triangle-cube-and-derivative-square-reembedding
  - THM-3064-pointed-cubic-norm-keller-decoder-and-inverse-different-boundary
script: 04-computation/jc_quartic_centered_norm_collision_covariant_thm3273.py
output: 05-knowledge/results/jc_quartic_centered_norm_collision_covariant_thm3273.out
script_sha256: f7e9af4c666716d7f1d9bd1db0b84762fcf2c9f8659c36c5acdf7e742d0b7785
output_sha256: ac5023086e2b2fd68fdfb847029a6c8ef483d1771badc1f7bf971f849ec8cad3
hash_basis: LF-normalized bytes
---

# THM-3273 -- one covariant measures every new centered-packet collision

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Packet map in depressed coordinates

Let `K` be a characteristic-zero field and take a monic depressed quartic

```text
f(X)=X^4+pX^2+qX+r.                                      (1)
```

The universal packet from
[THM-3271](THM-3271-universal-quartic-centered-norm-packet-and-local-singleton-projector.md)
is the cubic Tschirnhaus map

```text
N(X)=-(20X^3+18pX+27q)/27.                               (2)
```

If `xi=X mod f`, put

```text
P_N(Z)=Norm_(K[X]/(f)/K)(Z-N(xi))
      =Res_X(f(X),Z-N(X)).                                (3)
```

This is the scalar characteristic polynomial of the four-component
algebra-valued packet.  With the monic discriminant convention

```text
disc(product_i(Z-z_i))=product_(i<j)(z_i-z_j)^2,          (4)
```

the question is exactly when two packet values collide.

## 2. Pairwise factor and the collision covariant

For any two variables `u,v`, equation `(2)` gives

```text
N(u)-N(v)
 =(u-v) Q_p(u,v),                                        (5)

Q_p(u,v)=-(20(u^2+uv+v^2)+18p)/27.                       (6)
```

Let `alpha_1,...,alpha_4` be the roots of `f`.  Define

```text
C(p,q,r)=9p^6-1900p^4r+1000p^3q^2
         +110000p^2r^2-1000000r^3.                       (7)
```

Exact symmetric reduction of the six pair factors gives

```text
product_(i<j) Q_p(alpha_i,alpha_j)
   =-(2^6/3^18) C(p,q,r).                                (8)
```

This identity can be checked without choosing radicals: impose
`alpha_1+...+alpha_4=0`, eliminate `alpha_4`, substitute the elementary
symmetric functions for `p,q,r`, and expand.  It is also forced up to sign
by the discriminant calculation below; one labelled root control fixes the
minus sign.

Equation `(8)` identifies the precise information loss.  The ordinary root
difference product records collisions already present in the quartic.  The
new factor `C` records collisions introduced by applying the cubic packet
map before forgetting sheet labels.

## 3. Exact discriminant factorization

Multiplying the squares of `(5)` over all six pairs yields

```text
disc_Z(P_N)
 = (2^12/3^36) C(p,q,r)^2 disc_X(f).                     (9)
```

For reference,

```text
disc_X(f)=16p^4r-4p^3q^2-128p^2r^2
          +144pq^2r-27q^4+256r^3.                        (10)
```

Direct resultant expansion independently reproduces `(9)`.  In an
integer-coefficient normalization,

```text
3^12 P_N(Z)=
  531441 Z^4
 +944784 q Z^3
 +(2916p^3+174960pr+520506q^2) Z^2
 +(3888p^3q+25920pqr+116424q^3) Z
 +1296p^4r+972p^3q^2-28800p^2r^2
 -5040pq^2r+9261q^4+160000r^3.                           (11)
```

If `f` is separable, `(10)` is nonzero.  Therefore `(9)` proves the exact
equivalence

```text
two distinct quartic sheets have the same centered-norm value
  iff C(p,q,r)=0.                                         (12)
```

Thus `C=0` is not merely a sufficient warning or a discriminant shadow: it
is the complete new collision divisor in the separable depressed chart.

## 4. Affine covariance, not absolute invariance

For a general monic quartic, first translate by one quarter of the cubic
coefficient to obtain `(1)`.  This depression is canonical in
characteristic zero.  Under an affine root change

```text
X |-> uX+v,                       u!=0,                    (13)
```

the depressed coefficients transform as

```text
(p,q,r) |-> (u^2p,u^3q,u^4r).                             (14)
```

Every term in `(7)` has weight twelve, so

```text
C(u^2p,u^3q,u^4r)=u^12 C(p,q,r).                          (15)
```

Meanwhile quartic root differences scale by `u`, giving weight twelve for
`disc(f)`, and packet values scale by `u^3`, giving weight thirty-six for
`disc(P_N)`.  The weights in `(9)` are therefore exactly compatible:

```text
36=2*12+12.                                               (16)
```

The phrase **collision covariant** is intentional.  Formula `(7)` is a
relative invariant after the distinguished affine depression, not an
absolute binary-quartic invariant polynomial in arbitrary coordinates.

## 5. Tame pure-`C3` special fibre

Consider a special fibre whose roots are

```text
(a,m,m,m),                         d=a-m.                 (17)
```

Translate `m` to zero and depress the quartic `X^3(X-d)`.  Its coefficients
are

```text
p=-3d^2/8,              q=-d^3/8,              r=-3d^4/256. (18)
```

Substitution in `(7)` gives

```text
C=27d^12/8.                                               (19)
```

Hence `C` is a unit whenever `2,3,d` are units.  The quartic discriminant
vanishes on `(17)` because the three moving roots coalesce, but `C` does
not: the packet has exactly the inherited moving triple collision and no
additional fixed/moving collision.  This is the covariant explanation of
[THM-3272](THM-3272-tame-pure-c3-centered-norm-packet-integral-fixed-projector.md):
the fixed packet value remains separated, so its spectral denominator is a
unit.

The unit conclusion fails precisely on the already isolated boundaries:
residue characteristic two or three, or `d` nonunit.

## 6. Sharp scalar-packet hostile

Take

```text
f(X)=(X^2-9)(X^2-1)=X^4-10X^2+9.                         (20)
```

It is separable, but

```text
C(-10,0,9)=0.                                             (21)
```

The packet polynomial is

```text
P_N(Z)=Z^2(27Z-160)(27Z+160)/729,                         (22)
```

so the sheets `3` and `-3` collide at packet value zero.  This is the
minimal warning from THM-3271, now recognized as a point of the exact
covariant divisor.

For a noncollision control, the depressed quartic with roots
`0,1,2,-3` has `(p,q,r)=(-7,6,0)` and

```text
C=-11289159 !=0.                                          (23)
```

Both its quartic and packet discriminants are nonzero.

## 7. What this changes and what remains

The scalarization ledger is now exact:

```text
quartic algebra with tautological root
  --N_f--> algebra-valued centered packet
  --Norm characteristic polynomial--> scalar packet multiset

inherited branch collision: disc(f)=0;
new value collision:         C=0.                         (24)
```

At a good tame pure-`C3` place, `C` is a unit and scalarization does not lose
the fixed component.  Away from that lane, retaining only `P_N` requires an
explicit `C!=0` check; the full algebra-valued packet remains the safer
carrier.

The covariant contains no primitive-element chain-rule cofactor and imposes
no forbidden sign on a Keller map.  It does not prove that `C` is globally
nonzero on an arbitrary graph quartic, that a fixed boundary trace lies in
the affine source, or that the recovered THM-3230 cubeclass is trivial.

## 8. Exact companion

Run

```bash
python3 04-computation/jc_quartic_centered_norm_collision_covariant_thm3273.py
python3 -O 04-computation/jc_quartic_centered_norm_collision_covariant_thm3273.py
```

Both modes byte-match the stored transcript.  The companion derives the
resultant and both discriminants symbolically, checks all six pair quotient
identities and their product `(8)`, verifies affine weight twelve, proves
the special-fibre value `(19)`, and replays the hostile and positive controls
`(20)--(23)`.  Every truth-bearing check uses an explicit exception and
remains live under optimized execution.

No quartic Keller realization, cofactor/inverse-different compatibility,
pure-`C3` exclusion, `A4/S4` exclusion, `JC(2)`, or `DC(2)` follows.

**QED.**
