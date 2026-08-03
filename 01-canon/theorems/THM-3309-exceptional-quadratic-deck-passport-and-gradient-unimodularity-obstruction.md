---
id: THM-3309
title: "Exceptional quadratic deck passport and gradient-unimodularity obstruction"
status: >
  PROVED + VERIFIED-EXACT on the declared fixed slice.  For `C=c+x,d=k=1`
  in each of the two THM-3212 accessory fields, pull THM-3306's transverse
  degree-36 linear-row base locus to its nonsplit exceptional quadratic field
  `B_i=A_i[t]/(P2 t^2-Q2 t+R2)`.  The pointed deck restores a tautological
  relative root label, an anti-invariant inverse-different unit, and
  deck-equivariant first-normal motion.  The elimination pair `(U,W)` is a
  non-descending branch passport, with exact defect
  `U sigma(W)-sigma(U)W=-16 P2^2 F0prime`.  It is not a gradient Bezout row:
  direct pullback gives `P_x=P_z=0` in the nonzero field `B_i`, so the gradient
  ideal is proper, every `Jac(P,Q)` pulls back to zero, and the mate pipeline
  fails at gradient unimodularity before `mu(P)` is defined.  The result has
  relative degree two and total degree 72; it constructs no Keller cofactor,
  mate, inverse map, JC(2), or DC(2) consequence.
audit: >
  The exact companion internally reconstructs both accessory response fields,
  the degree-36 `linear_a`, the distinct degree-32 pre-reduction pivot `P2`,
  the quadratic deck arithmetic, the localized gradients, elimination pair,
  inverse different, and first-normal packet.  It pins THM-3306 and the two
  prior exact sidecars by hash but imports or executes none of them.  The true
  derivatives `P_x,P_z` are evaluated directly at `z=-t/V` independently of
  the triangular identities to `R1,R2`; both routes give zero.  Normal and
  optimized modes byte-match the frozen transcript, and the source has zero
  assertion nodes, zero floating literals, and no imported execution path.
source: root/creative-synthesis-recover/2026-08-03
depends_on:
  - THM-3212-centered-heptic-source-morse-obstruction-and-offcenter-clutch
  - THM-3306-affine-c-critical-section-square-discriminant-and-transverse-base-locus
related:
  - THM-3064-pointed-cubic-norm-keller-decoder-and-inverse-different-boundary
  - THM-3066-k4-initial-face-product-quotient-blind-to-keller-sheetwise-cofactor
script: 04-computation/jc_exceptional_quadratic_deck_cofactor_integrability_scout_20260803.py
output: 05-knowledge/results/jc_exceptional_quadratic_deck_cofactor_integrability_scout_20260803.out
script_sha256: c9809b5cafc5defca5aeda49bdb321f9c2d20b57405020c76bebc78e4e5dd2c6
output_sha256: 5026f0f17d0ff110d7f916a5b12858c04e09e2949209c7a75f6ab3e9cdc37b58
hash_basis: LF-normalized bytes
---

# THM-3309 -- exceptional quadratic deck passport and gradient-unimodularity obstruction

**PROVED + VERIFIED-EXACT on the declared fixed slice.**

The theorem distinguishes three notions that scalar elimination had put too
close together: elimination coefficients, an inverse-different generator,
and a primitive Keller cofactor.  Only the first two exist on this deck.  The
true gradient test fails before a Keller mate invariant can be formed.

## 1. Exceptional quadratic field

Work on `C=c+x,d=k=1` in either THM-3212 accessory field `K_i`.  THM-3306's
primitive linear PRS row has coefficient `linear_a(x)` of degree `36`.  At
its transverse base locus set

```text
A_i=K_i[x]/(linear_a).
```

The preceding quadratic PRS row has the distinct pre-reduction leading
coefficient

```text
P2=2V'+8V^2,                 deg(P2)=32,                     (1)
```

and exceptional polynomial

```text
F_0(t)=P2 t^2-Q2 t+R2.                                      (2)
```

The pinned exact nonsquare certificates for

```text
delta=Q2^2-4P2R2                                             (3)
```

make `(2)` irreducible and separable over `A_i`.  Hence

```text
B_i=A_i[t]/(F_0)
```

is a field of relative degree `2` and total degree `72` over `K_i`.  Put

```text
T=Q2/P2,               N=R2/P2,               sigma(t)=T-t. (4)
```

Then

```text
t+sigma(t)=T,              t sigma(t)=N,                    (5)
F_0(s)=P2(s-t)(s-sigma(t)).                                 (6)
```

Thus `(B_i,t)` restores a pointed relative root label, while forgetting the
pointing leaves one connected degree-two closed point.  After each of the
`36` geometric base embeddings there are two conjugate directions; there are
`72` geometric points over an algebraic closure of `K_i`, grouped into `36`
pairs.

Define

```text
F0prime=2P2t-Q2.                                             (7)
```

Exact reduction gives

```text
F0prime^2=delta,                sigma(F0prime)=-F0prime.     (8)
```

In particular, `F0prime` is a nonzero anti-invariant unit.

## 2. Elimination pair as branch passport

Let `R_1,R_2` be the localized gradient cubics, reserving `R2` in `(2)` for
the quadratic-row constant coefficient.  Their fraction-free linear
elimination identity uses

```text
Qelim=4P2 y+(6P2-4Q2),
U=P2^2-V'Qelim,
W=-4Qelim,                                                  (9)
```

and satisfies

```text
ell=U R_1+W R_2,
U-(V'/4)W=P2^2.                                             (10)
```

On the exceptional root `y=-t`, the invariant/anti-invariant decomposition
is

```text
Qelim=-24V^2-2F0prime,
W=96V^2+8F0prime,
U=P2^2+24V'V^2+2V'F0prime.                                 (11)
```

All three elements are nonzero and hence units in both fields `B_i`.  Their
traces and norms descend:

```text
Tr(Qelim)=-48V^2,
Tr(W)=192V^2,
Tr(U)=2P2^2+48V'V^2,                                       (12)

N(Qelim)=4(144V^4-delta),
N(W)=16N(Qelim),
N(U)=(P2^2+24V'V^2)^2-4(V')^2delta.                        (13)
```

The projective pair itself does not descend, because

```text
U sigma(W)-sigma(U)W=-16P2^2F0prime != 0.                  (14)
```

Thus `(U,W)` separates the two conjugate branches.  It is a branch passport,
not an invariant cofactor direction.

## 3. Inverse different is not a Keller cofactor

For the monic polynomial

```text
f_0(s)=s^2-Ts+N,
f_0'(t)=F0prime/P2,                                         (15)
```

one has

```text
Tr(f_0')=0,             N(f_0')=-delta/P2^2.                (16)
```

The inverse different therefore has the explicit unit generator

```text
eta=1/f_0'(t)=P2/F0prime,
sigma(eta)=-eta,        N(eta)=-P2^2/delta.                 (17)
```

This is ordinary etale primitive-element data.  A primitive Keller cofactor
would be chain-rule data of the form `q/f_0'`, with a numerator supplied by an
inverse or mate equation.  No such `q`, second coordinate, or inverse map is
present.  Replacing it by `1` is not a valid inference.

## 4. True gradient obstruction

For

```text
P(x,z)=(V(x)z^2+z+C(x))^2+A(x)z+E(x),                      (18)
```

put `y=Vz`.  The exact triangular change of gradient coordinates is

```text
R_1=V P_z,
R_2=V^3 P_x-(V'y/2)R_1.                                   (19)
```

Direct evaluation at `z=-t/V`, independent of `(19)`, gives

```text
P_x=P_z=0 in B_i.                                          (20)
```

The localized cubics also satisfy

```text
R_1(-t)=R_2(-t)=0,                                         (21)
```

and the unit `V` makes `(19)--(21)` a second exact verification of `(20)`.

Because `B_i` is a nonzero field, the evaluation homomorphism witnesses

```text
1 notin (P_x,P_z).                                         (22)
```

For every polynomial `Q`, regardless of its coefficients,

```text
Jac(P,Q)=P_xQ_z-P_zQ_x=0 in B_i.                           (23)
```

Equation `(10)` does not contradict this.  It says that the *coefficients*
`U,W` generate the unit ideal after `P2` is inverted; on the common root its
first line is only `0=0`.  It is not an identity of the form
`A P_x+B P_z=1`.

The mate-integrability class

```text
mu(P)=[A_x+B_z] in coker(P_x partial_z-P_z partial_x)
```

presupposes a gradient Bezout row `A P_x+B P_z=1`.  By `(22)` no such row
exists on this slice.  Therefore `mu(P)` is **undefined here**, rather than
uncomputed or proved nonzero.  Gradient unimodularity is the first failed
mate implication.

## 5. First-normal deck compatibility

In the fixed blow-up normalization write

```text
F(u,t)=F_0(t)+u(F12 t^2+F11 t+F10)+O(u^2),
dot(t)=-F_1(t)/F0prime.                                    (24)
```

If `T(u)=Q2(u)/P2(u)` denotes the trace of the monic quadratic, then

```text
dot(T)=(-P2 F11-Q2 F12)/P2^2.                              (25)
```

The exact quotient calculation in both fields gives

```text
Tr_B/A(dot(t))=dot(T),
sigma(dot(t))=dot(T)-dot(t).                               (26)
```

This is the derivative of `sigma_u(t(u))=T(u)-t(u)`: the two first-normal
velocities remain conjugate and distinct.  It is a local first-order result,
not global monodromy or a higher-normal classification.

## 6. Boundary

The proved source-to-target map is

```text
source:      transverse linear-row base locus plus quadratic PRS row;
target:      pointed etale deck (B_i,t);
map:         A_i -> B_i followed by y=-t;
preserved:   both critical roots, deck exchange, first-normal motion;
restored:    one tautological root coordinate after adjoining the sidecar;
destroyed:   descent of an individual branch and projective (U,W);
absent:      chain-rule numerator, gradient Bezout row, polynomial mate.    (27)
```

Only the two fixed accessory fields and `C=c+x,d=k=1` are covered.  The
degree-36 `linear_a` and degree-32 pre-reduction `P2` remain distinct, and the
field degrees are `2/72` as required by MISTAKE-362.  The theorem classifies
neither the degree-119 residual component nor deformations in `d,k,V,A`.  It
constructs no primitive Keller cofactor, mate, inverse map, JC(2), or DC(2)
consequence.
