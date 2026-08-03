---
id: THM-3279
title: "Affine transverse clutch critical no-go"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  For either THM-3212 cubic accessory response pair, every affine transverse
  clutch C_0 in the displayed B=1, constant-E_0 family leaves a critical
  point.  A nonconstant clutch meeting the owner divisor gives an explicit
  point; otherwise at least 50 units of saturated critical-resultant
  multiplicity remain off that divisor.
source: root/creative-synthesis/2026-08-03
audit: >
  The primary exact companion derives the localization, universal resultant,
  affine degree ledger, four T rows, two S walls and both sharp field
  controls.  Normal, optimized and stored outputs agree.  A fresh independent
  hostile audit uses a literal 6-by-6 Sylvester determinant instead of
  sympy.resultant, solves the response jet recursively rather than importing
  it, and implements exact Q[u]/(q) arithmetic from three rational
  coordinates.  It independently reproduces both characteristic-zero monic
  residual digests, the unit controls, degrees 96/52, exact S order two and
  no T escape.  Its normal, optimized and stored outputs agree; both sources
  have zero assertion nodes and zero floating literals.  The audit pins the
  primary script/output and inherited THM-3212 artifacts but deliberately
  does not pin this theorem file, avoiding a promotion hash cycle.
depends_on:
  - THM-3212-centered-heptic-source-morse-obstruction-and-offcenter-clutch
related:
  - THM-3225-affine-jacobian-clutch-resultant-and-two-boundary-no-escape
  - THM-3276-degree-at-most-eight-polynomial-clutch-critical-no-go
script: 04-computation/jc_affine_transverse_c0_clutch_no_go_thm3279.py
output: 05-knowledge/results/jc_affine_transverse_c0_clutch_no_go_thm3279.out
script_sha256: 06820b2476fc0f2cefe3982d054a7db09bb88b4892503580550a1b154564508a
output_sha256: 4a88b5ab31eed4c9a5f90f814a6a24a614db0afecc4cf1cab7fa32dae7c991c4
audit_script: 04-computation/jc_affine_transverse_c0_clutch_no_go_thm3279_independent_audit.py
audit_output: 05-knowledge/results/jc_affine_transverse_c0_clutch_no_go_thm3279_independent_audit.out
audit_script_sha256: dc345ef5f02fc922aea19931da005556f52a850638cd964ddd6b51854e5fd621
audit_output_sha256: ee554ed2f61e2bb5e5abe02db53ee029047982ec2f021dfb0e8f139da1edb4b4
hash_basis: LF-normalized bytes
---

# THM-3279 -- affine transverse clutch critical no-go

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The primary and independent companions complete normal/optimized/stored
audits through separate elimination and cubic-field implementations.

## 1. Statement

Let `K_i` be either cubic accessory field of
[THM-3212](THM-3212-centered-heptic-source-morse-obstruction-and-offcenter-clutch.md),
embedded in an algebraically closed characteristic-zero field `K_0`.  Retain
the response pair and owner divisor

```text
V=4SDT^2/Gamma^2,       A=2SET/Gamma,
g=ST,                   2VA'-AV'=2V,
deg V=16,               deg A=8.                         (1)
```

For every affine `C in K_0[x]` and constant `e in K_0`, put

```text
P_C(x,z)=(V(x)z^2+z+C(x))^2+A(x)z+e.                    (2)
```

Then `P_C` has a critical point.  More precisely, if `C` is nonconstant,
either:

1. `gcd(C,g)!=1`, and `(2)` has an explicit critical point over the common
   root; or
2. `gcd(C,g)=1`, and at least `50` units of saturated critical-resultant
   multiplicity lie at roots where `V!=0`.

The constant-`C` lane is already the `b=1,k=0` constant clutch excluded by
THM-3212.  At any critical point of `P_C`, `Jac(P_C,Q)=0` for every polynomial
`Q`.  Thus no member of `(2)` is one coordinate of a polynomial Keller pair.

## 2. Transverse gradient reduction

On `V!=0`, set

```text
y=Vz,                 L=y^2+y+CV.                       (3)
```

Multiplying the `z` derivative by `V` gives

```text
R_1=2L(2y+1)+VA.                                         (4)
```

Multiplying the fixed-`z` `x` derivative by `V^3`, then using `(1)` and
subtracting `(V'y/2)R_1`, gives

```text
R_2=V^2y+L(-V'y+2V^2C').                                (5)
```

The sign in this reduction is checked symbolically: raw minus `R_2` is
`(V'y/2)R_1`.  Thus `(4),(5)` generate the localized gradient ideal.

There is a useful universal calculation.  For

```text
L=y^2+By+CV,
R_1=2L(2y+B)+VA,
R_2=V^3e+V^2y+L(Jy+K),                                  (6)
```

exact elimination gives

```text
Res_y(R_1,R_2)=V^2 F,                                   (7)
```

where `F` has `40` monomials.  Substituting

```text
B=1,        J=-V',        K=2V^2C',        e=0          (8)
```

introduces one further factor:

```text
Res_y(R_1,R_2)=V^3 K_C.                                 (9)
```

Writing `d=C'`, the exact twenty-term factor is

```text
K_C=
 A^3(V')^3+8A^2CV^2(V')^2d+2A^2C(V')^3+32A^2V^5d^3
 +24A^2V^3V'd^2+24A^2V^3V'd+4A^2V(V')^2d
 +128ACV^5d^2+32ACV^3V'd+16ACV^3V'
 -32AV^4d^2-48AV^4d-16AV^4-8AV^2V'd-4AV^2V'
 +128C^2V^5d+32C^2V^3V'-32CV^4d-32CV^4-8CV^2V'.       (10)
```

## 3. Finite owner collision and infinity degree

Suppose first that nonconstant `C` meets `g` at `alpha`.  THM-3212 gives

```text
V(alpha)=A(alpha)=0,       C(alpha)=0.                  (11)
```

At `(alpha,0)`, both derivatives of `(2)` vanish.  Hence every such collision
is an explicit critical point.

Now assume `gcd(C,g)=1`, so `d!=0`.  In `(10)`, affine-degree bookkeeping
shows that

```text
32 A^2 V^5 d^3                                             (12)
```

is the unique degree-`96` term.  It cannot cancel.  Therefore

```text
deg K_C=96.                                               (13)
```

For the two passports define

```text
boundary_(4111)=S^3T^8x^9,
boundary_(3211)=S^3T^8x^6(x-1)^3.                       (14)
```

Both boundaries have degree `44`.

## 4. No escape through the four T branches

At a root of `T`, write

```text
V=v t^m+...,       A=2t/(2-m)+...,       C=c+dt,
m in {3,4,5,6}.                                         (15)
```

Substitution in `(10)` gives

```text
ord_t K_C=3m-1,
[t^(3m-1)]K_C=
  32m(m-1)/(m-2)^2 c v^3.                               (16)
```

The exponents are respectively `8,11,14,17`, exactly those encoded by
`(14)`.  Since `C` is a unit on `g`, `c!=0`; no further `T` multiplicity is
available.

At the simple `S` root the response identity supplies the universal `S^3`
factor.  Consequently

```text
K_C=boundary_i H_C,              deg H_C=52.             (17)
```

## 5. The two S walls and the untunable fifth jet

Use `t=S` and write

```text
V=v_1t+v_2t^2+v_3t^3+v_4t^4+...,
C=c+dt,                    c!=0.                         (18)
```

The response equation determines the jet of `A`.  Direct expansion of
`(10)` first gives

```text
[t^3]K_C=(32/3)c v_1^2(3c v_1^2+2v_2).                 (19)
```

Thus extra contact requires

```text
c=-2v_2/(3v_1^2).                                       (20)
```

On `(20)`, the next coefficient is

```text
[t^4]K_C=
 -64v_2/(135v_1)
 (90dv_1^3-45v_1^3+54v_1v_3-28v_2^2).                 (21)
```

The only second wall is therefore

```text
d=(45v_1^3-54v_1v_3+28v_2^2)/(90v_1^3).               (22)
```

After imposing both walls,

```text
[t^5]K_C=-32v_2 R_S/(315v_1^2),                         (23)

R_S=315v_1^3v_2+360v_1^2v_4-414v_1v_2v_3+148v_2^3.    (24)
```

Exact arithmetic in both cubic accessory fields proves

```text
v_1 !=0,       v_2 !=0,       R_S !=0.                  (25)
```

Hence `ord_S K_C<=5`, or equivalently

```text
ord_S H_C<=2.                                           (26)
```

Combining `(16),(17),(26)`, at least

```text
52-2=50                                                  (27)
```

units of `H_C` remain away from `g`.  Since `(4)` has constant leading
`y`-coefficient `4`, each off-owner resultant root supports an affine common
zero of `(4),(5)`, hence a critical point of `(2)`.  The count in `(27)` is
with multiplicity; only nonemptiness is used.

The bound is sharp as an owner-contact invoice.  In both fields the affine
control `C=c_wall+d_wall S` is a unit modulo `ST`, has

```text
(deg K_C,deg H_C)=(96,52),       ord_S H_C=2,
gcd(H_C,T)=1.                                             (28)
```

The monic residual digests are

```text
4111  c07d66a389368e57ff32cdcc1c10f134951a0f2008ed65b707033d8ea8844e8e
3211  41228f21741350a603d05efd59a11ec6157d5d0fb17b6b77f03f360cd81e78ad
```

## 6. Failure boundary and scope

The first failed repair is precise:

```text
two affine C-jets absorb all 52 residual units                 FALSE.
```

The value and slope of `C` buy exactly the walls `(20),(22)`; the field
invariant `(24)` is the first untunable coefficient.  This is transverse to
THM-3276, which varies `B` with `C=0`; it does not improve or weaken that
theorem.

The theorem concerns only the explicit family `(2)`.  It supplies no
polynomial second coordinate, marked inverse cover, branchwise cofactor,
Jelonek component, or classification of simultaneous `B,C,E_0`
deformations.  It proves neither `JC(2)` nor `DC(2)`.

## 7. Exact reproduction

Run

```text
python3 04-computation/jc_affine_transverse_c0_clutch_no_go_thm3279.py
python3 -O 04-computation/jc_affine_transverse_c0_clutch_no_go_thm3279.py
```

and compare LF-normalized bytes with the declared output.  The companion
pins four THM-3212 artifacts, derives the universal resultant and localization
identity, checks every local coefficient symbolically, rebuilds both cubic
fields, and verifies the sharp controls.  It has no assertion node, floating
literal, randomness, or fitted recurrence.

For the independent replay run

```text
python3 04-computation/jc_affine_transverse_c0_clutch_no_go_thm3279_independent_audit.py
python3 -O 04-computation/jc_affine_transverse_c0_clutch_no_go_thm3279_independent_audit.py
```

and compare with its declared audit output.  This companion derives the
resultant as a Sylvester determinant and uses its own exact cubic quotient
arithmetic; it imports neither the primary critical-factor expression nor
its algebraic-number-field objects.

QED.
