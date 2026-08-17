---
id: THM-3537
title: "Fixed Keller level-two old-L inertia and literal-x transverse index"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  At depth two, the old
  Jelonek prime L has coefficient four in the normalization discriminant and
  generic geometric inertia cycle type (4)(2)(1)^3.  On the canonical
  transverse DVR (a,b,c)=(2/27+t,1,1), the literal x_2 power order has
  discriminant exponent eight and exact local index length two.  This does
  not determine the generic literal-coordinate index along L or an all-level
  old-prime recurrence.
source: codex/old-L-inertia/2026-08-16
depends_on:
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
  - THM-2582-odd-block-discriminant-tower-and-composite-jelonek-square-class
  - THM-3529-fixed-keller-complete-packet-finite-sheet-unit
  - THM-3531-fixed-keller-intrinsic-all-level-discriminant-square-class
  - THM-3533-fixed-keller-newest-prime-reduced-different-and-index-square
related:
  - THM-3535-fixed-keller-full-wreath-and-all-level-linear-primitivity
  - THM-3536-berggren-angle-languages-signed-c4-and-harmonic-density
  - HYP-9033-discriminant-tower-and-genus-axis-of-the-keller-monoid
script: 04-computation/keller_level_two_old_L_newton_index_audit_20260816.py
output: 05-knowledge/results/keller_level_two_old_L_newton_index_audit_20260816.out
script_sha256: d2eb44fc2ad5c270ebddf7197ad2d3e3c905b7a9dad86fa80b83735551f55f5d
output_sha256: f159ba5453fc8e85c0f7d814fceae1a073ef41e6d4e0328c766c543b2c631fc0
semantic_sha256: 02cbf6a236a84ca153bea6fa60fb3c1267593b9e0694c0cddc5b1d7dc5bd8f59
hash_basis: LF-normalized bytes
---

# THM-3537 -- the first old Keller prime has a quartic-quadratic inertia lift

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Let `F` be the fixed sporadic Keller map and retain the inverse tower

```text
K_0 subset K_1 subset K_2,             [K_i:K_(i-1)]=3.       (1)
```

Let `B_i` be the normalization of `A=Q[a,b,c]` in `K_i`, and put

```text
P_0=L=27a^2c^2-18abc+16a+b^3c-b^2,
P_1=H.                                                       (2)
```

For a height-one prime `p`, write

```text
delta_i(p)=v_p(Disc(B_i/A)).                                (3)
```

## 1. The old-prime theorem

At the old prime `p_0=(L)`,

```text
delta_2(p_0)=4.                                             (4)
```

Over a strict henselization at the geometric generic point of `V(L)`, the
five inertia orbits on the nine depth-two leaves have lengths

```text
(4,2,1,1,1).                                               (5)
```

Equivalently, a tame inertia generator has cycle type

```text
(4)(2)(1)^3                                                (6)
```

and permutation different exponent

```text
(4-1)+(2-1)=4.                                             (7)
```

Thus the old component is not merely an even-parity passenger in the
degree-nine discriminant.  Its exact normalization multiplicity is four.

## 2. The parity squeeze proves the multiplicity without a degree-nine expansion

THM-3533 gives the one-step local decomposition above `p_0`.  There is one
regular finite branch `q_fin` with `e=1` and one escaping quadratic branch
`q_esc` with `e=2`; the latter has tame different exponent one.  Hence

```text
delta_1(p_0)=1.                                            (8)
```

The finite branch contributes no ramification in the next inverse step.
Indeed, THM-3529 identifies it with the explicit finite divisor and proves
that `L(q_fin)` is a unit.  Off `L=0`, the one-step inverse cover is finite
etale.  The independent audit also checked the canonical point directly:

```text
L(2,5/6,-7/8)=241465/1728 != 0.                          (8a)
```

Only the rank-three algebra above `q_esc` can contribute a new relative
different.  In residue characteristic zero it is tame, so its total
permutation different exponent `c` satisfies

```text
0 <= c <= 2.                                               (9)
```

Discriminant transitivity now reads

```text
delta_2(p_0)=3 delta_1(p_0)+c=3+c.                        (10)
```

THM-2582, equivalently the all-level class in THM-3531, says that `L` is an
even divisor of the depth-two discriminant: the unique odd class is the new
prime `H`.  Therefore `(10)` is even.  Among `3,4,5`, only `4` is even, so

```text
c=1,                 delta_2(p_0)=4.                     (11)
```

This small interval plus parity is the mechanism.  Square class alone would
not determine a positive multiplicity, and transitivity alone would leave
three possibilities.

Because a tame rank-three contribution of exponent one has type `(2,1)`,
the escaping `e=2` branch lifts to one branch of total ramification index
`2*2=4` and one of index `2*1=2`.  The finite branch has three unramified
children after strict henselization.  This gives exactly `(5)`.

## 3. Exact canonical transverse model

The cycle type has an independent literal-coordinate realization.  Use the
smooth transverse line

```text
(a,b,c)=(2/27+t,1,1),
L=t(27t+2),                 dL/dt at t=0 is 2.             (12)
```

The point avoids the new divisor:

```text
H(2/27,1,1)=61815040/27 != 0.                             (13)
```

Let `E_2(X)` be the saturated degree-nine eliminant for the literal first
coordinate `x_2` of a second preimage.  The exact companion reconstructs it
from the two inverse cubics and the rational inverse section; it does not
load a stored eliminant.  Its leading coefficient is a `t`-adic unit.

For the coefficient of `X^i`, record its `t`-valuation and initial
coefficient.  In increasing `i` the exact ledger is

```text
i :  0       1        2      3       4          5
v :  2       2        1      1       1          1
in: -1792  -58304   -3584   1664   248320/9  -9931808/27

i :  6        7        8       9
v :  0        0       +inf     0
in: -28672  -101376    0     -61815040/27.               (14)
```

The lower Newton hull is

```text
(0,2) -> (2,1) -> (6,0) -> (9,0),                       (15)
```

with slopes

```text
-1/2,              -1/4,              0.                 (16)
```

Writing `Y` for the residual variable, the three residual polynomials are

```text
-1792(2Y+1),
-3584(8Y+1),
-(256/27)(241465Y^3+10692Y+3024).                        (17)
```

All three are squarefree; the last has discriminant

```text
-3398871051182646293954560/27 != 0.                      (18)
```

Newton regularity therefore gives root-valuation/multiplicity pairs

```text
(1/2,2),              (1/4,4),              (0,3),       (19)
```

and, after geometric residue-field closure, the same cycles `(4,2,1,1,1)`
as the parity squeeze.  This is an exact second route to `(5)` at a clean
transverse point.

## 4. The literal `x_2` order pays an index tax of two

The exact eliminant discriminant has

```text
v_t(Disc(E_2))=8.                                        (20)
```

Let `R=Q[t]_(t)`, let `S` be the normalization of the rank-nine separable
`R`-algebra defined by `E_2`, and normalize the unit leading coefficient so
that `R[x_2]` is its literal power order.  The Newton data give

```text
v_t(Disc(S/R))=4.                                        (21)
```

The order-index identity is

```text
Disc(R[x_2]/R)=Disc(S/R) Index_R(S:R[x_2])^2.            (22)
```

Consequently

```text
length_R(S/R[x_2])=(8-4)/2=2.                            (23)
```

Thus this prescribed primitive observation is not locally maximal on the
canonical old-`L` transversal.  THM-3535's primitivity and THM-3533's
index-square warning are both sharp: full degree does not remove the local
order defect.

The extra four units in `(20)` have a visible source.  The quadratic and
quartic packets both reduce to `x_2=0`.  Their `4*2=8` cross-pairs have
difference valuation `min(1/4,1/2)=1/4`; after squaring, they contribute

```text
2 * 8 * (1/4)=4.                                         (24)
```

The internal quartic and quadratic differents contribute the other
`3+1=4`.  Hence

```text
normalization ramification 4 + cross-packet collision 4
 = power-order discriminant 8.                           (25)
```

This is the exact old-prime analogue of the repo's recurring distinction
between an intrinsic carrier and a coordinate realization.

## 5. Equality boundary and hostile controls

Every sidecar above is load-bearing.

1. **Parity without the degree-three ceiling is insufficient.**  It says
   only that `(10)` is even.  The tame rank-three bound shrinks the choices
   to one value.
2. **The finite-sheet unit is essential.**  Without THM-3529, a second
   relative contribution above `q_fin` could occur and the squeeze would
   fail.
3. **Normalization and power order differ.**  Their exponents are `4` and
   `8`; identifying them would erase the index length two.
4. **The transverse coordinate result is not a generic-index theorem.**
   Equation `(23)` is proved for the named DVR base change.  It does not by
   itself compute `length_(A_(L))(B_(2,(L))/A_(L)[x_2])` at the three-variable
   generic point.
5. **No all-level recurrence is proved.**  The separate numerical scout is
   stable through depth five and reports exponents `1,4,14,46,142`.  The
   tempting formula `2*3^(n-1)-2^(n-1)`, which matches the first four rows,
   predicts `146` and is therefore refuted at depth five.  The
   `VERIFIED-NUMERICAL` wreath-product profiles survive as discovery data;
   their cycle lengths already show that a naive pure-doubling story is
   false.

The theorem does not compute old-prime multiplicities at depth three or
beyond, the newest-prime index of `x_n`, any `y_n/z_n` coordinate index, or
the singularities of the full discriminant divisor.  It concerns one fixed
three-dimensional Keller map and has no arbitrary-map, planar Jacobian,
LRC(14), tournament, or physical-current consequence.

## 6. Exact companion

Reproduce with

```text
python -B 04-computation/keller_level_two_old_L_newton_index_audit_20260816.py
python -B -O 04-computation/keller_level_two_old_L_newton_index_audit_20260816.py
```

The normal and optimized transcripts byte-match the stored output.  The
companion verifies `(12)`--`(23)`, pins the raw `H` artifact before using
`(13)`, uses explicit failures rather than executable assertions, and has
semantic SHA-256

```text
02cbf6a236a84ca153bea6fa60fb3c1267593b9e0694c0cddc5b1d7dc5bd8f59.
```

The independent hostile audit rederived the tame degree-three ceiling and
finite-sheet unit at `(8a)`, reconstructed the cycle type from both the tower
and Newton routes, replayed the script normally and with `-O`, verified all
three hashes and the pinned `H` artifact, and found zero executable `assert`
nodes.  No mathematical or scope defect was found.

**QED.**
