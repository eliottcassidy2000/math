---
id: THM-3777
title: "Three-component tower reciprocal payment forces a horizontal divisor"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  In the all-degree tower A=1+xy, B=1+x^(2m+1)A^m,
  U=xAB, P=((m+1)B-1)/(mx), the reciprocal shear P_1=P-1/U pays the
  unique nonzero vertical residue and equals N_m/(AB).  The candidate proves
  that every nonconstant rational correction U+G(P_1) either retains a
  vertical pole or creates an explicit transverse horizontal pole
  N_m-rho AB=0.  If G is constant, a third shear can pay the two remaining
  vertical components only by restoring the original x-axis pole.  This
  closes that three-shear orientation uniformly in m; longer words remain
  open.
source: root + jc_sparse_direct_search / 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT by root.  The all-m formula for N_m, reduced
  denominator and vertical coefficients, exhaustive rational-G dichotomy,
  nonempty transverse horizontal divisor and genuine-pole argument, and
  compulsory restoration of the paid axis were rederived by hand.  Normal,
  optimized, and frozen exact transcripts agree byte for byte.
depends_on: []
related:
  - THM-3774-three-component-rational-keller-cover-tower
  - THM-3776-unequal-vertical-residue-three-target-shear-nonpolynomialization
script: 04-computation/jc2_three_component_horizontal_divisor_transfer_thm3777.py
output: 05-knowledge/results/jc2_three_component_horizontal_divisor_transfer_thm3777.out
script_sha256: 1825adbe74ec7f953064fed1f985088241242b61cc186d1440b2342e797b166a
output_sha256: 4da6b26146886b499e0579ea4b2fb783fe17c00227dd555f06902533cd9c5d86
semantic_sha256: b02db6cb6a9dd59f4b67313f2352f83bc1cca1d2a1392d4b42c025a69a4d16bb
hash_basis: raw LF bytes
---

# THM-3777 -- paying one vertical component creates a horizontal invoice

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  The proof is
self-contained and does not depend on THM-3774.

Let `k` be an algebraically closed field of characteristic zero and let
`m>=1` be an integer.  In `k[x,y]` define

```text
A=1+xy,
B=1+x^(2m+1)A^m,
U=xAB,
P=((m+1)B-1)/(mx).                                    (1)
```

Then `P in k(x,y)` and `U in k[x,y]` satisfy `J(P,U)=-1`.  The zero fibre of
`U` has the three pairwise coprime reduced factors `x,A,B`.  Its principal
coefficient packet for `P` is

```text
x=0: 1,                    A=0: 0,                    B=0: 0. (2)
```

Pay the unique nonzero coefficient by the rational target shear

```text
P_1=P-1/U.                                               (3)
```

For arbitrary `G,H in k(s)`, form the next two alternating target shears

```text
U_1=U+G(P_1),
P_2=P_1+H(U_1).                                         (4)
```

The theorem asserts:

1. if `G` is nonconstant, then `U_1` is not a polynomial.  More precisely,
   either it still has a pole on `A=0` and `B=0`, or it has a genuine pole on
   an explicit divisor transverse to `x=0`;
2. if `G` is constant, no choice of `H` makes `P_2` polynomial.

Consequently the opposite-start three-shear word `(3),(4)` cannot turn this
rational Keller seed into a polynomial pair for any `m`.  The conclusion is
scoped to this word and its prefixes.  Later shears might move the horizontal
invoice again and remain open.

## 1. The uniform rational Keller tower

Put

```text
t=x^2A.                                                  (5)
```

Then the definitions collapse to

```text
B=1+xt^m,
P=1/x+((m+1)/m)t^m,
U=t/x+t^(m+1).                                          (6)
```

Since `J_(x,y)(x,t)=x^3`, direct differentiation in `(x,t)` gives

```text
J_(x,t)(P,U)=-1/x^3,
J_(x,y)(P,U)=-1.                                        (7)
```

The same chart records the algebraic cover relation

```text
t^(m+1)-mPt+mU=0.                                       (8)
```

Neither `(7)` nor the divisor argument below assumes that this relation is
irreducible.

Let

```text
W=UP=A B((m+1)B-1)/m,                  R=W-1.           (9)
```

The restrictions of `W` to the three vertical factors are

```text
W|_(x=0)=1,                       W|_(A=0)=W|_(B=0)=0, (10)
```

which proves `(2)`.  The factors are pairwise coprime because `A=B=1` on
`x=0` and `B=1` on `A=0`.  They are reduced: `A` is irreducible, while any
common factor of `B` and

```text
B_y=m x^(2m+2)A^(m-1)
```

would have to divide `xA`, where `B` specializes to one.

## 2. The paid primitive has a reduced global denominator

Expanding `(9)` through `(6)` shows that `R=xN_m`, where

```text
N_m = y
    + ((2m+1)/m)x^(2m)A^(m+1)
    + ((m+1)/m)x^(4m+1)A^(2m+1).                       (11)
```

Therefore the paid primitive is exactly

```text
P_1=R/U=N_m/(AB).                                       (12)
```

This fraction is reduced.  Indeed, if an irreducible polynomial divided
both `N_m` and `AB`, it would divide `R=xN_m` and also `W`, because `AB`
divides `W`.  It would then divide `W-R=1`.  Thus

```text
gcd(N_m,AB)=1.                                          (13)
```

The restrictions needed later are equally exact:

```text
N_m|_(x=0)=y,
UP_1=R=-1 on A=0 and on B=0.                            (14)
```

Hence `P_1` is regular on the paid axis `x=0`, and it has equal coefficient
`-1` relative to `1/U` on every irreducible component above `A=0` and
`B=0`.

## 3. Every nonconstant middle shear has a global pole

Write `D=AB`, so `P_1=N_m/D` is reduced.  There are two exhaustive branches
for a nonconstant `G in k(s)`.

### 3a. A pole at infinity leaves the vertical debt

If `G` has pole order `d>=1` at infinity, then `(12)--(14)` show that
`G(P_1)` has a pole of positive order along every component of `D=0`.
Adding the polynomial `U` cannot cancel a negative valuation.  Thus `U_1`
is not polynomial.

### 3b. A finite pole creates a transverse horizontal divisor

Suppose `G` is regular at infinity.  Since it is nonconstant and `k` is
algebraically closed, it has a finite pole `rho in k`.  In reduced form its
denominator contains `(s-rho)^e`, while its numerator is nonzero at `rho`.
Pullback by `(12)` produces the divisor

```text
H_rho: N_m-rho AB=0.                                   (15)
```

It is a genuine new horizontal divisor.  Formula `(13)` gives

```text
gcd(N_m-rho AB,AB)=1,                                  (16)
```

so it contains no vertical component.  Even more concretely, on the paid
axis

```text
(N_m-rho AB)|_(x=0)=y-rho.                             (17)
```

Thus `H_rho` passes through `(0,rho)`, is smooth and transverse to the axis
there, and has `AB=1`.  At its generic point `P_1=rho`; the numerator of
`G` is a unit and the denominator has positive order.  Therefore
`G(P_1)`, and hence `U_1`, has a genuine pole on `H_rho`.

Every nonconstant rational function either has a pole at infinity or a
finite pole.  Sections 3a and 3b prove assertion 1, including rational
functions with both kinds of poles.

## 4. The constant middle branch restores the paid axis

It remains to take `G=g_0 in k`, so

```text
U_1=U+g_0.                                              (18)
```

Along `A=0` and `B=0`, equation `(14)` says

```text
P_1=-1/U+O(1).                                          (19)
```

For `P_2=P_1+H(U_1)` to be regular there, `H` must have a simple pole at
`g_0` with coefficient `+1`:

```text
H(s)=1/(s-g_0)+O(1).                                   (20)
```

A regular `H` cancels nothing, and a higher pole makes the valuation worse.
But `U_1-g_0=U` also vanishes simply on `x=0`, where `P_1` is regular by
`(14)`.  The compulsory term in `(20)` therefore creates an uncancelled
`1/U` pole on the paid axis.  The simplest choice makes the reversal literal:

```text
P_1+1/(U_1-g_0)=P.                                     (21)
```

Thus the constant branch restores the original axis debt, proving assertion
2 and the three-shear no-go.

## 5. Exact controls and next state

The companion named in the metadata checks `(6)--(14)` for `m=1,...,8`,
the symbolic all-`m` Jacobian and cover identities, reduced-denominator
certificates, the vertical packet, transverse divisors for several finite
poles, polynomial/pole-at-infinity controls, and the exact reversal `(21)`.
Normal and optimized executions must byte-match the frozen transcript before
promotion.

This mechanism is not THM-3776's unequal-residue obstruction.  Here the first
reciprocal shear turns `(1,0,0)` into `(0,-1,-1)`, so the remaining vertical
coefficients are locally equalizable.  Global polynomiality fails because a
nonconstant equalizer has a finite target pole, and its pullback is the new
curve `(15)`.  Longer alternating words must therefore carry a growing ledger
of horizontal fibres `H_rho`; this theorem identifies that state but does not
prove an all-word invariant.  **QED.**
