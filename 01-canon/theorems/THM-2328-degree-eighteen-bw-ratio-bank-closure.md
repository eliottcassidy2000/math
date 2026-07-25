---
id: THM-2328
title: "Degree-eighteen B--W ratio-bank closure"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. In the genuine
  nonsplit polynomial exact-square-prefix degree-eighteen branch of
  THM-2262/2297, all eight B--W ratio points in THM-2311's exactly
  two-sparse bank are empty. The two-point quadratic orbit has an
  absolutely irreducible genus-four trigonal spectrum with ten simple
  branches and one smooth totally ramified cubic fibre. The six-point
  sextic orbit has an absolutely irreducible genus-three spectrum with
  ten simple branches and one ordinary nonvertical node. Infinity is
  unramified throughout. Every rational Keller trajectory is therefore
  constant and gives the inherited nonsplit-deck contradiction. Together
  with THM-2314, THM-2316, THM-2320, and THM-2324, this closes all 31
  exactly-two-sparse ratios. This does not prove JC(2).
source: codex-2026-07-25-degree18-bw-bank
depends_on:
  - THM-2262-degree-eighteen-trigonal-spectral-discriminant-reduction
  - THM-2297-degree-eighteen-target-translation-normal-form
  - THM-2311-degree-eighteen-two-sparse-weighted-ratio-bank
related:
  - THM-2314-degree-eighteen-bd-linear-ratio-closure
  - THM-2316-degree-eighteen-cd-linear-ratio-closure
  - THM-2320-degree-eighteen-dw-ratio-bank-closure
  - THM-2324-degree-eighteen-bc-rational-ratio-closure
script: 04-computation/jc2_degree18_bw_ratio_bank_closure_thm2328.py
output: 05-knowledge/results/jc2_degree18_bw_ratio_bank_closure_thm2328.out
script_sha256: adfacdad45e3eed070e9f5d335411bcf88f71e00ec0fc773980949176efcfbdc
output_sha256: 598b1e6dc061a7a9e51851805db6fc30714e35cf44ea37fd09aaed1f4670608c
hash_basis: working-tree bytes (LF)
---

# THM-2328 -- the full B--W ratio bank has positive genus

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2311 reduces the exactly-two-sparse degree-eighteen branch to a finite
weighted-projective bank.  On the `B`--`W` line,

```text
t=W^2/B^5,
```

the complete off-axis bank is cut out by

```text
P_2(t)
 =5313800000
  +4508659468656t
  -136738899331083t^2,

P_6(t)
 =5511577600000000000000000000
  +4983290602536960000000000000000t
  -6564822237254419568640000000000t^2
  -3094052863483309848285092659200000t^3
  -81862566455344350924421142159812608t^4
  -744088924275617882256518828471658624t^5
  -2973811237322720333634598763466407943t^6.        (1)
```

Their degrees sum to `2+6=8`.  This theorem treats both complete Galois
orbits.

## 1. Weighted representatives and ratio fields

Use THM-2297's target-translation normal form

```text
(B,C,D,W),                         weights (2,3,4,5),

G_0(u,y;B,C,D,W)=0,               wt(u,y)=(2,1).    (2)
```

For a root `alpha` of either polynomial in (1), put

```text
(B,C,D,W)=(alpha,0,0,alpha^3).                     (3)
```

Then

```text
W^2/B^5=alpha^6/alpha^5=alpha.                     (4)
```

The constant terms in (1) are nonzero, so `alpha!=0` and (3) is
off-axis.

Modulo nineteen, the monic reduction of `P_2` is

```text
t^2+4t+5,                                          (5)
```

which has no root in `F_19`.  Thus `P_2` is irreducible over `Q`.

Modulo thirteen, the monic reduction of `P_6` is

```text
t^6+t^5-t^4-t^3-2t^2-t+4.                         (6)
```

Let this reduction be `p`.  Rabin's exact test gives

```text
gcd(p,t^(13^2)-t)=1,
gcd(p,t^(13^3)-t)=1,
p divides t^(13^6)-t.                              (7)
```

The prime divisors of six are two and three, so (7) proves that `p`, and
hence `P_6`, is irreducible.  Equations (5)--(7) therefore account for
exactly two and six distinct characteristic-zero ratios.

## 2. Uniform absolute irreducibility

There is a short coefficient proof that works for every off-axis
`B`--`W` specialization, not only the eight points in (1).

The spectral polynomial is cubic in `u` with constant leading coefficient
`-26040609`.  If it were reducible over an algebraic closure of the
constant field, it would have a root in the rational function field in
`y`.  After division by the leading coefficient, that root is integral
over the polynomial ring.  Since the latter is integrally closed, the
root is a polynomial in `y`.

If the root had degree `d>2`, the `u^3` term would have `y`-degree `3d`,
strictly larger than the possible competing degrees

```text
2d+2,                 d+4,                 6.       (8)
```

Hence every possible root has the form

```text
u=ay^2+by+c.                                        (9)
```

Define the common infinity cubic

```text
L(a)=1127-138915a+1607445a^2-26040609a^3.          (10)
```

After (9), the `y^6` and `y^5` coefficients are exactly

```text
L(a),                       b L'(a).                (11)
```

Moreover,

```text
Disc L=-153384762202971019112448 !=0.               (12)
```

Thus (11) forces `b=0`.  Once `b=0`, every remaining term is even in `y`
except the deepest odd term, and the `y^1` coefficient is exactly

```text
-5878656 W.                                        (13)
```

Equation (13) cannot vanish when `W!=0`.  This proves absolute
irreducibility uniformly.  As an independent machine control, adjoining
either `P_2(t)` or `P_6(t)` to the seven coefficient equations obtained
from (9) gives the reduced global Gröbner basis `{1}`.

Every specialized normalization below is consequently a connected
degree-three cover of the `y`-line.

## 3. The common branch pattern

Let `K=Q(alpha)` for either orbit and write

```text
G_alpha(u,y)=G_0(u,y;alpha,0,0,alpha^3),
Delta_alpha(y)=Disc_u G_alpha.                     (14)
```

Exact Euclidean arithmetic in `K[y]` gives one element `Y in K` and a
monic polynomial `h_10 in K[y]` such that

```text
Delta_alpha(y)
 =-153384762202971019112448
   *(y-Y)^2 h_10(y),                               (15)

deg h_10=10,
gcd(h_10,h'_10)=1,
gcd(h_10,y-Y)=1.                                   (16)
```

Thus there are ten distinct simple discriminant values away from `Y`.
Each contributes one ordinary ramification point of index two.  The local
meaning of the doubled value `Y` differs between the two ratio fields and
must be checked rather than inferred from (15).

For that check, write the exceptional cubic as

```text
a_3u^3+a_2u^2+a_1u+a_0
```

and use its standard triple-root covariants

```text
I=a_2^2-3a_3a_1,
J=2a_2^3-9a_3a_2a_1+27a_3^2a_0.                  (17)
```

All identities below are exact native arithmetic in the relevant power
basis, not floating-point root tests.

## 4. The quadratic orbit has genus four

For `P_2(alpha)=0`, evaluation at `y=Y` gives

```text
I=J=0,

G_alpha(u,Y)=-26040609(u-R)^3,

R=-a_2/(3a_3).                                     (18)
```

The point `(R,Y)` is smooth.  After clearing denominators and content, the
power-basis vector of `partial_y G_alpha(R,Y)`, in the highest-power-first
basis `(alpha,1)`, is

```text
(3060070329391879293, 3475058479525000),           (19)
```

which is nonzero.  Hence the exceptional fibre has one smooth point of
ramification index three and contributes two to the ramification divisor.
Together with the ten simple values in (16),

```text
R_total=10+2=12,

2g-2=3*(-2)+12=6,

g=4.                                                (20)
```

Every identity is an identity in the quadratic field, so both embeddings
have genus four.

## 5. The sextic orbit has genus three

For `P_6(alpha)=0`, the covariants in (17) are nonzero.  The repeated root
of the exceptional cubic is

```text
R=(9a_3a_0-a_2a_1)/(2I),                           (21)
```

and its third root is `S=-a_2/a_3-2R`.  Exact reduction gives

```text
G_alpha(u,Y)=-26040609(u-R)^2(u-S),

G_alpha(R,Y)
 =partial_u G_alpha(R,Y)
 =partial_y G_alpha(R,Y)=0.                         (22)
```

Let the tangent cone at `(R,Y)` be

```text
A U^2+B UV+C V^2.                                  (23)
```

After clearing denominators and content, the highest-power-first vector
of `A` in `(alpha^5,...,alpha,1)` is

```text
(31028486880250604115193956128190271227417532086873591,
 442483383802993521499669948929970393471906416638688,
 -295071948144665319720183819761328922699169146850304,
 -24345794298751232255212750527130807587448320000000,
 2017545934341343303155833912562276925440000000000,
 2784511768598360910718235115520000000000000000). (24)
```

The corresponding vector of `B^2-4AC` is

```text
(1950819895284140522250365395002160195454731331033763,
 236863098816083905325135771691922457285511608320992,
 20296498065720949838482428076714044255844497024000,
 533587777622434082055979241574918456384000000000,
 -47367739906029936163675236030935040000000000000,
 -53016344589311194767319040000000000000000000).  (25)
```

Both are nonzero.  Therefore `(R,Y)` is an ordinary node, and `A!=0`
makes both tangent branches nonvertical.  The two normalization branches
over the node, and the simple point `S`, are unramified over `y`.
Only the ten simple roots in (16) ramify, so

```text
R_total=10,

2g-2=3*(-2)+10=4,

g=3.                                                (26)
```

The calculation takes place in the full sextic field and therefore proves
the same local type and genus for all six embeddings at once.

## 6. Infinity is unramified

Use the weighted infinity chart

```text
z=1/y,                         v=u/y^2.              (27)
```

The `B` and `W` terms have lower weighted order, so at `z=0` every curve
has the same cubic

```text
L(v)=1127-138915v+1607445v^2-26040609v^3.          (28)
```

By (12), it has three distinct roots.  They give three distinct smooth
points above infinity, with `z` a local parameter at each.  Infinity is
unramified.  Equations (15)--(28) therefore include every ramification
contribution used in the genus counts.

## 7. Every B--W Keller trajectory is impossible

A putative survivor in the retained branch supplies

```text
(u,y) in C(x)^2,                  G_0(u,y)=0.         (29)
```

If the pair were nonconstant, it would extend to a nonconstant morphism
from `P^1` to the appropriate projective normalization.  Riemann--Hurwitz
forbids such a map to any of the genus-three or genus-four curves above.
Thus `u` and `y` are constant.

Undo the constant weighted curve isomorphism (3).  The wall `y=0` is
already empty by THM-2262.  The inherited first-flux identity

```text
Z=T^2=-2N_2/(5103y)                                (30)
```

then makes `Z`, nonzero `T`, and `q` constant.  The genuine nonsplit deck
fixes the algebraically closed constant field but sends `q` to `-q`,
contradicting `q!=0`.  All eight ratios in (1) are empty.

## 8. Exact gain and stopping boundary

The complete `B`--`W` bank of THM-2311 is empty:

```text
8 -> 0.                                             (31)
```

Combining the five audited plane closures gives

```text
BD:6,       CD:4,       DW:4,       BC:9,       BW:8,

6+4+4+9+8=31.                                      (32)
```

Thus no exactly-two-sparse point remains in the genuine nonsplit
polynomial exact-square-prefix degree-eighteen branch.

The connection ledger is

```text
source:
  the two irreducible factors in THM-2311's B--W bank;

map:
  a ratio root alpha goes to (B,W)=(alpha,alpha^3);

preserved:
  labelled weighted orbit, complete Galois orbit, spectral isomorphism
  class, normalization, local branch indices, and genus;

temporarily forgotten:
  the scaled third flux, Keller one-form, and whole-Faber sidecar;

why restoration is unnecessary:
  positive genus already makes every rational spectral trajectory
  constant;

hostile control:
  B=W=1 has squarefree degree-twelve branch discriminant, detecting an
  accidentally generic or identically repeated calculation.           (33)
```

This theorem does not classify any three- or four-sparse singular
trigonal stratum.  Split/even descent, other Newton edges, THM-2240's
`DC(2)` continuation axes, `JC(2)`, and `DC(2)` remain open.

## 9. Exact reproduction

Run

```bash
python3 04-computation/jc2_degree18_bw_ratio_bank_closure_thm2328.py
python3 -O 04-computation/jc2_degree18_bw_ratio_bank_closure_thm2328.py
```

Both executions are byte-identical to the stored output.  The companion
verifies the two finite-field irreducibility tests, the representative
identity (4), the coefficient proof (8)--(13), two independent global
Gröbner certificates, both complete number-field discriminant
factorizations (15)--(16), the smooth triple fibre (18)--(19), the node
factorization and tangent certificates (21)--(25), separable infinity,
both ramification counts, the bank arithmetic, and the squarefree hostile
control.  The Riemann--Hurwitz and deck steps are the mathematical proof
above, not computer assumptions. QED.

The independent hostile audit reconstructed the two discriminants by the
hand cubic formula over direct quotient rings `Q[a]/(P_2)` and
`Q[a]/(P_6)`, without using the companion's `CRootOf`, discriminant, or
orbit helpers. It reproduced the linear gcds, squarefree residuals,
triple and node factorizations, all three local power-basis certificates,
the genera, and the squarefree hostile control. A separate finite-field
Rabin implementation checked both ratio fields. Normal, optimized, and
stored transcripts and both hashes agree; the source-bank, infinity,
flux/deck, and `31=6+4+4+9+8` ledgers also passed.
