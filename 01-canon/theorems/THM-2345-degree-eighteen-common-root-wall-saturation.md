---
id: THM-2345
title: "Degree-eighteen common-root wall saturation"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. On the
  degree-eighteen wall 126D=25B^2, every B!=0 point whose structured
  Mordell polynomial has square-class degree at most four lies on the
  deeper wall 20BC+21W=0. Equivalently, in the normalized B=1 chart the
  degree-ten residual has gcd degree at least three exactly on
  20C+21W=0. THM-2297 closes the B=0 edge and THM-2342 closes the deep
  wall, so the complete wall 126D=25B^2 contains no Keller trajectory.
  Other parameter walls, the residual off-wall H_2/H_4 strata, and JC(2)
  remain open.
source: codex-2026-07-25-common-root-wall-saturation
depends_on:
  - THM-2297-degree-eighteen-target-translation-normal-form
  - THM-2332-degree-eighteen-genus-zero-square-class-and-dessin-trap
  - THM-2342-degree-eighteen-deep-wall-first-flux-cover-obstruction
related:
  - THM-2335-degree-eighteen-cyclic-square-class-stratum-empty
  - THM-2338-degree-eighteen-deep-common-root-wall-hurwitz-quartet
  - THM-2341-degree-eighteen-deep-wall-local-genus-split
script: 04-computation/jc2_degree18_common_root_wall_saturation_thm2345.py
output: 05-knowledge/results/jc2_degree18_common_root_wall_saturation_thm2345.out
script_sha256: b41537ed459f844cf963a1b66e720e1065b35e0525df49af365e4b31214da21a
output_sha256: 093b90634ef037d94dc29dc925152463768519733ef9eb88689d0361995708d9
hash_basis: working-tree bytes (LF)
---

# THM-2345 -- the common-root wall saturates to the deep wall

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2332 says that every degree-eighteen Keller survivor has structured
Mordell polynomial

```text
F=4P^3+49Q^2=H S^2,               deg(H) in {0,2,4}.       (1)
```

The previous deep-wall theorems closed the simultaneous locus

```text
126D=25B^2,                       20BC+21W=0.               (2)
```

The first equation in (2) looks much larger: after weighted-projective
normalization it is a two-parameter plane. This theorem proves that its
low-square-class part is nevertheless supported entirely on the second
equation. Thus the already-closed deep wall is not merely one useful
slice; it is the saturation of the whole common-root wall by the
square-class condition.

## 1. The wall exposes a residual decic

Retain THM-2332's exact covariants

```text
P
 =245y^4+1890By^2-24300B^2+122472D,

Q
 =539y^6+11340By^4+183708Cy^3
  +(72900B^2-367416D)y^2
  +(2361960BC+2480058W)y.                           (3)
```

On

```text
126D=25B^2,                                         (4)
```

direct expansion gives

```text
F=1323 y^2 G_10,                                    (5)
```

where

```text
G_10
 =55223y^10
  +1481760By^8+7334712Cy^7+12700800B^2y^6
  +(248618160BC+99018612W)y^5
  +(20412000B^3+1249949232C^2)y^4
  +(1984046400B^2C+2083248720BW)y^3
  +(32141551680BC^2+33748629264CW)y^2
  +206624260800B^2C^2
  +433910947680BCW+227803247532W^2.                (6)
```

The leading coefficient `55223` is constant and nonzero. The factor
`1323y^2` is a square up to a nonzero complex scalar, so `F` and `G_10`
have the same square class.

## 2. Low square class forces a triple gcd

Suppose (1) holds on (4). Absorb the square `y^2` from (5). Then

```text
G_10=H R^2
```

for some `R`, with

```text
10=deg(H)+2deg(R),                 deg(H)<=4,

deg(R)>=3.                                          (7)
```

Differentiating shows

```text
R divides gcd(G_10,dG_10/dy).                       (8)
```

Consequently every low-square-class point satisfies

```text
deg gcd(G_10,G_10')>=3.                             (9)
```

This implication is deliberately only necessary. It remains valid at
higher-multiplicity collisions and does not assume that the squarefree
part has generic root pattern.

## 3. Normalize the weighted chart

First suppose `B!=0`. The weighted action is

```text
(B,C,D,W,y)
 ->
(lambda^2B,lambda^3C,lambda^4D,lambda^5W,lambda y). (10)
```

Over `C`, choose `lambda` so that `lambda^2B=1`. Write the normalized
coordinates again as `(C,W)`. The invariant deep equation transforms as

```text
20C_normalized+21W_normalized
 =lambda^3 B^(-1)(20BC+21W).                        (11)
```

Hence its vanishing is independent of the choice of square root used in
the normalization. Put

```text
A=20C+21W                                             (B=1 chart).     (12)
```

On `A=0`, (6) recovers exactly the residual sextic of THM-2338:

```text
G_10=7y^4 G_6.                                      (13)
```

Thus `deg gcd(G_10,G_10')>=3` everywhere on the deep wall. The converse
is the substantive statement.

## 4. Exact subresultant saturation

Compute the subresultant sequence of the normalized degree-ten
polynomial and its derivative. Its degrees are

```text
10,9,8,7,6,5,4,3,2,1,0.                            (14)
```

Let

```text
Sres_2=a_2y^2+a_1y+a_0,

Sres_1=b_1y+b_0.                                    (15)
```

Because both leading coefficients are nonzero constants, specialization
commutes with this subresultant test. Condition (9) forces all five
coefficients in (15) to vanish.

Their exact deep factors are

```text
b_1=A^3 r_1,              b_0=A^4 r_0,

a_2=A s_2,                a_1=A^2 s_1,
                           a_0=A^2 s_0,              (16)
```

up to nonzero integer contents. Let

```text
I=(r_1,r_0,s_2,s_1,s_0) in Q[C,W].                 (17)
```

An exact grevlex Groebner basis, with variable order `(W,C)`, is

```text
[
 (15309C^2+250)
 (30618C^2+361)^3
 (199644669C^4+7654500C^2+62500),

 20C+21W
].                                                 (18)
```

In particular,

```text
A=20C+21W belongs to I.                             (19)
```

Suppose `A!=0` and (9) holds. Equations (15)--(16) then force every
generator of `I` to vanish. But (19) forces `A=0`, a contradiction.
Together with (13), this proves the exact equivalence

```text
deg gcd(G_10,G_10')>=3
 iff 20C+21W=0                         on B=1.       (20)
```

The first factor in (18) records the exceptional collision ratios seen
after approaching the deep wall. It is not an additional off-deep
component: the second basis element pins every common zero of the
stripped coefficient ideal back onto `A=0`.

## 5. Close the complete wall

Let a Keller trajectory lie on (4).

If `B=0`, then (4) also gives `D=0`. This is the complete plane
`B=D=0`, including its axes and origin, which THM-2297 already proves
empty.

If `B!=0`, THM-2332 gives (1), hence (7)--(9). After normalization,
(20) gives `20C+21W=0`, and (11) gives the invariant equation

```text
20BC+21W=0.                                        (21)
```

The trajectory therefore lies on the deep wall (2). THM-2342, using
THM-2338 and THM-2341, proves that complete wall empty. This is the
required contradiction.

Therefore

```text
no degree-eighteen Keller trajectory satisfies 126D=25B^2.             (22)
```

This closes a codimension-one weighted wall in the four-parameter
degree-eighteen cone. It does not close the off-wall `H_2` or `H_4`
strata, any other parameter wall, `JC(2)`, or `DC(2)`.

## 6. Exact reproduction

Run

```bash
python3 04-computation/jc2_degree18_common_root_wall_saturation_thm2345.py
python3 -O 04-computation/jc2_degree18_common_root_wall_saturation_thm2345.py
```

Both transcripts are byte-identical to the stored output. The companion
checks (5)--(6), the deep factorization (13), the complete degree profile
(14), every exact division in (16), the two-element basis (18), and the
zero remainder in (19). It also checks a generic deep point with gcd
degree three and two squarefree off-deep hostile controls. No executable
check uses Python `assert`.

The independent hostile audit rebuilt (5)--(6) from exact rational
samples, checked the square-class/gcd implication including higher
collisions, and recomputed the saturation with the reversed variable
order `(C,W)`. That independent basis again contains `20C+21W`; after
substituting `C=-21W/20`, its other basis element agrees with (18) up to
the nonzero scalar `512*10^12`. A separate hand-coded finite-field gcd
engine found exactly the deep line and no off-deep jump for every prime
from `29` through `97`. The audit also checked specialization, weighted
normalization, the `B=0` dependency, both execution modes, both hashes,
and the documentation gate. QED.
