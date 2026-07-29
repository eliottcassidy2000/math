---
id: THM-2949
title: "Fixed rank-thirty-five cofactor Newton atlas"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
  HOSTILE AUDIT.  For every translated four-slot support
  {n,n+a,n+b,n+M} with 0<a<b<M<=12, one fixed 35-by-35
  degree-seven Macaulay cofactor is nonzero at every integer n>=0.
  THM-2947's conjugate-pair corank gap therefore proves first-window
  SFC(4) throughout this 220-support atlas.  Width eleven is the first
  new boundary beyond the previous two-chart closure.  The exact
  companion verifies three outside-grid values per family, the first
  positive shifted Gregory--Newton base, and the complete ordinary
  sign-variation structure.  At width three all 216
  distinguished cofactors have degree 130 and share a degree-60
  selected-chart seam divisor; this is not a canonical cofactor or an
  arbitrary-width theorem.
source: codex-gmc-uniform-width-extension-2026-07-29
audit: Pending independent hostile audit.
depends_on:
  - THM-2947-conjugate-pair-corank-parity-and-one-minor-resultant-gate
  - THM-2943-width-seven-eight-two-chart-macaulay-resultant-closure
related:
  - THM-2940-consecutive-four-slot-macaulay-newton-short-closure
  - THM-2942-macaulay-extraneous-flag-factor-and-pluecker-mutation
  - THM-2944-width-nine-ten-two-chart-macaulay-resultant-closure
  - THM-2946-full-macaulay-maximal-minor-gcd-and-chart-free-resultant
script: 04-computation/gmc_fixed_rank_thirty_five_cofactor_newton_atlas_thm2949.py
output: 05-knowledge/results/gmc_fixed_rank_thirty_five_cofactor_newton_atlas_thm2949.out
script_sha256: 9a1c7068e079e232dc97fd6eb925621aa74b3d636380a85995f8e0db8b30aa54
output_sha256: 6bd2df044f1c07962fd1ff7c3b27a7dc20681f06d807bd4a45d40314027979ca
constructor_dependency_sha256: d2f8afeba7dd6c7950405a4845d7bf112b6c9872dd8161146446be8bbdaae0ba
parity_gate_dependency_sha256: d1bd09ff20925183f5488fcd8850469867f1dfad2bdb808504fc896708605744
hash_basis: LF-normalized bytes
---

# THM-2949 -- fixed rank-thirty-five cofactor Newton atlas

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
HOSTILE AUDIT.**

## 1. Finite-width statement

Let

```text
L:C[s] -> C,                         L(s^j)=j!.          (1)
```

For every

```text
n>=0,                    0<a<b<M<=12,                  (2)
```

and every nonzero polynomial

```text
H=c_0 s^n+c_1 s^(n+a)+c_2 s^(n+b)+c_3 s^(n+M),       (3)
```

at least one of

```text
L(H),                 L(H^2),                 L(H^3), L(H^4) (4)
```

is nonzero.  Thus first-window SFC(4) holds on all

```text
sum_(M=3)^12 C(M-1,2)=C(12,3)=220                    (5)
```

translated support types of width at most twelve.

Widths three through ten were already closed by the earlier fixed- or
two-chart theorems.  The `45` width-eleven supports are the first new
boundary here; width twelve adds another `55`.

## 2. The fixed cofactor

Use the mean-elimination coordinates and denominator-cleared ternary
forms

```text
Q,                         C,                         F              (6)
deg 2,                     deg 3,                     deg 4
```

from THM-2943.  In the degree-seven Macaulay matrix use the inherited
square chart

```text
selected global rows
  0,...,19; 21,...,29,35; 36,...,41.                 (7)
```

The last six are the selected quartic rows.  Delete selected-row
position `30`, namely global row `36`, which is the multiplier row
`x2^3 F`, and target column `0`, whose monomial is `x2^7`.  Thus this
is the corner diagonal adjugate coordinate in the inherited chart.
Call the resulting `35`-minor

```text
P_(M,a,b)(n).                                         (8)
```

It uses `20` quadratic rows, `10` cubic rows, and `5` quartic rows.
THM-2925's coefficient-degree law therefore gives

```text
deg P_(M,a,b)
 <=20(M-1)+10(2M-1)+5(3M-1)
 =55M-35.                                             (9)
```

The exact atlas attains equality in `(9)` for every one of the `220`
supports.

## 3. Shifted Gregory--Newton certificates

Normalize every polynomial in `(8)` to have positive leading
coefficient.  The companion evaluates it at every interpolation depth

```text
0,1,...,55M-35                                        (10)
```

and at the three independent outside-grid depths

```text
55M-34,                 55M-33,                 2(55M-35)+3. (11)
```

It then finds the smallest integer `B>=0` for which

```text
Delta^j P_(M,a,b)(B)>0          for 0<=j<=55M-35.     (12)
```

The finite prefix

```text
P_(M,a,b)(0),...,P_(M,a,b)(B-1)                       (13)
```

is checked directly and has no zero.  Equations `(12)--(13)` give

```text
P_(M,a,b)(n)!=0                    for every n>=0.     (14)
```

The stored transcript gives the complete exceptional-support/base
atlas, ordinary sign-variation census, and immutable coefficient and
Newton-vector digests for each width.

For auditability, here is the complete exceptional-base atlas in the
notation `a:b@B` for support `(0,a,b,M)` and first positive Newton base
`B`.  The `var` field is the ordinary coefficient-variation histogram.

```text
M=3;  var=0:1;                 maxB=0;  shifted=none
M=4;  var=0:3;                 maxB=0;  shifted=none
M=5;  var=0:6;                 maxB=0;  shifted=none
M=6;  var=0:9,1:1;             maxB=1;
      shifted=1:5@1
M=7;  var=0:12,1:3;            maxB=3;
      shifted=1:5@2,1:6@3,2:6@1
M=8;  var=0:14,1:6,2:1;        maxB=6;
      shifted=1:2@6,1:4@1,1:5@2,1:6@4,1:7@6,2:6@1,2:7@3
M=9;  var=0:15,1:10,2:3;       maxB=10;
      shifted=1:2@10,1:3@6,1:4@2,1:6@5,1:7@7,1:8@9,
              2:3@6,2:4@1,2:6@2,2:7@4,2:8@6,3:7@1,3:8@2
M=10; var=0:15,1:16,2:5;       maxB=14;
      shifted=1:2@14,1:3@10,1:4@5,1:5@2,1:6@5,1:7@9,
              1:8@12,1:9@14,2:3@10,2:4@5,2:5@2,2:7@5,
              2:8@8,2:9@10,3:4@5,3:5@1,3:7@1,3:8@4,
              3:9@6,4:5@1,4:9@2
M=11; var=0:18,1:20,2:7;       maxB=20;
      shifted=1:2@20,1:3@15,1:4@9,1:5@5,1:6@1,1:7@11,
              1:8@15,1:9@18,1:10@20,2:3@14,2:4@9,2:5@4,
              2:6@1,2:7@6,2:8@10,2:9@13,2:10@15,3:4@9,
              3:5@4,3:6@1,3:8@5,3:9@8,3:10@10,4:5@4,
              4:9@3,4:10@5,5:10@1
M=12; var=0:22,1:22,2:10,3:1; maxB=27;
      shifted=1:2@27,1:3@21,1:4@15,1:5@8,1:6@4,1:7@12,
              1:8@19,1:9@23,1:10@25,1:11@27,2:3@20,2:4@15,
              2:5@9,2:6@4,2:8@13,2:9@17,2:10@20,2:11@22,
              3:4@14,3:5@9,3:6@4,3:8@6,3:9@10,3:10@14,
              3:11@16,4:5@9,4:6@4,4:9@5,4:10@8,4:11@10,
              5:6@4,5:10@2,5:11@4
```

## 4. The ordinary sign hierarchy

After leading-sign normalization, the ordinary coefficient variations
give a useful but nonuniform secondary classification.  Two cases have
structural proofs:

```text
zero variations:   every coefficient is positive;
one variation:     negative ... negative | positive ... positive.
                                                               (15)
```

In the one-cut case,

```text
P(0)<0<P(B).                                           (16)
```

THM-2942's one-cut operator independently explains why all derivatives
and all positive-order forward differences are positive from `B`
onward.  Cases with two or more variations are instead discharged by
the direct shifted Newton certificate `(12)`.  The stored transcript
gives the complete variation histogram at every width.

The distinction is load-bearing.  For the first new-width support

```text
(0,1,2,11),                                            (17)
```

the normalized cofactor has two ordinary sign variations and

```text
P(11)>0,       P(12)<0,       P(19)<0,       P(20)>0. (18)
```

It therefore has exactly two positive real roots, one in `(11,12)` and
one in `(19,20)`, while it has no nonnegative integral root.

At width twelve, `(0,1,4,12)` already has three ordinary variations,
yet its exact prefix is nonzero and its positive Newton tail starts at
base `15`.  Thus the theorem is an integer-depth nonvanishing theorem,
not positivity of a resultant on the continuous real ray and not a
uniform one-cut theorem.

## 5. The complete width-three cofactor bank

At support `(0,1,2,3)`, the companion reconstructs all

```text
6*36=216                                                (19)
```

distinguished cofactors obtained by deleting one selected quartic row
and one target column.  Exact adjugate evaluation at `0,...,130`,
followed by the same three outside-grid controls, proves that all `216`
are nonzero polynomials of exact degree `130`.

Their primitive polynomial gcd is the degree-`60` polynomial

```text
G_60(n)
 =(n+1)^2 (2n+1)^5 (49n^2+99n+38)^5
  (n+2)^19 (n+3)^19
  (8788n^5+54873n^4+126718n^3
             +132729n^2+61288n+9732).                 (20)
```

Every factor in `(20)` is strictly positive for `n>=0`.  For the fixed
cofactor `(8)`, primitive normalization gives

```text
P_(3,1,2)/G_60=(n+1)^4 R_66(n),                       (21)
```

where all coefficients of the irreducible degree-`66` factor `R_66`
are strictly positive; its exact digest is stored by the companion.

Thus the single cofactor's positivity persists after removing a large
common selected-chart seam.  But `(20)` also shows why one cofactor is
not canonical: the whole `216`-bank carries shared chart data.  The
basis-relative fifth-compound energy from THM-2947 is the invariant
nonvanishing object.  No identification of `G_60` with a canonical
``half norm'' is asserted.

## 6. Projective consequence

For every physical depth, the quadratic in `(6)` is positive definite,
so `V_R(Q,C)` is empty.  The three-slot face restriction used in
THM-2943 gives `Q` not dividing `C`; hence `(Q,C)` is a regular
sequence.  THM-2947 applies.

Equation `(14)` supplies one nonzero `35`-minor.  The conjugate-pair
corank gap then forces full degree-seven Macaulay rank and

```text
Res(Q,C,F)!=0.                                         (22)
```

Thus `Q,C,F` have no common nonzero projective zero.  Mean elimination
turns this into `(4)`.

## 7. Scope

This is a finite exact theorem through `M=12`.  It does not prove:

```text
* the same cofactor works at width thirteen or at arbitrary width;
* positivity on every real depth;
* a canonical ordering of the three quotient conjugate pairs;
* that the selected-chart seam G_60 is an invariant norm; or
* a lower bound for the full fifth-compound energy.                  (23)
```

The exact width-thirteen continuation control

```text
(0,1,4,13)                                             (24)
```

still has a nonzero finite prefix and a positive shifted Newton tail
from base `21`; its ordinary coefficients have three sign variations.
This shows that the higher-cut regime persists, but one control is not
a width-thirteen closure.

The fixed cofactor is a cheap sufficient certificate because THM-2947
forbids corank one on the physical real locus.  Outside that locus a
single `35`-minor need not be a resultant gate.

In particular, no bare polynomial divisibility by the resultant is
being used: already at width three the fixed cofactor has degree `130`,
the resultant has degree `112`, and exact division has nonzero
remainder.  The implication `(14) => (22)` is the real-locus
conjugate-pair corank argument, not a factorization statement.

A separately reported Koszul exchange would replace the
`20Q+10C+5F` bank used here by a `21Q+9C+5F` bank and lower the formal
degree invoice from `55M-35` to `54M-35`.  That exchange is not used:
this theorem contains no unverified row replacement or basis-sign
identity, and its executable certificate remains exactly `(7)--(8)`.

## 8. Exact companion

Run

```text
python 04-computation/gmc_fixed_rank_thirty_five_cofactor_newton_atlas_thm2949.py
python -O 04-computation/gmc_fixed_rank_thirty_five_cofactor_newton_atlas_thm2949.py
```

The companion uses integer/rational FLINT arithmetic, exact
interpolation, exact adjugates, and explicit `require` gates.  Both
executions reproduce the stored transcript byte-for-byte.
