---
id: THM-3619
title: "Russell-cylinder even-fold higher-jet staircase"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE / PENDING INDEPENDENT AUDIT.
  On the doubly tuned THM-3612 even-fold locus, exact arbitrary-target-two-form
  collision matrices give four further necessary side-jet conditions through
  source order ten.  With q_m=Q^(m)(1), they are q_3=54, q_4=-270,
  q_5=1620, and q_6=-11340.  The displayed finite staircase is proved by
  exact jet algebra; its apparent all-order recurrence and rational limiting
  germ are CONJECTURAL beyond q_6.  No JC(2) counterexample is claimed.
source: root / audit_thm3612 higher-jet continuation, 2026-08-21
audit: PENDING -- provisional theorem and exact companion require hostile audit
depends_on:
  - THM-3612-russell-cylinder-even-fold-nongraph-collision-jet-rigidity
related:
  - THM-3605-russell-cylinder-graph-slice-puncture-no-filling
  - THM-3610-russell-cylinder-full-linear-projection-collision-rigidity
script: 04-computation/jc2_russell_cylinder_even_fold_higher_jet_staircase_thm3619.py
output: 05-knowledge/results/jc2_russell_cylinder_even_fold_higher_jet_staircase_thm3619.out
script_sha256: 7b5cc4f246c86e5cb4a061a23d371053ecb2fe9d54e0dc2d5e0b35569107cae7
output_sha256: d16a7df35ff58ca77708169da6d7aea8832c8f511cf0ae05135d6dca2f20873e
hash_basis: raw LF bytes
---

# THM-3619 -- Russell-cylinder even-fold higher-jet staircase

**RESERVED / PROVISIONAL PROOF CANDIDATE / PENDING INDEPENDENT AUDIT.**
This file continues THM-3612 only on its doubly tuned boundary.  It proves a
finite staircase through source order ten.  The tempting all-order pattern in
Section 6 is deliberately labelled **CONJECTURAL**.

All rings, formal germs, and differential forms are over `C`.

## 0. Setup, normalization, and statement

Retain the THM-3612 even fold

```text
E_Q(x,t)=(x,Q(x)+t^2,t),
Q even,                 Q(0)=-3/4,       Q(1)=Q(-1)=-3.      (1)
```

Write

```text
q_m=Q^(m)(1).                                               (2)
```

The three source points `p_-=(-1,0),p_0=(0,0),p_+=(1,0)` have
one common image `y_0` under the stabilized Russell map.  At the first-jet
escape of THM-3612, any hypothetical nonzero constant source Jacobian can be
normalized to the common value `12`.  THM-3612 proves the first two necessary
conditions

```text
q_1=9/2,                         q_2=-27/2.                 (3)
```

This theorem proves four more finite implications.  For every pair of regular
functions on the Russell target cylinder,

```text
(3) and q_3 != 54
    ==> source Jacobian is not a nonzero constant;

(3), q_3=54, and q_4 != -270
    ==> source Jacobian is not a nonzero constant;

(3), q_3=54, q_4=-270, and q_5 != 1620
    ==> source Jacobian is not a nonzero constant;

(3), q_3=54, q_4=-270, q_5=1620, and q_6 != -11340
    ==> source Jacobian is not a nonzero constant.          (4)
```

Thus every surviving even polynomial fold must satisfy

```text
(q_1,...,q_6)=(9/2,-27/2,54,-270,1620,-11340).             (5)
```

The locus `(5)` remains **OPEN**.  No assertion beyond `q_6` is part of the
proved statement.

## 1. Why arbitrary target two-forms suffice

Use the smooth target parameters of THM-3612

```text
c=C,                         y=Y,                 z=S+3/4.  (6)
```

Every target pair `F,G` gives the target two-form

```text
Omega=dF wedge dG
     =A(c,y,z) dc wedge dy
      +K(c,y,z) dc wedge dz
      +R(c,y,z) dy wedge dz.                              (7)
```

The proof below allows `A,K,R` to be completely arbitrary formal coefficient
functions.  This is a strict enlargement of the closed decomposable forms
coming from regular pairs.  Consequently, an obstruction for `(7)` is an
all-regular-target-pair obstruction and needs no linear-projection or
polynomial-in-ambient-coordinates hypothesis.

Let `phi_i` be the target germ of the fold at `p_i`, using local source
coordinates `xi=x-x_i,t`.  Define

```text
phi_i^*Omega=j_i(xi,t) dxi wedge dt.                       (8)
```

For a constant normalized Jacobian, every nonconstant coefficient of every
`j_i` vanishes and all three constant coefficients equal `12`.

## 2. The finite pullback universe and its unique vertical cokernel

For `N>=0`, put

```text
T_N=(C[c,y,z]_(degree <= N))^3,
S_N=direct_sum_(i=-,0,+) C[xi,t]_(degree <= N).            (9)
```

Substitution in `(7),(8)` gives an exact linear map

```text
P_N^Q:T_N -> S_N.                                         (10)
```

Its dimensions are

```text
dim T_N=3 binom(N+3,3),              dim S_N=3 binom(N+2,2). (11)
```

On the tuned tangent planes `(3)`, the linear target germs are

```text
             c              y                 z
p_-          3xi            -9xi+4t          -9xi/4
p_0          3xi            -9xi+4t           0
p_+          3xi            -9xi+4t           9xi/4.       (12)
```

For a homogeneous target coefficient of degree `N`, every source monomial
except the pure vertical `t^N` coefficient can be prescribed independently
on the three branches.  The three vertical values are affine in the three
equally spaced `z/c` slopes `-3/4,0,3/4`.  Hence the homogeneous symbol has

```text
rank 3(N+1)-1                                             (13)
```

and its one-dimensional cokernel is represented by

```text
Delta_N(j)=[t^N]j_- - 2[t^N]j_0 + [t^N]j_+.              (14)
```

This explains why each new even obstruction is one scalar after all lower
source rows have been killed.  Odd orders through nine introduce no new
constant-target obstruction in the exact matrices.

## 3. Exact higher-jet invoices

Direct expansion of the unrestricted map `(10)`, followed by exact row
reduction modulo all lower source rows, gives the following four identities.
The common constant has already been normalized to `12`.

```text
q_1=9/2, q_2=-27/2:
Delta_4 = -(32/9)(q_3-54);                               (15)

q_3=54:
Delta_6 = -(64/81)(q_4+270);                             (16)

q_4=-270:
Delta_8 = -(32/243)(q_5-1620);                           (17)

q_5=1620:
Delta_10= -(64/3645)(q_6+11340).                         (18)
```

Here “modulo all lower source rows” has a literal meaning: if the coefficients
of total source degree less than `N` agree with the constant `12`, then the
left side of the corresponding invoice is the displayed scalar.  Target
coefficient jets of degree `N` cancel from `Delta_N` by `(12)--(14)`; target
jets of higher degree cannot affect source order `N`.

The exact lower-row and homogeneous-symbol ranks are

| source invoice | lower map | lower rank | symbol rank | full rank |
|---|---:|---:|---:|---:|
| `Delta_4` | `30 x 60` | `26` | `14` | `40` |
| `Delta_6` | `63 x 168` | `57` | `20` | `77` |
| `Delta_8` | `108 x 360` | `100` | `26` | `126` |
| `Delta_10` | `165 x 660` | `155` | `32` | `187` |

On the hostile control immediately before each new tuning, adjoining the
constant target column raises the corresponding full rank by one.  On the
newly tuned value, it does not.  Thus `(15)--(18)` are necessary conditions,
not statistics of an intermediate projection.

## 4. `Q_dag` is closed at source order four

THM-3612 stopped at

```text
Q_dag=-3/4-27x^2/2+18x^4-27x^6/4.                       (19)
```

It has

```text
(q_1,q_2,q_3)=(9/2,-27/2,-378).                         (20)
```

Therefore `(15)` gives

```text
Delta_4=-(32/9)(-378-54)=1536.                          (21)
```

This contradicts a constant Jacobian.  In particular, the second-order
survivor displayed in THM-3612 cannot be extended through source order four,
even after allowing every higher target jet in both outputs.

One sparse representative of the exact row identity behind `(21)` is

```text
Delta_4 = 128[1]j_0
          +(2/3)([xi]j_- - [xi]j_+)
          -(4/9)([xi^2]j_- + [xi^2]j_+)
          +(2/3)([xi t^2]j_- - [xi t^2]j_+),            (22)
```

where `[1]j_0` denotes the constant coefficient common to the normalized
lower system.  Under `j_i=12`, `(22)` is exactly `1536`.

## 5. Explicit polynomial boundary controls

The following even polynomials retain `(1)` and successively meet the side
jet conditions.  They are hostile controls for the sharpness of the finite
invoices; they are not asserted to carry a Darboux pair.

```text
Q_3=-3/4-27x^2/2+9x^4+81x^6/4-27x^8+9x^10,             (23)

Q_4=-3/4-27x^2/2-45x^4/4+405x^6/4-297x^8/2
    +90x^10-81x^12/4,                                   (24)

Q_5=-3/4-27x^2/2-45x^4+270x^6-486x^8
    +855x^10/2-189x^12+135x^14/4,                       (25)

Q_6=-3/4-27x^2/2-189x^4/2+567x^6-2457x^8/2
    +2835x^10/2-1863x^12/2+1323x^14/4-99x^16/2.         (26)
```

Their tuned jets and next hostile invoices are

| control | matched side jets | next side jet | next invoice |
|---|---|---:|---:|
| `Q_3` | `q_1,q_2,q_3` | `q_4=7506` | `Delta_6=-6144` |
| `Q_4` | through `q_4` | `q_5=-127980` | `Delta_8=51200/3` |
| `Q_5` | through `q_5` | `q_6=2269620` | `Delta_10=-360448/9` |
| `Q_6` | through `q_6` | not tested here | **OPEN** |

They can be generated recursively from `Q_dag` by

```text
Q_3=Q_dag+9 x^4(x^2-1)^3,
Q_4=Q_3-(81/4)x^4(x^2-1)^4,
Q_5=Q_4+(135/4)x^4(x^2-1)^5,
Q_6=Q_5-(99/2)x^4(x^2-1)^6.                             (27)
```

The factor `x^4(x^2-1)^m` preserves the central value and every side jet of
order below `m`; its `m`-th derivative at `1` is `2^m m!`.  The independent
controls `x^2(x^2-1)^m` change central jets while producing the same new side
jet.  Exact paired probes with these two perturbations confirm that the four
quotient invoices depend on the displayed side jet, not on the central jet or
the next unused side jet.

## 6. Conjectural all-order pattern -- not proved

The four proved coefficients equal

```text
c_n=2^(n+3)/(3^(n-1)(n-1)!),                n=3,4,5,6, (28)
```

and `(15)--(18)` have the common form

```text
Delta_(2n-2)=-c_n (q_n+(n+1)q_(n-1)).                  (29)
```

It is natural to conjecture `(29)` for every `n>=3`.  If so, a formal survivor
would have to obey

```text
q_n=-(n+1)q_(n-1),
q_n=(-1)^(n-1) 9(n+1)!/4.                               (30)
```

Together with `Q(1)=-3`, the unique formal germ would be

```text
Q_infinity(x)=-3/4-9/(4x^2).                            (31)
```

The pole at `x=0` is incompatible with a polynomial even fold satisfying
`Q(0)=-3/4`.  This would close every polynomial even fold if the all-order
identity were proved.  It is **not** proved here: `(28)--(31)` are
`CONJECTURAL` beyond `q_6` and must not be used as an all-even-fold exclusion.

## 7. Scope and exact reproduction

The proved result closes `Q_dag` and every even polynomial fold failing one of
the finite side-jet conditions `(5)`.  It does not close:

- the finite tuned locus `(5)`, represented by `Q_6`;
- the conjectural recurrence beyond `q_6`;
- non-even folds;
- nonquadratic folds in the stable coordinate; or
- implicit source planes outside the THM-3612 family.

All these exits remain **OPEN**.

The deterministic companion constructs the exact truncated branch germs,
the unrestricted arbitrary-two-form pullback matrices, the sparse quotient
rows, every rank in the table, the controls `(19),(23)--(27)`, and the four
invoices `(15)--(18)`.  It also verifies the tangent-symbol cokernel and the
finite coefficient pattern `(28)`.

Reproduce with

```bash
python3 04-computation/jc2_russell_cylinder_even_fold_higher_jet_staircase_thm3619.py
python3 -O 04-computation/jc2_russell_cylinder_even_fold_higher_jet_staircase_thm3619.py
```

Both streams must be byte-identical to
`05-knowledge/results/jc2_russell_cylinder_even_fold_higher_jet_staircase_thm3619.out`.
