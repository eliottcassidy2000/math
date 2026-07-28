---
id: THM-2830
title: "Disjoint positive adjacent-cone factorial moment-three detection"
status: >
  PROOF-REDUCTION CANDIDATE + FINITE-EXACT; NOT PROVED.  The desired
  two-cone orientation is reduced coefficientwise to one symmetric
  four-index matching inequality.  Every one of its three complementary
  matchings implies the required quartic coefficient after averaging.
  Exact bounded and randomized audits are positive, but the universal
  matching inequality remains open.  Two tempting stronger routes
  (pairwise factorial dominance and coefficientwise Laguerre-ratio
  monotonicity) have exact counterexamples.
source: root/disjoint-adjacent-cone-factorial-orientation-2026-07-28
depends_on: []
related:
  - THM-2824-arbitrary-three-slot-factorial-moment-divisibility-and-atomic-orientation-boundary
  - THM-2828-lower-prefix-cone-factorial-moment-three-detection
script: 04-computation/gmc_disjoint_cone_matching_reduction_thm2830.py
output: 05-knowledge/results/gmc_disjoint_cone_matching_reduction_thm2830.out
script_sha256: a5c803e1117d9f0740aaa7e0ef91f00ff823f2a34777ffab57dcc53dc2b592d9
output_sha256: 74280f6e106718c4bed8eea08c566fd932e93919d368a9d8f176d7cc92780334
hash_basis: LF-normalized bytes
---

# THM-2830 -- two disjoint factorial prefix cones

**PROOF-REDUCTION CANDIDATE + FINITE-EXACT; NOT PROVED.**

Let

```text
L(s^n)=n!,                 f_n=s^n/n!,
d_i=f_(i+1)-f_i,

U=sum_(i<b) lambda_i(f_(i+1)-f_i),
V=sum_(j>=b) mu_j(f_(j+1)-f_j),
```

where both finite coefficient families are nonnegative and nonzero.  Prove
the mixed orientation

```text
2L(V^3)L(UV)-3L(UV^2)L(V^2)>=0,
```

with equality only on the adjacent singleton boundary.  Together with the
strict positive cubic tensor of THM-2828, this would make factorial moments
one through three detect every complex plane spanned by two such disjoint
positive cones, and give a many-versus-many two-charge moment-six theorem
when the constant slot is absent.

This file now gives the exact coefficient formula and reduces its positivity
to one four-index matching lemma.  That lemma has extensive exact evidence
but no universal proof.  Consequently this file still supplies no proved
dependency and no Gaussian-moment theorem.

## 1. Positive quadratic and cubic tensors

Put

```text
H(a,b)=L(d_a d_b)=binom(a+b,a),                       (1)
T(a,b,c)=L(d_a d_b d_c).                              (2)
```

THM-2828 proves `T(a,b,c)>0`.  Direct multiplication also gives the useful
two-atom identity

```text
d_p d_q
 =binom(p+q,p) f_(p+q)
  +binom(p+q+2,p+1)d_(p+q+1).                         (3)
```

Since `f_N=1+sum_(j<N)d_j`, the nonconstant adjacent-difference
coefficients of `(3)` are all nonnegative.  This explains the persistent
positive signal, but it does not by itself prove the quartic orientation.

For a fixed lower atom `d_i`, write

```text
D_i(V)=2L(V^3)L(d_iV)-3L(d_iV^2)L(V^2).               (4)
```

The desired orientation is

```text
D(U,V)=sum_(i<b)lambda_i D_i(V).                       (5)
```

Thus it is enough to prove `(4)` for `i<b`.

## 2. Exact quartic polarization

For four upper indices `a,b,c,d>i`, define

```text
P_i(a,b,c,d)
 = H(i,a)T(b,c,d)+H(i,b)T(a,c,d)
  +H(i,c)T(a,b,d)+H(i,d)T(a,b,c),                    (6)

N_i(a,b;c,d)
 = T(i,a,b)H(c,d)+T(i,c,d)H(a,b).                    (7)
```

There are three complementary matchings:

```text
(ab|cd),               (ac|bd),               (ad|bc). (8)
```

The completely polarized coefficient of `(4)`, divided by its positive
multinomial multiplicity, is

```text
C_i(a,b,c,d)
 =1/2 [
    P_i(a,b,c,d)
   -sum_(A subset {a,b,c,d}, |A|=2)
      T(i,A)H(A^c)
 ].                                                     (9)
```

Every unordered pair occurs once among the two sides of the three
matchings.  Therefore the single matching inequality

```text
P_i(a,b,c,d) >= 3 N_i(a,b;c,d)                         (10)
```

for each matching in `(8)` implies, after summing the three copies and
dividing by three,

```text
P_i(a,b,c,d)
 >=sum_(|A|=2)T(i,A)H(A^c),                            (11)
```

and hence `C_i>=0`.  Equality in `(10)` is expected only at

```text
a=b=c=d=i+1.                                           (12)
```

At `(12)`, all three matching inequalities and `(9)` are equalities.  If
`(10)--(12)` hold universally, every coefficient of `(4)` is nonnegative,
and the only zero of `(4)` is `V` proportional to `d_(i+1)`.  Equations
`(5)` then give exactly the adjacent-singleton equality claimed in the
target.

## 3. Factorial normalization of the live lemma

Put

```text
tau(a,b,c)
 =T(a,b,c)a!b!c!/(a+b+c)!
 =(a+b+c+1)
    [1/(a+1)+1/(b+1)+1/(c+1)]-1.                     (13)
```

Multiplying `(10)` by the positive common factor

```text
i!a!b!c!d!
```

turns it into

```text
sum_(x in {a,b,c,d})
  (i+x)!(a+b+c+d-x)! tau({a,b,c,d}\{x})

 >=3[
   (i+a+b)!(c+d)! tau(i,a,b)
  +(i+c+d)!(a+b)! tau(i,c,d)
 ].                                                     (14)
```

This is the current sharp target.  Exact computation says that the
difference in `(14)` is coordinatewise nondecreasing in each upper index
from the base point `(12)`.  A universal proof of that forward-difference
claim, or a direct factorial-word injection proving `(14)`, would finish
THM-2830.

## 4. Two false stronger routes

The first-moment comparison behind THM-2824 does extend below the upper
support.  But the tempting upgrade to all falling-factorial moments is
false.  In the exact pair block with

```text
i=0,              (y,z)=(3,44),              a=4,
```

the tilted fourth falling-factorial mean is

```text
Phi=1467522360/901,
[(y)_4+(z)_4]/2=1629012,

Phi-[(y)_4+(z)_4]/2=-217452/901<0.                    (15)
```

This is the first counterexample by total `y+z` in the audited universe.
For fixed `(i,y,a)` and `z` tending to infinity, the limiting sign is
controlled by

```text
2(i+1)/(i+a+1)+(i+1)/(y+1)-1,                         (16)
```

so `(15)` belongs to an infinite hostile family.  Pairwise increasing
convex order cannot prove `(14)`; the four-index compensation is
load-bearing.

There is a second exact boundary in the sign-Laguerre basis
`phi_k=(-1)^kL_k`.  Although

```text
d_j=sum_(k=1)^(j+1)binom(j,k-1)phi_k                  (17)
```

has positive coefficients, the coefficientwise quotient for `V^2/V`
need not increase.  For

```text
V=d_1+t d_3
```

the `k=3` adjacent determinant is

```text
2t(5640t^2+371t-3),                                   (18)
```

which is negative for sufficiently small positive `t`.  Pascal averaging
can repair this local failure, but raw Laguerre coefficient monotonicity
cannot be used as a theorem.

## 5. Exact evidence and scope

The companion checks:

1. `(1)--(3)` directly in the normalized factorial basis;
2. the polarization identity `(9)` against direct polynomial expansion;
3. all three matching inequalities in a bounded exact universe, including
   the equality boundary `(12)`;
4. coordinate forward differences of the normalized gap `(14)` in the
   same universe;
5. exact cone-level values of `(4)--(5)`; and
6. both hostile mechanisms `(15)` and `(18)`.

The older scout found `9,488` exhaustive cone pairs, `50,000` random
integer cone pairs, and `4,116` symmetrized quartic cells nonnegative.  The
new matching audit is stronger structurally, but remains finite evidence.

What is **proved in this file** is the algebraic reduction
`(1)--(14)` and the two stopping boundaries.  What is **open** is the
universal sign in `(10)` or `(14)`.  Therefore arbitrary disjoint positive
cones, the corresponding many-versus-many moment-six statement, general
HYP-8765, and SFC(3) remain open here.
