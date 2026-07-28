---
id: THM-2828
title: "Lower-prefix cone factorial moment-three detection"
status: >
  RESERVED / PROOF-COMPLETE THEOREM CANDIDATE + VERIFIED-EXACT /
  AWAITING INDEPENDENT HOSTILE AUDIT.  THM-2824 extends from a single lower
  interval to every nonzero nonnegative combination of lower adjacent
  factorial differences.  Equivalently, the lower coefficient vector may
  have arbitrarily many occupied slots provided all proper prefix sums are
  nonpositive.  An exact strictly positive cubic tensor formula supplies
  the missing orientation.  Every resulting binary plane is detected by
  factorial moments one through three.  When the constant slot is absent,
  the associated many-slot two-charge Gaussian envelope is detected by
  moment at most six.
source: root/audit-2809-2026-07-28
depends_on:
  - THM-2824-arbitrary-three-slot-factorial-moment-divisibility-and-atomic-orientation-boundary
related:
  - THM-2810-factorial-hankel-faithfulness-and-bounded-radial-carrier-no-go
  - THM-2815-optimal-finite-laguerre-carrier-and-radial-selector-access-boundary
  - HYP-8765-gmc2-radial-channel-return-tower
script: 04-computation/gmc_lower_prefix_cone_thm2828.py
output: 05-knowledge/results/gmc_lower_prefix_cone_thm2828.out
script_sha256: 2dcfc1112349335d43d96187af4ffcd3860c7c5761c6d93df8e887e6e4b440ad
output_sha256: 43acbc6a0d3e81cb431f4b2516b443a226ccf76630a7dc417ec98407541d996d
hash_basis: LF-normalized bytes
---

# THM-2828 -- lower-prefix cone factorial moment-three detection

**RESERVED / PROOF-COMPLETE THEOREM CANDIDATE + VERIFIED-EXACT /
AWAITING INDEPENDENT HOSTILE AUDIT.**

THM-2824 proves an atomic mixed orientation for one lower interval and one
upper interval.  Its lower occurrence is linear.  The only apparent
obstruction to replacing that lower interval by an arbitrary positive
combination is the strict cubic orientation used in the binary
quadratic-divisibility argument.

That obstruction has a closed positive answer: the cubic tensor of adjacent
factorial differences is entrywise strictly positive.  The result below
therefore allows arbitrarily many lower slots and arbitrary coefficients
inside a full prefix-majorization cone.

## 1. The lower prefix cone

Let

```text
L(s^n)=n!,                      f_n=s^n/n!,

d_i=f_(i+1)-f_i.                                      (1)
```

Fix integers `1<=b<c` and real numbers

```text
lambda_i>=0,                    0<=i<b,                (2)
```

not all zero.  Put

```text
U=sum_(i=0)^(b-1) lambda_i d_i,
V=f_c-f_b.                                             (3)
```

Both directions have mean zero:

```text
L(U)=L(V)=0.                                           (4)
```

If

```text
U=sum_(n=0)^b u_n f_n,                                 (5)
```

then coefficient comparison gives

```text
sum_(n=0)^i u_n=-lambda_i,          0<=i<b,
sum_(n=0)^b u_n=0.                                    (6)
```

Thus `(2)` is equivalent to the intrinsic coefficient condition

```text
sum_(n=0)^i u_n<=0 for every proper prefix i<b,
sum_(n=0)^b u_n=0,                                    (7)
```

with at least one strict proper-prefix inequality.  Raw coefficients
`u_n` may change sign many times; it is their cumulative profile that is
one-sided.

## 2. The adjacent-difference cubic tensor is strictly positive

For nonnegative integers `i,j,k`, write

```text
S=i+j+k,             E2=ij+ik+jk,          E3=ijk,

M_ijk=(i+j+k)!/(i!j!k!).                                (8)
```

Expanding the eight terms of `d_i d_j d_k` and dividing by `M_ijk`
gives

```text
L(d_i d_j d_k)
 =M_ijk [2(S+1)^2+S E2-E3]/
          ((i+1)(j+1)(k+1)).                           (9)
```

The bracket is strictly positive.  Indeed,

```text
S E2-E3
 =sum_sym i^2 j+2ijk>=0.                               (10)
```

Consequently

```text
L(d_i d_j d_k)>0                    for all i,j,k>=0.  (11)
```

For comparison, the corresponding quadratic tensor is the Pascal fixed
point

```text
L(d_i d_j)=binom(i+j,i)>0.                             (12)
```

Equations `(9)--(11)` imply

```text
t111=L(U^3)
 =sum_(i,j,k<b) lambda_i lambda_j lambda_k
    L(d_i d_j d_k)>0.                                  (13)
```

Strictness uses only that some `lambda_i` is positive.

## 3. THM-2824 is linear in the lower direction

Put

```text
g11=L(U^2),          g12=L(UV),          g22=L(V^2),

t111=L(U^3),         t122=L(UV^2),       t222=L(V^3).
                                                               (14)
```

Define the same division-free orientation as in THM-2824:

```text
D(U,V)=2t222 g12-3t122 g22.                            (15)
```

The lower direction occurs linearly in `(15)`.  Therefore

```text
D(U,V)=sum_(i=0)^(b-1) lambda_i D_i(b,c),              (16)
```

where `D_i(b,c)` is exactly the THM-2824 atom.  Its universal theorem and
`(2)` give

```text
D(U,V)>=0.                                             (17)
```

No claim about coefficientwise positivity of a second arbitrary cone
direction is used here.

## 4. Binary divisibility closes the plane

For

```text
H=alpha U+beta V,                    alpha,beta in C,   (18)
```

the quadratic form

```text
Q(alpha,beta)=L(H^2)                                  (19)
```

is positive definite over `R`: `U,V` are nonzero real polynomials and are
linearly independent because only `V` contains `f_c`.  Let

```text
C(alpha,beta)=L(H^3).                                  (20)
```

As in THM-2824, a common complex projective zero of `Q,C` exists exactly
when the real quadratic `Q` divides the real cubic `C`.  One of the two
division-free remainder invariants is

```text
I2
 =3t122 g11g22-2t222 g12g11-t111 g22^2
 =-g11 D(U,V)-t111 g22^2.                              (21)
```

Equations `(13),(17)` and Gram positivity give

```text
I2<0.                                                   (22)
```

Hence `Q` and `C` have no common projective zero.  Therefore

```text
L(H)=L(H^2)=L(H^3)=0                  implies H=0       (23)
```

for every plane `(18)` with lower direction in the prefix cone `(7)`.

## 5. Many-slot two-charge Gaussian consequence

If `lambda_0=0`, then both `U` and `V` are divisible by `s`.  Put

```text
h=H/s,
P=W+Z h(ZW),                         W=conj(Z),         (24)
```

for a standard complex Gaussian `Z`.  Charge balance gives

```text
E[P^(2j)]=binom(2j,j)L(H^j),          E[P^(2j-1)]=0,
                                                     j=1,2,3. (25)
```

Equation `(23)` now proves detection of every such many-slot two-charge
envelope by moment at most six.  Unlike the three-slot application of
THM-2824, `H` and `P` may have arbitrarily many monomials.  This still
occupies one prefix cone; it is not arbitrary radial-coefficient GMC2.

## 6. Exact companion and scope

The exact companion independently expands `(9)`, audits `(12)`, checks the
prefix-sum equivalence `(6)`, and tests the full invariant chain
`(13)--(22)` on a deterministic bank of nonnegative rational cone
directions.  It also checks every underlying THM-2824 atom in the bounded
bank used by the cone tests.  It uses only integers and fractions, explicit
exception gates, and no truth-bearing Python assertions.

This theorem is a genuine many-slot enlargement of THM-2824 and a partial
answer to the arbitrary-radial-coefficient problem.  It does not prove
entrywise positivity of the quartic polarization when both directions
vary, general Strong Factorial, general HYP-8765, or GMC2.

**Proof complete; awaiting independent hostile audit before promotion.**
