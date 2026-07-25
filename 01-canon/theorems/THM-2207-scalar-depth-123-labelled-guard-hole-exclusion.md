---
id: THM-2207
title: "Scalar depth-(1,2,3) labelled guard-hole exclusion"
status: >
  PROVED + VERIFIED-EXACT. In the scalar five-unit/three-blocker branch, the
  actual blocker valuation profile (1,2,3) is empty. On the primitive
  13^4 guard-safe annulus, the depth-one and depth-two masks are pullbacks
  from 1,014 and 78 sign classes. For every one of their 79,092 typed
  pairs, an exact labelled root-fibre branch certificate determines the five
  largest conditional capacities among all 13,182 unit sign classes. The
  unique minimum full-annulus deficit is 1,608. Direct full-torsion
  enumeration independently reconstructs every unit capacity on the hostile
  row. Together with THM-2198, THM-2204, and THM-2205, this empties every
  scalar profile whose unique deepest blocker has valuation at most three.
  Thus a surviving scalar branch has deepest valuation at least four. In the
  actual dyadic-terminal lane THM-2203 also gives the upper bound nineteen,
  leaving 1,136 finite valuation profiles; none of those deeper profiles is
  settled here, and this is not a proof of LRC(14).
source: codex-klein-2026-07-24-scalar-depth123-labelled-capacity
depends_on:
  - THM-2192-scalar-five-plus-three-root-sheet-chord-invoice
  - THM-2198-scalar-five-plus-three-image-pump-and-first-depth-exclusion
related:
  - THM-2201-cyclic-root-fibre-hasse-jet-transition-carrier
  - THM-2203-fixed-dyadic-coordinate-section-and-covector-intersection
  - THM-2204-scalar-depth-223-thirteen-lift-capacity-law
  - THM-2205-scalar-depth-113-exact-lift-capacity-exclusion
script: 04-computation/lrc14_scalar_depth123_branch_certificate_thm2207.py
output: 05-knowledge/results/lrc14_scalar_depth123_branch_certificate_thm2207.out
script_sha256: 8e5bd4bb0f6e9072fd9f1992a0d2e454d63e92f476395ca021f773fcf811f2df
output_sha256: 6710488fdb3b9b1f0a31f8b8ed771f079e75b4e7dadfb79d8863a7d929df3890
hash_basis: working-tree bytes (LF)
---

# THM-2207 -- scalar depth-`(1,2,3)` exclusion

Put

```text
D_a={t in R/Z:||at||<1/14},
C_H={t in R/Z:||Ht||>1/7}.                           (1)
```

In the scalar `5+3` branch of THM-2192 and THM-2198, suppose

```text
C_H subset union_(i=1)^5 D_(q_i)
             union union_(j=1)^3 D_(c_j)             (2)
```

almost everywhere. The coefficients `H,q_1,...,q_5` are positive
thirteen-units, and the three actual blockers are positive multiples of
thirteen. After relabelling, THM-2192 gives

```text
1<=lambda_1<=lambda_2<lambda_3,
lambda_j=nu_13(c_j).                                 (3)
```

This theorem excludes

```text
(lambda_1,lambda_2,lambda_3)=(1,2,3).                (4)
```

The finite object is a mixed-depth labelled root-fibre matrix. The proof
keeps the correlations between each unit label, each guard hole, and the two
different blocker scales. It does not replace a thirteen-lift family by its
average.

## 1. The mixed primitive layer

Assume (4) and put

```text
N=13^4=28561,                    Q=13^3=2197.         (5)
```

Multiplication of primitive numerators by `H modulo N` normalizes the guard
to one and replaces every terminal coefficient by its product with
`H^(-1) modulo N`. This is a bijection and preserves all three
thirteen-valuations. Define

```text
U_N={z mod N:13 does not divide z and 7||z||_N>N},
||z||_N=min(z mod N,-z mod N).                       (6)
```

The exact cardinality is

```text
|U_N|=18830.                                         (7)
```

The depth-three blocker is safe everywhere on `U_N`. If its normalized
coefficient is `13^3w`, then at a primitive numerator

```text
13^3w*z/13^4=wz/13
```

is a nonzero thirteenth root, whose norm is at least
`1/13>1/14`.

Every primitive phase `r modulo Q` has thirteen primitive roots

```text
z=r+kQ modulo N,                    k in F_13.        (8)
```

There are `phi(Q)=2028` such phases. Let

```text
h(r)=#{k:z in U_N}.                                  (9)
```

Exact root counting gives

```text
h(r) in {9,10},
#{r:h(r)=9}=1450,               #{r:h(r)=10}=578,
sum_r h(r)=18830.                                    (10)
```

A depth-one coefficient has the form `13u`. Its activity is constant on
the roots (8) and is the phase bit

```text
A_u(r)=1_[14||ur||_Q<Q].                             (11)
```

The unit part `u modulo Q`, up to sign, ranges through

```text
phi(Q)/2=1014                                        (12)
```

labels. A depth-two coefficient has the form `13^2v`; its phase bit is

```text
B_v(r)=1_[14||13vr||_Q<Q],                           (13)
```

where `v modulo 13^2`, up to sign, ranges through

```text
phi(13^2)/2=78                                       (14)
```

labels. The two blockers have different types, so their complete labelled
pair universe has

```text
1014*78=79092                                        (15)
```

rows. This includes every possible coincidence of their effective masks.

For a unit sign class `q modulo N`, define its guard-surviving root count
above a phase by

```text
w_q(r)=#{k:z=r+kQ belongs to U_N and
              14||qz||_N<N}.                        (16)
```

The root-window law of THM-2198 gives

```text
w_q(r) in {0,1,2}.                                   (17)
```

There are

```text
phi(N)/2=13182                                       (18)
```

unit sign classes.

For a mixed blocker pair `(u,v)`, put

```text
P_(u,v)={r:A_u(r)=B_v(r)=0}.                         (19)
```

The full residual size and the conditional capacity of a unit class are
exactly

```text
R(u,v)=sum_(r in P_(u,v)) h(r),
C_q(u,v)=sum_(r in P_(u,v)) w_q(r).                  (20)
```

If (2) held, the deepest blocker would contribute nothing on `U_N`, while
the five unit masks would have to cover this residual. Duplicate unit sign
classes are redundant in a union. Therefore the union has size at most the
sum of the five largest capacities among the distinct classes in (18).
Consequently a strict inequality

```text
R(u,v)>sum_(i=1)^5 C_(i)(u,v)                        (21)
```

for every pair contradicts (2).

## 2. The exact branch certificate

Equation (17) permits a short exact certificate for all the order
statistics in (21). For a fixed unit class put

```text
F_q=sum_r w_q(r),
X_q(A)=sum_(r in A) w_q(r),
X_q(B)=sum_(r in B) w_q(r).                          (22)
```

With `A={r:A_u(r)=1}` and `B={r:B_v(r)=1}`,
inclusion--exclusion gives

```text
C_q(u,v)
 =F_q-X_q(A)-X_q(B)+X_q(A intersection B).           (23)
```

By (17),

```text
X_q(A intersection B)<=2|A intersection B|.         (24)
```

The exact capacity is also no larger than either one-blocker capacity.
Hence

```text
C_q(u,v)<=U_q(u,v),

U_q(u,v)=min(
 F_q-X_q(A),
 F_q-X_q(B),
 F_q-X_q(A)-X_q(B)+2|A intersection B|
).                                                   (25)
```

For each pair the certificate visits unit classes in decreasing order of
the integer upper bound (25), computes their exact capacities from two
bitsets representing `w_q>=1` and `w_q>=2`, and maintains the five largest
exact pairs `(capacity,-label)`. It stops only when the next upper bound is
**strictly** smaller than the fifth retained capacity. Every unvisited exact
capacity is then strictly smaller as well. In particular, an equality at
the cutoff is evaluated, so the canonical smaller-label tie break cannot be
lost.

This certifies all 79,092 rows after

```text
940857 exact candidate evaluations,
average 11.895729024427 per row,
maximum 27 in one row.                               (26)
```

No pruning decision uses floating point. All stored counts are exact
integers within their audited dtypes.

## 3. Exact deficit and hostile row

The exhaustive result is

```text
R(u,v)-sum_(i=1)^5 C_(i)(u,v)>=1608                  (27)
```

on the full primitive annulus, for every one of the 79,092 typed pairs.
The minimum is unique in the least-positive sign-representative convention.
It occurs at

```text
(u,v)=(799,46),
R(u,v)=13526.                                        (28)
```

The exact five largest `(capacity,unit label)` pairs there are

```text
((2604,5193),(2472,10386),(2292,7773),
 (2288,10388),(2262,7775)).                          (29)
```

Their capacities sum to `11918`, and

```text
13526-11918=1608.                                    (30)
```

Two frozen digests protect the full scan:

```text
pair-table:
79b9b75f3732e47b43c2bba726906250ce0c1069b90d8952c33caa8a8364570f

branch trace:
2699c79c62d4c0e5805529b26dc4ce7ae5ed1f0146be25714b898ea42176802a.
                                                               (31)
```

The companion performs the following hostile controls.

1. For all 1,014 depth-one classes and all 78 depth-two classes, direct
   evaluation on every root in (8) agrees with the reduced phase bits
   (11) and (13).
2. On the unique hostile row, an independent direct enumeration of all
   18,830 full torsion points reconstructs all 13,182 unit capacities and
   reproduces (28)--(30).
3. All six depth-three sign classes modulo thirteen are checked directly
   to be inactive on `U_N`.
4. The unit classes split into 1,014 coefficient-lift families of size
   thirteen. On the hostile residual, every family satisfies THM-2204's
   exact sum law. The five hostile unit labels lie over base classes
   `(799,599,1015,597,1013)`; their complete labelled family rows are
   retained in the output.

The source handoff and independent audit map are recorded in
`agents/broadcast/MSG-2828-from-klein-2026-07-24-thm-2207-frozen-artifacts-and.md`.
Fresh ordinary and `python3 -O` runs are byte-identical to the canonical
output. Neither `assert` nor optimization-sensitive control flow is used.

Because `N` and `Q` are powers of thirteen, they are coprime to seven and
fourteen. Thus neither a guard equality nor a danger equality can occur on
these torsion universes. The strict predicates in (6), (11), (13), and
(16) therefore agree with the almost-everywhere cover without an endpoint
convention.

Equations (21) and (27) contradict (2), proving that the actual blocker
profile `(1,2,3)` is empty. ∎

## 4. Consequence and structural boundary

The unique-deepest ordering (3) leaves only these profiles through depth
three:

```text
depth 2: (1,1,2);
depth 3: (1,1,3), (1,2,3), (2,2,3).                 (32)
```

THM-2198 excludes the first, THM-2205 excludes `(1,1,3)`, this theorem
excludes `(1,2,3)`, and THM-2204 excludes `(2,2,3)`. Therefore

```text
every surviving scalar 5+3 branch has lambda_3>=4.   (33)
```

In the actual dyadic-terminal LRC lane, THM-2203 independently gives
`lambda_3<=19` and counts 1,140 profiles before these four exclusions.
Thus (32) leaves 1,136 profiles. This is a finite ledger, not a feasible
flat enumeration: the depth-twenty torsion layer and its labelled
owner-current sidecars remain enormous.

The mechanism here also identifies what a deeper proof must preserve. The
scalar family sum from THM-2204 is insufficient: the hostile top five live
in five different lift families, and their order is controlled by the
phasewise alignment of guard holes with labelled unit endpoints. THM-2201's
Hasse-jet carrier or an equivalent sparse correlation state is the natural
recursive object. A future depth transition must compress that state
without destroying the conditional top-five order statistic.
