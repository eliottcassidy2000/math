---
id: THM-2223
title: "Four-checkpoint Selberg relation-and-carry gate"
status: >
  PROVED + VERIFIED-EXACT. For four consecutive 169-checkpoints of a
  three-comb danger union, absence of a one-frequency-per-checkpoint
  relation of coefficient height at most sixteen bounds the common
  survivor mass by 26,873,856/200,533,921, strictly below the scalar
  residual floor 961/6930. Hence every putative lambda_1>=6 survivor from
  THM-2222 has a cross-checkpoint relation
  sum_r n_r 169^r d_(i_r)=0 with |n_r|<=16. Regrouping by blocker gives a
  two- or three-support coefficient vector of height at most 77,688,640
  whose balanced base-169 digit supports are disjoint. There are exactly
  44,264,640 raw nonzero checkpoint templates modulo simultaneous sign.
  The 13-adic first-digit invoice further forces two distinct blocker
  depths to differ by at most seven, unconditionally excluding the 21
  profiles (a,a,c) with 6<=a<c<=19 and c-a>=8. This reduces the remaining
  open cubic triple box to an explicit union of relation-and-carry planes;
  it does not prove their extremal bound or LRC(14).
source: codex-2026-07-24-four-checkpoint-selberg-carry
depends_on:
  - THM-2085-explicit-height-57-rank-seven-selberg-gate
  - THM-2222-scalar-transfer-parity-tower-and-four-checkpoint-survivor-reduction
related:
  - THM-2144-anisotropic-selberg-kraft-relation-box
  - THM-2163-radix-relation-carry-descent
  - THM-2210-nested-binomial-minorant-and-adaptive-moment-lp-hierarchy
script: 04-computation/lrc14_four_checkpoint_selberg_carry_gate_thm2223.py
output: 05-knowledge/results/lrc14_four_checkpoint_selberg_carry_gate_thm2223.out
script_sha256: cc4f9e4b1af41061061915a250851344993ce9f0ed714211948b17e575614a3a
output_sha256: e86de5b6493950b44651334f4912a07f6c78e8b344020f0116d9377f15ec0977
hash_basis: working-tree bytes (LF)
---

# THM-2223 -- four checkpoints force a bounded carry

This theorem applies the relation-and-carry spectrum to the finite extremal
problem isolated by THM-2222. It is a reduction of that problem, not its
solution.

## 1. Checkpoint event and relation packet

On `T=R/Z`, put

```text
D_a={x:||ax||<1/14},                  A=169.          (1)
```

For three distinct positive integers `d_1,d_2,d_3`, define

```text
U_r(d)=union_(j=1)^3 D_(A^r d_j),     0<=r<=3,
K_4(d)=intersection_(r=0)^3 U_r(d).                   (2)
```

For an integer `H>=1`, a one-frequency-per-checkpoint relation of height
`H` consists of choices

```text
i_r in {1,2,3},
n_r in [-H,H] intersection Z,         0<=r<=3,       (3)
```

with the `n_r` not all zero and

```text
sum_(r=0)^3 n_r A^r d_(i_r)=0.                        (4)
```

The indices attached to zero coefficients are immaterial. Notice that (4)
is a relation between checkpoint layers. No hypothesis about relations
among the three frequencies inside one fixed layer is needed.

## 2. Positive Selberg upper product

THM-2085, Section 2 records the periodic Vaaler interval sandwich. For the
basic interval `D_1` and degree `H`, let `V_H` be its upper polynomial.
Away from the two endpoints,

```text
1_(D_1)<=V_H,                 V_H>=0,
FourierSupport(V_H) subset [-H,H],
integral V_H=1/7+1/(H+1).                             (5)
```

Therefore, pointwise away from a finite null set,

```text
1_(K_4(d))(x)
 <=product_(r=0)^3
      sum_(j=1)^3 V_H(A^r d_j x).                    (6)
```

The nonnegativity in (5) is load-bearing: it permits multiplication of the
four upper bounds. This is a positive upper product for an intersection of
unions, whereas THM-2085 uses a signed lower tensor for a box.

Expand the right side of (6). There are `3^4=81` terms, one for each
`i=(i_0,i_1,i_2,i_3)`. A Fourier monomial in the corresponding term has
one frequency per checkpoint:

```text
sum_(r=0)^3 n_r A^r d_(i_r),          |n_r|<=H.      (7)
```

If (4) has no solution, the only zero frequency in every one of the
eighty-one products is `n_0=...=n_3=0`. One-dimensional Fourier
orthogonality then gives the exact constant term

```text
integral product_r V_H(A^r d_(i_r)x)
  =(1/7+1/(H+1))^4.                                  (8)
```

Integrating (6) proves the relation-free estimate

```text
measure K_4(d)
 <=81(1/7+1/(H+1))^4.                                (9)
```

Repeated frequencies cause no exception: any nonconstant zero-frequency
choice created by a repetition is exactly a relation (4).

## 3. Degree sixteen crosses the scalar floor

Set `H=16`. Then

```text
1/7+1/17=24/119,

81(24/119)^4
 =26873856/200533921,                                (10)

961/6930-26873856/200533921
 =925325143/198528581790
 >0.                                                 (11)
```

Consequently

```text
measure K_4(d)>=961/6930

 implies

there exist i_r,n_r satisfying (3)--(4) with H=16.   (12)
```

Degree sixteen is the first degree for this exact eighty-one-term
certificate. At degree fifteen its right side is

```text
81(1/7+1/16)^4
 =22667121/157351936
 >961/6930.                                         (13)
```

This is certificate minimality only; it is not an assertion that sixteen is
the smallest possible true relation height.

## 4. Regrouping exposes the carry sidecar

Given (4), collect the digits assigned to each blocker:

```text
a_j=sum_(r:i_r=j) n_r A^r,             1<=j<=3.     (14)
```

Then

```text
a_1d_1+a_2d_2+a_3d_3=0,                              (15)

|a_j|
 <=16(1+169+169^2+169^3)
 =77688640.                                         (16)
```

The vector `a` is nonzero. Indeed, if

```text
sum_(r=0)^R b_r169^r=0,       |b_r|<=16,            (17)
```

and `b_R` is the highest nonzero digit, then

```text
169^R
 <=|b_R|169^R
 <=16 sum_(r<R)169^r
 <169^R,                                             (18)
```

a contradiction. Thus balanced base-169 expansions with digits in
`[-16,16]` have no carry collision. Since the `d_j` are positive, (15) has
at least two nonzero coefficients and both signs occur.

The extra structure in (14) is stronger than a generic coefficient box:
at every digit position `r`, at most one of `a_1,a_2,a_3` has a nonzero
digit. The checkpoint layer owns that digit. This disjoint digit support is
the continuation sidecar that an ordinary height bound forgets.

There is one canonical zero option at each checkpoint and `3*32` nonzero
labelled-digit options. Hence the number of nonzero checkpoint words is

```text
(1+3*32)^4-1=97^4-1=88529280.                        (19)
```

Balanced uniqueness makes the map from these words to ordered coefficient
triples injective. Simultaneous sign has no fixed nonzero word, so there are
exactly

```text
(97^4-1)/2=44264640                                  (20)
```

raw unoriented templates before discarding same-sign planes or quotienting
blocker permutations.

## 5. Consequence for the THM-2222 branch

THM-2222 proves that every scalar survivor with `lambda_1>=6` supplies a
primitive sorted triple

```text
1<=d_1<d_2<d_3<=B_6,
B_6=3877322523365316,                                (21)
```

with

```text
measure K_4(d)>=961/6930.                            (22)
```

By (12)--(20), every such triple lies on one of the explicit bounded
relation-and-carry planes (15). Therefore THM-2222's cubic triple box can
be replaced by a finite union of two-dimensional integer slices, each
carrying its four checkpoint digits.

There is also an unconditional profile exclusion before any enumeration.
Use the unnormalized triple

```text
d_j=c_j/13^6,             nu_13(d_j)=lambda_j-6.     (22a)
```

THM-2222 may divide this triple by its common gcd before applying the
extremal problem. The relation supplied for the primitive triple pulls back
with the same coefficients after multiplication by that gcd. Equivalently,
the primitive valuations are `lambda_j-lambda_1`; the common shift cancels
from every comparison below.

For every active coefficient `a_j` in (14), let `r_j` be its least
checkpoint with a nonzero owned digit and put

```text
epsilon_j=nu_13(n_(r_j)) in {0,1}.                  (22b)
```

The last inclusion holds because a nonzero integer of absolute value at
most sixteen is either prime to thirteen or is `+/-13`. Later checkpoint
digits are divisible by two additional powers of thirteen, so

```text
nu_13(a_j)=2r_j+epsilon_j.                           (22c)
```

The active `r_j` are distinct: a checkpoint digit is assigned to only one
blocker. In a vanishing 13-adic sum the smallest term valuation occurs at
least twice. Thus (15) forces distinct active indices `j,k` with

```text
lambda_j-6+2r_j+epsilon_j
 =lambda_k-6+2r_k+epsilon_k.                         (22d)
```

Since `r_j!=r_k`, this implies

```text
1<=|lambda_j-lambda_k|<=7.                           (22e)
```

The deeper of the two blockers owns the earlier checkpoint. For a legal
profile

```text
lambda_1<=lambda_2<lambda_3<=19,     lambda_1>=6,    (22f)
```

condition (22e) is automatic when all three depths are distinct: three
distinct integers in `[6,19]` contain two within distance at most six. If
`lambda_1=lambda_2=a`, however, the only distinct depth is `lambda_3=c`.
Therefore the profiles

```text
(a,a,c),       6<=a<c<=19,       c-a>=8             (22g)
```

are impossible. Their number is

```text
sum_(a=6)^11 (12-a)=6+5+4+3+2+1=21.                (22h)
```

This conclusion uses only THM-2222's checkpoint lower bound and the
degree-sixteen relation gate, not the open extremal inequality `S_4(B_6)`.
It reduces THM-2222's `lambda_1>=6` packet from 455 profiles to 434.

This is the precise nonmixing branch of the relation-and-carry dichotomy:

```text
no bounded cross-checkpoint resonance
  => Selberg mixing already falls below the target;

survival
  => an explicit bounded relation with a disjoint base-169 carry word. (23)
```

The remaining problem is to bound `measure K_4(d)` uniformly on the
surviving relation planes, or to find a hostile plane and retain the
coordinate it exposes. Finiteness alone does not make that task
computationally small. Apart from the 21 profiles in (22g), no closure of
the `lambda_1>=6` branch or of LRC(14) is claimed.

## 6. Exact arithmetic referee

Run

```bash
python3 04-computation/lrc14_four_checkpoint_selberg_carry_gate_thm2223.py
python3 -O 04-computation/lrc14_four_checkpoint_selberg_carry_gate_thm2223.py
```

The companion verifies the first successful degree, every rational value in
(10)--(13), the coefficient-height invoice, the no-carry inequalities, and
the template and profile counts using exact integer and rational arithmetic.
The proof of (9) is the Fourier argument above.

QED.
