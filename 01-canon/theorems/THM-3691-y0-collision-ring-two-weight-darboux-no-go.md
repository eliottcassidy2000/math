---
id: THM-3691
title: "Two-weight Darboux no-go in the y=0 collision ring"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  In the y=0
  collision ring R of THM-3686, grade
  by wt(x)=1, wt(z)=-2.  If each of P,Q has at most two active nonconstant
  weight components, then {P,Q} is not a nonzero constant.  The proof is
  degree-unbounded: support-sum geometry leaves only crossed equal-gap
  weights; opposite-sign and weight-zero edges remove all nonpositive
  parameters; endpoint valuations force one negative weight to be -2; and
  the surviving even case factors through the homogeneous weight (-2,+1)
  Wronskian already excluded by THM-3686.  Consequently any Darboux pair in
  R, hence any planar Jacobian counterexample in this collision lane, must
  use at least three active grading weights in one output.  Three-weight
  pairs, JC(2), and arbitrary quartic C3 data remain OPEN.
source: jc-quartic-c3-construction / 2026-08-22
audit: >
  PASS -- root independently checked the crossed-support parameterization,
  all a/b nonpositive and weight-zero edges, endpoint valuations, UFD power
  exponents in the odd/even split, and the exact factorization through the
  actual homogeneous weight (-2,+1) bracket.
depends_on:
  - THM-3686-y0-collision-normalization-and-bracket-anatomy
related:
  - THM-3689-fully-transverse-two-by-two-sparse-support-gate
  - THM-3690-complete-normalized-two-by-two-sparse-support-closure
script: 04-computation/jacobian_y0_two_weight_darboux_no_go_thm3691.py
output: 05-knowledge/results/jacobian_y0_two_weight_darboux_no_go_thm3691.out
script_sha256: d49f33a8e8b58c3f997bb86790a377f3f82c9d7d699a6bd9ed0f85130786ec5c
output_sha256: d86d1a52346cb662924a91532bf84afbd91bb9073a8582a35461d1f6fe7cc164
hash_basis: LF-normalized bytes
---

# THM-3691 -- two active weights per output cannot carry the retained collision

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Let `R` be the collision-retaining target subalgebra from THM-3686:

```text
R=C[A,B,C] subset C[x,z],
A=3z(2-x^2z),
B=2xz(2-x^2z),
C=x(1-x^2z).                                           (1)
```

Give the source the grading

```text
wt(x)=1,                    wt(z)=-2.                   (2)
```

We prove that if `P,Q in R` each have at most two nonconstant homogeneous
components, then

```text
{P,Q}:=P_xQ_z-P_zQ_x                                   (3)
```

is not a nonzero constant.  Additive scalar components are bracket-invisible
and are never counted as active weights.

## 1. Weight modules and the commutation lemma

Put

```text
u=1-x^2z,                    h=1-u^2.                   (4)
```

Then

```text
A=3x^-2h,                    B=2x^-1h,        C=xu.    (5)
```

Every homogeneous element of weight `r` is therefore

```text
x^r p(u).                                                (6)
```

For two such elements, direct differentiation gives

```text
{x^r p(u),x^s q(u)}
 =x^(r+s+1)[s p'(u)q(u)-r p(u)q'(u)].                 (7)
```

The target monomials in `(5)` impose two endpoint divisibilities:

```text
r=-R<0  =>  h^ceil(R/2) divides p,
r= R>0  =>  u^R divides p.                             (8)
```

For the first implication, a monomial `A^iB^jC^k` of weight `-R` has
`2i+j=R+k>=R`, hence `i+j>=ceil(R/2)`.  For the second, weight `R` gives
`k=R+2i+j>=R`.

We will also use the following exact commutation rule.  If the bracket in
`(7)` vanishes, then

```text
s p'/p = r q'/q.                                       (9)
```

Consequently:

1. if `r,s` are strictly opposite in sign, `(9)` makes a positive product
   of powers of `p,q` constant; both polynomials would be constant, contrary
   to `(8)`;
2. if exactly one weight is zero and the other is nonzero, the weight-zero
   coefficient polynomial must be constant; and
3. if `r,s` have the same nonzero sign, unique factorization makes `p,q`
   powers of one common polynomial, with exponents determined by
   `gcd(|r|,|s|)`.

These statements allow arbitrary nonzero scalar factors, which will be
absorbed below.

## 2. The support-sum grid leaves one crossed pattern

One-weight against one-weight is the all-degree homogeneous no-go of
THM-3686.  The same argument handles one weight against two: whichever
bracket contributes the constant is a unique homogeneous constant bracket.
Thus suppose both active supports have size two:

```text
supp(P)={r_0<r_1},              supp(Q)={s_0<s_1}.      (10)
```

Their four bracket weights are `r_i+s_j+1`.  The lower-left and upper-right
weights are unique extrema.  If the constant occurs at either extremum, or
at only one cross, it comes from one homogeneous bracket and is impossible.
The only survivor has both crosses in the constant bucket:

```text
r_0+s_1=-1,                    r_1+s_0=-1.              (11)
```

In particular the two support gaps agree.  Write the weights as

```text
supp(P)={-a,b-1},               supp(Q)={-b,a-1},
a+b>1.                                                   (12)
```

The unique extreme brackets must vanish:

```text
{P_-a,Q_-b}=0,                 {P_(b-1),Q_(a-1)}=0,    (13)
```

and the crossed constant equation is

```text
{P_-a,Q_(a-1)}+{P_(b-1),Q_-b}=1.                       (14)
```

If `a<=0`, then `b>=2`, and the second pair in `(13)` has strictly opposite
nonzero weights `b-1>0` and `a-1<0`.  This contradicts the commutation lemma.
The case `b<=0` is symmetric.  Hence

```text
a,b>=1.                                                  (15)
```

## 3. The weight-zero edges

Suppose `a=1<b`.  The upper extreme in `(13)` pairs weight `b-1` with
weight zero.  The commutation lemma says that the weight-zero component is a
scalar, so it is inactive and may be removed.  Equation `(14)` then becomes
one homogeneous constant bracket, impossible by THM-3686.  The case
`b=1<a` is symmetric.

If `a=b=1`, write

```text
P=x^-1p(u)+r(u),                 Q=x^-1q(u)+s(u).       (16)
```

The lower extreme equation makes `p,q` proportional, while the two
weight-zero components commute automatically.  Equation `(14)` is a sum of
terms of the form

```text
p s' - r' q.                                            (17)
```

Both `p` and `q` are divisible by `h` by `(8)`.  Thus `(17)` is divisible by
the nonunit `h` and cannot equal one.  We have reduced to

```text
a,b>=2.                                                  (18)
```

## 4. Endpoint valuation forces one negative weight to be `-2`

Write the four coefficient polynomials as

```text
P_-a=x^-a p(u),             P_(b-1)=x^(b-1) r(u),
Q_-b=x^-b q(u),             Q_(a-1)=x^(a-1) s(u).      (19)
```

If `a>=3`, `(8)` makes `p` divisible by `h^2`.  Therefore both `p` and `p'`
vanish at `u=+1,-1`, so the first crossed bracket in `(14)` vanishes at both
endpoints.  If `b>=3`, the same holds for the second crossed bracket through
`q`.  Hence `a,b>=3` would make the left side of `(14)` vanish where the
right side is one.  Therefore

```text
min(a,b)=2.                                              (20)
```

Swap the two outputs, and change one sign to preserve `(14)`, if necessary.
We may assume

```text
a=2,                         b>=2.                       (21)
```

## 5. Odd `b` dies at the same endpoint

The lower commutation equation in `(13)` and unique factorization give

```text
p=alpha H^(2/d),              q=beta H^(b/d),
d=gcd(2,b).                                                (22)
```

If `b` is odd, then `d=1` and `p=alpha H^2`.  Because `h|p` and `h` is
squarefree, `h|H`; in fact `h^2|p`.  Both the `p`-cross and the `q=beta H^b`
cross in `(14)` vanish at `u=+1,-1`.  This is again impossible.  Thus

```text
b=2m is even.                                            (23)
```

## 6. The even case factors through the homogeneous no-go

When `b=2m`, equations `(13)` give polynomial cores, up to nonzero scalar
factors,

```text
p proportional H,             q proportional H^m,
r proportional K^(b-1),       s proportional K.         (24)
```

Here `h|H` and `u|K` by `(8)`.  The first crossed bracket is a scalar
multiple of

```text
D(H,K)=H'K+2HK'.                                        (25)
```

For the normalized cores, the second crossed bracket factors exactly as

```text
-b(r')q-(b-1)r(q')
=-(b-1)m H^(m-1)K^(b-2) D(H,K).                         (26)
```

Therefore the complete constant equation `(14)` has the form

```text
D(H,K)[c_1+c_2 H^(m-1)K^(b-2)]=1                       (27)
```

for constants `c_1,c_2` determined by the four nonzero scalar factors.

But `D(H,K)` is precisely the homogeneous bracket

```text
{x^-2H(u),xK(u)}.                                       (28)
```

Both entries lie in `R` by their origin in `(19)`.  THM-3686 proves in all
degrees that `(28)` is never a nonzero constant.  If it is zero, `(27)` is
zero; if it is nonconstant, it is a nonunit of `C[u]` and cannot divide one.
This contradiction closes the last support pattern.

We have proved

```text
if |supp(P)|<=2 and |supp(Q)|<=2, then {P,Q} notin C*.  (29)
```

## 7. Exact frontier

Every element of `R` takes the same value on the collision pair of THM-3686.
Thus `{P,Q}=1` in `R` would still be a planar Jacobian counterexample.  The
present theorem raises the first admissible grading complexity to

```text
at least three active weights in P or at least three active weights in Q. (30)
```

This is an all-degree statement; it is independent of ordinary target
degree and complements the finite target-degree-three gate in THM-3686 and
the source-monomial support gates THM-3689/3690.  It does not exclude a
`2x3`, `3x2`, or larger grading support, prove `JC(2)`, or settle arbitrary
quartic `C3` data.

## 8. Exact reproduction

Run

```bash
python3 -B 04-computation/jacobian_y0_two_weight_darboux_no_go_thm3691.py
python3 -B -O 04-computation/jacobian_y0_two_weight_darboux_no_go_thm3691.py
```

Both modes byte-match the stored transcript.  The companion checks the
Laurent bracket law, monomial divisibilities through target degree ten, the
crossed equal-gap support grammar on a hostile integer window, odd endpoint
controls, and the even factorization `(26)`.  The proof of `(29)` is the
degree-unbounded argument above.

**QED.**
