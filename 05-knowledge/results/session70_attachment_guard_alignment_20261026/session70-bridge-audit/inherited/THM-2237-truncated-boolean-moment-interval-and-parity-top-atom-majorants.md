---
id: THM-2237
title: "Truncated Boolean moment interval and parity top-atom majorants"
status: >
  PROVED. For every random variable K supported on {0,...,p}, its binomial
  moments through degree p-1 leave exactly one affine degree of freedom.
  The complete feasible interval for the top atom P(K=p) is given by
  alternating binomial inversion. Its exact upper endpoint is the minimum
  of ceil(p/2) explicit pointwise majorants, each supported only at K=p
  and one count of the opposite parity. The even- and odd-parity laws on
  the Boolean p-cube have identical proper marginals and realize top atoms
  2^(-(p-1)) and 0, so no statistic of moments through degree p-1 can
  recover the missing top Walsh character. For p=4 the two majorants are
  exactly THM-2222's adaptive cubic pair. This is a finite moment theorem,
  not the open four-checkpoint comb estimate or a proof of LRC(14).
source: codex-2026-07-24-parity-top-atom
related:
  - THM-2210-nested-binomial-minorant-and-adaptive-moment-lp-hierarchy
  - THM-2222-scalar-transfer-parity-tower-and-four-checkpoint-survivor-reduction
  - THM-2225-dyadic-critical-run-extractors-and-cyclic-checksum-shell-bisection
  - THM-2236-pointwise-nested-binomial-minorants-and-cubic-vertex-fan
---

# THM-2237 -- truncated Boolean moments and the top atom

THM-2222 found two sharp cubic majorants for a four-checkpoint all-hit
event. The pair is the first nontrivial member of a complete
parity-indexed family.

## 1. The one-dimensional truncated moment fibre

Fix an integer `p>=1`. Let `K` be supported on `{0,...,p}`, and write

```text
pi_k=P(K=k),
M_r=E binom(K,r),                    0<=r<=p-1.       (1)
```

For `0<=k<=p-1`, define

```text
A_k
 =sum_(r=k)^(p-1)(-1)^(r-k)binom(r,k)M_r.            (2)
```

If the missing top atom is denoted by `t=pi_p`, binomial inversion gives

```text
pi_k(t)=A_k+(-1)^(p-k)binom(p,k)t,
                                             0<=k<p. (3)
```

Indeed, full binomial inversion is

```text
pi_k=sum_(r=k)^p(-1)^(r-k)binom(r,k)M_r,
```

and `M_p=pi_p=t`. Thus every law with the prescribed truncated packet lies
on the affine line (3), and every nonnegative point on that line is such a
law.

Put

```text
L=max(
    0,
    max_(0<=k<p, p-k even)
      [-A_k/binom(p,k)]
  ),

U=min_(0<=k<p, p-k odd)
      [A_k/binom(p,k)].                              (4)
```

The maximum over an empty set in `L` is ignored. A packet is realizable if
and only if

```text
L<=U.                                                (5)
```

When it is realizable, its complete top-atom range is exactly

```text
L<=P(K=p)<=U.                                        (6)
```

To prove this, note that the inequalities in (3) with `p-k` even are
precisely the lower bounds in (4), while those with `p-k` odd are precisely
the upper bounds. The constraint `pi_p=t>=0` supplies the remaining lower
bound. Equation (3) preserves total mass because `M_0=1`, so there are no
other conditions. This proves (5)--(6).

## 2. The exact parity-indexed upper majorants

For every `k<p` of parity opposite to `p`, define

```text
Q_(p,k)(x)
 =1/binom(p,k)
   sum_(r=k)^(p-1)
     (-1)^(r-k)binom(r,k)binom(x,r).                 (7)
```

On the integer support `0<=x<=p`, this polynomial has the sparse value
table

```text
Q_(p,k)(x)
 =1/binom(p,k)                 if x=k,
 =1                            if x=p,
 =0                            otherwise.            (8)
```

For `x<k` the sum is empty. For `k<=x<p`, use

```text
binom(r,k)binom(x,r)
 =binom(x,k)binom(x-k,r-k)
```

and the alternating binomial sum. It is one when `x=k` and zero when
`k<x<p`. At `x=p`, the full sum through `r=p` vanishes; deleting its last
term leaves

```text
-(-1)^(p-k)binom(p,k)=binom(p,k),
```

because `p-k` is odd. This proves (8).

In particular, every `Q_(p,k)` is a pointwise majorant of `1_(x=p)`.
Equations (2), (7), and (8) give

```text
E Q_(p,k)(K)
 =A_k/binom(p,k)
 =pi_p+pi_k/binom(p,k).                              (9)
```

Combining (4), (6), and (9) proves the exact adaptive ceiling

```text
max{P(K=p): the moments M_0,...,M_(p-1) are fixed}

 =min_(0<=k<p, p-k odd) E Q_(p,k)(K).                (10)
```

Thus the upper endpoint is not merely LP-computable: it is the minimum of
exactly `ceil(p/2)` explicit sparse majorants. One of them binds precisely
when the corresponding opposite-parity atom vanishes.

For `p=4`, the two choices `k=3,1` are

```text
Q_(4,3)(K)=binom(K,3)/4,

Q_(4,1)(K)
 =[K-2binom(K,2)+3binom(K,3)]/4.                    (11)
```

These are exactly the two majorants in THM-2222.

## 3. The missing coordinate is top parity

Let `X=(X_1,...,X_p)` be uniform on either parity coset of the Boolean cube:

```text
X_1+...+X_p = epsilon mod 2.                         (12)
```

Every proper coordinate marginal is uniform. Therefore, for
`K=sum_i X_i`, both parity laws have the same truncated packet

```text
M_r=binom(p,r)2^(-r),                 0<=r<p.         (13)
```

But the all-one word belongs only to the coset with
`epsilon=p mod 2`. The two top atoms are consequently

```text
P(K=p)=2^(-(p-1))             and             0.     (14)
```

By (6), this common packet has the full feasible interval

```text
[0,2^(-(p-1))].                                      (15)
```

Thus even the entire collection of proper Boolean marginals, and a
fortiori the moments through degree `p-1`, need not determine the top
Walsh character

```text
(-1)^(X_1+...+X_p).                                  (16)
```

Indeed,

```text
E(-1)^K=sum_(r=0)^p(-2)^r M_r,
```

so its only coordinate absent from the packet is
`M_p=P(K=p)`.

This is an information-theoretic obstruction, not a weakness of a
particular polynomial basis.

## 4. Consequence for the four-checkpoint frontier

At `p=4`, (15) is `[0,1/8]`. Since

```text
1/8<961/6930,
```

the generic parity ambiguity does not itself defeat THM-2222's target.
It does prove that no further manipulation of the same moments
`M_0,...,M_3` can improve the exact adaptive ceiling (11). Any strict
advance must retain information absent from that packet: special comb
dynamics, cross-checkpoint arithmetic, or the top character together with
its carry/owner realization. THM-2223 and THM-2224 remain only reserved
namespaces for possible implementations of such a sidecar. QED.
