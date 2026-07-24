---
id: THM-2159
title: "Shunia finite integer-root extraction by cyclic spectral mixing"
status: >
  PROVED. Shunia's conjectured finite modular formula for the integer n-th
  root holds for every stated admissible pair. Multiplication by 1+y modulo
  y^n-a becomes a lazy directed walk on the n-cycle after the coefficients
  are scaled by rho^r. Its nonconstant Fourier modes contract far enough by
  K=2a^n to put the ratio of the last two quotient-ring coefficients strictly
  between floor(rho) and floor(rho)+1. Evaluation at X=a^K is lossless, gives
  the least modular remainders, and changes the consecutive ratio by less
  than half the algebraic distance to the adjacent integers. Hence natural
  division extracts floor(rho)+1 and subtracting one gives Nat.nthRoot.
  Independently verified by exact integer arithmetic on every admissible pair
  3<=a<=12 (21 cases). The LRC comparison is methodological only: OPEN-Q-108
  lacks the uniform arithmetic quantum used here.
source: opus-2026-07-24-puzzle-atlas
depends_on: []
related:
  - THM-2051
  - THM-2054
script: 04-computation/shunia_finite_root_spectral_extraction_opus_20260724.py
output: 05-knowledge/results/shunia_finite_root_spectral_extraction_opus_20260724.out
script_sha256: 1407c83122d2ba122b3f6935aef41240dc30ec50ed0fcd2304c4ffcacb0a472f
output_sha256: 337ceb56c1ae5963cc726ef0f6b96775b5d2ad0d312b684ce7c4b2e4f557c238
external:
  - "Joseph M. Shunia, Polynomial quotient rings and Kronecker substitution for deriving combinatorial identities, arXiv:2404.00332v6"
---

# THM-2159 -- Shunia finite-root extraction

## 1. Statement

Let `a,n` be natural numbers such that

```text
a>2,  1<n<=floor(log_2(a))+1,
```

and suppose that `a` is not a perfect `n`-th power. Put

```text
K=2a^n,  X=a^K,  M=X^n-a,  R_t=(X+1)^t mod M.
```

Here `mod` denotes the least nonnegative remainder and `div` natural-number
division. Then

```text
floor(a^(1/n)) = (R_(K+1) div R_K)-1.
```

This is exactly Shunia's Conjecture 6.1 in natural-number arithmetic.

## 2. Quotient-ring coefficients

Write `rho=a^(1/n)` and `d=floor(rho)`. For `t>=0`, define the unique
nonnegative coefficient vector

```text
F_t(y)=(1+y)^t mod (y^n-a)
      =sum_(r=0)^(n-1) c_(t,r)y^r.                    (1)
```

Explicitly,

```text
c_(t,r)=sum_(0<=j<=t, j=r mod n) binom(t,j)a^((j-r)/n).  (2)
```

Multiplication by `1+y` gives

```text
c_(t+1,0)=c_(t,0)+a c_(t,n-1),
c_(t+1,r)=c_(t,r)+c_(t,r-1),  1<=r<n.                (3)
```

Scale the coefficients by

```text
b_(t,r)=rho^r c_(t,r).
```

After division by `(1+rho)^t`, recurrence (3) is the transition rule of the
lazy directed random walk on `Z/nZ` which stays put with probability
`1/(1+rho)` and advances one step with probability `rho/(1+rho)`. If

```text
p_(t,r)=b_(t,r)/(1+rho)^t,
zeta=exp(2 pi i/n),
```

then the roots-of-unity filter gives

```text
p_(t,r)=1/n sum_(j=0)^(n-1)
  ((1+rho zeta^j)/(1+rho))^t zeta^(-jr).              (4)
```

Let

```text
theta=max_(1<=j<n) |1+rho zeta^j|/(1+rho),
eta_t=(n-1)theta^t.
```

Equation (4) yields the pointwise mixing estimate

```text
|n p_(t,r)-1|<=eta_t.                                 (5)
```

This is the spectral mechanism hidden by the modular formula.

## 3. Why the explicit exponent `K=2a^n` suffices

We record the elementary estimate used below:

```text
eta_K < 1/(8 n a^n).                                  (6)
```

Here are full details of the coarse bound. For `j!=0`,

```text
|1+rho zeta^j|^2/(1+rho)^2
 <=1-2rho(1-cos(2pi/n))/(1+rho)^2
 <=1-4/(a n^2).                                       (7)
```

The last inequality uses

```text
1-cos(2pi/n)=2sin^2(pi/n)>=8/n^2,
2rho/(1+rho)^2>=1/(2rho)>=1/(2a).
```

Thus

```text
theta<=exp(-2/(a n^2)),
eta_K<=(n-1)exp(-4a^(n-1)/n^2).                       (8)
```

The hypothesis on `n` says `a>=2^(n-1)`. Comparing logarithms in (8) proves
(6):

```text
4a^(n-1)/n^2 > log(8n(n-1)a^n).                      (9)
```

For completeness, the only small pairs not covered immediately by the
monotone right-minus-left comparison in (9) are:

```text
(n,a)=(2,3),(2,5),(2,6),(3,4).
```

For `n=2`, `theta=(rho-1)/(rho+1)`; it is `<1/3` at `a=3`, `<1/2` at
`a=5`, and the `a=6` inequality is direct. For `(n,a)=(3,4)`,
`theta<2/3`. These bounds give (6) by raising to `K=18,50,72,128`,
respectively. For fixed `n`, the difference between the left and right
sides of (9) has derivative

```text
(4(n-1)a^(n-1)/n^2-n)/a,
```

which is positive at, and hence beyond, each lower endpoint used below.
For `n=2,a>=7` and `n=3,a>=5`, (9) is already strict at that endpoint.
For `n>=4`, it is therefore enough to put `a=2^(n-1)`. At this value the
left side is

```text
4*2^((n-1)^2)/n^2,
```

whereas the right side is

```text
log(8n(n-1))+n(n-1)log(2)<(3/2)n^2.
```

The displayed left side exceeds `(3/2)n^2` at `n=4` (`128>24`), and
the inequality persists by induction because its numerator is multiplied
by `2^(2n-1)` when `n` increases by one, while the denominator changes
only from `n^2` to `(n+1)^2`. This proves (6) without an asymptotic
assumption.

Because `a` lies strictly between the consecutive powers `d^n` and
`(d+1)^n`,

```text
delta=min(rho-d,d+1-rho)>=1/(n a^(n-1)).              (10)
```

Indeed, factor the two differences of `n`-th powers; each resulting sum has
`n` terms and every term is at most `a^(n-1)`. By (5)--(6),

```text
c_(K,n-2)/c_(K,n-1)
 =rho p_(K,n-2)/p_(K,n-1),

|c_(K,n-2)/c_(K,n-1)-rho|
 <=2rho eta_K/(1-eta_K)
 <1/(2n a^(n-1))<=delta/2.                            (11)
```

Consequently,

```text
d < c_(K,n-2)/c_(K,n-1) < d+1.                       (12)
```

Using (3) in the last coordinate,

```text
d+1 < c_(K+1,n-1)/c_(K,n-1) < d+2,                  (13)
```

and the ratio in (13) stays at distance at least `delta/2` from either
endpoint.

## 4. Lossless Kronecker evaluation

Let `S_t=sum_r c_(t,r)`. Formula (2) gives

```text
S_t<=sum_(j=0)^t binom(t,j)rho^j=(1+rho)^t.           (14)
```

For `a>=3`, `rho<=sqrt(a)` and

```text
(1+rho)/a <=(1+sqrt(a))/a <11/12.
```

Since `K>=2a^2`,

```text
(1+rho)^(K+1)/a^K
 <=a(11/12)^(2a^2)<1.                                (15)
```

The last function is decreasing for `a>=3`, and its value at `a=3` is
`3(11/12)^18<1`. Hence `S_K,S_(K+1)<X`.

Evaluate (1) at `X=a^K`. Since `X^n=a mod M`,

```text
F_t(X)=(X+1)^t mod M.                                 (16)
```

Moreover, using (14)--(15),

```text
0<F_t(X)
 <=S_t X^(n-1)
 <=(X-1)X^(n-1)
 <X^n-a=M,                                            (17)
```

because `X^(n-1)>a`. Thus (16) is not merely a congruence:

```text
R_t=F_t(X),  t=K,K+1.                                 (18)
```

There is no hidden wrap or coefficient carry.

## 5. Lower digits cannot change the extracted integer

From (5)--(6),

```text
c_(t,n-1)
 >=(1+rho)^t/(2n rho^(n-1)),  t=K,K+1.               (19)
```

Therefore

```text
F_t(X)=c_(t,n-1)X^(n-1)(1+e_t),
0<=e_t<=2na/X.                                        (20)
```

Combine (13), (18), and (20):

```text
R_(K+1)/R_K
 =c_(K+1,n-1)/c_(K,n-1) * (1+e_(K+1))/(1+e_K).       (21)
```

The perturbation from the second factor is at most

```text
8na^2/X < 1/(2na^(n-1)) <= delta/2.                  (22)
```

The strict inequality is immediate from `X=a^(2a^n)`; its smallest
admissible case is `a=3,n=2`. More uniformly,
`n<=a` gives `K-n-1>=2a^2-a-1>=14`, so
`X/a^(n+1)>=a^14>16a^2>=16n^2`. Equations (13) and (22) now imply

```text
d+1 < R_(K+1)/R_K < d+2.
```

Taking natural-number division and subtracting one proves

```text
R_(K+1) div R_K - 1=d=floor(a^(1/n)).
```

QED.

## 6. Exact companion and scope

The companion

```text
04-computation/shunia_finite_root_spectral_extraction_opus_20260724.py
```

checks every admissible pair `3<=a<=12` (21 pairs) using integers only. It
independently constructs the quotient-ring coefficients, verifies both least
remainders against modular exponentiation, checks the last-coordinate ratio,
and checks the final natural-number quotient. Stored transcript:

```text
05-knowledge/results/shunia_finite_root_spectral_extraction_opus_20260724.out
```

All rows pass.

The reusable pattern is:

```text
dominant mode + explicit contraction + arithmetic quantum
  => a finite exact extraction from an asymptotic ratio.
```

THM-2051 and THM-2054 already provide LRC analogues of the first two pieces.
OPEN-Q-108 has no uniform denominator for `meas(G_C)`, so the third piece is
absent; this theorem does not prove LRC(14) by analogy.
