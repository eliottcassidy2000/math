---
id: THM-2165
title: "Shunia finite root extraction at the linear exponent K=a+1"
status: >
  PROVED. Shunia's finite modular root formula remains valid after replacing
  the conjectured exponent K=2a^n by the exponentially smaller K=a+1.
  For a>=501 a sharpened cyclic spectral gap beats the exact algebraic
  distance to the adjacent n-th powers, global Kronecker evaluation stays
  below the modulus, and lower radix positions cannot change the quotient.
  Every remaining admissible pair 3<=a<=500 is checked by exact integer
  arithmetic (3,462 pairs), including 61 independent modular-power controls.
  The single coefficient-sum carry at (a,n,t)=(3,2,K+1) is harmless: the
  evaluated quotient polynomial is still below the modulus. No pointwise
  minimality is claimed. Within the additive family K=a+c with c a
  nonnegative integer, c=1 is optimal: c=0 fails at (8,4) and (14,4).
source: opus-2026-07-24-puzzle-atlas
depends_on:
  - THM-2159-shunia-finite-root-spectral-extraction
related:
  - THM-2051
  - THM-2054
script: 04-computation/shunia_linear_exponent_exact_cutoff_opus_20260724.py
output: 05-knowledge/results/shunia_linear_exponent_exact_cutoff_opus_20260724.out
script_sha256: 5add81557bae0de7b490c6b8c8e336919c4560b05ceea8fdbc1b0bef1b5c5096
output_sha256: 3805b9c25ad3a5728380f12048f8b105b6556a0c4891706b8a95623390faa6ce
hash_basis: working-tree bytes (LF)
external:
  - "Joseph M. Shunia, Polynomial quotient rings and Kronecker substitution for deriving combinatorial identities, arXiv:2404.00332v6"
---

# THM-2165 -- Shunia extraction at a linear exponent

## 1. Statement

Let `a,n` be natural numbers such that

```text
a>2,  1<n<=floor(log_2(a))+1,
```

and suppose that `a` is not a perfect `n`-th power. Put

```text
K=a+1,  X=a^K,  M=X^n-a,  R_t=(X+1)^t mod M,
```

where `mod` is the least nonnegative remainder. Then

```text
floor(a^(1/n)) = (R_(K+1) div R_K)-1.                (1)
```

Thus the exponent `2a^n` in Shunia's Conjecture 6.1 can be replaced
uniformly by the linear exponent `a+1`.

No pointwise-minimality claim is made. In the uniform additive family

```text
K=a+c,             c in Z_(>=0),
```

however, the constant `c=1` is optimal: exact computation gives failures
for `c=0` at `(a,n)=(8,4)` and `(14,4)`. In both rows the true root is `1`
and the formula with `K=a` extracts `2`.

## 2. Inherited cyclic walk

We use the notation and the elementary quotient-ring identities of THM-2159.
Put

```text
rho=a^(1/n),       d=floor(rho),

F_t(y)=(1+y)^t mod (y^n-a)
      =sum_(r=0)^(n-1) c_(t,r)y^r.                   (2)
```

All `c_(t,r)` are nonnegative integers. After scaling by `rho^r`, their
normalized vector is the lazy directed walk on `Z/nZ`:

```text
p_(t,r)=rho^r c_(t,r)/(1+rho)^t,

|n p_(t,r)-1|<=eta_t,
eta_t=(n-1)theta^t,                                  (3)

theta=max_(1<=j<n)|1+rho exp(2 pi i j/n)|/(1+rho).
```

The exact spectral gap is

```text
theta^2
 <=1-[2rho(1-cos(2pi/n))/(1+rho)^2].                 (4)
```

The distance from `rho` to the adjacent integers has the stronger
factorization bound

```text
delta=min(rho-d,d+1-rho)
 >=1/[n(d+1)^(n-1)]
 >=1/[n(2rho)^(n-1)].                                (5)
```

Indeed both positive integers `a-d^n` and `(d+1)^n-a` are at least one;
factoring the corresponding differences of powers gives a sum of `n` terms,
each at most `(d+1)^(n-1)`.

## 3. Uniform spectral tail for `a>=501`

We prove

```text
eta_K < 1/[8n 2^(n-1)a] <= delta/(8rho).             (6)
```

All elementary logarithmic comparisons below can be made rational: use

```text
log 2<7/10,
```

which follows by truncating the positive exponential series at `7/10`.
The displayed target integers are bounded above by the indicated powers of
two; no floating-point estimate is a proof dependency.

For `n>=6`, the chord bound for sine on `[0,pi/6]` gives

```text
1-cos(2pi/n)=2sin^2(pi/n)>=18/n^2.                   (7)
```

We split at the integer-root boundary. This is the sharpening which makes a
linear exponent sufficient.

### The `d=1` branch

Here `1<rho<2`, so

```text
rho/(1+rho)^2>=2/9.
```

Since `a>=501`, this branch has `n>=9`. Equations (4), (7), and
`sqrt(1-x)<=exp(-x/2)` give

```text
theta^K<=exp(-4a/n^2).                                (8)
```

For `n=9`, at the lower endpoint `a=501`,

```text
4a/n^2=2004/81>24
  >log(8n(n-1)2^(n-1)a).
```

For the last inequality, the logarithm's argument is below `2^27`, so its
logarithm is below `189/10`.

For `n>=10`, the lower endpoint is `a=2^(n-1)`. At `n=10`,

```text
4*2^9/10^2=20.48
 >log(8*10*9*2^9*2^9).
```

The logarithm's argument is below `2^28`, so the right side is below `19.6`.
The left side is multiplied by more than `3/2` when `n` increases. The
right side increases by

```text
log((n+1)/(n-1))+2log 2<2
```

(`(n+1)/(n-1)<=11/9` and the exponential series gives
`log(11/9)<1/4`). This proves the inequality by induction. For fixed `n`,
`4a/n^2-log(Ca)` is increasing throughout this range. Hence (6) holds on
the `d=1` branch.

### The `d>=2` branch

Now `rho>=2`, and therefore

```text
[rho/(1+rho)]^2>=4/9.
```

Using (4) before discarding this factor gives

```text
theta^K
 <=exp[-K rho(1-cos(2pi/n))/(1+rho)^2]
 <=exp[-(4/9)(1-cos(2pi/n))rho^(n-1)].               (9)
```

For `n>=6`, (7) makes the last exponent at least

```text
8rho^(n-1)/n^2.                                      (10)
```

At `a=501`, the cases `n=6,7,8` have respectively

```text
rho<3,  rho<5/2,  rho<11/5.
```

Thus the quantities in (10) are greater than

```text
4008/108,       4008/(245/2),       20040/704,
```

all greater than `28`; their target logarithms are less than `18`.

For `n>=9`, this branch has `a>2^n`, so (10) is greater than

```text
8*2^(n-1)/n^2.
```

At `n=9` this exceeds `25`, while the target logarithm at the lower endpoint
is less than `20` (use the upper endpoint expression
`8n(n-1)2^(2n)<2^28`). The same multiplicative-versus-additive induction proves
the result for every larger `n`.

The remaining cycle sizes are stronger. At `a=501`:

```text
n=5: rho<4,  1-cos(2pi/5)>2/3,
     exponent in (9)>1002/27>37;

n=4: rho<5,  1-cos(pi/2)=1,
     exponent in (9)>2004/45>44;

n=3: exponent >=(2/3)rho^2>1002/24>41;

n=2: theta=(rho-1)/(rho+1)<exp(-2/rho),
     theta^K<exp(-2rho), with 2rho>44.
```

The respective target logarithms are below `15,15,12,10`. For each fixed
`n`, the power lower bound minus the logarithm increases from the displayed
endpoint. This completes the proof of (6).

From (3), (5), and (6),

```text
|c_(K,n-2)/c_(K,n-1)-rho|
 <=2rho eta_K/(1-eta_K)
 <delta/2.                                           (11)
```

The last-coordinate recurrence

```text
c_(K+1,n-1)=c_(K,n-1)+c_(K,n-2)
```

therefore yields

```text
d+1+delta/2
 <c_(K+1,n-1)/c_(K,n-1)
 <d+2-delta/2.                                       (12)
```

## 4. Global no-wrap and the lower radix positions

Let `S_t=sum_r c_(t,r)`. As in THM-2159,

```text
S_t<=(1+rho)^t.
```

For `t=K,K+1` and `a>=501`, `rho<=sqrt(a)` gives

```text
S_t<=(2sqrt(a))^(a+2)<a^(a+1)=X.                     (13)
```

The last inequality follows already for `a>=16` from
`2^(a+2)<4^a<=a^(a/2)`. Hence

```text
0<F_t(X)
 <=(X-1)X^(n-1)
 <X^n-a=M.                                           (14)
```

Since `X^n=a mod M`, (2) and (14) imply

```text
R_t=F_t(X),            t=K,K+1.                      (15)
```

This is a global no-wrap statement. It is stronger than needed to demand
that each coefficient be a base-`X` digit.

To control the lower positions, write

```text
F_t(X)=c_(t,n-1)X^(n-1)(1+e_t).
```

Equation (3), with (6), gives

```text
0<=e_t
 <=sum_(j=1)^(n-1) 2(rho/X)^j
 <=4rho/X
 <=2na/X.                                            (16)
```

Combining (12), (15), and (16), the difference between
`R_(K+1)/R_K` and the leading-coefficient ratio is at most

```text
8na^2/X.                                             (17)
```

The hypothesis on `n` gives `2^(n-1)<=a`; also `rho^(n-1)<a`. Thus (5)
implies

```text
delta/2>1/(2na^2).
```

Since `n<=a` and `X=a^(a+1)`,

```text
8na^2/X<1/(2na^2)                                   (18)
```

for every `a>=501`. Equations (12), (17), and (18) now give

```text
d+1<R_(K+1)/R_K<d+2.
```

Natural-number division and subtraction of one prove (1) throughout the
infinite tail.

## 5. Exact finite cutoff

The companion checks every admissible pair

```text
3<=a<=500,
2<=n<=floor(log_2(a))+1
```

using Python integers only. There are exactly `3,462` pairs. For each pair
it:

1. constructs both quotient-ring coefficient vectors by recurrence;
2. verifies the two strict leading-coordinate inequalities in (12);
3. evaluates both quotient polynomials at `X=a^(a+1)`;
4. verifies both values lie strictly between zero and `M`;
5. verifies the final natural-number quotient; and
6. on 61 hostile and boundary pairs, independently compares both values with
   Python's modular exponentiation.

It also independently replays the two exact `K=a` failures, establishing the
stated optimality inside the nonnegative additive family.

Normal and optimized runs reproduce

```text
row digest sha256:
1a99b0f13dcf475205035e6d2aaab441a596f55f6b895c5ba1e6124b2bd458c1
```

and the stored transcript exactly.

There is one instructive boundary:

```text
(a,n,t)=(3,2,K+1)
```

has `sum_r c_(t,r)>=X`. This refutes the stronger digitwise-no-carry claim,
but direct exact evaluation still gives `0<F_t(X)<M`. The proof uses the
correct global no-wrap predicate, not the false digitwise surrogate.

The finite cutoff plus the infinite tail prove (1) for every admissible pair.
QED.
