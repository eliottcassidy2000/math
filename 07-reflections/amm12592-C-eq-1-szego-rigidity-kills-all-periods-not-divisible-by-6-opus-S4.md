---
source: opus-2026-07-31-S4 (AMM 12592 minimal-C frontier; executing HYP-9061 sec 3's gamma=0 test)
status: >
  SUPERSEDED + CORRECTED by THM-3342 and MISTAKE-368.  The bounded-additive
  impossibility survives and is strengthened to T(n)-n != o(n).  The former
  inference C*>1 was invalid because endpoint nonattainment does not separate
  an infimum.  In the bounded-alphabet proof, apply Szego to (1-p)g (analytic
  across p=1), then recover rationality/eventual periodicity of g; a bare
  meromorphic pole is not itself the clean natural-boundary hypothesis.
tags: [amm12592, minimal-C, biased-coin, extractor, szego, carlson, rigidity, deficit-flow, HYP-9061, lower-bound]
related: [HYP-9061, THM-2160, THM-2225]
---

# AMM 12592, C = 1: Szego rigidity kills every period 6 does not divide (opus-S4)

> **Current correction.** THM-3342 canonizes the stronger theorem that no
> fixed fair extractor has `T(n)=n+o(n)`.  It does **not** imply `C*>1`: a
> family of different extractors with slopes `1+epsilon` could still have
> infimum one.  MISTAKE-368 records the quantifier failure.  The historical
> bounded-depth argument below should be read through the repaired analytic
> object `F=(1-p)g`, or simply as the bounded special case of THM-3342.

## Setup (HYP-9061 sec 1)

A deadline-`(Cn+D)` exactly-fair extractor exists iff there are decided-tree polynomials `W_m, V_m`
(integer Bernstein coefficients `w_{m,k} in {0,...,C(d_m,k)}`, values in `[0,1]` on `(0,1)`) of degree
`<= d_m = (C-1)m + D - 1` with

```
   sum_{m>=1} p^m q W_m(p)  +  sum_{m>=1} q^m p V_m(p)  =  1/2      for all p in (0,1),   q = 1-p.   (S)
```

`C* = 1 + gamma*`, `gamma* =` minimal degree-growth rate. **This note treats `C = 1`, i.e. `gamma = 0`**:
then `d_m = D-1` is CONSTANT, so each `W_m, V_m` is one of FINITELY MANY polynomials of degree `<= D-1`.

## Step 1 — the spine generating function has finite-set integer coefficients

Let `g(p) := sum_{m>=1} p^m W_m(p) = sum_N c_N p^N`. Writing `W_m(p) = sum_{k} a_{m,k} p^k` (ordinary
coefficients; integers, since the Bernstein coefficients are integers and `q^{d-k}=(1-p)^{d-k}` expands
integrally), `c_N = sum_{k=0}^{D-1} a_{N-k,k}` is a bounded sum of integers each from a FINITE set. Hence
`{c_N}` takes finitely many values.

## Step 2 — `g` is rational (repaired Szego step)

Near `p=1`, `S_2(p):=sum_{m>=1}q^m pV_m(p)` is analytic in
`|1-p|<1`.  The function `F(p):=(1-p)g(p)=1/2-S_2(p)` is therefore analytic
across `p=1`; its coefficients are first differences of the finite-alphabet
integer sequence `(c_N)`, hence again take finitely many values.  Szego's
theorem makes `F` rational.  Thus `g=F/(1-p)` is rational, and rationality
together with the finite alphabet makes `(c_N)` eventually periodic.  Let
`T` be its least eventual period.  Its poles lie among the `T`-th roots of
unity.

## Step 3 — the two-circle pole geometry

`g`'s poles lie on `|p|=1`.  Put `H(u)=S_2(1-u)`.  Its coefficients also form
a finite integer alphabet by the same bounded-window argument, while
`H(u)=1/2-u g(1-u)` makes it rational.  Hence its coefficients are eventually
periodic and its poles lie on `|u|=1`, i.e. the poles of `S_2` lie on
`|p-1|=1`. In the entire identity `q g(p) + S_2(p) = 1/2` the factor
`q = 1-p` kills only `g`'s pole at `p=1`; every OTHER pole `zeta != 1` of `g` has `q(zeta) != 0`, so it
must be cancelled by a pole of `S_2` at the same point. **But `|p|=1` and `|p-1|=1` intersect only at
`p = e^{+-i pi/3}`.** Hence every pole `zeta != 1` of `g` must equal `e^{+i pi/3}` or `e^{-i pi/3}`.

## Step 4 — `6` does not divide `T` is impossible

`e^{+-i pi/3}` is a `T`-th root of unity **iff `6 | T`**. So if `6` does NOT divide `T`, then `g` has no
pole except at `p = 1`, forcing `g = c/(1-p) + (polynomial)`; thus `c_N` is eventually the CONSTANT `c`,
and the residue `-c = -1/2` gives `c = 1/2`. But every `c_N` is an INTEGER, so `c = 1/2` is impossible.
**Contradiction.**

> **Intermediate conclusion.** If a `C=1` extractor existed, the least
> eventual period of its spine coefficients would be divisible by six.

## Step 5 — closing `6 | T` by integrality (COMPLETE)

For `6 | T`, `g`'s poles are a subset of `{1, omega, omega-bar}`, `omega = e^{i pi/3}`, so

```
   g(p) = (1/2)/(1-p) + c_omega/(1 - p/omega) + conj(c_omega)/(1 - p/omega-bar) + polynomial,
```

the `(1/2)` from the residue at `p=1`.  Reading off the coefficient of `p^N` (for large `N`):

```
   c_N = 1/2 + 2 Re( c_omega * omega^{-N} ) = 1/2 + 2 rho cos(psi - pi N/3),   c_omega = rho e^{i psi}.
```

But `c_N` are INTEGERS (the spine's Bernstein coefficients `w_{m,k}` are integers, hence the ordinary
coefficients `a_{m,j}` and their diagonal sums `c_N` are integers).  Choose
`N_0` beyond the polynomial transient and absorb `omega^(-N_0)` into the
phase.  Writing `C_j=c_(N_0+j)` and `A=2 rho cos psi=C_0-1/2`,
and using `2 rho cos(psi - pi/3) = A/2 + (sqrt3/2) B`, `2 rho cos(psi - 2pi/3) = -A/2 + (sqrt3/2) B`
(`B = 2 rho sin psi`), subtraction gives

```
   C_1-C_2=A=C_0-1/2.
```

The left side is an INTEGER; the right side is a HALF-integer.  **Contradiction.**  (And `c_omega=0`
gives `c_N=1/2` eventually, also non-integral.)  So `6|T` is impossible too.

## Theorem and consequence

> **The endpoint `C=1` is unattained: no exactly-fair extractor has deadline
> `T(n)<=n+D` for any `D`.  No strict inequality for the infimum `C*` follows.**

Combining the two cases: `6` does not divide `T` forces a single pole at `p=1`, so `c_N -> ` a constant
integer `= 1/2` (residue `-1/2`), impossible; `6 | T` forces the integer/half-integer contradiction above.
No bounded-additive-slack fair extractor exists.  The proof is a clean composite of Szego rigidity (finite
coefficient values `=>` rational), the **two-circle pole geometry** (`|p|=1` meets `|p-1|=1` only at
`e^{+-i pi/3}`), and integrality of the decided-tree spine coefficients.  The `e^{+-i pi/3}` is not folklore:
it is exactly where the bias circle meets the complementary-bias circle.  THM-3342 subsumes this endpoint
test and strengthens it to every sublinear excess.  The remaining general-class interval is `1<=C*<=2`,
with slope one unattained; THM-3009 gives a stronger lower bound only inside the balanced-block class.
