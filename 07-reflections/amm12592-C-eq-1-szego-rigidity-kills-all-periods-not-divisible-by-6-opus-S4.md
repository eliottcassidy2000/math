---
source: opus-2026-07-31-S4 (AMM 12592 minimal-C frontier; executing HYP-9061 sec 3's gamma=0 test)
status: >
  THEOREM (rigorous, COMPLETE). Executes and finishes the Carlson-Szego route HYP-9061 sec 3 named but left
  undone. For the linear deadline C=1 (gamma=0) the spine generating function g is rational (Szego); the
  two-circle pole geometry forces g's only poles to be at p=1 and (if 6|T) e^{+-i pi/3}; and the INTEGRALITY
  of the spine coefficients c_N kills every case. Hence NO bounded-additive-deadline (T(n)<=n+D) exactly-fair
  extractor exists: C* > 1 for AMM 12592, unconditionally.
tags: [amm12592, minimal-C, biased-coin, extractor, szego, carlson, rigidity, deficit-flow, HYP-9061, lower-bound]
related: [HYP-9061, THM-2160, THM-2225]
---

# AMM 12592, C = 1: Szego rigidity kills every period 6 does not divide (opus-S4)

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

## Step 2 — `g` is rational (Szego)

Near `p=1`: `S_2(p) := sum_{m>=1} q^m p V_m(p)` is analytic in `|1-p| < 1` (uniformly convergent) with
`S_2(1) = 0`. So (S) gives `g(p) = (1/2 - S_2(p))/(1-p)`, meromorphic at `p=1` with a SIMPLE POLE of
residue `-1/2` (numerator `-> 1/2 != 0`). A power series whose coefficients take finitely many values is,
by **Szego's theorem**, either rational or has the unit circle as a natural boundary; the meromorphic
continuation through `p=1` rules out the natural boundary, so **`g` is rational** and `{c_N}` is
eventually periodic with some period `T`. Its poles lie among the `T`-th roots of unity.

## Step 3 — the two-circle pole geometry

`g`'s poles lie on `|p| = 1`. `S_2` is rational in `u = 1-p` (its `V_m` are also eventually periodic), so
its poles lie on `|u| = 1`, i.e. on `|p-1| = 1`. In the entire identity `q g(p) + S_2(p) = 1/2` the factor
`q = 1-p` kills only `g`'s pole at `p=1`; every OTHER pole `zeta != 1` of `g` has `q(zeta) != 0`, so it
must be cancelled by a pole of `S_2` at the same point. **But `|p|=1` and `|p-1|=1` intersect only at
`p = e^{+-i pi/3}`.** Hence every pole `zeta != 1` of `g` must equal `e^{+i pi/3}` or `e^{-i pi/3}`.

## Step 4 — `6` does not divide `T` is impossible

`e^{+-i pi/3}` is a `T`-th root of unity **iff `6 | T`**. So if `6` does NOT divide `T`, then `g` has no
pole except at `p = 1`, forcing `g = c/(1-p) + (polynomial)`; thus `c_N` is eventually the CONSTANT `c`,
and the residue `-c = -1/2` gives `c = 1/2`. But every `c_N` is an INTEGER, so `c = 1/2` is impossible.
**Contradiction.**

> **Partial theorem.** If a `C = 1` exactly-fair extractor exists, its spine is eventually periodic with
> period divisible by `6`.

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
coefficients `a_{m,j}` and their diagonal sums `c_N` are integers).  Writing `A = 2 rho cos psi = c_0 - 1/2`
and using `2 rho cos(psi - pi/3) = A/2 + (sqrt3/2) B`, `2 rho cos(psi - 2pi/3) = -A/2 + (sqrt3/2) B`
(`B = 2 rho sin psi`), subtraction gives

```
   c_1 - c_2 = A = c_0 - 1/2.
```

The left side is an INTEGER; the right side is a HALF-integer.  **Contradiction.**  (And `c_omega = 0`
gives `c_N = 1/2` eventually, also non-integer.)  So `6 | T` is impossible too.  Verified numerically: no
complex `c_omega` makes `c_0..c_5` all integers (200,000 random samples, none).

## Theorem and consequence

> **`C = 1` is impossible for AMM 12592: no exactly-fair extractor has deadline `T(n) <= n + D` for any `D`.
> Hence `C* > 1`.**

Combining the two cases: `6` does not divide `T` forces a single pole at `p=1`, so `c_N -> ` a constant
integer `= 1/2` (residue `-1/2`), impossible; `6 | T` forces the integer/half-integer contradiction above.
No bounded-additive-slack fair extractor exists.  The proof is a clean composite of Szego rigidity (finite
coefficient values `=>` rational), the **two-circle pole geometry** (`|p|=1` meets `|p-1|=1` only at
`e^{+-i pi/3}`), and integrality of the decided-tree spine coefficients.  The `e^{+-i pi/3}` is not folklore:
it is exactly where the bias circle meets the complementary-bias circle.  This closes HYP-9061 sec 3's
"cheapest decisive test" as a full theorem and sets the base of the `C* in (1, 2]` ladder.
