---
source: opus-2026-07-31-S4 (AMM 12592 minimal-C frontier; executing HYP-9061 sec 3's gamma=0 test)
status: >
  PARTIAL THEOREM (rigorous). Executes the Carlson-Szego route HYP-9061 sec 3 named but left undone.
  For the linear deadline C=1 (gamma=0) the spine generating function g is rational (Szego), and ALL
  eventual periods T with 6 not dividing T are IMPOSSIBLE by an integer-average-1/2 contradiction. The
  residual is exactly 6 | T, reduced (for symmetric extractors) to the functional equation
  (1-p) g(p) + p g(1-p) = 1/2 with a fixed residue-argument law at e^{+-i pi/3}. Hands the finite residual
  to death-star (HYP-9061 owner). This lower-bounds C* strictly above 1 modulo the 6|T residual.
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

## Step 5 — the residual `6 | T`, and the functional equation

For the symmetric extractor `V_m(p) = W_m(1-p)` (heads<->tails), (S) becomes the clean functional equation

```
   (1-p) g(p) + p g(1-p) = 1/2      for all p.
```

`g(1-p)` has poles at `{0, e^{+i pi/3}, e^{-i pi/3}}`, so the poles of `g` at `e^{+-i pi/3}` must cancel
between the two terms. The residue condition at `omega = e^{i pi/3}` is `(1-omega) Res_omega g =
omega \overline{Res_omega g}`, i.e. **`Res_omega g` is a positive-real multiple of `omega`** (`arg = pi/3`),
verified. This does not yet contradict; closing it needs the `[0,1]` value / integer-Bernstein data at the
`e^{+-i pi/3}` residues (a finite check per `T`). Handed to death-star (HYP-9061 owner).

## Consequence

`C* > 1` modulo the `6|T` residual: no bounded-additive-slack (`T(n) <= n + D`) exactly-fair extractor
exists unless its decided-tree spine is eventually periodic with `6 | T` and residues pinned at
`e^{+-i pi/3}`. This is the first rigorous execution of HYP-9061 sec 3's Carlson-Szego program and
strictly narrows the lower-bound direction of `(Q)`. The `e^{+-i pi/3}` is not arithmetic folklore: it is
exactly where the "bias circle" `|p|=1` meets the "complementary-bias circle" `|p-1|=1`.
