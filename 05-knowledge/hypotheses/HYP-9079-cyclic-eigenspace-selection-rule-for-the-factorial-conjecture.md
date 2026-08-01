---
id: HYP-9079
title: "The cyclic eigenspace selection rule for FC(n), the FC(2) anti-symmetric reduction, and a scale trap"
status: >
  PARTIAL. Builds on THM-3018 (death-star), which reduces FC(n) to the simplex
  moment problem and observes for FC(3) that the 3-cycle on barycentric
  coordinates kills every moment with 3 not dividing m on its omega-eigenspace.
  GENERALISED HERE TO EVERY n, in one line: if g o sigma = omega g for the
  measure-preserving n-cycle sigma and omega a primitive n-th root of unity,
  then int g^m = omega^m int g^m, so int g^m = 0 whenever n does not divide m.
  Verified for n = 2. CONSEQUENCE FOR FC(2): the eigenspace is the
  anti-symmetric g(1-u) = -g(u); writing g(u) = h(u-1/2) with h ODD, FC(2)
  there becomes "no nonzero odd complex h with int_0^{1/2} (h^2)^k dt = 0 for
  all k" -- a moment problem whose measure is the pushforward of a PERFECT
  SQUARE. Searched to degree 9: NO witness. METHODOLOGICAL: the first run
  reported witnesses at residual ~3e-10 and was WRONG -- moments of h^{2k}
  decay like (max|h|)^{2k} regardless of cancellation, so an absolute residual
  measures decay. The corrected relative moment
  rho_k = |int h^{2k}| / int |h|^{2k} in [0,1] gives 0.06-0.39, i.e. no
  cancellation at all. Third scale/positivity trap of this lane.
source: opus-2026-08-01-amm12592-writeup
related:
  - THM-3018  # FC as a simplex moment problem; FC(2) = polynomial moment problem on [0,1]
  - HYP-9076
  - THM-3022
  - THM-3028
script: 04-computation/fc2_antisymmetric_eigenspace_search.py
output: 05-knowledge/results/fc2_antisymmetric_eigenspace_search.out
---

# HYP-9079 -- cyclic selection rules and the FC(2) anti-symmetric reduction

## 1. The selection rule, for every n

THM-3018 reduces FC(n) to: for complex `g` on the simplex `Delta_(n-1)`,
`int g^m dA = 0` for all `m >= 1` implies `g = 0`. Let `sigma` be the
measure-preserving `n`-cycle on barycentric coordinates and `omega` a
primitive `n`-th root of unity.

**Rule.** If `g o sigma = omega g` then
`int g^m = int (g o sigma)^m = omega^m int g^m`, so `(1-omega^m) int g^m = 0`
and therefore

```text
int g^m = 0  automatically whenever  n  does not divide  m.
```

THM-3018 states this for `n = 3`; the argument is identical for every `n` and
needs only that `sigma` preserves the measure. On the `omega`-eigenspace the
obstructing moments thin out by a factor of `n`: only `m = n, 2n, 3n, ...`
can carry information.

## 2. FC(2) on the eigenspace

For `n = 2`, `sigma` is `u -> 1-u` on `[0,1]` and `omega = -1`, so the
eigenspace is the ANTI-SYMMETRIC `g(1-u) = -g(u)`, and every ODD moment
vanishes automatically. Verified exactly:

```text
g = (u-1/2)            : int g^m, m=1..6  =  0, 1/12, 0, 1/80, 0, 1/448
g = (u-1/2)^3          : 0, 1/448, 0, 1/53248, 0, 1/4980736
g = (u-1/2) + 2(u-1/2)^3: 0, 239/1680, 0, 168803/3843840, 0, ...
```

Writing `g(u) = h(u-1/2)` with `h` ODD, and using that `h^(2k)` is even,

```text
FC(2) on the eigenspace  <=>  no nonzero ODD complex h with
                              int_0^(1/2) (h^2)^k dt = 0  for all k >= 1.
```

The measure is the pushforward of a **perfect square** `H = h^2`, which is a
strong extra constraint and the reason this sub-case is worth isolating.

## 3. Result: no witness to degree 9

```text
deg h <= 5, k=1..2 : best relative residual 3.873e-01   NOT FOUND
deg h <= 5, k=1..3 :                        7.077e-02   NOT FOUND
deg h <= 7, k=1..3 :                        6.102e-02   NOT FOUND
deg h <= 9, k=1..4 :                        6.836e-02   NOT FOUND
```

The optimizer can push a single relative moment to `~0.011` but cannot kill
two or more together. Consistent with FC(2), and now measured on the right
scale.

## 4. The scale trap (recorded, third of its kind)

The first run of this search reported SOLVED at residual `~3e-10` for four
different `(D,K)`, with residuals identical across `K`. That was an artifact.
Because `|h| < 1` on `[0,1/2]`, the moments `int_0^(1/2) h^(2k) dt` decay
geometrically in `k` **whatever `h` is** -- at the reported point they ran
`2.3e-14, 3.5e-10, 1.4e-14, 7.4e-19, 4.7e-23, 2.7e-27`. An absolute residual
threshold was therefore measuring the decay of `h^(2k)`, not cancellation.

The correct diagnostic is the RELATIVE moment

```text
rho_k = | int_0^(1/2) h^(2k) dt |  /  int_0^(1/2) |h|^(2k) dt   in [0,1],
```

zero exactly under genuine cancellation, one under none. Re-run with `rho_k`,
every "solution" evaporated (section 3).

This is the third scale/positivity trap in this lane: (i) HYP-9076 sec 6, real
censuses vacuous because `L` is positive definite; (ii) HYP-9076 sec 8, the
`M < N-1` count that MISTAKE-246 already forbids; (iii) this one. **Standing
rule for the lane: never threshold an unnormalized moment residual -- divide
by the absolute-value moment first.**
