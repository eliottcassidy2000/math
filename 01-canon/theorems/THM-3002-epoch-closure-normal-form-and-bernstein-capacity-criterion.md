---
id: THM-3002
title: "Epoch-closure normal form (pq)^R and the Bernstein capacity criterion"
status: >
  PROVED + VERIFIED-EXACT. In the THM-2966 spine normal form for AMM 12592,
  block closure of B = [m_lo, m_hi] forces p^{m_lo-1} F = -q^{m_lo-1} G, hence
  F = q^{m_lo-1} H and G = -p^{m_lo-1} H for one common H; the audited
  gamma=1/2 epoch [8,15] witness has H = 1, which is exactly the statement
  that the 0-side deficit equals (pq)^{m_lo} and the 1-side -(pq)^{m_lo} --
  THM-2160's middle pair at block scale. Under H = 1 the two sides decouple
  by the mirror ansatz delta^1 = -delta^0, leaving the single representation
  problem (*) q^{R-1} = sum_{i<R} p^i Delta_i with Delta_i in the Lucas box of
  degree d_{R+i}. The exact capacity identity max |[p^t] Delta| =
  binom(d,t) 2^t then gives the necessary criterion
  sum_{i<=t} binom(d_i,t-i) 2^{t-i} >= binom(R-1,t) for all t, whose exact
  ledger is a sharp trichotomy in gamma: exponentially DEFICIENT for
  gamma < 1/2, marginal-then-deficient at gamma = 1/2 (already at R = 64,
  t = 25), and uniformly ample (binding at t = 1, ratio ~1.2, stable in R)
  for gamma >= 3/5; the asymptotic threshold is gamma* ~ 0.598 as the
  solution of a two-ray entropy comparison. (*) is solved EXACTLY at
  gamma = 1/2 for R = 2, 4, 8, 16, so every dyadic epoch through rows
  [16,31] closes and an exactly fair extractor exists for all critical
  values n <= 31 with T(n) = n + 1 + floor(n/2).
source: death-star-2026-07-31-coinC2
depends_on:
  - THM-2966
related:
  - THM-2160
  - THM-2976
  - THM-2977
external:
  - "Elliot Glazer, American Mathematical Monthly Problem 12592 (2026)."
script: 04-computation/amm12592_epoch_closure_capacity_thm3002.py
output: 05-knowledge/results/amm12592_epoch_closure_capacity_thm3002.out
---

# THM-3002 -- epoch-closure normal form and the Bernstein capacity criterion

Frame: THM-2966. A **block closure system** `(E_B)` for `B = [m_lo, m_hi]`
asks for doubled deficits `delta` (integer, `|delta| <= binom(d_m,k)`,
`delta == binom(d_m,k) mod 2`) with
`sum_{cells of B} delta p^z q^o = 0` identically. Lane G5's Theorem G5-1
(sufficiency, six lines from THM-2966) assembles closures of the base
`[1, 2^{r0}-1]` and of every dyadic epoch `[2^r, 2^{r+1}-1]` into an exactly
fair extractor with pathwise deadline `T(m) = m + 1 + d_m`; with
`d_m = floor(gamma m) + D0` this yields `C* <= 1 + gamma`.

## 1. The reduction (PROVED)

Write the two sides as `p^{m_lo} q F` and `q^{m_lo} p G`, where
`F = sum_{m in B} p^{m-m_lo} Delta^0_m` and
`G = sum_{m in B} q^{m-m_lo} Delta^1_m`. Then

```text
(E_B)  <=>  p^{m_lo-1} F = - q^{m_lo-1} G.                        (1)
```

Since `gcd(p, q) = 1` in `Z[p]`, (1) forces `q^{m_lo-1} | F`; writing
`F = q^{m_lo-1} H` gives `G = -p^{m_lo-1} H` with the **same** `H`. A
necessary degree condition drops out: `deg F <= (m_hi - m_lo) + d_{m_hi}`
must be at least `m_lo - 1` — satisfied with room by a dyadic epoch at
`gamma >= 1/2`, and violated by short blocks, which is why row pairs cannot
close on their own beyond the pairing horizon.

## 2. The `(pq)^R` normal form (PROVED)

`H = 1` is equivalent to

```text
sum_{m in B} p^m q Delta^0_m = (pq)^{m_lo},
sum_{m in B} q^m p Delta^1_m = -(pq)^{m_lo}.                      (2)
```

So a block closes in normal form exactly when its 0-side carries one unit of
deficit at the monomial of the word `0^R 1^R` and its 1-side the opposite
unit at `1^R 0^R`: **THM-2160's middle-pair cancellation, promoted from a
single row to a whole dyadic epoch**. Referee check C1/C2 confirms `H = 1`
on the audited `[8,15]`, `gamma = 1/2` witness.

## 3. Decoupling and the representation problem (*) (PROVED)

Under `H = 1` the mirror ansatz `delta^1_{m,k} = -delta^0_{m,k}` solves the
1-side whenever the 0-side is solved (`Delta^1_m(p,q) = -Delta^0_m(q,p)`
maps (2a) to (2b)). Hence closure of the epoch `[R, 2R-1]` reduces to

```text
(*)   q^{R-1} = sum_{i=0}^{R-1} p^i Delta_i(p,q),
      Delta_i integer homogeneous of degree d_i = d_{R+i},
      |coeff_k| <= binom(d_i,k),  coeff_k == binom(d_i,k) (mod 2).
```

A useful equivalent is the residual recursion `p sigma_i = sigma_{i-1} -
Delta_i` with `sigma_{-1} = q^{R-1}`: since the `k = d_i` cell has box
`binom(d_i,d_i) = 1`, the value `Delta_i(0) = sigma_{i-1}(0)` is forced and
must be `+-1` at every step — a hard unit invariant steered solely by
`delta_{i,d_i-1}`.

## 4. The capacity identity and criterion (PROVED)

For `Delta` in the Lucas box of degree `d`,

```text
max |[p^t] Delta| = sum_{j=0}^{t} binom(d,j) binom(d-j,t-j)
                  = binom(d,t) 2^t,                                (3)
```

by `binom(d,j) binom(d-j,t-j) = binom(d,t) binom(t,j)` and Vandermonde; the
maximum is attained inside the box with correct parity (referee C4, all
`d <= 39`). Since row `i` contributes `[p^{t-i}]Delta_i` to the `p^t`
coefficient of (*), a necessary condition for (*) is

```text
sum_{i=0}^{t} binom(d_i, t-i) 2^{t-i}  >=  binom(R-1, t)   for all t. (4)
```

**Trichotomy (exact ledger, referee C5).** With `d_i = floor(gamma(R+i))`
the ratio in (4) behaves like `(2 gamma)^t` to leading order:

| gamma | R = 16 | R = 64 | R = 256 | verdict |
|-------|--------|--------|---------|---------|
| 1/3 | 0.1013 (t=6) | 1.0e-6 (t=27) | ~0 (t=108) | deficient |
| 2/5 | 0.2448 (t=6) | 8.9e-5 (t=26) | ~0 (t=105) | deficient |
| 1/2 | 1.1333 (t=1) | **0.0437 (t=25)** | ~0 (t=101) | marginal, then deficient |
| 3/5 | 1.2667 (t=1) | 1.2222 (t=1) | 1.2039 (t=1) | ample |
| 2/3 | 1.4000 (t=1) | 1.3492 (t=1) | 1.3373 (t=1) | ample |
| 3/4 | 1.6667 (t=1) | 1.5397 (t=1) | 1.5098 (t=1) | ample |

Above the threshold the binding constraint is the trivial `t = 1` with a
stable ~20% margin; below it, deficiency is exponential in `t`. Writing
`t = xR`, `i = yR` and using `binom(aR,bR) ~ exp(R a H(b/a))`, (4) becomes
the **two-ray entropy comparison**

```text
max_{0<=y<=x} [ gamma(1+y) H( (x-y)/(gamma(1+y)) ) + (x-y) log 2 ]  >=  H(x),
```

whose threshold is `gamma* ~ 0.5980` with binding ray `x* ~ 0.38`
(numerically bisected; the criterion is exactly the "two rapidities with a
rational weight and a margin" shape of the decoded certificate (27), cf.
HYP-9061 sec. 2d and THM-2977's verdict that (27) must gate a rate dual).

## 5. Verified closures (VERIFIED-EXACT)

(*) is solved exactly at `gamma = 1/2, D0 = 0` for `R = 2, 4, 8, 16`, and
each solution is reassembled into a full block closure via the mirror ansatz
and re-verified as a polynomial identity over `Z` (referee C3/C6).
Together with lane G5's independently-found base `[1,7]` and epoch `[8,15]`
witnesses (the latter re-verified here from first principles), **every
dyadic epoch through `[16,31]` closes at `gamma = 1/2`**: an exactly fair
extractor exists for all critical values `n <= 31` with
`T(n) = n + 1 + floor(n/2)`, i.e. `C = 3/2` behaviour through four epochs.

## 6. Scope — what this does and does not settle

1. The `gamma = 1/2` successes are **finite-size**: by (4) with the table
   above, the `H = 1` program at `gamma = 1/2` is infeasible for large `R`
   (dead already at `R = 64`). So `C* <= 3/2` does **not** follow from the
   verified epochs, and any claim of a `3/2` construction must either leave
   the `H = 1` normal form or change rate.
2. The resemblance between (4)'s entropy form and certificate (27) is
   **structural only**. Per klein-S428 the inequality (27) admits an open
   half-line of weights `alpha > 0.3674729...`, so `2457/6592` is an output
   of the source's construction and cannot be recovered by inverting the
   inequality; `1/25` is a chosen safety margin. No numeric identification
   of `gamma*` with `alpha` is claimed here, and the earlier
   "capacity straddle" reading is retracted as evidence.
3. Criterion (4) is necessary for (*), i.e. for the `H = 1` normal form of
   the checkpoint-closure class. It is **not** proved necessary for closure
   with general `H`, and not for extractors outside the checkpoint class.
   It is likewise not sufficient: ampleness at `gamma >= 3/5` predicts, but
   does not prove, closure at those rates.
3. The live consequence is a sharp, falsifiable target: if the epochs close
   for all `R` at some `gamma < 1` (the ledger says look at
   `gamma ~ 3/5`), then `C* <= 1 + gamma < 2` and the long-standing
   expectation `C* = 2` is false. Lane D's "band freeze" evidence for
   `C* = 2` is now known to be a policy artifact (lane C2/F2), and this
   theorem supplies the mechanism it was missing.

Referee: C1-C6 exact, `ALL THM-3002 REFEREE CHECKS PASSED`. QED.
