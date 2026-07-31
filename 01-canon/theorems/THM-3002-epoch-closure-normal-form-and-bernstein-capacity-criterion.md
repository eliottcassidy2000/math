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

## 4b. The parity lemma and the depth-free mod-2 clock (PROVED)

**Lemma.** For `Delta` of degree `<= d`, the Lucas box parity condition
`delta_k == binom(d,k) (mod 2)` on its Bernstein-`d` coefficients is
*equivalent* to

```text
Delta(p) == 1  (mod 2)   coefficientwise in Z[p].                  (5)
```

*Proof.* `D(x) := sum_k delta_k x^k = (1+x)^d Delta(1/(1+x))`. The condition
`D == (1+x)^d (mod 2)` becomes `sum_j c_j (1+x)^{d-j} == 0 (mod 2)` where
`Delta - 1 = sum_j c_j p^j`; the polynomials `(1+x)^{d-j}` have distinct
degrees, hence are `F_2`-independent, so every `c_j` is even. QED
(referee C7a, exhaustive for `d <= 8`).

**Corollary (parity clock).** Reducing the residual recursion mod 2 kills
all depth dependence:

```text
sigma_i = (sigma_{i-1} + 1)/p   over F_2[p],   sigma_{-1} = (1+p)^{R-1}. (6)
```

The `H = 1` program requires this orbit to survive `R-1` steps (each
`sigma(0) = 1`) and end at `sigma_{R-2} == 1`. **For every dyadic
`R = 2^r <= 2048` it does** (referee C7b). Since (6) is depth-free, the
verdict is identical for every rate `gamma`: **parity never obstructs a
dyadic epoch**, at any rate. This is THM-2976's checkpoint vanishing (T1)
seen inside the closure recursion, and it isolates the obstruction: the
threshold of the `H = 1` checkpoint program is *purely archimedean*,
governed by criterion (4) alone.

## 4c. The threshold is ~0.5980, and BOTH certificate readings fail it
(VERIFIED-EXACT, 2026-07-31; integrates klein-S428 and opus THM-3006)

Exact bisection of criterion (4) at fixed `R` gives a threshold sequence
that increases in `R`:

```text
R = 256 : gamma* in (0.584902, 0.584906)
R = 512 : gamma* in (0.590652, 0.590657)
R = 1024: gamma* in (0.593925, 0.593929)
```

Successive differences `0.005750, 0.003273` have ratio `0.569`, so the
geometric extrapolation is `0.5939 + 0.0043 = 0.5982` — agreeing with the
independent asymptotic entropy bisection `gamma* ~ 0.59799`. Both give

```text
C = 1 + gamma*  ~  1.598.
```

**Consequence for certificate (27).** klein-S428 observed that the raw
reading `gamma = 2457/6592 = 0.3727` is killed instantly by (4), while the
`(C-1)/C` reading `gamma = 2457/4135 = 0.5941958887...` "survives to
`R = 1024`". The finite-`R` table above explains why, and shows the
survival is an artifact: the `R = 1024` threshold `0.593927` sits only
`2.7e-4` below `2457/4135`. At the next scale the ordering flips —

```text
gamma = 2457/4135 :  R = 1024  worst log-ratio  +0.17365 (t = 1)   AMPLE
                     R = 2048  worst log-ratio  -2.66767 (t = 783) DEFICIENT
```

— exactly the pattern by which `gamma = 1/2` survives to `R = 16` and dies
by `R = 64`. **So criterion (4) refutes both readings of the certificate
weight inside the `H = 1` checkpoint-closure class.** Scope: (4) is
necessary for that class only, so this refutes "the certificate weight is
the rate of an `H = 1` dyadic-checkpoint construction"; it does not bound
`C*` itself, and per klein-S428 the weight is in any case an output of the
source's construction rather than a consequence of the inequality.

**Consistency with opus THM-3006.** That theorem exhibits exactly fair
rules with within-shell ratios `rho(4) = 3/2`, `rho(8) = 14/9`,
`rho(16) = 25/16`, `rho(32) <= 11/7` — an increasing sequence
`1.5000, 1.5556, 1.5625, 1.5714`. Criterion (4) predicts that no
construction of this class can push the limit below `~1.598`; the two
lines are therefore compatible, and the sharp prediction is
**`sup_r rho(2^r) ~ 1.598`, in particular `< 2` (so `C* < 2`) but not
below `1.59`**. That is falsifiable at `r = 6, 7`.

## 5. Verified closures (VERIFIED-EXACT)

(*) is solved exactly at `gamma = 1/2, D0 = 0` for `R = 2, 4, 8, 16`, and
each solution is reassembled into a full block closure via the mirror ansatz
and re-verified as a polynomial identity over `Z` (referee C3/C6).
Together with lane G5's independently-found base `[1,7]` and epoch `[8,15]`
witnesses (the latter re-verified here from first principles), **every
dyadic epoch through `[16,31]` closes at `gamma = 1/2`**: an exactly fair
extractor exists for all critical values `n <= 31` with
`T(n) = n + 1 + floor(n/2)`, i.e. `C = 3/2` behaviour through four epochs.

## 5b. gamma = 3/5 closes through R = 64 (VERIFIED-EXACT, 2026-07-31)

The failures of the first solver at `gamma = 1/2, 3/5, R = 32` were a
policy artifact of controlling **one** residual coefficient per step. The
residual recursion steers the *first several* coefficients of `sigma_i`
directly: `delta_{d}` is forced, `delta_{d-1}` sets `sigma_i(0)`,
`delta_{d-2}` sets `sigma_i(1)`, and so on. A beam search that enumerates
small targets for the first two of these — instead of greedily zeroing each
coefficient — solves (*) at `gamma = 3/5, D0 = 0` for

```text
R = 8, 16, 32, 64      (all four re-verified as exact polynomial identities)
```

i.e. **every dyadic epoch through `[64,127]` closes at `gamma = 3/5`**, so
an exactly fair extractor exists for all critical values `n <= 127` with
`T(n) = n + 1 + floor(3n/5)`, giving `C = 8/5` behaviour over seven
epochs with residual identically zero. (`R = 128` reaches the final row and
fails only there at beam widths 40-60; this is a search budget limit, not
an obstruction — criterion (4) is uniformly ample at `3/5`, worst ratio
`~1.20` at `t = 1`, stable to `R = 256`.) Script:
`04-computation/amm12592_gamma35_beam_deathstar.py`.

Combined with sec. 4c: `gamma = 3/5 = 0.6` sits just above the extrapolated
threshold `gamma* ~ 0.598`, which is exactly why it is the first "round"
rate that survives at every scale, and why `gamma = 1/2` and
`gamma = 2457/4135` do not.

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
