---
id: THM-3028
title: "AMM 12592: a certified sup-norm bound on admissible blocks, and a propagated tail bound on every epoch residual"
status: >
  PROVED + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE AUDIT. Two necessary
  conditions for the gamma-epoch recursion, both proved from capacity+parity
  alone, plus a numerical warning that invalidates float diagnostics in this
  lane.
  (SUP) EVERY admissible block satisfies |Delta(t)| <= 1 for all t in [0,1].
  Proof: parity gives delta_k = binom(d,k) - 2 m_k with 0 <= m_k <= binom(d,k);
  x + (1-x) = 1 gives sum_k binom(d,k) B_{d,k}(x) = 1 identically, where
  B_{d,k}(x) = x^{d-k}(1-x)^k >= 0 on [0,1]; so Delta = 1 - 2 sum_k m_k B_{d,k}
  lies in [1-2, 1] = [-1,1]. TIGHT: sup-norm exactly 1 is attained.
  (PROP) Unrolling p sigma_i = sigma_{i-1} - Delta_i gives the tail identity
  sigma_{j-1} = sum_{i>=j} x^{i-j} Delta_i, so with (SUP)
      |sigma_{j-1}(t)| <= (1 - t^{R-j})/(1-t)     for all t in [0,1),
  a necessary condition on EVERY intermediate residual, not just the last.
  VERIFIED TIGHT: the exact violation ratio along the true R = 16, 32, 64
  solution paths is exactly 1.000000, and used as a search filter it still
  lets R = 16, 32, 64 solve -- so it is sound and saturated.
  (FLOAT) WARNING: epoch residuals have alternating coefficients reaching
  ~1e54, so floating-point evaluation is meaningless. Float reported a
  violation ratio of 3.198872 on the true R = 64 path where the exact value is
  1.000000. All (SUP)/(PROP) tests must use exact rational/integer arithmetic.
  DIAGNOSTIC CONSEQUENCE, refining the THM-3002 5b correction: under the
  certified filter, R = 128 fails at a row that MOVES STRONGLY with the
  constant slack D0 -- row 42 (D0=0), 51 (D0=1), 97 (D0=2), 90 (D0=4), 80
  (D0=8) out of 128 -- with only MARGINAL misses (best achievable ratio 1.04
  to 1.25, not orders of magnitude) and a starved beam (2-16 states passing)
  immediately before death. A structural wall would not move 55 rows under a
  constant that leaves C = 8/5 unchanged. So the R = 128 failure is evidence
  of SEARCH LIMITATION; the question remains formally open, but the balance of
  evidence now leans away from an obstruction.
source: death-star-2026-08-01-coinC2
depends_on:
  - THM-3002
related:
  - THM-3026
  - THM-3024
  - THM-3007
  - HYP-9061
script: 04-computation/amm12592_certified_tail_bounds_thm3028.py
output: 05-knowledge/results/amm12592_certified_tail_bounds_thm3028.out
---

# THM-3028 -- certified bounds on epoch residuals

## 1. Setup

The `gamma`-epoch recursion is `sigma_{-1} = q^{R-1}` (`q = 1-x`) and
`x sigma_i = sigma_{i-1} - Delta_i`, where each block is written in the basis

```text
B_{d,k}(x) = x^{d-k}(1-x)^k,      Delta_i = sum_k delta_{i,k} B_{d_i,k},
|delta_{i,k}| <= binom(d_i,k)   (capacity),
 delta_{i,k} = binom(d_i,k) mod 2  (parity).
```

Call such a block **admissible at `d`** (cf. THM-3026).

## 2. (SUP) Admissible blocks are bounded by 1 on the unit interval

Parity says `delta_k` and `binom(d,k)` have the same parity, so

```text
delta_k = binom(d,k) - 2 m_k,     and capacity gives   0 <= m_k <= binom(d,k).
```

Since `x + (1-x) = 1`,

```text
sum_k binom(d,k) B_{d,k}(x) = (x + (1-x))^d = 1     identically,
```

and `B_{d,k}(t) >= 0` for `t` in `[0,1]`. Therefore

```text
Delta(t) = 1 - 2 sum_k m_k B_{d,k}(t),      0 <= sum_k m_k B_{d,k}(t) <= 1,
```

so

```text
|Delta(t)| <= 1      for every t in [0,1].                              (SUP)
```

The bound is **attained** (take all `m_k = binom(d,k)`, giving `Delta = -1`).
Verified on 2000 random admissible blocks with `d <= 12`: maximum sup-norm
exactly `1.000000`.

## 3. (PROP) The bound propagates to every residual

Unrolling `sigma_{i-1} = Delta_i + x sigma_i` gives the tail identity

```text
sigma_{j-1}(x) = sum_{i=j}^{R-1} x^{i-j} Delta_i(x),
```

so by (SUP), for every `t` in `[0,1)`,

```text
|sigma_{j-1}(t)| <= sum_{i=j}^{R-1} t^{i-j} = (1 - t^{R-j})/(1-t).      (PROP)
```

This is a **necessary condition on every intermediate residual**, derived from
capacity and parity alone -- not a heuristic. It is also **saturated**: the
exact maximum violation ratio along the true `R = 16, 32, 64` solution paths is
exactly `1.000000`. Used as a search filter it is sound: `R = 16, 32, 64` all
still solve through it.

**Ordering caveat.** Because true paths *saturate* (PROP), sorting a beam by
ascending violation ratio is exactly backwards -- it discards the true path
first. (PROP) must be used as a **filter**, with the ordering left to a
separate objective.

## 4. (FLOAT) Floating point is worthless in this lane

Epoch residuals have alternating coefficients reaching `~1e54` at `R = 128`.
Evaluating them in double precision is catastrophic cancellation:

```text
true R = 64 solution path:   exact max violation ratio = 1.000000
                             float max violation ratio = 3.198872
```

The float figure suggested (PROP) was false. It is not; the arithmetic was.
**All (SUP)/(PROP) testing must be exact** -- evaluate `b^D sigma(a/b) =
sum_k c_k a^k b^{D-k}` in integers and cross-multiply against the bound.

## 5. What this says about R = 128

Applying (PROP) as a certified filter (ordering by residual length, which is
what closes `R <= 64`), the `R = 128` search dies at a row that **moves
strongly with the constant slack `D0`**:

```text
D0      death row (of 128)     best achievable ratio at death   states passing just before
 0            42                        1.252                            9
 1            51                        1.169                            7
 2            97                        1.183                            9
 4            90                        1.164                            2
 8            80                        1.040                           10
```

`D0` is a *constant* additive slack: `T(n) = n + 1 + floor(3n/5) + D0` gives
`C = 1 + 3/5 = 8/5` for every fixed `D0`, while enlarging each capacity to
`binom(d_i + D0, k)`. **A structural wall would not move 55 rows under a
constant that leaves the asymptotic constant unchanged**, and the misses are
marginal (`1.04` to `1.25`) rather than gross, with the beam starved to `2-16`
states immediately before death.

This **refines** the scope correction recorded in THM-3002 section 5b. That
correction withdrew the unsupported claim that `R = 128` was "a search budget
limit, not an obstruction"; the withdrawal stands, because nothing there
established it. What the certified filter adds is genuine evidence, and it
points back toward search limitation: the failure is `D0`-mobile and marginal.
The honest status is unchanged in form -- **`R = 128` is formally open, and
`C <= 8/5` is verified for `n <= 127` only** -- but the balance of evidence now
leans away from an obstruction.

## 6. Scope

(SUP) and (PROP) are proofs, valid for any `gamma` and any `D0`, and depend
only on capacity and parity. Section 5 is a finite search diagnostic and proves
nothing about the existence of an `R = 128` solution in either direction. The
tightness observation (true paths saturate (PROP)) is a verified fact about
three solution paths, not a theorem about all of them.

Referee: `amm12592_certified_tail_bounds_thm3028.py` -- 2000 random blocks for
(SUP); exact-arithmetic verification that the true `R = 16,32,64` paths
saturate (PROP) at ratio exactly `1`; the float-vs-exact discrepancy on
`R = 64`; soundness of the filter (`R = 16,32,64` still solve); and the `D0`
death-row table.
