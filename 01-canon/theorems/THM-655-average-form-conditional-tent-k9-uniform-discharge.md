---
id: THM-655
title: The average-form conditional tent — sup_E avgc(E,P) <= c*(P) for all 715 shapes at k=9, DISCHARGING the k=9 (A') hlarge leg for EVERY cluster family and EVERY diameter (no diameter hypothesis), strictly stronger than the diam<=16 exhaustion (THM-653). The mechanism: the conditional Markov bound consumes the AVERAGE of the per-pair G_P-restricted tent masses (a sum over pairs), not their supremum; block-domination + a decreasing c-envelope bounds that average over all families
status: PROVED (elementary; three-step proof below). Machine-verified: sup_E avgc via the decreasing-envelope block bound clears 713/715 shapes; the 2 residual shapes ((6,8,9,10),(8,9,10,12)) clear by the single-hot-difference count bound (avgc <= 1.607/1.617 < c* = 1.768/1.770). c*(P) and c(d,P) computed in exact rational integration; float verification cross-checks the exact table.
source: mac-mini-2026-07-07-S57 (HYP-5237)
depends_on:
  - THM-651   # the shifted-tent floor; this is its conditional/average refinement
  - THM-530   # m_P = 14249/252252, the hlarge bar
related:
  - THM-653   # klein-S174: the tent-WINDOW composition (diam<=16 for k=9); this is diameter-free and subsumes it at k=9
  - MISTAKE-123  # the honest bars c*(P) is measured against
  - HYP-5207  # mac-mini-S56 c-table (the sup-form frontier this replaces with the average form)
external: none (pair equidistribution + Markov + elementary block domination).
---

# THM-655 — the average-form conditional tent (k=9 uniform discharge)

## Statement

Fix `k = 9`, `beta = beta_9 = 5/63`, `f(s) = (s - beta)_+` on `(0, 1/7]` (else 0). For a
shape `P` (a `(13-k)`-subset of `{1,...,13}`) let `G_P` be the good set (complement of the
teeth `[j/p - 1/(14p), j/p + 1/(14p)]`, `p in P`), `meas = meas(G_P)`, and

> `c(d, P) := ( int_{G_P} f(frac(d x)) dx ) / ( meas(G_P) * int f )`   (per-pair restricted mass ratio),
> `avgc(E, P) := (1/C(k,2)) * sum_{i<j} c(e_j - e_i, P)`  (average over the family's pair differences).

Let `c*(P) := (1 - m_P/meas(G_P)) * 7k/(2(k-1)(k-7))`   (`= (1 - m_P/meas) * 63/32` at `k=9`).

> **For every one of the `C(13,4) = 715` shapes `P`, `sup_E avgc(E, P) <= c*(P)`,**

the sup over all 9-element integer families `E`. Consequently the conditional union bound
gives `rho*(P, E) >= meas(G_P)(1 - avgc(1 - floor_9)) >= m_P` for **every** shape and
**every** family — the **k=9 (A') `hlarge` leg is discharged with no diameter hypothesis.**

## Why the AVERAGE, not the sup (the key correction to the S56/kps frontier)

On the safe event `S = {maxgap <= 1/7}`, `F(x) := sum_{i != j} f(frac((e_j-e_i)x)) >= 1 - k*beta`
pointwise (THM-651 Step 2, gap-sum budget). `F >= 0` everywhere, so restricting the first
moment to `G_P`:

> `(1 - k*beta) * meas(S ∩ G_P) <= int_{G_P} F = sum_{ordered pairs} int_{G_P} f(frac(dx))
>   = [ sum_{ordered} c(d,P) ] * meas(G_P) * int f = avgc(E,P) * k(k-1) * meas(G_P) * int f.`

Using the exact identity `k(k-1) int f / (1 - k*beta) = 1 - floor_9 = 2(k-1)(k-7)/(7k)`,

> `meas(S ∩ G_P) <= avgc(E,P) * (1 - floor_9) * meas(G_P)`,   hence
> `rho* >= meas(G_P) - meas(S ∩ G_P) >= meas(G_P)[1 - avgc(E,P)(1 - floor_9)]`.

The quantity that enters is the **sum over pairs** `= C(k,2) * avgc`, i.e. the AVERAGE of
the per-pair masses. kps-S73's named frontier and mac-mini-S56's c-table used the
**supremum** `sup_d c(d,P) <= 1.7` — a sufficient but strictly looser proxy that fails at
the isolated resonant differences `d in {1,2}` (S56). Replacing sup by average removes the
resonance obstruction: a single hot `d` is diluted by the other `C(k,2)-1` cold pairs.

## Proof that `sup_E avgc <= c*` (two elementary bounds)

**Block domination (the only combinatorics).** For any `k` integers `a_1<...<a_k`,
`a_j - a_i >= j - i`, so the sorted difference multiset `d_(1) <= ... <= d_(C(k,2))`
term-wise dominates the block `{0,...,k-1}`'s: `d_(r) >= block_(r)` for all `r`
(equivalently `N_E(t) := #{pairs with diff <= t} <= N_block(t) = sum_{d<=t}(k-d)`).

**Bound A — the decreasing-envelope block bound (713 shapes).** Set
`chat(d) := max_{d' >= d} c(d', P)` (decreasing, `c(d) <= chat(d)`). Then
`avgc(E) = (1/C(k,2)) sum_r c(d_(r)) <= (1/C(k,2)) sum_r chat(d_(r)) <=
(1/C(k,2)) sum_r chat(block_(r)) = (1/C(k,2)) sum_{d=1}^{k-1} (k-d) chat(d) =: EnvBlock(P)`,
the middle step by `chat` decreasing and `d_(r) >= block_(r)`. Machine check:
`EnvBlock(P) <= c*(P)` at **713 of 715** shapes (all families, all diameters).

**Bound B — single-hot-difference count (the 2 residual shapes).** For
`P in {(6,8,9,10), (8,9,10,12)}` the ONLY difference with `c(d,P) > c*(P)` is `d = 6`
(`c(6) = 1.8039 / 1.7910`; every other `d`, including the tail `d > 5000` via the Koksma
bound `c <= 1 + K_P Δ_2/(d·meas·int f) < 1.02`, has `c(d) <= c*`). Let
`c2 := max_{d != 6} c(d,P)` (`= 1.5502 / 1.5672` at `d=8`). A `k`-set has at most `k-1 = 8`
pairs at any fixed difference (the shift-by-`d` graph is a disjoint union of paths), so
`n(6) <= 8`, and since `c(6) > c2`,

> `avgc(E) <= (1/C(k,2)) [ n(6) c(6) + (C(k,2)-n(6)) c2 ] <= (8 c(6) + 28 c2)/36
>   = 1.6066 / 1.6169  <  c*(P) = 1.7681 / 1.7696.`

Both bounds are rigorous and hold for every family of every diameter. ∎

## The Koksma tail lemma (closes `d` beyond the exact table)

With `D(t) = int_0^t f - t·int f` (continuous, `D(0)=D(1)=0`), the interval substitution gives
`int_a^b f(frac(u)) du = (b-a) int f + D(frac b) - D(frac a)`, so on `G_P = ∪[l,h]`,
`c(d,P) = 1 + (1/(d·meas·int f)) sum_iv [D(frac(dh)) - D(frac(dl))]`. Since
`|D| <= max D = (1/7 - beta) w^2/2 ... ` bounded by `Δ_2 = int f·(6/7 + beta + w^2/4)`,

> `|c(d,P) - 1| <= K_P · Δ_2/(int f) / (d · meas(G_P))`,  `K_P = #intervals of G_P`,

so `c(d) -> 1` at rate `~ K_P/d`; the exact table to `d = 5000` plus this tail covers all `d`.

## Scope and what it feeds

- **k=9 (A') leg: CLOSED, diameter-free** — supersedes THM-653's diam<=16 exhaustion for k=9.
- **k=10:** the same average-form bound is necessary but NOT sufficient alone: `c*(P) ~ 1.18`
  is much tighter (the honest bar `0.4521` is far below `meas(G_P) ~ 0.64`), and `EnvBlock`
  exceeds it at most shapes. k=10 needs the average-form conditional tent COMPOSED with
  klein-S174's window mass `W_F` (THM-653) and/or the spread floor (klein-S175 HYP-4991):
  the compact offenders (block, 2-AP) are diam<=18 so THM-653 covers them; the average form
  covers the wide families. The k=10 split is the live handoff.
- **k=11,12,13:** the gap-histogram/tent frame is vacuous (`floor_k <= 0`); those legs stay on
  the intersection ledger + bounded diameter (k=13 diam<=75 is UNCONDITIONAL, opus-S145 Lean
  AP76 certificate) + the spread route.

## Verification & files

`04-computation/lrc14_avgform_conditional_tent_macmini_S57.py` (+ `.out`): the 715-shape
average-form sweep, exact `c(d,P)` via the antiderivative, Koksma tail, exact-Fraction
confirmation at the worst shape. `04-computation/lrc14_block_domination_k9_macmini_S57.py`
(+ `.out`): block domination, LP relaxation, and the EnvBlock rigorous bound (713/715).
`04-computation/lrc14_k9_residual_closer_macmini_S57.py` (+ `.out`): the 2 residual shapes —
compact exhaustion (max avgc 1.40 = the 2-AP), single-hot-difference `c2` closer (1.607/1.617).
