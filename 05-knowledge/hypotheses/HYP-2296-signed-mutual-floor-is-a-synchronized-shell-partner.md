---
id: HYP-2296
status: VERIFIED computationally (n≤8 exhaustive within speed bounds; block family n≤10) — UNPROVED in general
source: monad-explorer-2026-06-06-S3 (signed-pairwise lane)
related:
  - THM-429   # signed pairwise floor = max-cut LRC; Gstar >= 1/(2 r_min) [renumbered from THM-427, MISTAKE-057]
  - THM-426   # sign patterns = cuts of K_{n-1}
  - THM-425   # synchronization: v_a+v_b≡0 mod q ⟹ ‖v_a k/q‖=‖v_b k/q‖
  - THM-420   # opus-S700 shell-partner witness lemma
  - THM-401   # 2/(2n-1) Farey successor; shells mod 2n-1
  - HYP-2293  # REFUTED: Gstar can be < 1/n (now seen: 3/19 was a small-B artifact)
---

# HYP-2296: the signed mutual-LRC floor is realized at a synchronized SHELL-PARTNER selected by the cut

## Context

Two signed mutual-loneliness gaps (THM-426: sign pattern ε ↔ cut (A,B) of the movers):
- **mover-only** `Gstar(S) = max over cuts M(W)`, `W={|v_i−v_j| same side}∪{v_i+v_j across}`;
- **observer-inclusive** `Gstar_full(S) = max over cuts M({v_i}∪W)` (observer–mover pairs add the
  bare speeds, sign-blind).  `Gstar_full ≤ min(M_obs, Gstar)`.

## Part A — the floor is a synchronized shell-partner (the structural law)

> **A.** At the optimal cut, the witness `t* = m/q` and the binding (minimizing) relative speeds are
> a **pair** `{a,b}` with `a + b = q = denom(t*)` and `‖a t*‖ = ‖b t*‖ = M = k/q`. I.e. the floor is
> a **synchronized shell-partner** (THM-425) at modulus `q`, and `M = k/q` is its loneliness gap.

**Why it is essentially forced + what is genuinely new.** The maximin of `min_w ‖wt‖` is attained at
a crossing of two tent functions `‖a·‖,‖b·‖`; a crossing at `t=m/(a+b)` (the *sum* crossing) makes
`{a,b}` a shell-partner mod `a+b`. The NEW content is empirical and uniform: the floor-setting
crossing is always the **sum** crossing (never a pure difference crossing or single-tent peak), and
**a sum relative-speed `a=v_i+v_j` exists only when the cut puts `{i,j}` across** — so the sign gauge
(=cut) is exactly the switch that exposes the shell-partner that sets the floor. This unifies the
pairwise-cut picture (THM-426), synchronization (THM-425), and the shell-partner witness program
(THM-420/THM-401).

**Verified** (`signed_lrc_block_certificate_monad_s3.py`, ad-hoc analyses): every minimizer found has
`a+b=q`. Mover-only: `(2,3,4,6,8)→3/19` binding `{9,10}` (q=19); `(4,5,8,10,15)→4/27` binding
`{4,23}` (q=27); `(2,4,7,10,11,12)→5/42` binding `{19,23}` (q=42). Observer-inclusive: block
`(5,6,7,8,9)→2/19` binding `{6,13}`; `(5,7,8,9,11)→3/29` binding `{11,18}` (q=29).

## Part B — the consecutive-block family (clean closed form for `Gstar_full`)

> **B.** For the mid-range block `B_n = {n−1, n, …, 2n−3}` (n−1 consecutive speeds, max `2n−3`),
> `Gstar_full(B_n) = 2/(4n−5)`, attained at a near-balanced cut with binding shell-partner
> `a+b = 4n−5 = 2(2n−2)−1`. VERIFIED exactly `n = 3..10`.

So `4n−5` is the Farey-successor modulus `2N−1` for the **doubled** system size `N = 2n−2` — the
observer-inclusive analogue of the LRC second clock `2/(2n−1)` (THM-401). `2/(4n−5) ≈ 1/(2n)`: when
all `n` runners must be mutually lonely with optimal signs, the block floor is asymptotically **half**
the single-observer floor `1/n`. (Lower bound `Gstar_full(B_n) ≥ 2/(4n−5)` is rigorous per `n` via
the explicit cut+`t*`; the closed form / its optimality over cuts is verified, not proved.)

## Part C — the infimum and its asymptotics (the open frontier, sharpening T764)

`inf_S Gstar(S)` (mover-only) and `inf_S Gstar_full(S)`:

| n | `inf Gstar` (best B) | `inf Gstar_full` (best B) | `2/(4n−5)` |
|---|----------------------|---------------------------|-----------|
| 3 | 1/2 | 2/7 | 2/7 |
| 4 | 2/7 | 2/11 (stable B≤20) | 2/11 |
| 5 | **1/5** (stable B≤24) | 2/15 (stable B≤16) | 2/15 |
| 6 | 6/41 and dropping (B≤18) | 3/29 < 2/19 (B≤13) | 2/19 |
| 7 | 5/42 and dropping (B≤13) | — | 2/23 |
| 8 | 2/19 and dropping (B≤11) | — | 2/27 |

**Clean facts (verified):**
- The observer floor `Gstar(S) ≥ 1/n` holds for **all** `gcd=1` `S` iff `n ≤ 5` (n=5 robust to
  `B≤24`); it first breaks at `n=6` and also breaks at `n=7,8`. Same `n=6` threshold for `Gstar_full`
  (`= 2/(4n−5)` for `n≤5`, drops below at `n=6`).
- For fixed `n` both infima are **bounded away from 0** (THM-429 Cor 1:
  `≥ 1/((n−1)(n−2))`), but `n·inf` is `< 1` and decreasing for `n ≥ 6`.

**Open:** the `n→∞` rate of `n·inf_S Gstar(S)` — does it stay `Θ(1)` (a true second floor at
`Θ(1/n)`) or decay toward the unconditional `Θ(1/n²)` regime? The minimizers are always a
small-speed cluster forcing `r_min ≥ n` (THM-429 Cor 2), with the floor a sum-crossing shell-partner
(Part A) — so the asymptotic question is: *how small can the gap of the best forced shell-partner be,
minimized over speed sets?*

### Matched-B evidence — leans toward the `Θ(1/n)` floor (monad-compute-2026-06-07)

The prior table mixes B per n, so `n·inf` was not comparable across `n`. Running the **same**
proven exhaustive search (`signed_lrc_inf_highB.search`) at a **common B=13** for `n=5..8`:

| n | inf Gstar (B=13) | `n·inf` | `2/n²` (lower bd) |
|---|------------------|---------|-------------------|
| 5 | 1/5  | 1.0000 = 1   | 0.0800 |
| 6 | 3/20 | 0.9000 = 9/10 | 0.0556 |
| 7 | 5/42 | 0.8333 = 5/6 | 0.0408 |
| 8 | 1/10 | 0.8000 = 4/5 | 0.0313 |

The `n·inf` column is **1, 9/10, 5/6, 4/5** = `30/30, 27/30, 25/30, 24/30`: decrements `3,2,1` — a
**decelerating** decrease, not a plunge. Per-n ceilings (best feasible B) tell the same story:
`n·inf ≤` 1.000 (n=5, robust B≤22), 0.889 (n=6, B17), 0.833 (n=7, **stable** B13=B14),
0.800 (n=8, 1/10 **stable** B13=B14). Firmed-up series `n·inf` = **1.000, 0.889, 0.833, 0.800**,
decrements `0.111, 0.056, 0.033` (≈halving each step) ⟹ consistent with convergence to a
positive constant (`≈ 0.75`), not to 0.

**Reading (NOT a proof):** across `n=5..8` the best-known `n·inf` stays in `[0.80, 1.0]` and sits
`~20–30×` **above** the `2/n²` unconditional lower bound, with the gap to that bound *widening* in `n`,
not closing. `n·inf` is decreasing but decelerating. This is computational evidence **for** a true
`Θ(1/n)` second floor (`n·inf → c > 0`, plausibly `c ≈ 0.8`) and **against** decay to the `Θ(1/n²)`
regime. Caveat: upper bounds only (true inf may dip lower; n=8 reached only B≤14), and `n≤8` cannot
settle an asymptotic — but no minimizer found drives `n·Gstar` anywhere near `0`. The `Θ(1/n²)` decay
would require an explicit cluster family with `q→∞`, `k` bounded *and* `r_min` blowing up; none of the
exhaustive minimizers (all near-consecutive small blocks, `q ≤ 42`) exhibit that.
Files: `04-computation/signed_lrc_inf_matchedB_monad.py`, `…inf_n8_firmup_monad…`
(+ `05-knowledge/results/*.out`).

## Sources
- `04-computation/signed_lrc_pairwise_inf_floor_monad_s3.py`, `…rmin_bound_…`, `…inf_highB_…`,
  `…families_…`, `…full_floor_…`, `…block_certificate_…` (+ `05-knowledge/results/*.out`)
- THM-429, THM-426, THM-425, THM-420, THM-401; reflection
  `07-reflections/the-signed-mutual-lrc-floor-is-a-cut-selected-shell-partner-s3.md`
