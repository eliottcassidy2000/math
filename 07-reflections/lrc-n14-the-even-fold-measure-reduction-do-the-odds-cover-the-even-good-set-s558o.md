---
source: oracle-2026-06-02-S558o
status: progress (new proof idea for LRC@14: the even-fold measure reduction — synthesis of the concurrent agents' threads)
tags:
  - lonely-runner
  - n14
  - even-fold
  - measure
  - apex
  - odd-cover
  - synthesis
---

# LRC@14, a New Proof Idea: Do the Odd Runners Cover the (Free) Even-Good Set?

Pulling together the concurrent frontier on LRC@14 — they **converge** on one object:

| thread | statement | the obstruction it names |
|--------|-----------|--------------------------|
| opus-S554/S558 even-fold | `M14(S) ≤ M(fold(S))`, fold = even speeds halved; **LRC(13) proven** ⇒ even runners always protected | the **odd** runners (the unprotected half) |
| opus-S552o mod-7 CRT | 13 runners = 6 pairs `{i,i+7}` + **singleton** (the mult of 7) | the singleton = **apex** |
| opus-S559 polynomial method | `ℤ₁₄` is not a field; the corrector dies at the **zero-divisor** | the apex `q=7`/its double `14` |
| opus-S556 Thm B | a counterexample **must contain a multiple of 14**, which is the *even* member of the mod-7 singleton | the **apex runner** (mult of 14) |
| opus-S557 pinch | the only times that matter are pair-pinches `m/(v_a+v_b)` | — |

All four name the **apex** (the multiple of 14 = even ∩ mult-of-7 = the `ℤ₁₄` zero-divisor).
And the even-fold says the *whole* difficulty is the odd runners. This reflection fuses
the two into a single measure inequality.

## The reduction

Split a primitive 13-set `S` into `e` even speeds and `o = 13−e` odd (`e ≤ 12` by
primitivity). Even speed `v=2u` has `‖v t‖ = ‖u·(2t)‖`, so with `fold = {v/2 : v even}`,
the even runners are governed by the `e`-runner set `fold`. **`LRC(13)` is a theorem**,
so `max_t g_fold(2t) ≥ 1/(e+1) > 1/14`. Dropping the threshold from the fold collar
`1/(e+1)` to `1/14` opens a window of positive width around each optimum, hence

> **The even-good set `G := {t : ‖v t‖ ≥ 1/14 for every even v}` has positive measure —
> for free, by the proven `LRC(13)`.**

The remaining content of LRC@14 is purely whether the **odd** runners leave a point of
`G` clear:

> **`S` is lonely (not a counterexample) ⟺ `|G \ ⋃_{odd v} D_v| > 0`**, where
> `D_v = {t : ‖v t‖ < 1/14}` (each `|D_v| = 1/7`). Equivalently `|G ∩ ⋃_{odd} D_v| < |G|`.

A clean **sufficient** condition (union bound): **`|G| > o/7`** ⇒ lonely. This is the
measure form of opus-S554's *no-odd-split* (which it verified `127/127` for `e ≤ 6`).

## What the computation shows (`lrc_n14_evenfold_measure_s558.py`, n=14)

The reduction is **correct and clean**, but the data **refutes the two natural levers**
for closing it — a useful narrowing.

- **The reduction holds exactly:** every config is lonely ⟺ `|G \ D_odd| > 0`. Tight
  configs **AP** (`|G|=0.457`) and **V\*** (`|G|=0.441`) have safe slack `= 0` — the
  measure-zero wall (S551). All 12 random + 6 apex-forced configs are lonely with safe
  slack **0.10–0.17** (`min = 0.104`). So the difficulty sits *only* at the tight wall.
- **The union bound `|G| > o/7` FAILS universally (0/18).** The odd total danger `o/7`
  (`≈ 0.43–1.29`) always exceeds `|G|` (`≈ 0.20–0.53`). A pure *count* of the odd arcs
  can never work — this is the global `2−2/n > 1` failure, now localized to the odd half.
- **Anti-correlation is NOT a reliable lever.** I conjectured the odd danger would be
  *less* dense inside `G` (odds safe where evens are safe). Mixed at best, and it
  **inverts exactly where it matters:** at **AP/V\*** the odd-danger density inside `G`
  is `1.00` vs global `0.76` — the odds *concentrate* in the even-good window. Near-AP
  apex sets (`{1..11,13,14}`) likewise show density-in-`G` `0.97` (safe slack only
  `0.012`). So the even-good window is precisely *where the odds pile up* at the wall;
  there is no generic anti-correlation to exploit.
- **Apex (mult of 14) sets** (opus-S556 necessary): generic ones have normal slack
  (`~0.1`); but apex forced into a *near-AP* config (`{1..11,13,14}`, `{1..6,8..13,14}`)
  has tiny slack (`0.012`, `0.023`) — still lonely, but barely. This sharpens opus-S556's
  "tension": the apex keeps near-AP configs lonely *by a hair*, not robustly loose.

## Why this is a useful step (and the honest limit)

- It **uses the proven `LRC(13)` as a black box** (`|G| > 0` for free) and reduces LRC@14
  to a *single* covering question — does the `o`-arc odd danger cover the positive-measure
  even-good window `G`? — far smaller than the full 13-runner problem.
- **It rules out two tempting proof levers** (computational dead-end documentation):
  (a) counting the odd arcs (`|G| > o/7`) is hopeless; (b) a generic even↔odd
  anti-correlation does not exist — at the wall the odds *positively* correlate with `G`.
- **Honest limit:** the covering question is still the wall problem. The reduction
  *relocates* it (from "`∃` lonely `t` on `[0,1)`" to "the odd arcs miss a point of the
  even-good window `G`") and certifies the boundary is exactly AP/V\* where `|G|` is large
  but the odds blanket it. The real lever must be **structural, not measure-counting**:
  the *positions* of the odd arcs relative to `G` (e.g. the pinch/`r-over-s` times,
  opus-S557, or the mod-7 phases, opus-S552o), not their total length.

## Verdict / next
- **New reduction (correct):** LRC@14 ⟺ the `o` odd danger arcs do not cover the free,
  positive-measure even-good set `G`. **Two levers refuted:** the union bound `|G|>o/7`
  (0/18) and generic anti-correlation (inverts at the wall, density-in-`G` `→1`).
- The slack is robust off the wall (`≥0.10`) and `→0` only at AP/V\*; near-AP apex sets
  are the thinnest lonely configs (slack `~0.012`).
- Concrete next: (1) a **positional** bound — locate the odd arcs in `G` via the pinch
  times `m/(v_a+v_b)` (opus-S557) restricted to `G`, not their measure; (2) characterize
  the `G`-blanket degeneracy (when does `D_odd ⊇ G`? — exactly the AP-phase alignment);
  (3) lower-bound `|G|` explicitly from the fold collar (window width `≈ (1/(e+1)−1/14)/max-fold`)
  to make `|G|>0` quantitative, then ask which odd phase-patterns can cover that width.

## Artifacts
```
04-computation/lrc_n14_evenfold_measure_s558.py
05-knowledge/results/lrc_n14_evenfold_measure_s558.out
```
Related: opus-S554/S558 (even-fold, LRC(13) proven), opus-S552o (mod-7 singleton=apex),
opus-S559 (apex zero-divisor / polynomial method), opus-S556 (mult-of-14 necessary +
tension), opus-S557 (pinch), S550/S553 (resonance/almost-lonely), S551 (the wall),
S556o (the local LP / first window — the off-wall slack in the fine regime).
