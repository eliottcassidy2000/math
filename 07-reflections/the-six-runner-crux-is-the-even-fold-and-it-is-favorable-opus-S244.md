---
source: opus-2026-07-11-S244
status: A creative synthesis (past-work mining) — the ~6-runner crux IS the oracle's even-fold (S558o),
  iterated over all small primes to the ≤6 coprime-to-30030 core. For the residual (divisor-complete)
  families the core SPREADS and misses the structured good set (danger density 0.23 << 1 ⟹ lonely); the
  "blanket" obstruction (density → 1) is exactly the AP/V* wall = bucket A, dispatched by t=1/14. Unifies
  five threads.
tags:
  - lrc14
  - six-runner
  - even-fold
  - coprime-core
  - favorable
  - synthesis
  - past-work
---

# The ~6-runner crux is the even-fold — and it is favorable for the residual

**opus-2026-07-11-S244.** Owner: work the ~6-runner problem creatively, look back at past work. Mining the
repo turned up the exact object — the **oracle's even-fold** (S554–S558o) — and connecting it to my recent
work (coprime core, S243) with S239 reframes the crux as *favorable*.

## The past work: the even-fold (oracle S558o)

`M14(S) ≤ M(fold(S))`, `fold` = even speeds halved. Since **LRC(13) is proven**, the even runners are
protected: the **even-good set** `G = {t : ‖v t‖ ≥ 1/14 for every even v}` has **positive measure for free**.
The whole of LRC(14) reduces to: *do the odd runners leave a point of `G` clear?* — i.e.
`LRC(14) ⟺ |G \ ⋃_{odd} D_v| > 0`. S558o refuted two levers (the union bound `|G| > o/7`, and generic
even↔odd anti-correlation) and found the obstruction is **positional**: at **AP/V*** the odd arcs *blanket*
`G` (density → 1.00).

## The synthesis: iterate the fold to the coprime core

The even-fold folds out `div-by-2`. **Iterate over all small primes**: fold out `div-by-3, 5, 7, 11, 13` —
each protected by LRC(≤13). What survives is the **coprime-to-30030 core** (`≤ 6`, opus-S243). So the
even-fold's "odd runners" sharpen to the ≤6 coprime core, and:

> **`LRC(14) ⟺` the ≤6 coprime core does not blanket the structured good set `G' = {t : all non-core safe}`.**

This is the *same object* five independent threads reach: the even-fold (S558o), the coprime core (S243),
`spread = bad coverer` (S239), klein's ~6-odd shrink (S263), and mac-mini's ≤6 decorrelation lifts (cont.49).

## The favorable reframing (verified)

For **divisor-complete (residual)** families, the core-danger density in `G'` is **0.23 (mean), always < 1**
— the ≤6 core **spreads and misses** `G'`, so a point of `G'` is clear ⟹ **lonely**. (Core size ≤ 4 in the
sample.) By contrast the **AP `{1..13}` and V*** have density = **1.000** — the *blanket*. But those are
**bucket A** (no multiple of 14), dispatched by the elementary `t = 1/14` witness, and are **not**
divisor-complete.

So S558o's blanket obstruction is the **AP/V* wall = bucket A**, handled elementarily. For the residual, the
core is *spread* (S239 bad coverer), its danger-density in `G'` is far below 1 (0.23), and it misses `G'`.
**The ~6-runner crux is the favorable direction of the anti-concentration** — quantitatively robust (density
0.23 ≪ 1, huge margin), not marginal.

## Honest status

This is a **synthesis + favorable reframing**, not a proof. The remaining content is "spread core ⟹
danger-density in `G' < 1`," which is the anti-concentration in its favorable direction — S238 showed it has
no bounded-window shortcut, but here it is far from marginal (0.23). The value: the ~6-runner crux, past work
(even-fold), and my recent bricks (auto-safe, coprime core, spread reframing) are **one object**, and for the
residual it is favorable — the difficulty (blanket, density → 1) lives only at the AP/V* wall, which is
already dispatched. So the residual's ~6-runner problem is not the hard blanket case; it is the spread case,
where the core robustly misses `G'`.

→ oracle S554–S558o (the even-fold, the reduction), opus-S243 (coprime core = the iterated fold), opus-S239
(spread = bad coverer), opus-S241/S242 (auto-safe + pigeonhole), klein S263 (~6-odd shrink), mac-mini cont.49
(≤6 lifts), the AP/V* wall = bucket A (S233). Files: `lrc14_sixrunner_evenfold_synthesis_opus_S244.py`
(+`.out`).
