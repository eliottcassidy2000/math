# The LRC source-menu is exactly half of A000568 — until n=7, where it collapses

oracle-2026-06-01-S520

Building on the S512 capstone (HYP-1987): LRC, in the observer-marked
reformulation, asks whether a primitive speed set's clock-walk reaches a SOURCE
cell, where the runner sub-tournament is a half-turn tournament of `n-1` points
confined to the safe arc `[1/n, 1-1/n]`. S512 computed the *reachable* menu of
such source classes as `1, 2, 6, 6` for `n=4..7`, each from a single speed box.

This session asked two follow-ups S512 left open: **is each count converged
(or just a lower bound from one box), and what is n=8?**

## What the saturation study found

Re-running `analyze_targets` over growing boxes:

```
n=5: 2 stable through max_speed 26 (14161 sets)   -> converged
n=6: 6 stable through max_speed 15 ( 2981 sets)   -> converged
n=7: 6 stable through max_speed 14 ( 2996 sets)   -> converged
n=8: 11 at boxes 9,11  ->  12 at box 12 (H=24 appears) -> NOT converged, >=12
```

So the menu sizes for `n<=7` are box-stable (strong evidence they are the true
values), but **n=8 is only a lower bound `>=12`** — a new class surfaced the
moment the box grew. The honest extended sequence is `1, 2, 6, 6, >=12`.
LRC is re-verified on the tournament side for n=8 as well: **0 source-avoidance
failures** in every box tested.

## The pattern worth flagging

Line the menu up against `A000568(n-1)` (the ambient iso-class count):

```
n         4     5     6     7      8
menu      1     2     6     6     >=12
A000568   2     4    12    56     456
ratio   1/2   1/2   1/2  0.107  >=0.026
```

For `n = 4, 5, 6` the reachable menu is **exactly `A000568(n-1)/2`** — the LRC
win-set is precisely half of all iso classes. Then at `n=7` it **collapses**:
6 instead of 28. At `n=8` it collapses further.

The "exactly half" coincidence is almost certainly the complement involution:
for small `n` the arc is permissive enough that a class is reachable iff its
complement is not (or some Z_2 pairing), so the menu picks out one of each
complement pair — the same `/2` that runs through the merged metagraph
`V_merged = (A000568 + SC)/2` (CLAUDE.md). What is striking is not that it holds
but that it **breaks sharply at n=7**, the same threshold where so much else in
this project first goes non-generic (width formula fails, level/decreasing
metagraph edges appear, odd holes in E_7). The arc length `L = 1-2/n` crosses a
structural value where most of `A000568(n-1)` becomes *unreachable* faster than
it grows.

## Why this matters for the proof program

HYP-1987's hope is that the win-set is a **vanishing fraction** of an enormous
multiplicatively-shaped target, so any LRC counterexample must thread a huge
complement. This data supports that hope *and* sharpens it: the vanishing does
not start gradually — it switches on at `n=7`. The pre-n7 regime (`menu =
A000568(n-1)/2`) is exactly the "easy" range where LRC is folklore-verifiable;
the collapse marks the onset of the hard regime. A proof mechanism should
explain the `/2 -> collapse` transition, not just the asymptotic smallness.

## Methodological caveat (for HYP-1987 and future menu counts)

A reachable-menu count from a single speed box is a **lower bound**, not the
exact menu, until box-saturation is demonstrated. n=8 is the cautionary case:
the single-box snapshot (11) undercounts the saturating value (>=12, true value
still open). Always run the box-growth check before quoting a menu size as exact.

See: HYP-1987, `07-reflections/lrc-tournament-analysis-capstone-s512.md`,
`05-knowledge/results/lrc_source_menu_saturation_s520.out`,
`04-computation/lrc_source_reachability_n8_s520.py`.
