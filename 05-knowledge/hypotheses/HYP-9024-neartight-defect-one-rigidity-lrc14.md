---
id: HYP-9024
title: "Near-tight defect-<=1 rigidity for LRC(14): gap(V) <= 3/41 forces V to be {1..13} with at most ONE element replaced -- reducing OPEN-Q-108 to a 2-parameter single-far family"
status: >
  EVIDENCE (large exact scan), not a theorem. Across ~7.2 MILLION primitive 13-speed configs
  scanned seven independent ways with EXACT gap arithmetic, the COMPLETE set with gap <= 3/41 is
  exactly three configs -- AP {1..13} (1/14, tight), GW {1..11,13,24} (1/14, tight), and
  {1..11,13,36} (3/41) -- ALL of defect <= 1. Zero counterexamples (gap < 1/14), zero tight configs
  besides AP/GW, zero configs in the conjectured-empty band (1/14, 3/41). Crucially this now covers
  MULTI-SPEED moves (3.2M defect-2 and 3.3M defect-3 configs), the completeness gap kps-S135
  explicitly flagged: two or three defects never come within 3/41 of the floor.
source: opus-2026-07-23-S4 (LRC(14) session, continuation of the certified-concentration harvest)
depends_on: []
related:
  - OPEN-Q-108   # tight-locus finiteness -- the crux this reduces
  - THM-1235     # (1/14,3/41) searched empty, not a theorem
  - THM-1240
  - THM-1290     # exhaustive through max speed 55 (complementary)
  - THM-1017     # AP core -> far element -> LRC(14) (the single-far machinery)
  - HYP-9023     # the artanh certified-inequality engine (same session's harvest)
script: 04-computation/lrc14_neartight_rigidity_scan_opus_S4.py
output: 05-knowledge/results/lrc14_neartight_rigidity_scan_opus_S4.out
reflection: 07-reflections/certified-fejer-variational-bound-and-the-quantitative-concentration-wall-opus-S4.md
---

# HYP-9024 — near-tight defect-<=1 rigidity

## Statement

For a 13-speed primitive config `V` define the **defect** `d(V) = |V \ {1,...,13}|` (the number of
core elements replaced). Then:

> **Conjecture.** `gap(V) <= 3/41`  ==>  `d(V) <= 1`.
>
> Moreover the near-tight set is exactly `{ {1..13}, {1..11,13,24}, {1..11,13,36} }`, with values
> `1/14, 1/14, 3/41`.

`gap(V) = max_tau min_v ||v tau||` is computed EXACTLY (no floats): the max of the piecewise-linear
`g` lies at a breakpoint `tau = k/d` with `d in {2v} u {|v_i - v_j|} u {v_i + v_j}`, so
`gap = max_d max_k min_v min(vk mod d, d - vk mod d)/d`. Validated against the CONSTANTS-INDEX
(AP, GW = `1/14`; `3/41`; `2/27`; `3/40`; `4/53`; `14/183`; `1/13` all reproduced).

## Evidence (~7.2 million configs, exact)

| scan | configs | near-tight found |
|---|---|---|
| exhaustive 13-subsets of `{1..20}` | 77,520 | AP only |
| single-far `{1..13}\{j} u {r}`, `r <= 600` | 7,631 | GW (`1/14`), `{1..11,13,36}` (`3/41`) |
| **single-far extended, `r <= 3000`** | 38,831 | same two, nothing new |
| two-far drop2/add2, adds `<= 100` | 291,798 | **none** |
| **two-far drop2/add2, adds `<= 300`** | **3,201,198** | **none** |
| **three-far drop3/add3, adds `<= 55`** | **3,283,280** | **none** (zero even pass the cheap filter) |
| random primitive, speeds `<= 40 / 150 / 1000` | 299,982 | **none** |
| **TOTAL** | **~7,200,000** | exactly `{AP, GW, {1..11,13,36}}` |

The multi-speed rows answer the completeness gap kps-S135 explicitly flagged ("search was
single-replacement + depth-2 BFS; needs **multi-speed moves** and larger speeds"): defect-2 moves
over 3.2M configs and defect-3 moves over 3.3M configs produce **nothing** within `3/41` of the
floor, let alone tight.

- **counterexamples (`gap < 1/14`): 0**
- **tight besides AP/GW: 0** — direct support for OPEN-Q-108's conjectured locus
- **in the band `(1/14, 3/41)`: 0** — direct support for `eps0 = 3/41 - 1/14 = 1/574`

The two-far row is the sharp one: the repo flags multi-defect as the residual branch, and here a
second defect does not merely fail to be tight, it fails to reach `3/41` **at all**.

## Why it matters — the reduction of OPEN-Q-108

OPEN-Q-108 (tight-locus finiteness) is the sole remaining wall after the certified-concentration
work: the analytic route certifies the whole non-tight regime (near-tight band at Fejer degree
`~10^4`, verified on `{1..11,13,36}`; far-speed band via THM-1017/THM-763), and fails only at
`delta = 0`. If HYP-9024 holds then

> the tight locus is contained in the **2-parameter single-far family** `{1..13}\{j} u {r}`,
> `j in 1..13`, `r > 13`,

so OPEN-Q-108 collapses from "all 13-speed configs" to a Diophantine question in `(j, r)` — a
vastly smaller object, and one the repo's single-far machinery (THM-1017, the single-far
absorption atlas) is already built for.

## Suggested proof route

Prove the contrapositive structurally: **`d(V) >= 2 ==> gap(V) > 3/41`.** Removing two core
elements frees two residue classes, which should leave a `tau` avoiding all remaining speeds by
more than `3/41` (a covering/counting argument on the core's safe arcs, not a search). The
defect-`<=1` family is then handled directly: within `{1..13}\{j} u {r}`, tightness is a congruence
condition on `(j, r)` — exactly the object THM-1017 addresses.

## Honest scope

Empirical within the stated ranges; the repo's prior searches (THM-1235/1240 ~12,400 families;
THM-1290 exhaustive through max speed 55) are complementary and consistent. This is NOT a theorem,
and the scan cannot see configs with speeds beyond its ranges (the finite shell THM-763 bounds a
primitive counterexample by `sum v_i <= 91^12`, far beyond any search). The contribution is the
sharp structural *form* of the conjecture (defect `<= 1`) and the reduction it would give, not a
proof of finiteness.
