---
id: HYP-9024
title: "Near-tight defect-<=1 rigidity for LRC(14): gap(V) <= 3/41 forces V to be {1..13} with at most ONE element replaced -- reducing OPEN-Q-108 to a 2-parameter single-far family"
status: >
  PARTLY PROVED (defects 2 and 3 are now THEOREMS; see the closure section), plus large-scan
  evidence elsewhere. Across ~7.2 MILLION primitive 13-speed configs
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

## PROVED: defects 2 and 3 are CLOSED

The suggested contrapositive route is now realized for `d = 2, 3`, combining klein-S415's covering
lemma with one further elementary step.

**Band-width criterion (new, elementary, exact).** Write `V = W u {r}` with `r` the largest far
speed. `gap(V) <= h`  iff  `Lon_h(W)` is covered by `D_r = {tau : ||r tau|| <= h}`. But `D_r` is a
union of bands of width `2h/r` separated by gaps of width `(1-2h)/r > 0`, so an interval longer than
`2h/r` cannot lie inside a single band and must meet a gap. Hence every arc of `Lon_h(W)` has length
`<= 2h/r`, giving

```
    r  <=  2h / L_max(W).
```

This bounds the LARGEST far speed directly, and is sharper than bounding `sum 1/r` because it uses
the actual lonely structure of `W`.

**Defect 2 (THEOREM).** klein's lemma (`k=2`, factor `29/6`) gives `min(far) <= 70`. For each of the
78 cores and each `s in 14..70` (4,446 pairs, exact rationals) the criterion gives
`r <= 2h/L_max(C u {s})`, with **max `r_max` = 73** (worst: `drop(6,10)`, `s=40`,
`L_max = 0.002003`). So both far speeds are `<= 73`, a finite region contained in the exhaustive
adds-`<=300` scan (3,201,198 configs, zero hits). **No defect-2 config has `gap <= 3/41`.**
(`lrc14_defect2_closure_opus_S4.py`)

**Defect 3 (THEOREM).** Recursion: klein's lemma (`k=3`, factor `23/6`) gives `s1 <= 112`; the lemma
again with core `C u {s1}` (`k=2`) gives `s2 <= 141`; the band criterion gives `s3 <= 82`. Since
`s1 <= s2 <= s3 <= 82`, ALL far speeds are `<= 82`. Exhaustive scan of that proved region — all 286
ten-cores x all triples from `14..82` = **14,984,684 configs** — finds **zero** near-tight.
**No defect-3 config has `gap <= 3/41`.** (`lrc14_defect3_{bounds,closure_scan}_opus_S4.py`)

**Remaining.** `d = 1` (the 2-parameter family `{1..13}\{j} u {r}`, scanned to `r <= 3000`: only GW
tight and `{1..11,13,36}` at `3/41`) and `d >= 4`.

*Corrected note (tested, `lrc14_defectk_quicktest_opus_S4.py`).* I first expected the recursion to get
EASIER as `d` grows (smaller core => larger `L_max(C)`). **That is false.** The minimum over cores of
`B_k(C) = L_max(C)(1-2kh)/(2h)` does NOT grow with `k` —

```
k = 2:0.0285  3:0.0266  4:0.0298  5:0.0253  6:0.0131      (factor shrinks 29/6, 23/6, 17/6, 11/6, 5/6)
```

— because the worst core is adversarial and its `L_max` stays `~0.005-0.01` while the factor decays.
**Recorded negative:** the immediate test "`min_C B_k(C) > k/14`" (all far speeds `>= 14` force
`sum 1/s_i <= k/14`) FAILS at every `k` (`0.03` vs `0.14-0.43`), so the lemma alone never closes a
defect level — `d=2,3` genuinely needed the band criterion PLUS an exhaustive scan of the finite
region. For `d=4` that region is `~618M` configs (715 cores x C(69,4)), so brute force is infeasible;
`d>=4` needs better pruning or a structurally different argument, and `d>=7` fails klein's hypothesis
`h < 1/(2k)` outright. Consequently the reduction of OPEN-Q-108 to a defect-1 question is conditional
on `d >= 4` as well.

## Honest scope

Empirical within the stated ranges; the repo's prior searches (THM-1235/1240 ~12,400 families;
THM-1290 exhaustive through max speed 55) are complementary and consistent. This is NOT a theorem,
and the scan cannot see configs with speeds beyond its ranges (the finite shell THM-763 bounds a
primitive counterexample by `sum v_i <= 91^12`, far beyond any search). The contribution is the
sharp structural *form* of the conjecture (defect `<= 1`) and the reduction it would give, not a
proof of finiteness.
