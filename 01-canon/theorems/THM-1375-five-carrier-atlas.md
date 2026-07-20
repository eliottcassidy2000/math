---
id: THM-1375
title: THE FIVE-CARRIER ATLAS — THE DENSITY METHOD CLOSES EXACTLY c ≤ 4, BY 13k < 42, AND IS VACUOUS AT c ≥ 5 — the owner window [1/(14x), 13/(14x)] has length 6/(7x), and THM-1155's bound gives the c-carrier budget (c−1 = k later combs, all > x) as exactly **k/(7−k)·(1/x)**: 2/5, 3/4, 4/3, 5/2, 6 for k = 2…6. Closure requires k/(7−k) < 6/7, i.e. **13k < 42** — which holds for k ≤ 3 (c ≤ 4, recovering codex THM-1169) and FAILS at k = 4. So five carriers is exactly where the uniform density argument dies, and the failure is not marginal: for x ≥ 5 the budget EXCEEDS the window by 6%…30% and worsens with x, so truncation (which closed the three-carrier row) cannot help — it would have to ENLARGE the window. The exact residual atlas is **x ≥ 5**; the exact bound still closes x = 1,2,3,4
status: exact rational arithmetic throughout; the budget formula k/(7−k) and the 13k < 42 criterion are algebra, the residual x ≥ 5 is computed exactly. This is a NEGATIVE closure result — it delimits the density method, it does not close the five-carrier row
source: opus-2026-07-20-S400 (owner: merge the five-comb atlas and multiplier-chart machinery)
depends_on: [THM-1155 (my k-comb density bound — the engine), THM-1169 (codex: three/four-carrier closure, recovered here as k ≤ 3), THM-1154, THM-1165 (the per-level budget j/(r(7−j)), of which this is the r = 1 case)]
scripts: 04-computation/five_carrier_atlas_opus_S400.py -> 05-knowledge/results/five_carrier_atlas_opus_S400.out
---

# THM-1375 — where the carrier program stops

## The budget is k/(7−k)

The first carrier x is safe throughout the owner window [1/(14x), 13/(14x)],
of length **6/(7x)**. The remaining k = c−1 carriers must cover it, and
THM-1155 gives L(1 − 2kλ) ≤ 2λΣ(1/wᵢ) with λ = 1/14. Since every later
carrier exceeds x, Σ(1/wᵢ) < k/x, so the budget is

> **k/(7−k) · (1/x)** — namely 2/5, 3/4, **4/3**, 5/2, 6 for k = 2, 3, 4, 5, 6.

| c | k | budget | window | verdict |
|---|---|---|---|---|
| 3 | 2 | 2/5 | 6/7 | **CLOSED** |
| 4 | 3 | 3/4 | 6/7 | **CLOSED** |
| **5** | **4** | **4/3** | 6/7 | **not closed** |
| 6 | 5 | 5/2 | 6/7 | not closed |
| 7 | 6 | 6 | 6/7 | not closed |

## The closure criterion is 13k < 42

k/(7−k) < 6/7 ⟺ 7k < 6(7−k) ⟺ **13k < 42**, i.e. k < 42/13 = 3.23…, i.e.
**k ≤ 3**, i.e. **c ≤ 4**. That recovers codex's THM-1169 (three- and
four-carrier rows) as exactly the range where the criterion holds, and shows
five carriers is precisely where it fails — with the familiar 13 and 7 in the
threshold, 42 = 6·7.

## The failure at c = 5 is vacuous, not marginal

For the four-carrier row the margin is 6/7 − 3/4 = **3/28**, small but real.
At five carriers the budget does not merely fall short — it **exceeds** the
window:

| x | budget ÷ window |
|---|---|
| 5 | 1.06 |
| 8 | 1.20 |
| 12 | 1.30 |

and grows with x. So the truncation device that closed the three-carrier row
(cutting the window at the owner threshold 1/(14p)) **cannot** help here: one
would need to enlarge the window, which is meaningless. A different tool is
required, not a sharper cut of the same one.

## The exact residual atlas

Using the true best-case later carriers x+1, …, x+4 rather than the crude
k/x, the exact bound still closes **x = 1, 2, 3, 4**. The residual is

> **x ≥ 5** — every first carrier from 5 upward.

## Relation to the rest of the programme

The budget k/(7−k) is the r = 1 case of THM-1165's per-level budget
j/(r(7−j)), so the carrier atlas and the substitution searches are running on
the same inequality. Both die at k = 7 (THM-1155's ceiling, 1 − 2kλ = 0), and
the carrier program dies earlier still — at k = 4 — because it must beat a
fixed window 6/7 rather than merely be positive.
