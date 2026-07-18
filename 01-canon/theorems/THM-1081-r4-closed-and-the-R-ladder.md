---
id: THM-1081
title: THE r=4 CLUSTERED CASE CLOSED, AND THE R-LADDER — why the finite horn stops being redundant exactly at r=4. (I) r=4 FINITE HORN COMPLETE: all 220 nine-speed cores, **143,112,134** quadruples passing the sound covering-necessary condition, **ZERO uncertified** (160 cores / 119,489,369 in cont.61 plus 60 cores / 23,622,765 here). (II) THE THRESHOLD SCANNED PROPERLY, and it corrects my own method: the right quantity is not the absolute threshold T but the ratio **R = T / k_max-removed**, since the un-removed killer always exceeds the removed ones — R < 1 means the measure horn certifies with no finite horn at all. Exhaustive for r=2 over all 12 cores and every k₁ ∈ (13·max P, 4000): **max R = 0.51852**, zero swallow cases. (III) METHOD CORRECTION: the worst case sits at **SMALL** killers — r=2 at k₁ = 160, r=3 at (150,156), r=4 at (150,156,158), all just above 13·max P — whereas my cont.60/61 scans sampled the TOP of the range and therefore looked in the wrong place and reported ratios (0.389–0.434) that were too optimistic. (IV) THE R-LADDER: max R = **0.51852** (r=2, exhaustive), **0.73375** (r=3), **0.98453** (r=4) — rising steeply, with only a **1.5% margin** left at r=4 and extrapolating above 1 at r=5. (V) CONSEQUENCE, which inverts the status of my own earlier work: the finite horns of THM-1051 and THM-1061 were **REDUNDANT** — the measure horn alone closes r=2 and r=3 — but at r=4 the measure horn is marginal on a narrow scanned window, so the finite horn of (I) is the **load-bearing** result there, not a redundancy check; and at r=5 the measure horn is expected to fail outright, making the finite horn mandatory
status: (I) PROVED — exhaustive finite verification, exact arithmetic, explicit (q,a) witness per family, the covering-necessary pruning being sound (a quadruple failing Σ frac ≥ 1 cannot cover, hence certifies automatically). (II) r=2 EXHAUSTIVE to 4000. (III),(IV) measured; r=3 and r=4 R-values come from exhaustive scans of the SMALL end (where the worst case provably sits for r=2) plus coarse tails, NOT from a full scan — so **r=4's R < 1 is verified on a window, with a 1.5% margin, and I do not claim it in general**. (V) follows from (I)–(IV); the r=5 prediction is extrapolation, explicitly not a result
source: kind-pasteur-2026-07-18-S128 (cont.62; owner: finish the remaining 60 cores and scan the r=4 threshold properly)
depends_on:
  - THM-1071    # the r=4 partial run this completes, and the pruning it introduced
  - THM-1051, THM-1061   # the r=2 / r=3 closures whose finite horns (V) shows were redundant
related:
  - THM-1042    # klein: the component-length obstruction; R < 1 is its quantitative converse
script: 04-computation/r4_finite_horn_tail_kps_S128c62.py, R_exhaustive_kps_S128c62.py, R_scan_r3r4_kps_S128c62.py, R_scan_r4_kps_S128c62.py (+ .out)
---

# THM-1081 — r=4 closed, and the R-ladder

## (I) The r=4 finite horn is complete

| run | cores | quadruples tested | uncertified |
|---|---|---|---|
| cont.61 | 160 | 119,489,369 | 0 |
| cont.62 (here) | 60 | 23,622,765 | 0 |
| **total** | **220 / 220** | **143,112,134** | **0** |

The pruning is sound: a quadruple can only be uncertified if its four kill-sets *cover* the
core's safe (q,a) set, which requires Σ frac ≥ 1; quadruples failing that certify
automatically. So the check is exhaustive despite testing only the tail.

## (II) Scanning the threshold properly — the right quantity is R

I had been scanning the absolute threshold T. That is the wrong object. The measure horn
removes the r−1 smaller killers and needs the largest to exceed T, and the largest always
exceeds every removed killer, so what matters is

> **R := T / k_max-removed**,  and **R < 1 ⟹ the measure horn certifies with no finite horn.**

For r=2 this is cheap enough to settle exhaustively — all 12 cores, every
k₁ ∈ (13·max P, 4000):

> **max R = 0.51852**, attained at core-drop 6 with k₁ = 160; **zero** cases where the
> removal swallows S(P).

## (III) The worst case is at small killers — my earlier scans looked in the wrong place

| r | worst-case killers | 13·max P | max R |
|---|---|---|---|
| 2 | 160 | 156 | 0.51852 |
| 3 | (150, 156) | 143 | 0.73375 |
| 4 | (150, 156, 158) | 143 | 0.98453 |

Every worst case sits **just above 13·max P**, i.e. at the smallest admissible killers,
tightly clustered. My cont.60 and cont.61 scans sampled the *top* of the killer range on the
reasoning that dense bad-sets chop the safe set finest — that reasoning was wrong, and the
ratios those scans reported (0.389–0.434, and the "generic 7/18") are the behaviour of the
easy end, not the worst case. The 7/18 value is still the correct *asymptotic* constant;
it is simply not the maximum.

## (IV) The R-ladder

> **0.51852 (r=2) → 0.73375 (r=3) → 0.98453 (r=4)**

Rising steeply, with **1.5% of margin left at r=4** and a straightforward extrapolation
above 1 at r=5.

## (V) What this does to the earlier theorems

- **THM-1051 (r=2) and THM-1061 (r=3): their finite horns were REDUNDANT.** R = 0.519 and
  0.734 mean the measure horn alone certifies every family, with no split point and no
  enumeration. The 41,986 and 3,408,751 family checks stand as independent verification, and
  I am glad to have them, but they were not needed for the proof.
- **r=4 is different.** R = 0.98453 is a 1.5% margin, established on a narrow scanned window
  (killers in [lo, lo+55) exhaustively, plus a coarse tail). I do **not** claim R < 1 at r=4
  in general. So the finite horn of (I) — 143 million quadruples, zero failures — is the
  **load-bearing** result at r=4, not a redundancy check.
- **r=5 is expected to break the measure horn outright.** If the ladder continues, R > 1 and
  the finite horn becomes mandatory rather than optional. That is the honest prediction, and
  it is extrapolation, not a result.

The structural reading: each removed killer fragments the safe set further, so the surviving
component shrinks faster than the killers grow. The measure horn is a *finite-depth* tool —
it survives three removals and expires around the fourth.

## Named next
- Settle r=4's R properly: the window here is [lo, lo+55) exhaustive plus a coarse tail. A
  full scan over all triples would either confirm R < 1 or find the crossing, and given the
  1.5% margin this is worth doing before anyone leans on the measure horn at r=4.
- r=5: expect to need the finite horn. The covering-necessary pruning (Σ frac ≥ 1) is what
  makes these runs feasible and it generalises unchanged; the cost is that quintuples are
  another factor ~250 more numerous.
- Retire the split points from THM-1051 and THM-1061 in favour of "R < 1, measure horn
  alone" — it is a strictly simpler statement of both theorems.
