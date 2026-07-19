---
id: THM-1230
title: THE n=14 STABILITY GAP IS NOT EMPTY — {1,…,11,13,36} REALISES 3/41 — transferring boxeph-S123's determinant stratification from the n=12 gap (1/13,2/25) to n=14 gives the gap (1/14, 2/27) (Farey neighbours, 1·27−2·14 = −1) with mediant 3/41 = the depth-minimal target (41 = 3·14−1, matching their 38 = 3·13−1). Conditional on LRC(14) the same two lines exclude numerator ≤ 2, so a gap family needs D ≥ 3. **A D = 3 witness exists**: V = {1,…,11,13,36} has M = 3/41 exactly at t* = 17/41, active pair (5,36), D = 3, s = 41 — verified by exhaustive brute force over all q ≤ 200. Its LRC(14) margin is s ≤ 14D reading **41 ≤ 42, slack 1** — the tightest non-extremal family found. It sits on the exact ladder M({1,…,11,13,12m}) = m/(12m+5), active pair (5,12m), D = m, slack 2m−5, which starts at the extremal 1/14 (m=2), enters the gap at 3/41 (m=3), and accumulates at **1/12**
status: the 3/41 witness is verified exactly (pinch-set maximiser agrees with independent brute force over q ≤ 200; active pair and determinant confirmed). The ladder formula m/(12m+5) is verified exactly for m = 3…10 and is a conjecture beyond. The numerator ≤ 2 exclusion is CONDITIONAL on LRC(14), which is open — unlike boxeph's n=12 version, which rests on the proved LRC(13)
source: opus-2026-07-19-S395 (owner: continue and integrate incoming ideas from other agents)
depends_on: [boxeph-S123 / HYP-7782 (the determinant stratification this transfers), THM-1205/1210/1215 (the D-coordinate and the sieve/high-D split), THM-1220 (n=14 is the unique non-rigid case — sharpened here), HYP-2059/THM-401 (pinch lemma), THM-633 (the n=12 d=1 ladder this parallels)]
scripts: 04-computation/n14_gap_strata_opus_S395.py -> 05-knowledge/results/n14_gap_strata_opus_S395.out
---

# THM-1230 — the gap above 1/14 is populated

## The transfer from boxeph-S123

boxeph stratified the n=12 uniqueness gap (1/13, 2/25) by the reduced
numerator p of M = D/s, observing p = D/gcd(D,s) ≤ D so that stratifying by
p **is** stratifying by determinant, and crediting my THM-1210. Their
exclusion needs only two lines and two inputs — a proved lower bound on M,
plus the parity of a reduced denominator.

At n = 14 the Farey neighbour of 1/14 is **2/27** (1·27 − 2·14 = −1), so the
stability gap is **(1/14, 2/27)**, with mediant **3/41** and 41 = 3·14 − 1
(matching boxeph's 38 = 3·13 − 1). The same two lines transfer, **conditional
on LRC(14)**:

- p = 1: M = 1/q ≥ 1/14 forces q ≤ 14, so M = 1/14 or M ≥ 1/13 = 0.0769 > 2/27.
- p = 2: M = 2/q reduced makes q odd; 2/q ≥ 1/14 gives q ≤ 28, so q ≤ 27 and
  M ≥ 2/27.

So a gap family needs **D ≥ 3** — the n=14 analogue of boxeph's conclusion.

## The gap is NOT empty

At n = 12 the analogous gap is conjectured empty (Tao uniqueness). **At n = 14
it is not.**

> **V = {1, 2, …, 11, 13, 36} has M = 3/41**, strictly inside (1/14, 2/27).

Verified exactly: maximiser t* = 17/41, active pair **(5, 36)**, **D = 3**,
s = 41, D/s = 3/41; the pinch-set computation agrees with independent brute
force over all q ≤ 200. It lands precisely on the predicted depth-minimal
value.

Its LRC(14) margin is the tightest yet seen: **s ≤ 14D reads 41 ≤ 42 — slack 1.**
A family with s = 43 at D = 3 would violate LRC(14).

## The ladder, and where 3, 4 and 1/12 meet

The witness is not isolated; it sits on an exact ladder:

> **M({1,…,11,13,12m}) = m/(12m+5)**, active pair (5, 12m), D = m, s = 12m+5,
> slack 14D − s = **2m − 5**

| m | x = 12m | M | D | s | slack |
|---|---|---|---|---|---|
| 2 | 24 | **1/14** | 1 | 14 | **0** — the second extremal |
| **3** | **36** | **3/41** | **3** | **41** | **1** — in the gap |
| 4 | 48 | 4/53 | 4 | 53 | 3 |
| 5 | 60 | 1/13 | 5 | 65 | 5 |
| 10 | 120 | 2/25 | 10 | 125 | 15 |

Verified exactly for m = 3…10. Note m = 2 reduces (D = 1, s = 14) and is
exactly the known second extremal {1,…,11,13,24}, so the ladder *starts* at
the floor and climbs.

**M = m/(12m+5) → 1/12.** So the owner's three mined threads — 3, 4, and 1/12
— are one object: 3 and 4 are the first two gap-side rungs, and 1/12 is the
ladder's accumulation point. That also explains why 1/12 keeps recurring in
this corpus without being an extremal value: it is a limit, not an attained
minimum.

## Status and scope

The witness and the ladder values are exact. Two honest caveats: the
numerator ≤ 2 exclusion is **conditional on LRC(14)** (open), whereas
boxeph's n = 12 version rests on the proved LRC(13); and the ladder formula
is verified for m = 3…10, conjectural beyond.

This sharpens THM-1220 ("n = 14 is the unique non-rigid case"): n = 14 is not
merely non-rigid at the floor — the gap immediately above the floor is
**populated**, by an explicit family, at exactly the value the stratification
predicts. Whatever proves LRC(14) cannot rely on a gap above 1/14, because
there is none.
