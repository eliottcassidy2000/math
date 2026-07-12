---
id: THM-709
title: The doubling-tight locus is a SINGLETON — among all doubling rearrangements AP[e→2e] of {1..13}, ONLY e = 12 (THM-708's point {1..11,13,24}) has M = 1/14; every other e ∈ {7..11,13} escapes (M = 1/11, 2/23, 2/23, 2/27, 2/25, 2/27, exact q<500 scans), all double-doublings escape, and e = 7 flips to COVERING (14|14) with a strict witness — so the non-AP tight boundary that the clean-ruler supply must exclude is (within this family) exactly ONE point
status: VERIFIED-EXACT (per-family exact M over all p/q, q < 500; escapes attained at explicit small-q times). MECHANISM (noted, not proved): 2e = 24 is the unique double adjacent to BOTH flanking escape rulers (24 = 23+1 = 25−1), aliasing onto ±(runner 1)'s position at each — the freed 12-slot buys no new clearance; for other e the doubled runner lands fresh at some flanking ruler and the freed slot is exploitable. CONSEQUENCE: THM-708's tight point is ISOLATED; the tight locus relevant to kps THM-707's clean-ruler supply = the AP line plus isolated aliasing points, not a positive-dimensional family — the supply's excluded boundary is thin, as the census (72/72) suggested.
source: mac-mini-2026-07-09-S65 (cont.33, 2026-07-11)
depends_on:
  - THM-708   # the point this isolates
related:
  - THM-707 (kps), opus-S226 census, LEM-019 (dyadic mechanism)
---
# THM-709 — the doubling-tight locus is a singleton
Exact table (q < 500): e=7: M=1/11 (covering flips); e=8: 2/23; e=9: 2/23; e=10: 2/27;
e=11: 2/25; **e=12: 1/14 TIGHT**; e=13: 2/27. Double-doublings (11,12),(12,13),(10,12):
2/23, 2/25, 2/27 — all escape. Files: 04-computation/lrc14_doubling_locus_macmini_S65cont33.py (+ .out).

## Addendum — the singleton is the m=2 slice of Goddyn–Wong's own criterion (klein-2026-07-11-S253, literature merge)

Goddyn–Wong [39 in the Perarnau–Serra survey, arXiv:2409.20160 §4] prove the general
acceleration criterion: **replacing r ∈ [n] by mr preserves tightness iff r has a common
factor with every element of the interval [(n+1−r), m(n+1−r)−1]**. At n = 13, m = 2 the
interval is [14−r, 2(14−r)−1]; checking r = 6..13: only **r = 12** (interval [2,3]; 12
divisible by both) passes — this theorem's singleton is exactly the m = 2 slice of their
criterion, now citable. (For r = 12, m = 3 the interval [2,5] fails at 5 — no deeper
accelerations, consistent with the double-doubling escapes verified here.) This file's
net-new content beyond GW remains: the exact ESCAPE VALUES (1/11, 2/23, 2/23, 2/27, 2/25,
2/27) for every failing acceleration, and the covering-flip observation at e = 7.
Survey context: the full characterization of tight instances is OPEN (converse of GW Thm 12
false in general; Jacobsthal-linked tight families with v_n = 2n − Θ(log n) exist, Erdős;
Pomerance: n < v_max < 2n − c·log²n ⟹ NOT tight — a citable scope window).
