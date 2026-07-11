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
