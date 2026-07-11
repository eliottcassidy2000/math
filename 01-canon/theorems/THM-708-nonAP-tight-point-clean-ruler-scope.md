---
id: THM-708
title: The M-tight locus contains a non-AP point — {1,…,11,13,24} has M = 1/14 EXACTLY (attained at t = 1/14; exact scan all p/q, q < 700) — resolving opus-S226's clean-ruler edge case structurally (tight families admit no clean ruler and need none: the family is non-covering, 14 ∤ every element, sieve-dispatched with equality); the mechanism is the DOUBLING rearrangement 12 → 24 = 2·12 (‖24t‖ = ‖2·(12t)‖ keeps the removed element's lattice constraint), and it is RARE: sibling rearrangements 11→25, 13→15, 6→20 all jump to M ∈ {1/11, 1/13, 2/23}
status: VERIFIED-EXACT (M = 1/14 as max over all rationals q < 700 — the attained-max of a piecewise-linear function; the doubling mechanism is the noted structure, its generality untested beyond the probe). CONSEQUENCES: (i) opus's near-AP clean-ruler census is 72/72 on the actual residual class (the lone exception is not residual: not covering AND tight); (ii) kps THM-707's clean-ruler supply is correctly scoped to residual/covering families — the tight locus (AP + at least this doubling point) is exactly where clean rulers must and do fail; (iii) any extremal-uniqueness claim "M = 1/14 ⟹ AP/dilate" is FALSE as stated — the correct locus includes doubling rearrangements (flag for THM-530-adjacent write-ups).
source: mac-mini-2026-07-09-S65 (cont.32, 2026-07-11)
depends_on:
  - THM-707 (kps)  # the clean-ruler reduction this scopes
related:
  - opus-S226 (the census + edge flag), LEM-019 (the dyadic mechanism), THM-530 (extremal claims to scope)
---
# THM-708 — the non-AP tight point and the clean-ruler scope
v = {1,…,11,13,24}: not covering (14 divides nothing) ⟹ t = 1/14 gives min‖v/14‖ = 1/14 (runners
1, 13 on the boundary; 24 ≡ 10 doubles position 10, 12 empty) — LONELY with equality. Exact sweep of
all p/q, q < 700: best clearance = 1/14 exactly ⟹ M = 1/14, tight, NON-AP. The doubling 12→24 is the
mechanism; siblings without the doubling relation escape to M > 1/14. Files:
04-computation/lrc14_tight_point_macmini_S65cont32.py (inline in session commands; probe table in log).
