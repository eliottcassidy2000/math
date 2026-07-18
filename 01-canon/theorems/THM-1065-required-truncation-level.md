---
id: THM-1065
title: THE REQUIRED TRUNCATION LEVEL, MEASURED — on the THM-1060 primitive family class the first certifying odd Bonferroni level is 11, 9, 9, 7 at min speeds 56, 248, 808, 2648: the level DECREASES as the family grows, level 5 provably NEVER suffices (BONF5 → ≈ −0.16), and level 13 is exact by construction; the per-family cost of full inclusion–exclusion is only 8191 subset-measures (seconds), so the obstruction is NOT computational per family — it is that infinitely many families demand a UNIFORM argument, which extends THM-1026's ledger from five slots to seven (and to eleven in the correlated corner)
status: measured exactly at four scales (all S_k for k = 1..13 by subset DFS; truth verified independently by exact uncovered measure); the monotone level/scale relation is a measured trend across four points, not a proved law
source: opus-2026-07-17-S366 (owner: work the stronger tool, higher level truncation)
depends_on: [THM-1060 (the infinite primitive failure class this measures), THM-1026 (the five-slot ledger this extends), THM-1050/1055 (why census is unavailable)]
scripts: 04-computation/higher_truncation_opus_S366.py -> 05-knowledge/results/higher_truncation_opus_S366.out
---

# THM-1065 — how much stronger a tool is needed

> **ROUTE REFUTED (opus-S367), see THM-1070.** The closing proposal below -- that the same sawtooth/containment technique is the natural candidate for each remaining slot -- is WRONG. Both new B7 slots admit a VALID bound and neither is usable: the S6 containment floor is tight only at k=2 and loses ~5x per additional speed (exact/floor 3.5, 24.5, 114, 200, 2101 at k=2..6), and the iterated-fragmentation S7 upper bound is ~1190x too loose. The technique filled the k=2 slot precisely because k=2 is where it is sharp. The ledger needs a bound using the JOINT alignment of the k combs -- a k-fold analogue of THM-965's folded identity -- not another iteration of a per-step worst case.

THM-1060 showed BONF5 fails on an infinite primitive family class while
the truth stays ≈ +0.115. Since full inclusion–exclusion at level 13 is
exact, *some* odd truncation certifies. This file measures which.

| L | min speed | BONF5 | BONF7 | BONF9 | BONF11 | truth | first positive |
|---|---|---|---|---|---|---|---|
| 7 | 56 | −1.1411 | −1.0553 | −0.2025 | +0.1032 | 0.1207 | **11** |
| 31 | 248 | −0.4074 | −0.1839 | +0.0468 | +0.1164 | 0.1204 | **9** |
| 101 | 808 | −0.2368 | −0.0030 | +0.0908 | +0.1130 | 0.1142 | **9** |
| 331 | 2648 | −0.1975 | +0.0444 | +0.1054 | +0.1136 | 0.1140 | **7** |

**Three readings.**

1. **The required level decreases as the family grows** — 11, 9, 9, 7.
   Smaller (more strongly correlated) families need deeper truncation.
   The trend is measured across four points, not proved.

2. **Level 5 is provably insufficient**, not merely inadequate here:
   THM-1060 established BONF5 → ≈ −0.16 along this class, so no amount of
   scale rescues it. Level 7 appears to suffice for the bulk, with the
   correlated corner reaching 11.

3. **Per-family cost is negligible.** All 8191 subset-measures ran in
   seconds; empty-intersection pruning bought nothing (these combs are
   dense enough that intersections rarely vanish), so full
   inclusion–exclusion is simply cheap at n = 13. **Deciding any single
   family exactly is easy.**

**What this locates.** The obstruction is not computational per family
and not a missing region of family space — it is that infinitely many
families require a UNIFORM argument. Concretely, THM-1026's ledger

> B₅ = 1 − S₁ + S₂ − S₃ + S₄ − S₅, needing lower bounds on S₂, S₄ and
> upper bounds on S₃, S₅

must extend to **B₇** (adding a lower bound on S₆ and an upper bound on
S₇) and, for the correlated corner, to **B₁₁**. My S₂ floor (THM-1025)
fills one slot of five; the extended ledger has seven to eleven, and the
same sawtooth/containment technique is the natural candidate for each —
S₄ already has its separated-regime floor from S360.

**The honest shape of the frontier.** LRC(14) in this program now needs:
a uniform level-7 (or level-11) certificate, built slot by slot the way
S₂ was. That is a longer ledger than the level-5 one, but it is the same
kind of object, and every slot filled so far was filled by containment
and counting rather than by new analytic machinery.
