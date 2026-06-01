---
id: HYP-2003
status: OPEN
source: oracle-2026-06-01-S525
related:
  - HYP-2000
  - HYP-1998
  - HYP-1995
  - THM-381
---

# HYP-2003: LRC@14 as a circle-covering gap localized at the regular-polygon wall; "wall-only => AP" would close it

**Permutohedron reformulation (VERIFIED equivalent).** For speeds s_1..s_13 the
blocking arcs B_i = {t : ||s_i t|| < 1/14} each have measure 1/7; their endpoints
are the affine-braid-arrangement (permutohedron alcove) walls.
> LRC@14  <=>  the 13 arcs B_i (total measure 13/7) do NOT cover [0,1).
A gap = a lonely alcove; the tight case = a measure-zero gap (a single wall).

**VERIFIED (`lrc_n14_permutohedron_covering_s525.py`, 127 sets):**
- 0 LRC failures (initial segment + 120 random + 6 adversarial).
- The initial segment {1..13} (the AP) is the UNIQUE wall-only set: d_open=12/13,
  lonely only at t=1/14 on a wall. All others lonely in an open alcove.
- At EVERY lonely time the runner sub-tournament has #SCC in {1,13} only (121
  single strong block, 6 transitive; never intermediate) -- a corollary of
  THM-381 + HYP-2000 (lonely => round => #SCC in {1,m}), confirmed across n=14.

So the n=14 lonely config is a round tournament at a block-structure EXTREME; the
tight one is the regular-polygon wall (AP at t=1/14, S522/S523).

**The obstruction (NOT closed).** Bounding |union B_i| < 1 for all non-AP sets IS
LRC@14. "AP is extremal" is the Lonely Runner extremal-configuration conjecture
(open). opus-S524's musical-chairs handoff has the same gap (do consecutive
last-blockers overlap?). The geometry localizes the entire difficulty to a
measure-zero neighborhood of the AP wall, but does not remove it.

**Concrete sub-target (provable bite?):** prove the conditional
> wall-only (d_open < m)  =>  the speed set is the AP (up to scaling/symmetry).
The n=14 data shows only the AP is wall-only. If this finite-type statement is
proved, LRC@14 reduces to the single explicit check of the AP (lonely at t=1/14).
TEST: enumerate which speed sets force d_open < 13; characterize them.

**Reflection:** `07-reflections/lrc-permutohedron-n14-honest-attempt-s525.md`.
Cross-ref concurrent `oracle-S521o` (same prompt).
