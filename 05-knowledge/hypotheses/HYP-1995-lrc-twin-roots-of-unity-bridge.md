---
id: HYP-1995
status: OPEN
source: oracle-2026-06-01-S522
related:
  - HYP-1987
  - HYP-1994
  - HYP-1981
  - THM-369
---

# HYP-1995: the regular polygon (roots of unity) is the shared extremal target of LRC and twin-Goldbach

**Frame.** Model the half-turn runner clock as points on the unit circle; the
tournament `i->j iff frac(θ_j-θ_i)∈(0,1/2)`. Equally spaced points = regular
`m`-gon = `m`-th roots of unity.

**VERIFIED (`tournaments_as_polygons_lrc_s522.py`, m=3..7):**
1. Regular `m`-gon half-turn tournament = rotational circulant `R_m` (connection
   `{1..(m-1)/2}`), for m=3,5,7.
2. Regular **even**-gon = degenerate WALL (antipodal pairs at frac=1/2 are ties).
3. `R_7` (regular heptagon, H=175) **is a reachable LRC source class** in the
   S520 menu; Paley `P_7` (H=189) is **not**. The menu selects the geometric
   (evenly spaced) regular polygon, not the arithmetic (QR/Paley) one.

**OPEN claims to test:**
- **(A) R_m always reachable, unique, geometric.** For every odd `m`, `R_m` is a
  reachable LRC source class, is the unique all-regular class in the menu, and is
  `≠ Paley_m` (m≡3 mod4). TEST: extend the S520 source-menu computation to
  m=9,11 (n=10,12) and check `canon(R_m)`/`canon(P_m)` membership + H.
- **(B) Even/odd polygon parity drives the menu collapse.** The S520 menu equals
  `A000568(n-1)/2` for n≤6 then collapses at n=7 (=m=6, the hexagon, an even-gon
  wall). Conjecture: the clean `/2` is the complement (= antipodal map
  `θ↦θ+1/2`) pairing, and the collapse is antipodal degeneracy acquiring fixed
  walls at the even-gon. TEST: count menu classes fixed vs paired by `T↔T^op`
  across n; correlate the `/2`→collapse with even-gon walls.
- **(C) Shared template.** Both LRC and twin-Goldbach are "an additive orbit must
  meet an arithmetically-thinned, rotation-symmetric root-of-unity target on a
  circle." twin-Goldbach = the finite/solved instance on the mod-6 hexagon
  (K=`{6k±1 prime}`, HYP-1994); LRC = the open instance on the real track
  (regular-polygon menu, HYP-1987). Seek an LRC "necklace reduction" analogous to
  twin-Goldbach's reduction to binary Goldbach on `K`: quotient the non-regular
  menu (the "thickening") to leave a covering statement on the polygon skeleton
  `{R_m}`.

**Why it matters.** It identifies the *extremal vertex* of the HYP-1987 target as
the roots of unity, and pins that LRC reachability is decided by circle geometry
(arc) rather than field structure (Paley/QR) — sharpening where a proof should
look. Reflection: `07-reflections/tournaments-as-polygons-lrc-twin-bridge-s522.md`.
