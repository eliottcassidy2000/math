---
id: THM-747
title: The vertices-as-runners phase sum S(W), triangulated -- by-vertex pair-sector collapse (per-end 516 -> 98.1, with an exact-zero vertex), exact periodicity (Q = 97020 / 8820; period max |S| = 71.23 / 62.93, a LOOKUP), and the tower verdict (integer time pairs mirror phases EQUALLY -- level-1 symmetry gives no level-2 cancellation; the observed 10-100x is phase equidistribution, resolved exactly by the scan); the tail lane's asymptotic coefficient is now a finite exact object
status: PROVED (the by-vertex identity and periodicity) + COMPUTED-EXACT (full-period scans, exact confirmation at argmax) + HONEST (standing W0 = 339/513 unchanged; the scan form applies for W >= Wz = 924/588; the deliverable is the asymptotic coefficient and the terminal map)
source: opus-2026-07-14-S282 (owner prompt: bound S(W) with the perspective frame; understand the critical remaining unknowns; triangulate the structure)
depends_on:
  - THM-746 (defines S(W), T(W), Z(W)); THM-745 (the exact first-order content); THM-742/743 (assembly)
related:
  - the S276 delta-census (the by-vertex couplings collapse to pair-difference form); klein S294-S301
  - NB two THM-744s exist in canon (opus F-telescoping, klein shadow-gap) -- ID housekeeping flagged
---

# THM-747 -- the phase sum triangulated

## LENS 1: by-vertex (the pair sector, one level down)

Regrouping the 344 / 308 ends by their 57 / 54 distinct vertices (integer time W: theta1 =
1 - {u W}, theta2 = {u W}):

>  S(W) = D0 + Sum_v b_v {u_v W},   D0 = 33.6429 / 39.2143 (exact deterministic mass),
>  |S(W)| <= |S_bar| + B1/2 = 98.14 / 93.29     [per-end S1 was 516.35 / 357.70 -- 5.3x / 3.8x]

The couplings collapse in the pair sector: at shared (swap) vertices b_v carries the
difference structure, and at u = 1/6 (shape 1) TWENTY-TWO ends sum to b_v = 0 EXACTLY --
total cancellation at a single vertex.  43/57 vertices are shared (multi-end).

## LENS 2: the runner lens -- exact periodicity (the promised lookup)

{u_v W} has period den(u_v) in the integer time W, so S(W) is periodic mod
Q = lcm(vertex denominators) = 97020 / 8820.  Full-period float scan + exact Fraction
confirmation at the argmax:

>  **max_W |S(W)| = 71.2259 (at W = 50754) / 62.9254 (at W = 5519)**;  max|T(W)| = 2636 / 1459.

The triangulation: per-end 516 > by-vertex 98.1 > EXACT 71.2 (measured typical 2.6-49).
The level-2 bound is no longer an estimate: it is a finite exact computation.

## LENS 3: the tower verdict

Integer time makes mirrored phases EQUAL ({(1-u)W} = {-u W}), not opposite: the level-1
mirror symmetry is EXHAUSTED at level 2 (verified: only 7/57 mirror pairs have equal b_v --
the coupling multisets differ because start/end types swap).  The typical smallness of S(W)
(2.6-49 vs max 71 vs bound 98) is genuine PHASE EQUIDISTRIBUTION across the 57 vertex-runners
-- and at level 2, unlike level 1, it can be resolved EXACTLY by the finite scan.  The tower
alternates: constrained level (heights fenced by the origin band, mirror cancels) -> free
level (phases unconstrained, periodicity resolves).  No level-3 is needed: the scan closes it.

## Assembly (honest)

The scan-based form applies for W >= Wz = 2/min(du) = 924 / 588 (below Wz some segments have
no crossings and Z(W) != 0; THM-746's charge covers that range).  For W >= Wz: C1 = 2.112 /
8.357, C2 = C2(743) + 2 S_max + T_max/2 = 5158 / 2832, and the bound sits ~9x below Area --
so the STANDING W0 = 339/513 (min-form with THM-743) is unchanged; what is nailed is the
ASYMPTOTIC error coefficient of the tail family:

>  limsup_W W |L - Area| <= 2(#comp + |Xi_sv|) = 2.112 / 8.357,  with the 1/W correction
>  bounded exactly by S_max + ... -- every constant a rational or a finite lookup.

## The critical remaining unknowns (named)

U1. **The 743 pots** (vertex capping, one-signed, Sum delta = 2172/1116; quadratic wedge
    remainders) now dominate C2 at mid-W.  Both are ALSO periodic mod Q: their exact
    period-extremes are the same kind of lookup (unscanned).  With them scanned, the tail
    lane's entire error budget becomes a finite exact object -- NO analytic unknown remains
    in this lane.
U2. **The compact core** (bounded-Vmax bodies): outside this program; kps's exact-certificate
    territory.
U3. **The general multi-speed equidistribution** (klein-S300's capstone): the fleet-level
    residual.  This lane contributes exact structure (strand identity, pairing theorem,
    periodicity) but no bridge to it -- the honest boundary.

## Files

04-computation/lrc14_phase_sum_triangulation_thm747_opus_S282.py (+.out): by-vertex regroup +
identity, pair-collapse tables, mirror census, full-period scans with exact argmax
confirmation, assembly, unknowns.
