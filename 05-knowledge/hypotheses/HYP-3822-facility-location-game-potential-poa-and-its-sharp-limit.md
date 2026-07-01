---
id: HYP-3822
title: THE COVERING-MIN AS AN ADVERSARIAL FACILITY-LOCATION GAME, AND THE SHARP LIMIT OF THE POTENTIAL / PRICE-OF-ANARCHY IMPORT. The covering-min is a defender-attacker facility-location game on the circle: runners = facilities x_i(t)=v_i t mod 1, the observer = the client at 0, loneliness(t)=min_i||v_i t||, M(S)=max_t loneliness = the defender's minimax objective, covering-min = min over covering sets S of M(S). The lonely set L=(t : min_i||v_i t|| >= r) and inf meas = inf_S |L| = the FLOOR. POTENTIAL (Rosenthal congestion): coverage C(t)=#{i:||v_i t||<r}; A=int C=(n-1)*2r (1st moment); Phi=int C(C-1)=sum_{i!=j}|D_i cap D_j| (pairwise congestion); C_max=max C<=n-1. THEOREM |L| >= 1 - A + Phi/C_max (valid, beats the union bound 1-A). SHARP LIMIT: at the critical radius r=1/14 (A~1.86~2) the bound goes NEGATIVE, and the MOMENT LP (min m_0 given A,Phi,support) = EXACTLY 0 -- a same-moment coverage distribution with ZERO loneliness exists -- so the GLOBAL potential/PoA CANNOT force inf meas>0; the floor is a LOCAL/ARITHMETIC fact (HYP-3571 Gamma_0(N) congruence moment), not a global-moment one. This SHARPENS HYP-3817: reach for the ARITHMETIC/local moment; the global empirical moment is as blind to the atom as a transform. Blaschke(arXiv:2604.16750): facility dynamics = the circle maps B_a, involution I(z)=1/zbar = the fold iota, Arnold tongues T_{p/q} = the rational resonances (the local a/q piles), Herman ring = the Diophantine deep well. Kaczynski boundary functions: C(t) is a TAME boundary function (Baire-1, finite arc-union); the floor = its local continuity/tameness.
status: MIXED -- GAME FRAME + potential bound CONFIRMED (exact, computed n=14); the IMPOSSIBILITY (global moments can't force the floor: LP min m_0=0) CONFIRMED; the connections (Blaschke Arnold-tongue, Kaczynski boundary-function tameness) are structural analogies, exact where arithmetic, analogical where dynamical.
source: mac-mini-2026-07-01-S91
related:
  - HYP-3817   # fixed-point instruments: covering/moment not transform -- THIS sharpens it (global moment insufficient)
  - HYP-3571   # the 2nd-moment floor: the ARITHMETIC (Gamma_0(N) congruence) moment that DOES work
  - HYP-3785   # the LRC spectral gap: Fourier/Delsarte transforms blind to the atom (same wall)
  - HYP-3796   # S80 Kaczmarz/Christoffel/Blaschke (the same Blaschke/circle-map dynamics, POCS witness)
  - HYP-3797   # S79 loop-function dictionary (AGL(1) runner=dilation; the game's moves)
results:
  - 04-computation/facility_location_potential_poa_macmini_20260701.py
  - 04-computation/facility_location_moment_lp_macmini_20260701.py
  - 05-knowledge/results/facility_location_potential_poa_macmini_20260701.out
  - 05-knowledge/results/facility_location_moment_lp_macmini_20260701.out
---

# HYP-3822 -- the facility-location game, the potential import, and its sharp limit

Owner's directive: the covering-min is an adversarial facility-location game on the circle (runners =
facilities, the lonely observer = the point farthest from all, covering-min = the adversary's
min-over-configs value); import a potential-function / price-of-anarchy argument for `inf meas`.

## 1. The game (attacker-defender facility location)
Circle `R/Z`. Runners are **facilities** at `x_i(t)=v_i t mod 1`. The observer is the **client at 0**.
- **loneliness**`(t) = d(0,t) = min_i ||v_i t||` (distance to the nearest facility).
- `M(S) = max_t loneliness` = the worst client-service distance = the covering objective.
- **covering-min** `= min over covering sets S of M(S)` `= 14/183` at `n=14` (the DEFENDER minimizes the
  worst loneliness; constrained to covering sets = integer speeds meeting the AP/`phi(14)` structure).
- **lonely set** `L_S(r) = {t : min_i ||v_i t|| >= r}`; **inf meas** `= inf_S |L_S(r)|` = the FLOOR.

Attacker-defender: DEFENDER picks `S` to cover the observer (minimize loneliness); ATTACKER (nature)
picks the time `t` to maximize it. Value `= min_S max_t` loneliness `=` covering-min (a minimax).

## 2. The potential import (Rosenthal congestion) -- a valid bound
Coverage count `C(t) = #{i : ||v_i t|| < r}` (facilities within `r` of the observer). Two moments:
- `A = int C dt = sum_i |D_i| = (n-1)*2r` (each danger set `D_i={t:||v_i t||<r}` has measure `2r`).
- `Phi = int C(C-1) dt = sum_{i!=j} |D_i cap D_j|` = the **Rosenthal / congestion potential** (total
  pairwise overlap).

**THEOREM (potential lower bound on the lonely measure).** `|L| >= 1 - A + Phi / C_max`, `C_max=max C<=n-1`.
*Proof.* With `m_k=|{C=k}|`: `A=sum k m_k`, `|danger|=sum_{k>=1}m_k`, `Phi=sum k(k-1)m_k`. For `k>=2`,
`k(k-1)<=C_max(k-1)`, so `Phi<=C_max sum(k-1)m_k=C_max(A-|danger|)`, giving `|danger|<=A-Phi/C_max`, hence
`|L|=1-|danger| >= 1-A+Phi/C_max`. QED

This **beats the union bound** `|L|>=1-A` (overlap = wasted coverage = gaps elsewhere) and is **live in the
sub-critical regime** (computed, `S={1..12,182}`, `n=14`):

| `r` | `A` | union `1-A` | potential `1-A+Phi/C_max` | actual `|L|` |
|---|---|---|---|---|
| 0.0357 | 0.928 | 0.072 | **0.225** | 0.432 |
| 0.0400 | 1.040 | -0.040 (dead) | **0.133** | 0.373 |
| 0.0500 | 1.300 | -0.300 | -0.079 | 0.232 |
| **1/14=0.0714** | 1.857 | -0.857 | **-0.517 (dead)** | 0.024 |

So the congestion potential **extends** the union bound's range (union dies at `A=1`, `r=0.0385`; the
potential stays positive to `r~0.047`) -- but **dies before the critical radius** `r=1/14`.

## 3. The sharp limit -- the global potential CANNOT force the floor (moment LP)
Is the potential just loose, or is the **global moment** genuinely insufficient? Solve the moment LP:
`min m_0` s.t. `sum m_k=1, sum k m_k=A, sum k(k-1)m_k=Phi, m_k>=0` on support `{0..n-1}`.

**Result: `min m_0 = 0.000000`** at BOTH `r=1/14` and `r=14/183`. The witness distribution
`{C=1: 0.39, C=2: 0.59, C=12: 0.025}` has the SAME `(A,Phi)` as the actual coverage but **zero lonely
measure.** So **the global congestion moments provably cannot force `inf meas > 0`** at the critical radius
-- a valid coverage distribution with the same first two moments is never lonely. The floor is **not** a
global-moment fact; it is **local/arithmetic** (which `t` carries which coverage).

## 4. The sharpened principle (this is the point)
This is the **PoA-lens sharpening of HYP-3817**. That principle said: for a fixed-point/atom extremum, reach
for a covering or a **moment**, not a transform. This session shows the moment must be the **right** moment:
- the **global empirical** 2nd moment `Phi` is **insufficient** (LP min `= 0`) -- it is an averaging
  instrument, blind to the atom, exactly like the Fourier/Delsarte transforms hit the spectral gap (HYP-3785);
- the moment that **works** is the **arithmetic / congruence** 2nd moment (HYP-3571, `Gamma_0(N)`,
  set-independent `CV(N_R)^2`) -- it carries the LOCAL rational structure the global `Phi` averages over.

So: **"reach for a moment" -> "reach for the ARITHMETIC/local moment."** The global moment is as blind to
the pointwise atom as a transform; only the congruence moment (or a covering) sees it.

## 5. The analytic frame: Blaschke dynamics + Kaczynski boundary functions
- **Blaschke / Riemann sphere (arXiv:2604.16750):** the facility dynamics `t -> {v_i t}` are the linear case
  of the circle maps `B_a(z)=e^{2pi i t} z^{d+1}((z-c)/(1-cbar z))^d`; the involution `I(z)=1/zbar` **is the
  complement/fold `iota`** (the involution atlas, S88); the **Arnold tongues** `T_{p/q}` are the rational
  resonances = the LOCAL `a/q` congestion piles the global moment averages away; the **Herman-ring /
  Diophantine** regime (absent for `|a|<=2d+1`) is the deep well where the atom lives (`t*=n/Phi6`, bounded
  partial quotients = the S79/Herman rigidity). The floor lives at a Diophantine rotation number inside an
  Arnold tongue -- a LOCAL object, invisible to the endomorphism-regime average.
- **Kaczynski boundary functions (PhD 1967, radial/curvilinear limits on the disk):** the coverage `C(t)`
  and the lonely indicator are **tame boundary functions** -- Baire-1, finite unions of arcs (each `D_i` is
  a finite arc-union). The floor is a fact about this boundary function's **zero set** (a union of intervals
  with measure controlled by the arc structure), a LOCAL continuity/tameness fact, not a moment. The
  **Bagemihl ambiguous points** (where the boundary function jumps) = the **endpoints of the lonely arcs**;
  their finiteness (tameness) is what makes `|L|` well-defined and drives the Fourier decay `N/(pi m)` (S75).

## Honest scope
The game frame and the potential bound are exact (theorem + computation, `n=14`). The impossibility (`min
m_0 = 0`, global moments insufficient) is a solid LP fact -- a genuine LIMITATION of the PoA/potential import,
and the session's real content: it explains WHY the naive congestion argument fails and points precisely at
the arithmetic moment (HYP-3571) as the fix. The Blaschke/Kaczynski connections are exact where arithmetic
(the involution, the tameness) and structural analogies where dynamical (Arnold tongue = resonance).
