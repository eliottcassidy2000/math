---
id: THM-530
title: The G_P-intersection uniform floor for LRC(14) S3 — the via-max ρ*_{2/7} target is REFUTED (exact zeros on admissible (P,E), all reconstructed S still lonely), and the CORRECT global-witness density ρ*_{1/7}=meas(G_P∩{maxgap>1/7}) admits a clean two-branch floor: k≤7 PROVED unconditional (pigeonhole μ_{1/7}=1 ⟹ ρ*=meas(G_P)≥m_P), k≥8 a union bound ρ*_{1/7}≥meas(G_P)+μ_{1/7}(E)−1≥1891/5880 contingent only on a WEAK three-distance floor μ_{1/7}(E)≥1−min_P meas(G_P) (≥0.32 slack)
status: MIXED. REFUTED (exact, with realizability): the dispatched via-max target "ρ*_{2/7}(P,E)=meas(G_P∩{maxgap>2/7})≥c₀>0 uniformly" is FALSE — exact ZEROS on an explicit admissible family E={0,2,…,k}, P={1,2,3,…}, k≥7 (e.g. k=7 P={1,2,3,6,12,13} E={0,2,3,4,5,6,8}: meas(G_P)=515/1092, meas(Good^{2/7})=13/35, ρ*_{2/7}=0 EXACT); BUT all reconstructed covering sets S=P∪{Vmax−e} are LONELY (M(S)>1/14, 8/8 Vmax-sweep, gcd=1), so LRC(14) is NOT threatened — the via-max criterion is merely not necessary. PROVED: the admissible-|P| floor m_P=min_{|P|≤10}meas(G_P)=14249/252252 (EXACT); the k≤7 branch (pigeonhole μ_{1/7}=1 ⟹ ρ*_{1/7}=meas(G_P)≥m_P, k≤6 unconditional + k=7 a.e., exhaustive 18564 shapes spread≤18). VERIFIED (exact, large search, ZERO failures): the k≥8 union bound ρ*_{1/7}≥meas(G_P)+μ_{1/7}(E)−1≥1891/5880≈0.322; consecutive minimizes μ_{1/7} at every k≥8 (bounded+large spread); the global-witness quasi-independence ratio R_{1/7}≥67053/84241≈0.796 (vs via-max 0.068); no ρ*_{1/7}=0 anywhere. OPEN: the k≥8 union-bound branch is contingent on a WEAK universal three-distance floor μ_{1/7}(E)≥1−min_P meas(G_P) (the 1/7 spread bound, slack ≥0.32) — strictly easier than B(k); and the upstream global-witness reformulation (gap>1/7 ⟹ M≥1/14) is THM-527/kps-S4 (assumed, consistency-checked here). LRC(14) NOT proved.
source: kind-pasteur-2026-06-18-S5 (angle: gp-intersection-uniform)
depends_on:
  - THM-527   # lonely-density reformulation; global-witness (gap>1/7) variant
  - THM-528   # the four-window G_P coupling (Angle B); this sharpens/corrects the threshold
  - THM-523   # covering-set reduction; meas(G_P) floor
related:
  - OPEN-Q-108
  - HYP-2586   # the universal three-distance floor μ_min(k) — here the 1/7 variant
  - HYP-2590   # four-window scale decoupling
external: Lonely Runner Conjecture (13 speeds = LRC(14), first open case). Steinhaus three-gap.
---

# THM-530 — The G_P-intersection uniform floor (global-witness threshold) and the via-max refutation

## Setup (THM-527 notation)
Admissible (P,E): `P ⊆ {1,…,13}`, `E` with `0∈E`, `k=|E|≥3`, `|P|+|E|=13` (so `|P|≤10`).
`G_P = {x: ‖px‖≥1/14 ∀p∈P}`. For a threshold `θ`, `Good_E^θ={x: maxgap{frac(e_i x)}>θ}`,
`ρ*_θ(P,E)=meas(G_P∩Good_E^θ)`, `μ_θ(E)=meas(Good_E^θ)`.

THM-527 part A uses the **via-max** criterion `θ=2/7` (anchored at `v=Vmax`); the kps-S4
convergence note records the **global-witness** criterion `θ=1/7` (`good x ⟺ x∈G_P` and the
cluster phases leave a circular gap `>1/7`), of which `>2/7` is only a sufficient special case.

## A. The admissible-|P| floor (PROVED, EXACT)
`m_P := min_{|P|≤10} meas(G_P) = 14249/252252 ≈ 0.0565 > 0`, attained at `|P|=10` (k=3),
`P={1,2,3,5,7,8,9,11,12,13}`. (Finite exact computation; replaces the loose `7/858`.)
So `meas(G_P) ≥ m_P` for **every** admissible `P`. The min over each size `|P|=psz` is
monotone decreasing in `psz`: `6/7, 66/91, 55/91, 1979/4004, 2243/5880, 3029/10780,
45107/229320, 2479/17640, 10601/114660, 14249/252252` for `psz=1..10`.

## B. REFUTATION of the via-max (θ=2/7) uniform floor (EXACT)
`ρ*_{2/7}(P,E)=0` exactly on the admissible perforated near-AP family
`E={0,2,3,…,k}` (drop the “1”), `P={1,2,3,…}`, `k≥7`. Witnesses:

| k | P | E | meas(G_P) | μ_{2/7}(E) | ρ*_{2/7} |
|---|---|---|-----------|-----------|----------|
| 7 | {1,2,3,6,12,13} | {0,2,3,4,5,6,8} | 515/1092 | 13/35 | **0** |
| 8 | {1,2,3,12,13} | {0,2,3,4,5,6,7,8} | 177/364·… | 0.336 | **0** |
| 9 | {1,2,3,13} | {0,…,9}∖{1} | 0.571 | 0.287 | **0** |
| 10| {1,2,3} | {0,…,10}∖{1} | 0.690 | 0.247 | **0** |

Both factors are bounded away from 0, yet the intersection is **empty** — the
anti-correlation pathology Angle B feared is REAL at `θ=2/7`. **Realizability:** every
covering 13-set `S=P∪{Vmax−e:e∈E}` reconstructed from these (8 values of `Vmax` each) is
**lonely** (`M(S)>1/14`, exact; `gcd=1`). So the via-max zeros are **not** LRC(14)
counterexamples — the via-max criterion is sufficient but **not necessary**, and "ρ*_{2/7}≥c₀"
as a target is simply false. (This corrects the dispatch premise and refines THM-528, whose
"consecutive case proved" is at θ=2/7 and whose exception list is a symptom of the same
non-necessity.)

## C. The correct object: global-witness ρ*_{1/7}, two-branch floor

**Branch k≤7 (|P|≥6): PROVED unconditional.** Pigeonhole on `1/7`-gaps: `k≤6` points have
`k` circular gaps summing to 1, so `maxgap > 1/7` **always** (6·(1/7)<1) ⟹ `μ_{1/7}(E)=1`
for all `x`, all `E`. `k=7`: `maxgap≤1/7` forces all 7 gaps `=1/7` exactly (the shifted
`1/7`-grid), a measure-zero set of `x` ⟹ `μ_{1/7}=1` a.e. (exhaustive: 18564 shapes
spread≤18, min `μ_{1/7}=1`). Hence `ρ*_{1/7}=meas(G_P)≥m_P` for `k≤7`.

**Branch k≥8 (|P|≤5): union bound.** Inclusion–exclusion gives
`ρ*_{1/7}(P,E) ≥ meas(G_P) + μ_{1/7}(E) − 1`. With `μ_{1/7}(consecutive_k)` exact
(`691/735, 247/294, 38/49, 1381/2205, 13823/24255, 477/1078` for `k=8..13`) and
`min_P meas(G_P)` exact, the per-pair minimum is
`min_{k≥8,P}[meas(G_P)+μ_{1/7}(E)−1] = 1891/5880 ≈ 0.322 > 0` (at `k=8`,
`P={1,5,7,8,9}`). **Required:** `μ_{1/7}(E) ≥ thr_k := 1 − min_P meas(G_P)`
(`thr_8=3637/5880≈0.619`, decreasing to `thr_13=0`); the searched `μ_{1/7}` exceeds `thr_k`
by `≥0.32` at every `k` (bounded **and** large spread; consecutive is the minimizer).

**Combined floor candidate** `c₀' = min( m_P, 1891/5880 ) = m_P = 14249/252252` (the binding
value is the k=3 end of Branch 1, where Good=whole circle and `ρ*=meas(G_P)`).

## D. Quasi-independence (VERIFIED) — the right threshold makes it benign
For `θ=1/7`, `R_{1/7}=ρ*_{1/7}/(meas(G_P)·μ_{1/7})` satisfies `R_{1/7}≥67053/84241≈0.796`
over all consecutive cases, and `R_{1/7}=1` for all `k≤7`. (Contrast the via-max
`R_{2/7}≥0.068`.) No `ρ*_{1/7}=0` over structured + 8000 random bounded-spread pairs.

## E. Honest remaining gap
Branch k≥8 is contingent on the **weak universal three-distance floor**
`μ_{1/7}(E) ≥ thr_k` (the `1/7` spread bound) — strictly easier than the `2/7` lemma B(k)
because of the `≥0.32` slack and because for `θ=1/7` "large spread RAISES μ" is even more
favourable (`μ_{1/7}→` the iid value, well above `thr_k`). Also relied upon (upstream,
assumed): the global-witness reformulation `gap>1/7 ⟹ M≥1/14` (THM-527/kps-S4); a direct
consistency check found **no** non-lonely reconstructed `S`. **LRC(14) remains open**, but
the G_P coupling is now established to be a non-obstruction at the *correct* threshold, with
the residual reduced to a weak `1/7`-floor with large margin.

## Scripts (all exact)
- `04-computation/lrc14_Bk_gp-intersection-uniform_kps-S5-wf.py` (main: STEP0–5)
- `04-computation/lrc14_Bk_gp_unionbound_kps-S5-wf.py` (k≥8 union floor)
- `04-computation/lrc14_Bk_gp_slack_kps-S5-wf.py` (thr_k and slack)
- `04-computation/lrc14_Bk_largespread_clean_kps-S5-wf.py` (μ_{1/7} at all spreads)
- `04-computation/lrc14_Bk_k7edge_kps-S5-wf.py` (k=7 pigeonhole edge, exhaustive)
- `04-computation/lrc14_Bk_mu17_isolated_kps-S5-wf.py` (isolated μ_{1/7}, anti-corruption)
- `04-computation/lrc14_Bk_globalwitness_valid_kps-S5-wf.py` (loneliness consistency)
- outputs in `05-knowledge/results/lrc14_Bk_*_kps-S5-wf.out`

## Main-loop confirmation + the FINAL-LEMMA framing (kind-pasteur-2026-06-18-S5, post-workflow)
Independently re-verified with an exact `mu_theta` engine (order-cell + gap=θ breakpoints): all anchors
reproduce (`μ_{2/7}`: consec₄=19/21, consec₆=4/7; `μ_{1/7}`: consec₇=1, consec₈=691/735, consec₁₃=477/1078);
`thr_8=3637/5880`; and **k=8 EXHAUSTIVE over all 8-subsets of {0..14}: consecutive {0..7} is the unique
minimizer (μ_{1/7}=691/735≈0.940), 0/3432 below thr_8** (margin 0.32). So THM-530's two-branch floor is
confirmed; the combined floor is **c₀'=m_P=14249/252252** (NOT the k≥8 piece 1891/5880).

> **THE STATE OF LRC(14): reduced to ONE lemma.** k≤7 PROVED; k≥8 needs only the **1/7-SPREAD BOUND**
> `μ_{1/7}(E) ≥ thr_k` (8≤k≤12), with consecutive minimizing and ≥0.32 slack (binding at k=8;
> thr_13=0 trivial). This is *strictly easier* than the 2/7 lemma (which is FALSE — μ_{2/7} has no floor),
> survived the same adversarial large-spread descent that crushed μ_{2/7}, and is the single open
> analytic gap (plus the upstream finite-Vmax / integer-vs-real glue, THM-527-A). → HYP-2600.

## CANON CORRECTIONS surfaced this session (VERIFIED exact)
1. **μ_{2/7} has NO uniform floor (corrects THM-527-C / my own kps-S5 2/7 work).** Exact k=13 witnesses
   `μ_{2/7}(E)<1/14`, e.g. `E=[0,12,17,20,24,26,27,28,32,36,37,47,60]` → `3736/85785≈0.0436`. So the
   2/7 bounded-spread "floor ~0.11, spread O(k)≤30" estimate is REFUTED; the 2/7 measure descends with
   spread without a useful floor. This is HARMLESS — μ_{2/7} (via-max criterion) is sufficient-not-
   necessary; the LRC object is μ_{1/7}, where consecutive minimizes and the floor is genuine.
2. **THM-527-C value fix:** consecutive `μ_{2/7}(consec_7)=83/210` (THM-527-F's 0.395), NOT 13/35 —
   13/35 is the PERFORATED minimizer `μ_{2/7}({0,2,3,4,5,6,8})`. (The two were conflated in THM-527-C.)
3. **THM-528 threshold (FLAGGED for mac-mini, not unilaterally changed):** THM-528's "G_P anti-correlation
   impossible / not the obstruction" holds at the 2/7 threshold where it is FALSE — exact admissible
   witness `P={1,2,3,6,12,13}, E={0,2,3,4,5,6,8}` gives `ρ*_{2/7}=0` (both factors >0, disjoint). The
   four-window LOWER bounds are correct as pure-μ statements but give no positive G_P-intersection floor
   at 2/7. RESOLUTION (coordinate): re-point THM-528 at θ=1/7, where the same witness gives ρ*_{1/7}>0 and
   R_{1/7}≥0.796. See message to mac-mini.
