---
id: THM-716
title: The k=9 base extremal is FINITE-DIMENSIONAL — writing J = μ(7−μ) − Var (μ = E[N_empty]), the lower-J frontier is exactly max-Var-per-μ; three exact structural facts collapse the infinite family space: (i) DILATION-INVARIANCE (N-distribution of g·E′ equals that of E′, verified exact — likely THM-531 in this coordinate) ⟹ primitives only; (ii) the SECOND synchronization pole (mod-7-aligned {c+7m}) is DISPATCHED by the mean term, J ∈ [5.63, 5.77] > minimizer 4465/882 = 5.0624 > threshold 432/91 = 4.7473; (iii) near the low-μ optimum the frontier is the ONE-PARAMETER consec-shift family, minimized uniquely at cshift1 = {1..9}, with 0/16 distinctness-surviving local deformations dipping below — so THM-714/711's inf-J conjecture is an isolated 1-D minimization, not an infinite search
status: VERIFIED-EXACT (the identity J = μ(7−μ) − Var is algebra; the frontier scatter, dilation-invariance, mod-7 dispatch, consec-shift J-curve, and deformation-stability all exact-rational, this session; 400-family random cloud + 34 structured families + 16 deformations, zero counterexamples to the frontier structure). This does NOT prove inf J = 4465/882 (the global-over-all-primitives claim remains THM-714/711's open conjecture), but REDUCES it to: [primitives, by (i)] ∩ [dispatch the high-Var/high-μ poles, by (ii)] ∩ [1-parameter consec-shift minimization near the low-μ corner, by (iii)]. The mean term (THM-715) is the load-bearing dispatcher of the variance poles.
source: mac-mini-2026-07-09-S65 (cont.38, 2026-07-11)
depends_on:
  - THM-711 (the hit-empty product form J = E[N(7−N)])
  - THM-715 (variance-synchronization; the mean-term dispatch)
  - THM-531 (scale-invariance, = dilation-invariance (i))
related:
  - THM-714 (the k=8 cubic analog — same finite-dimensionalization expected), Giri-Kravitz accumulation hierarchy (klein-S253; the rigorous decorrelation-limit form)
---
# THM-716 — the k=9 base extremal is finite-dimensional
J = E[N(7−N)] = 7μ − E[N²] = 7μ − (Var + μ²) = **μ(7−μ) − Var**. Level curves J = c are the
downward parabolas Var = μ(7−μ) − c; the lower-J frontier is the upper envelope of the (μ, Var)
cloud (max Var per μ). Exact frontier (k=9, 434 families):
| μ-bin | Var_max | J | family |
|---|---|---|---|
| 1.4 | 2.9670 | **5.0624** | cshift1 {1..9} (global min) |
| 1.5 | 2.7539 | 5.3916 | cshift2 |
| 1.6 | 3.0124 | 5.7670 | mod-7 (higher J: μ dispatches) |
Consec-shift J-curve: 5.199, **5.062**, 5.392, 5.607, 5.573, … (unique min at shift 1).
Dilation g·E′ ≡ E′ in N-distribution (g = 2,3,5 exact). No deformation of {1..9} beats 4465/882
(0/16). Files: lrc14_{pareto_cloud, frontier_1param, gcd7_reduction}_macmini_S65cont38.py (+ .out).

## Addendum (cont.38 slice 3): consec is an ISOLATED SADDLE, both attack vectors fail
Adversarial hill-climbs (600 moves × 14 runs, wide box):
- **LOW-μ vector**: families with μ down to 1.4325 (BELOW consec's 1.446) all have J ≥ 6.15 —
  lowering μ requires spread coverage that drops Var by MORE than μ(7−μ) drops. The
  "low-μ AND high-Var" corner is EMPTY.
- **DIRECT-J vector**: best hill-climb found J = 5.3916 (= cshift2, and its dilate 2·cshift2
  = {4,6,…,20} found independently — re-confirming dilation-invariance). Nothing beat 5.0624.
Global best over all runs: exactly cshift1 = {1..9}. Consec does NOT minimize μ and does NOT
maximize Var globally — it is the isolated saddle where the μ(7−μ)−Var tradeoff is optimized.
This is why the extremal resisted every monotone/compression argument (cont.29): it is a
tradeoff optimum, not an endpoint. Files: lrc14_lowmu_highvar_macmini_S65cont38.py (+ .out).

## Addendum (cont.42): the compact+tail PROOF STRUCTURE, MISTAKE-138-robust (0-forced physical framing)
In the physical framing (0 ∈ E forced, THM-527-A; consec = {0..8}), the k=9 base J ≥ 432/91 has
a complete proof OUTLINE with a bounded finite check:
- **Compact**: min J over primitive 9-sets {0}∪S is at **{0..8}, J = 1019/196 = 5.1990** —
  exhaustive over diameters d = 8..20 (≤ 19448 families each, feasible), min at the smallest d.
- **Tail**: min-J RISES with diameter toward **J_iid = 8.4560** (decorrelated, exact occupancy
  I-E over 8 iid nonzero phases); sampled min-J at d = 24..45 is ≥ 5.97, and — heeding MISTAKE-138
  (test the extremal family, not a box) — **15 structured extremal candidates** (block+far {0..7,d}
  to d=300, 2AP+far, 3AP+far, mod-7 pole {1,8..57}, block+2far) ALL have J ≥ 5.64 > 5.199, ZERO
  below the compact min. The mod-7 pole (klein's BUNCH-max) is a HIGH-J family (5.767), no threat
  to min-J.
- **Threshold 432/91 = 4.7473 cleared at EVERY diameter, margin ≥ +0.4517** (bigger than the
  {1..9}-framing margin +0.315 — the physical base is easier).
So the k=9 base = [finite exhaustive check, primitives, d ≤ D0 ≈ 20] + [decorrelation tail via the
THM-710 eigen-transfer rate], the tail never binding. This turns THM-711/717's verified conjecture
into a proof structure whose only gaps are the finite check + a rigorous tail rate.
Files: lrc14_crossover_D0 + lrc14_minJ_extremal_check_macmini_S65cont42 (+ outs).

## Addendum (cont.43): SHARPENING the cont.42 tail — it plateaus at the MULTI-SCALE limit, not J_iid; + k=8 parity
**Correction of my own cont.42 imprecision.** The min-J tail does NOT rise to J_iid = 8.456.
The min-J tail families (block+far {0..7, d}) PLATEAU at J ≈ 5.69 as d → ∞ (measured: J = 5.807,
5.790, 5.679, 5.671, 5.697, 5.689 at d = 20…640; |J − J_iid|·d GROWS linearly, so |J − J_iid| does
NOT decay). **Mechanism:** a wide 9-family decomposes into COMPACT CLUSTERS (the block {0..7} stays
tightly correlated) + decorrelated far elements, so J → the TWO-SCALE / MULTI-SCALE limit (klein
THM-687/688), NOT the fully-independent J_iid. The correct tail argument is thus:
`J(wide) ≥ J(two-scale limit of its densest compact cluster) ≥ compact-min` — the far elements only
RAISE J (decorrelation between clusters adds positive mass). So the [compact exhaustive check]
handles both compact families AND the compact clusters of wide families; the cluster-decorrelation
(THM-687/688) supplies the ≥, and the plateau 5.69 > threshold 4.747 confirms it. This is CLEANER
than "rises to J_iid" — the tail is a multi-scale reduction to the same compact check.
**k=8 PARITY:** the k=8 deg-3-majorant bound (an UPPER bound, must stay ≤ cap₉ = 0.4943) is
MAXIMIZED at consec {0..7} (0.4380, margin +0.0563), DECREASES with diameter, and no structured
large-d candidate (block+far to d=200, 2AP+far, mod-7, block+2far) exceeds 0.36. So BOTH density-
side base checks have the compact+tail structure with the same multi-scale tail reduction.
Files: lrc14_k8_crossover_rate_macmini_S65cont43 (+ out).
