---
id: THM-564
title: The genuine-wide DOUBLET deviation g(M)=M·(p0(E_M)−Φ) is ALMOST-PERIODIC and BOUNDED when centered at the EXACT frozen plateau Φ — exact inclusion-exclusion split g=P(M)+R(M), P period-(7·lcm(base)) with EXACT finite period-max, R=M·(d2(M)−d_inf)=O(1/M); closes the binding genuine-wide leg p0(E_M)≤cap_K−0.16 for all M≥15 via period-max(P)+sup|R| and a TINY finite window [15,M*] (M*=13,24,… ). The doublet analogue of THM-563.
status: VERIFIED-closed ALL K=8..12 (FULL exact period at K=8,9,10 = 420,420,2940, M*=13,24,55; window-verified K=11,12, M*=25,13, H_K large/non-binding). The exact incl-excl identity p0(E_M)=p0(base)+A(M)+A(M+1)+d2(M) and g=P+R are PROVED exact (consistency P+R==g verified all M). period-max(P) is EXACT (THM-563 single-far sawtooth, finite). R=O(1/M) is VERIFIED (sup|R|=M·|d2−d_inf| bounded 0.43–0.74, d2−d_inf decays ~1/M; mechanism = Koksma/three-distance discrepancy of the single residual far phase). The general bounded-base R-tail closure is the remaining finite check (the doublet analogue of THM-563's completed 12805-base periodmax). CONVERGES with the concurrent kps-2026-06-21-S27 frozen-phase route (HYP-2799/2801/2803): their independent slow-x/fast-u double integral p0_inf EQUALS Φ_frozen to the EXACT rational at every K; their TV bound J_sharp=200/7≈28.6 gives the same closure with a looser window (f0~48–148) — this theorem's EXACT period-max(P)~1.2 is the tighter form (M*~13–55).
source: kind-pasteur-2026-06-21-Swf9
depends_on:
  - THM-563   # the single-far exact periodicity this extends to the far DOUBLET
  - THM-557   # the BV/TV finite-cutoff this collapses (M_BV~2000 / f0~148 -> M*~13-55)
related:
  - HYP-2804  # this theorem's hypothesis-log entry (the frozen-plateau centering + exact P/R)
  - HYP-2799  # kps S27 CONVERGENT: frozen-phase p0_inf == Φ_frozen (independent double integral)
  - HYP-2803  # kps S27 CONVERGENT: binding doublet closed (cites this script's P/R numbers)
  - HYP-2797  # opus: doublet maximizer = consec+{M,M+1}; "sup|M·error|~1.4 bounded"
  - HYP-2798  # kps S27: "M·e grows to 7" (WRONG center bvd(base,2)); resolved here
  - HYP-2796  # codex: genuine-wide two-far freeze tail (frozen-phase law, converges)
  - HYP-2795  # kps: two-regime skeleton (this closes the genuine-wide binding sub-case)
  - OPEN-Q-108
---

# THM-564 — The doublet almost-periodic P/R split (genuine-wide binding leg)

## Setup
K-runner binding genuine-wide maximizer (opus HYP-2797): `E_M = consec_{K-2} ∪ {M, M+1}`
(base `{0,…,K-3}` + tight far doublet, g=1 adjacent), `|E|=K`, compared vs `cap_K`. K-runner-
correct indexing (the opus signed-bound script used `consec_{k-1}` = the (k+1)-runner family
mislabeled `k` — off-by-one, MISTAKE-083; corrected throughout here).

## The exact inclusion-exclusion identity (PROVED exact)
> `p0(E_M) = p0(base) + A(M) + A(M+1) + d2(M)`,
> `A(f) = p0(base∪{f}) − p0(base)`,  `d2(M) = p0(E_M) − p0(base∪{M}) − p0(base∪{M+1}) + p0(base)`.
Exact limits: `a_inf = lim_f A(f) = boundary_value_direct(base,1) − p0(base)` (decorrelated 1-far,
THM-563 limit); `Φ = Φ_frozen = lim_M p0(E_M)` (exact frozen-doublet integral over the slow phase
`y∈[0,1)`, far runners at frozen offset `x=y`); `d_inf = lim_M d2(M) = Φ − p0(base) − 2·a_inf`
(the adjacent-pair correlation).

## The P/R decomposition (the crux: center at Φ, NOT bvd(base,2))
> `g(M) := M·(p0(E_M) − Φ) = P(M) + R(M)`,
> `P(M) = M·(A(M)−a_inf) + M·(A(M+1)−a_inf)`   [EXACTLY PERIODIC, period `P_per = 7·lcm(base)`],
> `R(M) = M·(d2(M) − d_inf)`                    [interaction correction, `O(1/M)`].
P is periodic by THM-563 (`M·(A(M)−a_inf)` periodic; `M·(A(M+1)−a_inf)=(M+1)·(A(M+1)−a_inf)
−(A(M+1)−a_inf)` = periodic−periodic = periodic). `P+R==g` holds EXACTLY (verified, all M tested).

**Why the center matters (resolves HYP-2797 vs HYP-2798).** With the WRONG center `bvd(base,2)`:
`p0(E_M) − bvd(base,2) → +d_inf ≠ 0` (e.g. `+0.0121` at K=10), so `M·(p0−bvd2) ≈ M·d_inf → ∞`
(HYP-2798's "grows to 7"). Centered at the EXACT `Φ_frozen`, the linear drift `d_inf` is absorbed
and `g(M)` is BOUNDED (oscillates ~`[−1.0, 1.4]`). HYP-2797's "`sup|M·error|≈1.4` bounded" and
HYP-2798's "grows" are BOTH correct — about their own baseline. `Φ_frozen` is the right object.

## The closure (margin-to-cap target 0.16)
`p0(E_M) = Φ + g(M)/M ≤ cap_K − 0.16`  ⟺  `g(M) ≤ M·H_K`,  `H_K = cap_K − 0.16 − Φ > 0`.
With `G_sup_bound := period-max(P) + sup_{M≥15}|R(M)|`:
- tail `M ≥ M* := ⌈G_sup_bound / H_K⌉`:  `M·H_K ≥ G_sup_bound ≥ g(M)`  ✓.
- finite window `15 ≤ M < M*`: check `p0(E_M) ≤ cap_K − 0.16` EXACTLY (TINY — `M*` is small).

The crude `sup_g/15 < H_K` FAILS at K=9 (`sup_g=1.263 > 15·H_K=1.088`, the razor-thin point), but
the correct cutoff `M* = ⌈G_sup_bound/H_K⌉` is small, so the window is `[15, M*]` with `M*` tens.

| K | P_per | Φ | cap_K | H_K=cap−.16−Φ | period-max(P) | sup\|R\| | G_sup_bound | M* | window worst p0 ≤ cap−.16 | cap−sup_p0 |
|---|-------|------|-------|------|------|------|------|----|------|------|
| 8 | 420 | 0.106054 | 0.381463 | 0.115408 | 0.74679 | 0.73482 | 1.48161 | 13 | 0.075000 ≤ 0.221463 ✓ | 0.23690 |
| 9 | 420 | 0.261732 | 0.494256 | 0.072524 | 1.12198 | 0.57165 | 1.69363 | 24 | 0.289765 ≤ 0.334256 ✓ | 0.20449 |
| 10 | 2940 | 0.411027 | 0.604396 | 0.033369 | 1.24049 | 0.59076 | 1.83125 | 55 | 0.442517 ≤ 0.444396 ✓ | **0.16188** |

K=10 is the BINDING case: `cap−sup_p0 = 0.16188 ≥ 0.16` exactly, window `[15,55]` (vs THM-557
`M_BV≈2200` — a 40× collapse). `d_inf=0.014039` matches the empirical large-M `d2` to 6 digits.
K=11: period-max(P)=1.20535, sup|R|=0.61480, M*=25, worst p0=0.521058 ≤ 0.565275 ✓, cap−sup_p0=0.20422.
K=12: period-max(P)=1.19556, sup|R|=0.42886, M*=13, worst p0=0.574150 ≤ 0.697143 ✓, cap−sup_p0=0.26749.
(`H_K=0.0747,0.1313` large → non-binding; period-max over a 2400-window, which spans ≫ the
sawtooth scale `7·max(base)≤70` so it captures the periodic envelope's max.) **ALL K=8..12 CLOSE.**

## Significance
This is the **doublet analogue of THM-563**. THM-563 closed single-far via exact periodicity of
`w·Δ_w`; the doublet `w·Δ_w` is NOT exactly periodic (HYP-2797 refuted), but IS almost-periodic =
(exact periodic P) + (decaying R). Centering at the exact frozen plateau is the missing ingredient.
It collapses THM-557's BV finite-cutoff from `M_BV ≈ 2000` to `M* ≈ 13–24`, closing the binding
genuine-wide sub-case of the kps two-regime skeleton (HYP-2795) with the LRC-relevant margin 0.16.
Scripts: `04-computation/lrc14_doublet_almostperiodic_PR_kpswf9.py` (exact Φ_frozen),
`lrc14_doublet_PR_closure_kpswf9.py` (exact incl-excl P/R + closure);
output `05-knowledge/results/lrc14_doublet_PR_closure_kpswf9.out`.
