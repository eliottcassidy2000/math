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
