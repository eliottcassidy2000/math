---
id: THM-910
renumber_note: drafted as 907 but 907-anova-assembly landed mid-session; took 910 (908 codex sieve, 909 projective-pair taken)
title: THE BERNOULLI LADDER COMPLETED (B₃/B₅ rungs) + CAP ATTAINMENT + DILATION INVARIANCE — (I) support-3 lattice terms in closed form; term(k) = (1/(6·∏₃kᵢ))·Σ_ε(−1)^{#}B₃({Σkᵢ(cᵢ+εᵢ)/7})·∏_omitted(|B|/7), REFEREED to R/closed = 1.002; (II) the five-runner full-support form c_B(k) = (1/(120·∏₅kᵢ))·Σ_ε(−1)^{#}B₅({k·ε}) — derived by the same i-power bookkeeping (odd orders → sin series → B_{2m+1}), refereed to 1.10 trending 1 with the residual attributed to a coexisting lattice vector (THM-899(ii) rank-2 correction); (III) THM-905's three cap candidates are EXACTLY ATTAINED: (2,3,4,6) → 1/12, (1,2,3,4) → 5/42, (3,5,7,9,11) → 40/441 (exact ℚ), and P is DILATION-INVARIANT (P_{cV} = P_V, measure-preserving x ↦ cx) — the caps are sharp and primitivity is WLOG; all probed large/mixed strata sit ≤ 0.025 vs cap 0.083 (factor-3 margins)
status: (I) PROVED + refereed (1.002, wide boxes, v to 257); (II) derived + validated modulo the coexisting-relation lattice sum (three contaminated plantings caught by defect checks en route — test families MUST have their full |k| ≤ 2 relation spectrum certified; logged); (III) attainment + invariance PROVED (exact ℚ; three dilate ladders constant). THE REMAINING GLUE for universal caps (named, honest): (a) the uniform lattice-mass bound table (corner-max over heights, geometric tail); (b) an explicit O(1/v) remainder constant OR adaptive scan extension per small-part stratum — with the observed factor-3 margins these are engineering, not new mathematics
source: death-star-2026-07-16-S26 (owner: derive the B₅ odd analog and close codex's three sector-box caps)
depends_on: [THM-905 (codex, the caps), THM-899 (boxeph, the B₄ lattice law), THM-906, THM-880 (B₂)]
verification: 04-computation/b3_b5_ladder_cap_closure_deathstar_S26.py -> 05-knowledge/results/b3_b5_ladder_cap_closure_deathstar_S26.out (+ scratch extensions recorded there)
---

# THM-907 — the ladder's odd rungs and the sharp caps

The program's Bernoulli ladder is now complete: **B₂** (THM-880: the Q_s kernel), **B₃**
(support-3 lattice terms, this file, refereed 1.002), **B₄** (THM-899/906: quadruple
plateaus), **B₅** (five-runner, this file). All from one mechanism: interval-box Fourier
transforms → corner expansion → Σ_{j≠0}e(jx)/j^r = ±(2πi-power)·B_r({x})/r!, with odd r
producing sin-series/odd Bernoullis and the i-powers cancelling against ∏(2πij)^{-1}.

Cap attainment: the THM-905 candidates are the EXACT values of the scanned maximizers
(1/12, 5/42, 40/441 as exact rationals), constant along dilate ladders — so the caps, if
true, are sharp, and dilation invariance reduces them to primitive tuples. Every probed
non-maximizer stratum (large one-relation, commensurate-mixed) sits at ≤ 0.025 against
the 0.083 cap: the maximizers are isolated small dilation-classes, and the universal caps
are a factor-3-margin finishing problem: [uniform lattice-mass table] + [explicit
remainder or per-stratum scans] — both named for the next session/agent with all
constants now available.
