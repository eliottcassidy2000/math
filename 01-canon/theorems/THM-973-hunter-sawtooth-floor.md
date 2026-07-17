---
id: THM-973
title: THE HUNTER–SAWTOOTH FLOOR — the 7-wall's measure estimate on nearly-equal pairs, closed: Hunter's inequality is tree-feasible at exactly k = 7, giving μ(U₇) ≥ Σ path-tree consecutive-pair overlaps; the sawtooth pair floor ρ(r) (exact, λ = 1/14, tabulated) bounds each edge; minimizing over ∏r_i ≤ 13 gives the explicit global floor μ(U₇) ≥ 547/6006 = 0.091076 (sorted-ratio law: Sum_j rho(13^(1/j))); Markov over window positions converts it to a dead-position bound ≤ 0.908924 < 1 — the wall closes modulo window choice (BlockSplitLift's machinery)
status: sawtooth-sum formula verified vs direct intervals 40/40 exact; Hunter floor violated 0/30 on random 7-blocks (exact μ(U₇) range [0.3152, 0.3747]); floor table + minimization exact
source: opus-2026-07-17-S341 (owner: work the measure estimate on one pair of nearly-equal runners)
depends_on: [boxeph LRCTreeHunter (tree_hunter_add_le = link (H) kernel-pure + muNum = the sawtooth integer sum LEM-042, formalized CONCURRENTLY tonight), THM-856 (the sawtooth law), THM-928(B) (beat structure), the S340 crumb law, BlockSplitLift (the window-choice consumer)]
scripts: 04-computation/hunter_sawtooth_floor_opus_S341.py -> 05-knowledge/results/hunter_sawtooth_floor_opus_S341.out
---

# THM-973 — the Hunter–sawtooth floor

**The chain (each link exact).**
(H) Hunter's inequality — μ(∪A_i) ≤ ΣμA_i − Σ_tree μ(A_i ∩ A_j), any sets,
any spanning tree — applied to a 7-block with the SORTED PATH tree: since
Σμ(D_i) = 7·(1/7) = 1, it gives μ(U₇) ≥ Σ_{i=1}^{6} μ(D_i ∩ D_{i+1}):
k = 7 is exactly where the union bound dies AND exactly where the 6-edge
tree resurrects it — the 7.5-wall from the inside.
(S) The sawtooth pair floor: μ(D_a ∩ D_b) exact by the sawtooth sum
(verified 40/40 against direct intervals); the tabulated floor ρ(r) =
min over pairs with b/a ≤ r is decreasing in r (table in the out file).
(P) The consecutive ratios of a comparable 7-block satisfy ∏ r_i ≤ 13;
minimizing the 6-edge floor sum over allocations gives
**μ(U₇) ≥ 547/6006 = 0.091076 (sorted-ratio law: Sum_j rho(13^(1/j)))** (worst allocation identified in the out file).
(M) Markov over window positions x: E_x[covered-in-window] = L·μ(∪), so
the dead-position fraction is ≤ 1 − μ(U₇) ≤ **0.908924** < 1.

**What it closes.** The 7-wall's analytic content is DONE: every
7-comparable-block has a quantified global uncovered floor and hence a
quantified live-position set at every window scale; combined with any
formalized window-choice mechanism (the fleet's BlockSplitLift already
navigates window candidates), the wall is crossed. S340's empirics said
the dead set is ≤ 7%; this chain proves ≤ 0.908924 unconditionally — the gap
between them is the (unneeded) sharpness left in the sawtooth floors.

**The instrument note.** The estimate the wall needed was, in the end, the
project's OLDEST exact object — the codex/opus sawtooth pair law (THM-856,
S313) — composed with a 19th-century inequality and Fubini. The close pair
that S340 identified empirically (ratio ≈ 13^(1/6)) is where the path-tree
floor concentrates.
