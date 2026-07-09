---
id: THM-675
title: The (C1) signed-box run — the per-relation closed form is EXACT at first order (unit costs — Schur −0.00129, midpoint −0.00116, ratio +0.00120: ratio relations HELP) but FIRST-ORDER ADDITIVITY FAILS BOTH WAYS (destructive interference on generic sets: predicted deficit 3–20× above measured, H-divergent; constructive on coherent/dilated sets: 2× below measured — the deficit is governed by relation-lattice COHERENCE, not scalar budgets, converging with opus-S181); the B5 certificate's BOUNDARY is found (large-V near-dilations: 0/260 certified vs 244/260 exactly live) and DEPTH ESCALATION rescues it (B7 = 100 ≈ LM = 102 at q=280, nearly tight; B11: 206/260; B13 = exact); the branch is scale-persistent (avg-B5 = −4.88 at V=260), so (C2) is necessary — but its job is precisely the longest-AP ≥ k−6 family that LEM-012's ELEMENTARY Dirichlet clustering already handles on the τ-line: THE TWO ROUTES COMPOSE ALONG THE FLEET'S EXISTING DICHOTOMY
status: MIXED, decisive. PROVED: the first-order closed form Δ(j) = Π k(j_l)·(−1)^{|U|}·T_{|U|} with k(j) = sin(πj/7)/(πj) + O(j/q) (symmetric-window Dirichlet kernel; T_2 = +0.0700, T_3 = +0.4898, T_4 = −0.2857, T_5 = 1). REFUTED (load-bearing): first-order additivity of the box — the single-resonance expansion is non-additive at realistic relation densities in BOTH directions; scalar budgets (Σ|Δ|, E_H, E3) cannot decide the deficit (4th proxy failure; the governing structure is lattice coherence/Freiman dimension — opus-S181's law, now at the modular level). MEASURED: the B5 boundary and the depth-escalation ladder; the liveness-law amendment (16/260 dead supra-Vmax moduli on 40→41 — the first at large V; HYP-5731(d) now reads "almost all + branch exceptions"). ARCHITECTURE: the composed dichotomy (below) — the branch needs no new mathematics.
source: klein-2026-07-09-S212 (HYP-5766; owner-directed "run the signed low-relation box (C1)")
depends_on:
  - THM-671   # the B5 certificate whose boundary this maps
  - THM-673   # the aggregated skeleton; (C1)/(C2) named there
  - THM-604   # the truncation identities; B7(13) = 0.1346 > 0 iid
  - LEM-012   # longest-AP ≥ k−6 ⟹ elementary Dirichlet clustering (the τ-line leg of the composition)
related:
  - HYP-5766, HYP-5732, HYP-5731 (amended), HYP-5683/opus-S181 (coherence law)
---

# THM-675 — the signed-box run: interference, the B5 boundary, and the composed dichotomy

## 1. The closed form (proved; the unit-cost table)

Kill window symmetric ⟹ K̂_q(j) = sin(πj(2c−1)/q)/(q sin(πj/q)) — real, explicit,
= k(j) + O(|j|/q) with k(j) = sin(πj/7)/(πj) (signs period-14 in j; zero at 7 | j).
An exact relation j (support U, u = |U|) shifts avg_q B5/(q−1) at first order by

> Δ(j) = Π_{l∈U} k(j_l) · (−1)^u · T_u,  T_u = Σ_{e≤5−u}(−1)^e C(13−u,e) 7^{−e}.

Unit costs: Schur (1,1,−1): **−0.00129**; midpoint (1,−2,1): **−0.00116**;
ratio (2,−1): **+0.00120** (ratio relations RAISE the supply — new and unexpected).

## 2. The interference refutation (measured; load-bearing)

Predicted deficit −ΣΔ vs measured (V = 91..280): generic instances OVERSHOOT 3–20×
(pred +0.15..+0.58 at H=2..3, growing with H — the expansion does not converge
additively; measured +0.004..+0.072): DESTRUCTIVE interference. The V=260
near-dilations UNDERSHOOT (pred +2.7 vs measured +5.0): CONSTRUCTIVE. **No scalar
functional of the relation list decides the deficit; the relation-lattice COHERENCE
does** — the modular-frame instance of opus-S181's dimension/coherence law and the
4th failure of order-blind proxies (E2, R3, envelopes, now Σ|Δ|).

## 3. The B5 boundary and the escalation ladder (measured)

40→41 near-dilation (V=260, primitive, covering; M(S) = 2/17, witness 43/340;
exact LM: 244/260 moduli live, q=280 has LM=102):
**B5 certifies 0/260** — the first total B5 failure — while **B7 = 100 at q=280**
(3/260 certified, nearly tight), B9: 4/260, B11: 206/260, B13 = exact = hist[0].
The odd-depth Bonferroni ladder is monotone-improving and EXACT at D = #classes;
low depth is an ANALYSIS device, not a computational limit. Amendments:
- THM-671 part 5: "100%" was corpus-limited; the boundary is the coherent branch,
  rescued at depth 7.
- HYP-5731(d) liveness law: 16/260 supra-Vmax moduli on 40→41 are DEAD — the first
  at large V; the law is "almost all q ∈ (V, 2V], with a structured exception set
  on the near-dilated branch."

## 4. The composed dichotomy (the architecture; the session's synthesis)

The branch that breaks the modular certificates is the COHERENT family — and it is
exactly the **longest-AP ≥ k−6** family that the τ-line machinery already handles:
LEM-012's Dirichlet clustering is ELEMENTARY and PROVED there (klein-S196/S197;
near-dilations have longest-AP = 12). So:

> **[coherent / longest-AP ≥ k−6]: LEM-012, elementary, τ-line — PROVED.**
> **[generic / dissociated]: aggregated B5, verified +0.05..+0.12 margins — the
> a-priori bound is the remaining open item, now known to require an
> interference-aware (not scalar-budget) argument.**
> **[middle]: depth-7 escalation or LEM-013's census machinery.**

The two proof frames (τ-line and modular) are complementary along the SAME
dichotomy the fleet proved for good-period existence. No new mathematics is needed
for the branch; the open item is the generic-branch a-priori deficit bound, where
interference works FOR us (the truth is better than first order) but proving
cancellation is the hard direction (Mertens caveat) — candidate: exact S_d
evaluation via CRT/lattice geometry per (S, q) class, bypassing resonance
expansions entirely.

## Files

`04-computation/lrc14_signed_box_klein_S212.py` (+ `.out`; includes the corrected
T_u — an off-by-one in the truncation range was caught in-session),
`05-knowledge/results/lrc14_branch_exact_check_klein_S212.out` (exact LM + M(S) +
the escalation ladder).
