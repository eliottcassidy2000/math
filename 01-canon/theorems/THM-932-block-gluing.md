---
id: THM-932
title: THE LOCAL-DENSITY BLOCK-GLUING THEOREM — locally-certified blocks compose across scale gaps: μ(∩ W_{B_i}) ≥ m₁·∏η_i − Σ κ(V_{i−1})·ℓ_i with every input exactly computable; the cascade (THM-928A) is the all-singleton case (single-speed local density is EXACT: η_x(k/x) = 1 − 2λ); the LOCALIZATION SCALE of a certificate: the THM-928(C) packet has η_X = 0 below 1/862 but η_X > 0 from ℓ = 1/300 — the uncovered set is spread at the SLOWEST runner's resolution; the reduction corollary: any packet with a big enough scale gap splits into independently certifiable blocks, so the certificate program's remaining hard case is the GAP-FREE bounded-ratio core
status: lemmas G1/G2/G3 PROVED (proofs in-file) + verified exactly (40+40+12 configs, 0 violations); G1/cascade formalized sorry-free in TournamentH7/CascadeGluing.lean and the full recurrence in TournamentH7/LRCCascadeGluing.lean; dilation law exact; composition demos: v1 honest-negative at gaps 18–28 (diagnosed: crude block-κ overcounts 30×), sharp-input v2 at gap 300 — see results block
source: opus-2026-07-16-S333 (owner: prove the local-density block-gluing theorem; aim to close LRC(14) mathematics then formalization)
depends_on: [THM-928 (the two-scale certificate this composes), THM-883/FragmentationLemma.lean (the proved Lean rung below the cascade step)]
scripts: 04-computation/block_gluing_opus_S333.py, block_gluing2_opus_S333.py → 05-knowledge/results/*.out
formalization: 04-computation/lean/TournamentH7/TournamentH7/CascadeGluing.lean; 04-computation/lean/TournamentH7/TournamentH7/LRCCascadeGluing.lean
---

# THM-932 — the local-density block-gluing theorem

**Setting.** Circle 𝕋 = [0,1), gap parameter λ ∈ (0, ½) (LRC(14): λ = 1/14;
certificate convention: λ = 1/13). D_x = {t : ‖xt‖ < λ}; for a block of
speeds B, W_B = 𝕋 ∖ ∪_{x∈B} D_x. κ(·) = circle component count.

**Definition (local density floor).** η_B(ℓ) = min over circle intervals I
with |I| = ℓ of μ(I ∩ W_B)/ℓ. Computable exactly: φ(a) = μ([a,a+ℓ] ∩ W_B)
is piecewise linear in a, so the min is attained among the ≤ 4κ(W_B)
breakpoint candidates {endpoints, endpoints − ℓ}.

## The three lemmas (proved; verified 0 violations)

**G1 (sampling).** For any union V of κ(V) intervals, any block B, any ℓ > 0:
μ(V ∩ W_B) ≥ η_B(ℓ)·(μ(V) − κ(V)·ℓ).
*Proof.* Per component J, |J| = L: tile J from its left end by ⌊L/ℓ⌋
disjoint intervals of length exactly ℓ; each holds ≥ η·ℓ of W_B;
⌊L/ℓ⌋·ℓ ≥ L − ℓ (components with L < ℓ contribute ≥ 0 ≥ η(L − ℓ)). Sum. ∎

**G2 (components).** κ(U ∩ V) ≤ κ(U) + κ(V) (each intersection component's
left endpoint is a left endpoint of U or of V, injectively);
κ(W_B) ≤ Σ_{x∈B} x (complement of ≤ Σx arcs). ∎

**G3 (single-speed exactness).** D_x is (1/x)-periodic, so any interval of
length k/x (k ∈ ℕ) meets it in measure exactly 2λ·(k/x):
**η_{x}(k/x) = 1 − 2λ exactly.** ∎ (Verified 12/12 across x, k.)

## The theorem

**Block gluing.** Let B₁ < … < B_r be blocks (every speed of B_{i+1}
exceeds every speed of B_i), V_i = W_{B₁} ∩ … ∩ W_{B_i}. For any scales
ℓ_i > 0, with η_i = η_{B_i}(ℓ_i) and m₁ = μ(W_{B₁}):

> **μ(V_r) ≥ m₁·∏_{i=2}^r η_i − Σ_{i=2}^r κ(V_{i−1})·ℓ_i·∏_{j=i+1}^r η_j
> ≥ m₁·∏ η_i − Σ κ(V_{i−1})·ℓ_i.**

*Proof.* Induction: μ(V_i) ≥ η_i·(μ(V_{i−1}) − κ(V_{i−1})·ℓ_i) by G1;
unroll, using η ≤ 1 on the loss terms. ∎

**Sharp-input form (the certificate shape).** Every hypothesis quantity —
m₁, η_i (breakpoint sweep), κ(V_{i−1}) (exact interval arithmetic) — is
exactly computable from the blocks alone, in time O(Σ x·μ) (subtract combs
only inside the current uncovered set). The a-priori form
κ(V_{i−1}) ≤ Σ_{j<i} Σ_{x∈B_j} x (G2) is valid but overcounts ~30× in
practice (the v1 demo's honest negative); the sharp inputs are the point.

**The relative-loss law.** κ(V_i) ≈ (Σ_{x∈B_i} x)·μ(V_i)-scale while
ℓ_{i+1} = c/min(B_{i+1}), so loss_{i+1}/μ_i ≈ c·b_i/G_{i+1} (b_i = block
size, G = junction gap) — scale-free. Coercivity needs
∏η_j(c-scale) ≳ Σ c·b_i/G_i. Empirically (v2 demos) the coercive regime
starts at junction gaps ≈ 46–104 for 4–7-speed blocks with optimized
scales — the optimizer trades η(ℓ) against κ·ℓ per junction.

## Corollaries

**1 (cascade recovery).** All-singleton blocks, ℓ_i = 1/x_i: η_i = 1 − 2λ
exactly (G3), and the theorem reproduces THM-928(A)'s shape (with the
slightly cruder κ·ℓ loss; the (L1)-sharpened constant 2/R stands in
THM-928 for the singleton case).

**2 (the reduction).** If a 13-packet splits at a junction with gap G and
both sides carry local certificates with m₁∏η > Σκℓ, loneliness follows
with NO cross-block condition — no Sidon, no ratio-cleanliness, nothing,
across the gap. **The certificate program's remaining hard case is
exactly the gap-free bounded-ratio core** — matching the S297 frontier
map's capped-envelope regime from the other side.

## The localization scale (new invariant of a certificate)

For the THM-928(C) certified packet X = {300, …, 2208} (μ = 0.117744,
κ(W_X) = 2030 exact):

> η_X(ℓ) = 0 for ℓ ≤ 1/862; η_X(1/300) = 0.0073; η_X(2/300) = 0.0446;
> η_X(4/300) = 0.0765; η_X(16/300) = 0.1044 → μ = 0.1177.

**The uncovered set is spread at the resolution of the SLOWEST runner
(1/300) and has genuine deserts below 1/862.** A global certificate
(BONF5 > 0) says nothing local; the η-profile is the local upgrade — the
quantity a gluing partner actually needs. The dilation law
η_{cB}(ℓ) = η_B(cℓ) (verified exact) makes one profile serve all scales.

## Composition demos (exact)

**v1 (honest negative, kept as data):** 4+4+5 blocks at gaps 27.6/18.4
with the crude G2 block-κ: bound −0.39 (valid: exact 0.1162 ≥ bound;
not coercive). Diagnosis: crude κ overcounts (67,560 vs true need), and
the relative-loss law demands bigger gaps.

**v2 (sharp inputs):** two-block 7+6 at junction gap 240.8, and three-block
5+4+4 at junction gaps 104.4 / 46.0, composed from X-prefix bases —
RESULTS:

> two-block 7+6 at junction gap 240.8: bound +0.047612 (COERCIVE); exact mu = 0.115395, bound holds: True
>
> three-block 5+4+4 at junction gaps 104.4/46.0: bound +0.024575 (COERCIVE); exact mu(V3) = 0.114998, holds True, nonempty True

## Formalization (closed)

Klein S317's first-pushed `TournamentH7.CascadeGluing` closes the original
draft sorry-free:

- `cascade_step`: the complement lower bound from the proved fragmentation
  theorem, including the measurable lifted arc grid and exact ENNReal subtraction;
- `window_floor_sample`: the full `floor(L/ell)` half-open tiling argument,
  with endpoint-null conversion back to closed intervals;
- `union_floor_sample`: G1 over a finite disjoint interval union, discarding
  short components only after proving their real contributions nonpositive;

The independently developed `TournamentH7.LRCCascadeGluing` was reduced after
pulling that work and now contributes only the nonduplicate recurrence layer:

- `block_gluing_sharp_closed`: the fully unrolled sharp recurrence with
  `eta_i*kappa_i*ell_i` losses and suffix density products;
- `block_gluing_closed`: Opus's published coarser recurrence, derived formally
  from `0<=eta_i<=1` and nonnegative component losses.

The module has no `sorry` and no `native_decide`; every declaration's axiom
audit is exactly `[propext, Classical.choice, Quot.sound]`.
