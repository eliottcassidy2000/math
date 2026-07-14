---
id: THM-740
title: The two-cluster hierarchical-clocks lemma — V = B ∪ (W₁+J₁) ∪ (W₂+J₂) has L ≥ Area₂ − C_in/W₁ − C_out(W₁)/W₂ with the PRODUCT-AREA Area₂ = ∫_{G_B} A_{J₁}(u)·A_{J₂}(u) du exact; separated clusters decouple multiplicatively; closes the separated cone and reduces every fixed W₁ to a THM-739 instance
status: PROVED (outer ride = THM-739 instantiated; inner ride proved below) + VERIFIED-EXACT (product-area confirmed against exact L to ~1e-3 at moderate scales; all bound instances valid; both test shapes)
source: opus-2026-07-13-S274 (owner prompt: prove the two-cluster hierarchical clocks extension)
depends_on:
  - THM-739 (opus-S273)   # the single-cluster clock lemma; the outer ride is its direct instance at base B' = B ∪ (W₁+J₁)
related:
  - THM-737 (multiplicative pack clock), THM-668, THM-731/732, kps THM-735 (Bonferroni j≤6)
  - HYP-6545 (opus; the S273 backlog lead this discharges)
---

# THM-740 — the two-cluster hierarchical-clocks lemma

## Statement

Let `B` (base), `J₁, J₂` (offset patterns, `max Jᵢ = Jᵢₘ`), `|B| + |J₁| + |J₂| = 13`, and for
scales `W₁, W₂` set `V = B ∪ (W₁+J₁) ∪ (W₂+J₂)`. With `A_J`, `G_B` as in THM-739 and

- **`Area₂(B, J₁, J₂) = ∫_{G_B} A_{J₁}(u) · A_{J₂}(u) du`** (the **product-area**; exact rational,
  piecewise-quadratic Simpson integration),
- `C_in = V(A_{J₁}A_{J₂}) + V(A_{J₂}) + 2#comp(G_B) + 2Σ_B b + |J₁|J₁ₘ`,
- `C_out(W₁) = V(A_{J₂}) + 2#comp(G₁) + 2Σ_{B'} b + |J₂|J₂ₘ` where `B' = B ∪ (W₁+J₁)`,
  `G₁ = G_{B'}` (so `C_out ≤ α + 4|J₁|(W₁+J₁ₘ)` with `α = V(A_{J₂}) + 4Σ_B b + |J₂|J₂ₘ + 2`),

we have:

> **(i) outer ride (= THM-739 at base `B'`):** `L(V) ≥ ∫_{G₁} A_{J₂}(u) du − C_out(W₁)/W₂`;
> **(ii) inner ride (new):** `∫_{G₁} A_{J₂}(u) du ≥ Area₂ − C_in/W₁`;
> **(iii) combined:** `L(V) ≥ Area₂ − C_in/W₁ − C_out(W₁)/W₂`.

**Corollary (separated cone + ladder).** If `Area₂ > 0`: every `(W₁, W₂)` with
`W₁ ≥ 3C_in/Area₂` and `W₂ ≥ 3C_out(W₁)/Area₂` has `L ≥ Area₂/3 > 0`. For every *fixed* `W₁`
(any value), the `W₂`-family is a THM-739 instance with base `B'` and area
`∫_{G₁}A_{J₂}` — verified positive (≈ Area₂) at the tested `W₁` — so the two-parameter family
reduces to [the cone] ∪ [per-`W₁` single-cluster closures + finite checks].

## Proof of (ii) — the inner ride

Branch the `W₁`-clock: `u = (n+σ)/W₁`, `x_n = n/W₁`, `σ ∈ [0,1)`.

1. **Base freezes** (frozen fan): if `‖b x_n‖ ≥ 1/14 + b/W₁ ∀b ∈ B`, every `σ` is base-safe;
   the fattening costs measure `≤ 2Σ_B b/W₁` (as in THM-739).
2. **Cluster 1 = fattened AP in σ:** `(W₁+j)u ≡ σ + j·u = σ + j x_n + O(J₁ₘ/W₁)`; the safe-σ
   set has measure `≥ A_{J₁}(x_n) − |J₁|J₁ₘ/W₁`.
3. **The weight rides along:** on branch `n`, `A_{J₂}(u) = A_{J₂}(x_n + σ/W₁) ≥ A_{J₂}(x_n) −
   osc_n`, where `osc_n` is the oscillation of `A_{J₂}` on the branch interval; the branches
   partition `[0,1)`, so `Σ_n osc_n ≤ V(A_{J₂})`.
4. For `a, b ≤ 1`, `x, y ≥ 0`: `(a−x)(b−y) ≥ ab − x − y` (and when either factor is negative
   the branch integral is `≥ 0 ≥ ab − x − y` directly). Hence the branch-`n` contribution is
   `≥ A_{J₁}(x_n)A_{J₂}(x_n) − |J₁|J₁ₘ/W₁ − osc_n`.
5. Riemann-sum `f = A_{J₁}A_{J₂}` (BV, `sup ≤ 1`) over `{x_n ∈ G_B^{fat}}`: error
   `≤ (V(f) + 2#comp(G_B))/W₁` plus the step-1 measure loss. Collecting terms gives `C_in`. ∎

(Part (i) is THM-739 verbatim with base `B'`; part (iii) is (i)+(ii).)

## Exact verification (companion script/out)

| shape | `Area₂` (exact) | `C_in` | true `L` at (120, 12800) | `∫_{G₁}A_{J₂}` at `W₁`=30/60/120 |
|---|---|---|---|---|
| A: `{1}∪(W₁+{0..5})∪(W₂+{0..5})` | **30077/308700 = 0.097431** | 44.7 | 0.096424 | 0.100291 / 0.097927 / 0.096449 |
| B: `{1,2}∪(W₁+{0..4})∪(W₂+{0..5})` | **166417/1852200 = 0.089848** | 40.3 | 0.088831 | 0.094066 / 0.090981 / 0.088857 |

- The **product-area prediction is confirmed**: exact `L` hugs `Area₂` to ~1e-3 already at
  moderate scales (the S273 taste test's `L ≈ 0.10` was `Area₂ = 0.0974`); the per-`W₁` areas
  oscillate around `Area₂` (non-monotone convergence), all positive.
- All bound instances valid. HONEST LOOSENESS: as in THM-739 (whose crude `C` is ~3000×
  conservative), the rigorous cone here is astronomically separated (`W₁ ≳ 1.4k`, `W₂ ≳ 10⁶`
  for shape A); the *true* convergence is orders faster. The named sharpening (second-order
  treatment of the within-branch drift, which is equidistributed rather than adversarial)
  applies to both rides and would collapse the thresholds; not needed for the closure logic.

## Structure notes

1. **Decoupling = the product.** Separated clusters contribute *multiplicative* admissibility
   factors: cross-cluster pair-arcs (the T(n−2) sector between clusters) carry no constraint in
   the limit — they factorize. Within-cluster pairs = coherent (`A_J`); base–cluster pairs =
   frozen fan. Every pair-sector role of the owner's `(n−1)²` perspective decomposition now has
   a theorem.
2. **Separation helps loneliness:** `Area₂({1},{0..5},{0..5}) = 0.0974 >
   Area({1},{0..11}) = 0.0750` — splitting one 12-cluster into two separated 6-clusters raises
   the safe measure. Coherence is what makes covering efficient; the covering-min extremals are
   maximally coherent (THM-730's AP), consistent with the whole S270–S274 arc.
3. **k clusters (induction sketch):** `W₁ ≪ … ≪ W_k`: iterate the inner ride against the
   product weight `∏_{i>1} A_{J_i}`; the same four steps give
   `L ≥ ∫_{G_B} ∏ᵢ A_{J_i} − Σᵢ Cᵢ/Wᵢ` with `Cᵢ` depending on `W_{i−1}` (fully-separated cone).
   The k=2 case above is the induction step's engine; not separately filed.
4. **Remaining after this:** the comparable-scales strip (`W₂ ≍ W₁` — the clusters merge into
   one wide pattern; interpolates toward the genuinely-spread two-scale territory = the density
   route's home turf), and the constant sharpening. Backlog lead updated.

## Files

`04-computation/lrc14_two_cluster_hierarchical_thm740_opus_S274.py` (+ `.out`): exact
piecewise-quadratic Simpson integration of the product-area (two-half self-check per piece),
exact `V` with interior-vertex handling, inner-bound and combined-bound instances, both shapes.
