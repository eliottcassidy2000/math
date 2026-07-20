---
id: THM-1289
title: THE n=14 FLOOR IS ISOLATED FROM ABOVE, BY CITATION — there exists δ > 0 such that no 13-speed integer family has M ∈ (1/14, 1/14 + δ). More generally, NO real number is approached from above by 13-speed M-values (every attained value has empty spectrum-space immediately above it), and at 12 speeds (1/13, 1/13 + δ′) is likewise empty — the ineffective, all-heights sibling of HYP-7920's effective bounded-height micro-gap 113/1466. Consequence: THM-1268's second horn ("…or 1/14 is an accumulation point from above and there is no gap at all") is DEAD; the (1/14, 3/41) question is now "how big is δ, and what finite/conditional content sits in [1/14+δ, 3/41]"
status: CITATION-GRADE modulo one named residual check. Source theorem: Giri–Kravitz, "The structure of Lonely Runner spectra", arXiv:2304.01462 v4 (2025-03-30), PUBLISHED — Math. Proc. Camb. Phil. Soc., DOI 10.1017/S0305004125101497 — Theorem 1.4: "Let n ≥ 2 … 1 ≤ k < n … Then the sets S_k(n), S*_{k'}(n) have only upper accumulation points…" with their Definition: x is an UPPER accumulation point of S if S ∩ (x, x+ε) ≠ ∅ for all ε > 0 (lower defined analogously). Pinned this session across four independent extraction passes; internally cross-validated by the mirror Theorem 1.3 (Kravitz [21] Thm 6.5: lower accumulation points of the tilde/M-side spectrum S̃(n) contain S̃(n−1) — the rising s/(ns+1) ladders) and by their explicit S(n) = 1/2 − S̃(n). RESIDUAL CHECK (small): a final human/PDF glance at the printed Theorem 1.4 phrase "have only upper accumulation points" — everything else (definitions, version history, journal status, the withdrawn-v3 story, footnote-2 error mechanism) is triple-consistent. δ is INEFFECTIVE (their argument is qualitative degeneration/compactness); effectivizing it is the named follow-up (backlog S401/S402)
source: opus-2026-07-19-S402 (owner: do the G-K pinning session to promote HYP-7930)
depends_on: [G-K Theorem 1.4 (citation, pinned), the translation lemma below (proved in-file), THM-1050 (dilation invariance — spectrum over primitive = over all), LRC(≤13) NOT needed for this statement]
scripts: 04-computation/lrc_gk_pinning_controls_opus_S402.py -> 05-knowledge/results/lrc_gk_pinning_controls_opus_S402.out (Fan–Sun positive control 4/4 exact + the gridmax boundary value)
---

# THM-1289 — the floor is isolated from above (and every value is)

## The translation lemma (proved here)

For x ∈ [0,1), circle distance to 1/2 satisfies d(x, 1/2) = 1/2 − d(x, 0) (check both
halves). Hence for a 1-dim subtorus T_v = {tv mod 1} with v ∈ ℤ^n primitive and all
coordinates nonzero (⟺ T_v not contained in the union of coordinate hyperplanes):
D(T_v) = min_t max_i d(v_i t, 1/2) = 1/2 − max_t min_i ‖v_i t‖ = 1/2 − M(v),
both extrema attained (piecewise-linear continuous on a compact circle). Every
nondegenerate 1-dim subtorus arises this way, and M is dilation/sign/permutation
invariant, so G-K's S(n) is exactly {1/2 − M(v)} over n-speed families, n = #speeds.
The map x ↦ 1/2 − x swaps upper and lower accumulation.

## The statement and its one-line proof

G-K Theorem 1.4 (published v4): S(13) — their n = 13 speeds, k = 1 — has only upper
accumulation points; equivalently the M-side spectrum S̃(13) has only LOWER accumulation
points: **no real x has 13-speed M-values in (x, x+ε) for every ε**. Apply at x = 1/14:
if (1/14, 1/14+δ) met the spectrum for every δ, then 1/14 would be approached from above
in S̃(13), i.e. 1/2 − 1/14 would be a lower accumulation point of S(13) — contradicting
Theorem 1.4. ∎ Identically at 12 speeds above 1/13, and at any x whatsoever.

Mirror-consistency check: their Theorem 1.3 says each S̃(n−1)-value IS approached from
below (the classical rising ladders s/(ns+1) → 1/n) — from-below approach is the allowed
direction, from-above the forbidden one; the two theorems are two sides of one coin, which
is strong evidence the extraction of the 1.4 phrase is faithful.

## What this does and does not give

- KILLS: the accumulation horn of THM-1268's dichotomy; any strategy premised on values
  m ↓ 1/14 (slack ladders D/(14D−1) realize only finitely often near the floor... NO —
  careful: finitely-near-the-floor needs more; what dies is any INFINITE descent INTO
  (1/14, 1/14+δ); rungs above 1/14+δ are untouched by this theorem).
- DOES NOT give: finiteness of the whole window (1/14, 3/41]. The proven G-K chain
  terminates in S*₀(n−2) — D-values of finite proper subgroups = GRID-LONELINESS values
  of the 11-torus — which populate the window's closure (gridmax((1,…,11); 14) = 1/14,
  verified exactly). From-below accumulation onto interior gridmax values is not excluded
  by the published theorem; full finiteness needs their Conjecture 1.5
  (acc(S(n)) = S(n−1); their stated sufficient condition: S₂(n) = S(n−1)). See
  HYP-7930-UPDATE for the split and the new gridmax-window question this localizes to.
- δ ineffective; the effective companion below bounded height is HYP-7920's cage
  (12 speeds, (1/13, 113/1466) empty to height 258,276).

## Version archaeology (why the repo missed this for 13 days)

v1/v2 (2023) claimed acc(S(n)) = S(n−1) as the main theorem — the sentence the repo's
proof-map quoted in the MISTAKE-117 era. v3 (2023-12) is a WITHDRAWAL notice. v4
(2025-03, published) proves Theorem 1.4 and demotes the equality to Conjecture 1.5; their
footnote 2 names the error: disconnectedness when intersecting subtori with coordinate
subspaces — i.e. intersections are subGROUPS, not subtori — exactly the star-subtlety
that forces the proven chain into S*₀-land. MISTAKE-117 retracted the repo's SUP-misuse
of this paper; the present theorem is the paper's actual (published) content used for
what it governs. MISTAKE-190 records the S401 headline's propagated mis-chain.
