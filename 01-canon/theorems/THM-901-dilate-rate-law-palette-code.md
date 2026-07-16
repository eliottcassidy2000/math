---
id: THM-901
renumber_note: was THM-898 (S23); boxeph-S33 first-pushed 898 (four-runner box-hit) and mac-mini-S116 also used 898 (resistant-31) — my dilate law takes 901, leaving 899/900 as buffer for their two-way sort
title: THE DILATE RATE LAW — exact [n,k,d] of the movie palette code under dilation. For E = c·E₀ (c ≥ 2): merged events n = c·n₀, states V = V₀, hence k = (c−1)n₀ + k₀ and rate = 1 − (V₀−1)/(c·n₀); d = 2 (every base edge acquires c parallel copies — exactly n₀ parallel pairs); structurally palette_code(cE₀) = (c−1)·n₀ parallel-edge doublets ⊕ lift(palette_code(E₀)). For Hamiltonian bases (V₀ = n₀, k₀ = 1 — the generic small-core case per HYP-7036): rate = 1 − (n₀−1)/(c·n₀) ≈ 1 − 1/c, explaining the S22 census dilate row exactly (1009/1260 = 0.8008)
status: PROVED (dilate covariance: the cE₀-movie is the base movie traversed c times — states pull back under x ↦ cx mod 1; the collision/merge structure scales to the (1/(7c))-grid so merged events multiply exactly by c) + REFEREED 8/8 (two bases × c = 2,3,4,5: [n,V,k] all exact; parallel-pair count = n₀ = 252 exactly in all rows)
source: death-star-2026-07-16-S23 (owner directive: derive the exact rate formula for dilates); completes HYP-7036's named next step (i)
depends_on: [HYP-7036 (the palette-code census + structure theorem), HYP-7027, THM-889 (dilate covariance)]
verification: 04-computation/dilate_rate_law_and_diagpoly_n5_deathstar_S23.py -> 05-knowledge/results/dilate_rate_law_and_diagpoly_n5_deathstar_S23.out
---

# THM-898 — the dilate rate law

**Statement.** Let E₀ have merged-event count n₀, state count V₀, cycle rank k₀ = n₀ − V₀ + 1
in its wall movie. For E = c·E₀, c ≥ 2:

> n = c·n₀, V = V₀, k = (c−1)·n₀ + k₀, **rate = 1 − (V₀ − 1)/(c·n₀)**, d = 2,
> and palette_code(cE₀) ≅ (c−1)·n₀ parallel-edge doublets ⊕ lift(palette_code(E₀)).

**Proof.** σ_{cE₀}(x) = σ_{E₀}(cx mod 1) (sections of cf at x = sections of f at cx), so the
cE₀-movie is the E₀-movie traversed c times; each base event at y pulls back to the c
positions (y + j)/c, and the collision-merge structure is preserved (the shared-boundary
grid scales from (1/7)Z-type to (1/(7c))Z-type), so merged events multiply exactly by c
while the visited state set is unchanged. Each base transition edge (u,v) thus appears c
times: n₀ parallel families, each contributing (c−1) independent 2-cycles (the doublets,
weight-2 palette words pairing the same base wall's c dilate copies); modulo these, the
quotient is the base movie. k = n − V + 1 = (c−1)n₀ + k₀ counts them. d = 2 from any
parallel pair. ∎ (Referee: 8/8 exact, two bases.)

**Corollaries.** (i) Hamiltonian bases (V₀ = n₀ — generic small cores, HYP-7036): rate
≈ 1 − 1/c: the census's 0.801 at c = 5 is 1009/1260 exactly. (ii) rate → 1 as c → ∞: exact
dilates saturate the coding face of coherence — the maximal-coherence endpoint of the
HYP-7036 rate law, now with the constant. (iii) The doublet subcode (weight-2, one per
base wall) is the coding-theoretic fingerprint of exact dilation — a decidable dilate
detector (presence of n₀ disjoint weight-2 palette words).

**Scope note.** NEAR-dilates (consecutive cores) do NOT obey the law (V grows with c:
434/819/1029/1113 at c = 10/20/30/50 — the difference-system modulation splits states);
their V(c) asymptote is the named open follow-up.
