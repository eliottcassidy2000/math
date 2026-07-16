---
id: LEM-030
title: DISCRETE RIESZ REARRANGEMENT ON Z_n (the arc-Green lemma) — for every symmetric kernel g decreasing in circular distance, the contiguous arc maximizes E_g(A) = Σ_{{s,t}⊆A} g(s−t) among m-subsets of Z_n; proof by POLARIZATION (two-point symmetrization): the four-point exchange inequality collapses, via the reflection isometry d(x,y) = d(σx,σy), to "same-side pairs are closer" (d(x,y) ≤ d(x,σy) for x,y ∈ H), so each polarization weakly increases E_g; iteration terminates at the arc (odd n: every reflection has one fixed point; a set polarized for all reflections is an arc); strict kernels give uniqueness up to rotation/reflection
status: PROVED (the polarization step and termination as above) + REFEREED EXHAUSTIVELY (step monotonicity: ALL subsets × ALL reflections, n = 5, 7, 9, 11, 13 — zero failures; 300 random iteration trials all reach arcs; size-balance corollary exact to n = 41)
source: death-star-2026-07-16-S29 (owner: prove the arc-Green lemma); the cycle case of Baernstein–Taylor-style polarization, written self-contained for the repo
depends_on: [THM-913 (consumer), HYP-7097 (the Fourier reduction that named this lemma)]
verification: 04-computation/arc_green_polarization_deathstar_S29.py -> 05-knowledge/results/arc_green_polarization_deathstar_S29.out
---

# LEM-030 — discrete Riesz on Z_n by polarization

Let σ(x) = c − x be a reflection of Z_n (n odd), x₀ its fixed point, H the closed half
{x : d(x, x₀) ≤ (n−1)/2 on the x₀-side}. Polarize A: every element whose σ-partner is
absent moves to its H-representative; both-present or fixed elements stay. For x, y ∈ H:
d(x, y) ≤ d(x, σy) (same-side proximity) and d(x,y) = d(σx, σy) (isometry) — so in every
mixed orbit configuration the polarized pair-distances are coordinatewise ≤ the original,
and g decreasing gives E_g(A^σ) ≥ E_g(A). Iterating over the n reflections: the
center-of-mass potential strictly improves at every proper move, and a set unmoved by all
polarizations is an arc. ∎

Consequences: (i) with g = the cycle Green's function (Fourier 1/(4sin²(πk/n)), decreasing
✓), THM-913's coloring-universality holds for ALL odd n — the contiguous split is the
exact max-cut of the class circulant; (ii) equivalently the arc minimizes Σ_{pairs} d(n−d);
(iii) the lemma is a general repo tool (any circular symmetric-decreasing pair-energy —
the covering program's tent kernels qualify).
