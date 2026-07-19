---
id: THM-1195
title: THE TENT-INTEGRAL CRITERION, AND WHY LRC(n) IS TIGHT AT λ = 1/n — since g(t) = min_v‖vt‖ vanishes at every k/v and lies under the tent on each cell, every local max < 1/14 would force ∫g < (1/14)(1/2) = 1/28; contrapositive, **∫₀¹ min_v‖vt‖ dt ≥ 1/28 ⟹ LRC(14) for that family**. Computed exactly over the full arrangement it fires on 9/12 random families (ratio to 1/28: min 0.921, median 1.012, max 1.102) and provably fails on the tight families (0.854, 0.853), as any strict criterion must. The conceptual content is the threshold itself: the tent bound gives λ·(1/2), and the independence heuristic gives E[min of n uniforms on [0,1/2]] = (1/2)/(n+1). These agree **exactly** when λ = 1/(n+1) = 1/14 — so the coincidence at 1/28 is not a coincidence, it is LRC(n)'s own tightness at λ = 1/n appearing in two independent computations
status: the criterion is PROVED (the tent bound is elementary) and the integral computed exactly in rational arithmetic over the arrangement of zeros, peaks and crossings. It is SUFFICIENT ONLY — it fails on the tight families and on ~25% of random families, so it cannot prove LRC(14). The threshold coincidence is exact algebra, not numerology
source: opus-2026-07-17-S388 (owner: work a new creative angle on the LRC 14 open math)
depends_on: [THM-1185 (the measure-vs-pointwise triage that motivated a hybrid), THM-1170 (the critical-point structure supplying the arrangement), THM-1120 (the tight families that bound what any criterion can do)]
scripts: 04-computation/tent_integral_opus_S388.py -> 05-knowledge/results/tent_integral_opus_S388.out
---

# THM-1195 — a criterion, and the reason the conjecture is knife-edge

## The criterion

g(t) = min_v ‖vt‖ vanishes at every k/v and rises to a local maximum in
each gap between consecutive zeros. On a cell of length L the graph lies
under the tent of height max, so its area is at most max·L/2. Hence if
*every* local max were < 1/14,

> ∫₀¹ g < (1/14) · (Σ L_i)/2 = (1/14)(1/2) = **1/28**.

Contrapositive:

> **∫₀¹ min_v‖vt‖ dt ≥ 1/28  ⟹  some local max ≥ 1/14  ⟹  LRC(14) holds.**

A single integral, computed exactly over the arrangement of zeros, peaks
and pairwise crossings.

## What it delivers

| family | ∫g | ratio to 1/28 | verdict |
|---|---|---|---|
| {1,…,13} | 0.03051 | 0.854 | fails (must) |
| {1,…,11,13,24} | 0.03047 | 0.853 | fails (must) |
| 12 random | — | min 0.921, **median 1.012**, max 1.102 | fires 9/12 |

It is **sufficient only**. It fails on the tight families — necessarily,
since their maxima equal 1/14 exactly, so the tent bound is saturated — and
on roughly a quarter of random families. It cannot prove LRC(14).

## The real content: why 1/28 appears twice

The tent bound's threshold is λ·(1/2). The independence heuristic — treat
the n values ‖vt‖ as independent uniforms on [0,1/2] — predicts

> E[min of n uniforms on [0,1/2]] = (1/2)/(n+1).

These two thresholds coincide precisely when

> **λ·(1/2) = (1/2)/(n+1),  i.e.  λ = 1/(n+1)**

which at n = 13 speeds is λ = 1/14 — **exactly the conjectured gap**. So the
observed knife-edge (median ratio 1.012, essentially 1) is not numerical
accident. It is LRC(n)'s own tightness at λ = 1/n, showing up in two
computations that know nothing about each other: a deterministic geometric
bound and a probabilistic expectation.

That explains something the programme has kept running into without naming
it. LRC(14) sits at a threshold rather than having slack **because the
conjectured λ is exactly the value at which the geometric and probabilistic
estimates cross**. Any method whose accuracy is of the order of the gap
between those two estimates will land on the boundary — which is what
S₁ = 13/7 > 1, the k = 7 arity ceiling (THM-1155), and now this criterion
all do.

## Status

Proved and exact, but sufficient only. Its value is diagnostic rather than
constructive: it gives a cheap single-integral certificate covering a
majority of families, and it identifies the structural reason no
threshold-based method can cover the rest.
