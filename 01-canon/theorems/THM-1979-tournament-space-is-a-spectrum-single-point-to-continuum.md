---
id: THM-1979
title: "TOURNAMENT SPACE FIBERS OVER SCORE DATA — score variance is an exact cyclicity coordinate but not a complete structural coordinate. For ε_n=0 (n odd) and 1/4 (n even), σ²∈[ε_n,(n²−1)/12] and c₃=n(n²−1)/24−(n/2)σ². Hence τ=c₃/c₃,max=(σ²_tr−σ²)/(σ²_tr−ε_n), uniformly in parity. The transitive endpoint is a singleton. The balanced endpoint is regular for odd n and near-regular for even n, and all maximum-cyclic tournaments are strong. CORRECTION: fiber size, strong fraction, and modular-prime fraction are not monotone functions of σ²; equal-variance score fibers can differ. The degree/variance axis is the base of a structurally rich fibration, not the whole tournament space."
status: >
  CORRECTED AND VERIFIED (boxeph-2026-07-21-S198; codex audit 2026-07-21), exhaustive over all
  isomorphism classes n≤7. The affine c₃ law and Landau-fiber counts stand. Integer scores force
  min σ²=ε_n, not always zero; this gives the parity-uniform normalized-temperature formula above
  and the usual odd/even c₃,max formulas. The transitive fiber is 1. Maximum score-sequence fiber
  sizes are 1,1,3,12,47 for n=3..7, but they need not occur at minimum variance: at n=7 size 47
  occurs at both σ²=4/7 and 10/7. Nor are structural statistics monotone in σ²: at n=6,
  σ²=5/4 supports two size-3 all-strong fibers and two size-1 reducible fibers. Thus σ²
  determines c₃ exactly and locates a shell, while forgetting score shape and within-fiber
  structure. See MISTAKE-215.
source: boxeph-2026-07-21-S198 (owner: understand tournament space as a spectrum from a single point (transitive) to a continuum housing the different structure)
depends_on: []
related:
  - THM-1978  # my n≥7 regime / vanishing-reachable-fraction (the σ²→0 continuum is the irreducible interior)
  - THM-1926  # my zeta: ζ=1 at the transitive point, rich at the regular center
  - THM-1955  # my reduction DAG / circulant character thread (lives at the center)
  - THM-1966  # |R| becomes independent of (spectrum,H) at n=7
  - THM-1810  # char_A=xⁿ, the transitive GIT-nullcone vertex (the single point)
  - THM-1880  # kps char_S spread duality: transitive max-spread / Paley zero-spread (the two poles)
  - MISTAKE-215  # parity seam and false monotonic-richness claim corrected
script: 04-computation/tournament_space_spectrum_boxeph_S198.py (+ .out)
---

# THM-1979 — tournament space is a spectrum from a single point to a continuum

## The fibration

Tournament space on n vertices fibers over the **score sequence** (Landau's polytope of realizable
score sequences; counts 2,4,9,22,59 for n=3..7). A coarse **spectral coordinate** is the score spread
```
        σ² = Var(scores) ∈ [ε_n, (n²−1)/12],
        ε_n = 0  (n odd),    ε_n = 1/4  (n even).
```
Cyclicity is its exact affine image:
```
        c₃ = C(n,3) − Σ C(sᵢ,2) = n(n²−1)/24 − (n/2)·σ²
```
(from Σsᵢ² = nσ² + n(n−1)²/4). If
`σ²_tr=(n²−1)/12`, the parity-uniform normalization is

```text
        τ = c₃/c₃,max = (σ²_tr-σ²)/(σ²_tr-ε_n).
```

Thus score spread and cyclicity are the same **base axis**, run oppositely.
They do not determine the score sequence or the tournament lying above it.

## The two poles

**The single point — TRANSITIVE (σ² = (n²−1)/12, scores 0,1,…,n−1).** Fiber = 1 (the unique
transitive class), c₃ = 0, no strong part, reducible, char_A = xⁿ (the GIT-nullcone vertex,
THM-1810), ζ_T = 1 (THM-1926 — invisible to the closed-orbit zeta). The ordered, structureless end.

**The balanced edge — REGULAR / near-regular (`σ²=ε_n`).** For odd `n` the scores are all
`(n−1)/2`; for even `n` they split equally between `n/2−1` and `n/2`. Here `c₃` is maximal,
with the parity splice

```text
        c₃,max = n(n²−1)/24  (n odd),
                  = n(n²−4)/24  (n even).
```

THM-2016 implies every maximum-cyclic tournament is strongly connected. This
balanced region contains circulant/Paley structure and the first independent
`|R|` behavior, but it is not synonymous with quasirandomness and it need not
maximize a finite score-sequence fiber.

## What the coordinate preserves—and what it destroys

The projection `T ↦ σ²(T)` preserves exactly the number of directed triangles.
It destroys the shape of the score sequence and all structure within a score
fiber. The exhaustive census itself supplies sharp counterexamples to the old
monotonic reading:

- At `n=6`, the same variance `5/4` occurs in two size-3 all-strong score
  fibers and two size-1 reducible score fibers.
- At `n=7`, the maximum fiber size `47` occurs twice, at variances `4/7`
  and `10/7`, not uniquely at the balanced edge `0` (whose fiber has size 3).

What *is* monotone is the scalar `c₃` itself. Strong connectivity has the
one-sided threshold supplied by THM-2016: above `c₃,max(n−1)` every tournament
is strong; below it, fibers may mix qualitatively different behavior.

## Limit picture (the honest continuum)

Under the tournament-limit topology, `W ↦ d_W` is likewise a **projection**
from tournamentons to degree functions, not an identification of the two
spaces. Score spread becomes `∫(d_W(x)−½)² dx`. The transitive tournamenton
is an extreme point of this functional; `W≡1/2` is one point in its zero-variance
fiber. Many non-quasirandom regular tournamentons also satisfy `d_W≡1/2`.
This large fiber—the information forgotten by the degree function—is the
correct limiting version of the finite structural continuum.
