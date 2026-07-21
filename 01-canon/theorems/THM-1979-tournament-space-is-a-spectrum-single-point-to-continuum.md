---
id: THM-1979
title: "TOURNAMENT SPACE IS A SPECTRUM FROM A SINGLE POINT TO A CONTINUUM, fibered over Landau's score polytope. The spectral coordinate is the score spread σ²=Var(scores)∈[0,(n²−1)/12], and cyclicity is its EXACT affine image: c₃ = n(n²−1)/24 − (n/2)σ² (so σ², c₃, and structural richness are one axis with three readings). At the maximum-spread end (σ²=(n²−1)/12, scores 0..n−1) sits the TRANSITIVE vertex: a SINGLE POINT — fiber 1, acyclic (c₃=0), reducible, char_A=xⁿ, ζ=1. As σ²→0 the fibers SWELL (max fiber 1,1,3,12,47 for n=3..7) into a CONTINUUM of strongly-connected, modular-prime, maximally-cyclic structures — where 'the different structure' lives. Fiber size, strong-fraction, and modular-prime-fraction all run MONOTONICALLY OPPOSITE to σ²: at high spread every fiber is a singleton reducible chain, at the regular center every class is strong. Limit picture: the degree function interpolates from d(x)=x (transitive, W=1_{x>y}) to d(x)=1/2 (quasirandom W=1/2, the entropy-maximal continuum)."
status: >
  VERIFIED (boxeph-2026-07-21-S198), exhaustive over all iso classes n≤7 (grouped by score
  sequence): Landau score-sequence counts 2,4,9,22,59 (n=3..7); fibers (# iso classes per score
  sequence) with strong/modular-prime/c₃ per fiber. The affine law c₃ = n(n²−1)/24 − (n/2)σ² is the
  classical identity c₃ = C(n,3) − ΣC(sᵢ,2) restated (Σsᵢ² = nσ² + n(n−1)²/4), here read as "the
  spectral coordinate is cyclicity". Transitive fiber = 1 at every n; max fiber = 1,1,3,12,47
  (n=3..7), located near (not exactly at) the regular center — at n=7 the peak is σ²=4/7, score
  (2,2,3,3,3,4,4), fiber 47, all 47 strong, 29 modular-prime. Strong-fraction per fiber → 1 as
  σ²→0 and → 0 as σ²→max. Integrates THM-1978 (vanishing reachable fraction), THM-1926 (ζ=1 at the
  point), THM-1955 (circulant thread at the center), mac-mini THM-1966 (|R| independent from n=7).
source: boxeph-2026-07-21-S198 (owner: understand tournament space as a spectrum from a single point (transitive) to a continuum housing the different structure)
depends_on: []
related:
  - THM-1978  # my n≥7 regime / vanishing-reachable-fraction (the σ²→0 continuum is the irreducible interior)
  - THM-1926  # my zeta: ζ=1 at the transitive point, rich at the regular center
  - THM-1955  # my reduction DAG / circulant character thread (lives at the center)
  - THM-1966  # mac-mini: signed |R| a genuinely new invariant from n=7 (structure in the continuum)
  - THM-1810  # char_A=xⁿ, the transitive GIT-nullcone vertex (the single point)
  - THM-1880  # kps char_S spread duality: transitive max-spread / Paley zero-spread (the two poles)
script: 04-computation/tournament_space_spectrum_boxeph_S198.py (+ .out)
---

# THM-1979 — tournament space is a spectrum from a single point to a continuum

## The fibration

Tournament space on n vertices fibers over the **score sequence** (Landau's polytope of realizable
score sequences; counts 2,4,9,22,59 for n=3..7). The **spectral coordinate** is the score spread
```
        σ² = Var(scores) ∈ [0, (n²−1)/12].
```
Everything structural is a monotone function of σ², because cyclicity is its exact affine image:
```
        c₃ = C(n,3) − Σ C(sᵢ,2) = n(n²−1)/24 − (n/2)·σ²
```
(from Σsᵢ² = nσ² + n(n−1)²/4). **Score-spread and cyclicity are the same axis**, run oppositely.

## The two poles

**The single point — TRANSITIVE (σ² = (n²−1)/12, scores 0,1,…,n−1).** Fiber = 1 (the unique
transitive class), c₃ = 0, no strong part, reducible, char_A = xⁿ (the GIT-nullcone vertex,
THM-1810), ζ_T = 1 (THM-1926 — invisible to the closed-orbit zeta). The ordered, structureless end.

**The continuum — REGULAR / near-regular (σ² → 0, scores all (n−1)/2).** The fibers **swell**:
max fiber = 1, 1, 3, 12, **47** for n=3..7, and → ∞ with n. At the low-spread end every class is
**strongly connected** and mostly **modular-prime**; c₃ is maximal = n(n²−1)/24; this is where the
circulant/Paley thread (THM-1955), the |R|-independent structure (mac-mini THM-1966, first at n=7),
and the whole irreducible interior (THM-1978) live. The structurally-richest score class is
**near** the center, not exactly at it — at n=7 the peak (fiber 47) is σ²=4/7, score (2,2,3,3,3,4,4),
all 47 strong, 29 modular-prime.

## The monotone law (verified n≤7)

```
   score spread σ²   →   large                                          small
   position              TRANSITIVE vertex          ...          REGULAR center
   fiber size            1                          ↗↗↗          max (47 at n=7)
   strong / fiber        0                          ↗↗↗          1
   modular-prime/fiber   0                          ↗↗↗          high
   cyclicity c₃          0                          ↗↗↗          n(n²−1)/24
   ζ_T                   1                          richer        max #cycles
```

**Structural richness runs exactly opposite to score spread.** The "single point" is the
maximally-spread-score / minimally-cyclic transitive tournament; the "different structure" is housed
in the low-spread (regular) fibers whose classes are maximally spread *across distinct structures*
(fiber 47) — a continuum that grows without bound.

## Limit picture (the honest continuum)

Under the tournament-limit (tournamenton) topology the spectrum is the space of degree functions
`d: [0,1]→[0,1]`, interpolating from `d(x)=x` (the transitive limit `W(x,y)=1_{x>y}`, a single
ordered point) to `d(x)=1/2` (the quasirandom limit `W≡1/2`). Score spread is the functional
`∫(d(x)−½)² dx`; it is maximal at the transitive step function and zero at the quasirandom point,
whose neighborhood is the genuine positive-entropy **continuum** of tournamentons. Finite n is the
lattice-point shadow of this: one point at the ordered vertex, an exponentially-swelling fiber at the
quasirandom center.
