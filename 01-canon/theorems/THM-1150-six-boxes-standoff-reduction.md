---
id: THM-1150
title: THE MAXIMISER REDUCED TO ONE DIOPHANTINE STANDOFF STATEMENT — my own tube framing is refuted and replaced by six isolated boxes. (0) THE TUBE CLAIM OF THM-1149 IS WRONG: distance to the line L = {(s,2s,3s)} does not separate the bad region — points inside B reach sup-distance 0.228 while points outside dip to 0.026. B is not a tube. (I) WHAT IT ACTUALLY IS: badness depends only on g = (frac(−d₂u), frac(−d₃u), frac(−d₄u)), so B ⊂ T³ is a FIXED region independent of d, and **B is six isolated boxes around the six permutations of (1/4, 1/2, 3/4)** — all 260 sampled points of B lie within sup-distance **0.0917** of a centre, every one inside 1/8, and the measured volume |B| = 0.003367 gives box half-width ρ = (|B|/48)^{1/3} = **0.0412**, the same order. (II) THE FLOW PICTURE: u ↦ g is a CLOSED GEODESIC of direction (d₂,d₃,d₄), so the bad measure is its **sojourn time in B** — and with |B| = 0.0034 against a (1,2,3) sojourn of 2/21 = 0.0952, the geodesic is **28× concentrated**, nowhere near equidistributed. (III) THE CENTRE-HITTING CRITERION IS EXACT AND ALGEBRAIC: passing through (1/4,1/2,3/4) requires d₂u ≡ −1/4, d₃u ≡ −1/2, d₄u ≡ −3/4 (mod 1); the first gives u = (m − 1/4)/d₂ and substituting forces **d₃/d₂ = 2 and d₄/d₂ = 3** — exactly proportionality to (1,2,3). (IV) VERIFIED STARKLY: minimum sup-distance from the geodesic to a centre is **0.00000** for (1,2,3), (2,4,6), (3,6,9) — they hit exactly — and **0.046 to 0.107** for (1,2,4), (1,3,5), (2,3,4), (1,2,5), (3,5,7), (1,4,7), (2,4,7). Sojourn is correspondingly 0.0957 for the aligned family and 0.000000 for every other direction tested. (V) SO THE WHOLE MAXIMISER CLAIM NOW RESTS ON ONE STATEMENT: **every non-proportional integer direction keeps standoff > ρ ≈ 0.041 from all six centres.** The smallest observed standoff is 0.04583 at (2,4,7) — above ρ, but only just, which is why some triples show tiny nonzero sojourn
status: (0) REFUTED — my own THM-1149 tube framing. (I) verified by sampling plus the volume match; the exact box shape is not derived. (II) exact (the flow identification) with the concentration factor measured. (III) PROVED — elementary algebra. (IV) MEASURED on ten directions. (V) is the remaining gap and it is NOT proved: the standoff lower bound is verified on seven non-proportional triples, not established in general. **I was asked twice to prove the maximiser claim and have not done so** — what this delivers is a reduction of it to a single sharp Diophantine statement, with everything else exact. Uniform r=5 remains OPEN
source: kind-pasteur-2026-07-18-S128 (cont.78; owner: prove d ∝ (1,2,3) is the maximiser — second request)
depends_on:
  - THM-1149    # whose tube framing this refutes and whose sketch it sharpens
  - THM-1148, THM-1147
related: [MISTAKE-171]
script: 04-computation/tube_geometry_kps_S128c78.py, six_boxes_kps_S128c78.py (+ .out)
---

# THM-1150 — six boxes, and the standoff that is left

## (0) My tube framing was wrong

THM-1149 proposed that B is a tube around L = {(s,2s,3s)}. It is not. Sup-distance to L
does not separate:

| | min | median | max |
|---|---|---|---|
| points **inside** B | 0.0025 | 0.0875 | **0.228** |
| points **outside** B | **0.0255** | 0.1445 | 0.2415 |

Inside-points reach further from L than many outside-points. The tube picture is refuted.

## (I) B is six isolated boxes

Badness depends only on g = (frac(−d₂u), frac(−d₃u), frac(−d₄u)), so B ⊂ T³ is a **fixed**
region independent of d. And it is concentrated at the six permutations of **(1/4, 1/2, 3/4)**:

- all 260 sampled points of B lie within sup-distance **0.0917** of a centre, every one
  inside 1/8;
- the measured volume |B| = **0.003367** gives box half-width ρ = (|B|/48)^{1/3} = **0.0412**,
  the same order as the observed spread.

This is the correct replacement for the tube: not one connected axis, but six separated
lumps.

## (II) The flow picture

The map u ↦ g is a **closed geodesic of direction (d₂,d₃,d₄)** in T³, so

> **bad measure(d) = sojourn time of that geodesic in B.**

With |B| = 0.0034 against a (1,2,3) sojourn of 2/21 = 0.0952, the geodesic is **28×
concentrated** — far from equidistributed, which is exactly why the direction matters.

## (III) The centre-hitting criterion — proved

Passing through (1/4, 1/2, 3/4) requires

> d₂u ≡ −1/4,  d₃u ≡ −1/2,  d₄u ≡ −3/4  (mod 1).

The first gives u = (m − 1/4)/d₂; substituting into the second and third forces
**d₃/d₂ = 2 and d₄/d₂ = 3**. So a geodesic passes through a centre **iff** its direction is
proportional to (1,2,3). This step is elementary and complete.

## (IV) Verified starkly

| direction | hits a centre exactly | min sup-distance | sojourn |
|---|---|---|---|
| (1,2,3) | **YES** | **0.00000** | 0.0957 |
| (2,4,6) | **YES** | **0.00000** | 0.0962 |
| (3,6,9) | **YES** | **0.00000** | 0.0957 |
| (3,5,7) | no | 0.05000 | 0.000000 |
| (2,4,7) | no | **0.04583** | 0.000000 |
| (2,3,4) | no | 0.07143 | 0.000000 |
| (1,2,4), (1,3,5), (1,2,5), (1,4,7) | no | 0.068–0.107 | 0.000000 |

## (V) What is left — one statement

Everything reduces to:

> **Every non-proportional integer direction keeps standoff > ρ ≈ 0.041 from all six centres.**

If true, non-proportional geodesics miss B entirely and the maximiser claim follows. The
smallest observed standoff is **0.04583** at (2,4,7) — above ρ, but only just, which is
precisely why a few triples show tiny nonzero sojourn rather than exact zero. The margin is
thin enough that this needs proof, not sampling.

## Honest status

**I was asked twice to prove that d ∝ (1,2,3) is the maximiser, and I have not proved it.**
What this session delivers is a genuine reduction: my own tube framing refuted, the correct
six-box geometry identified, the centre-hitting criterion proved algebraically, and the whole
claim reduced to a single Diophantine standoff bound — with a thin margin (0.0458 against
0.0412) that makes the remaining step non-trivial rather than routine. **Uniform r=5 remains
open.**

## Named next
- Prove (V). For fixed d the standoff is min over u of the sup-distance from
  (frac(−d₂u), frac(−d₃u), frac(−d₄u)) to the six centres — a computable piecewise-linear
  minimisation. What is needed is a lower bound over all non-proportional d, which is a
  simultaneous-approximation statement: how well can (d₂,d₃,d₄)·u track (1,2,3)/4 mod 1
  without d being proportional to (1,2,3)?
- Derive the box shape exactly rather than sampling it, so ρ is a theorem rather than a
  volume estimate. The interior region is already explicit (sorted edges with spacing
  constraints); what is missing is the boundary/partial-tooth cases.
