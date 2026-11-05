---
id: HYP-9471
title: "Seam-gadget semigroup thickness by piece order. On the plane P = {S = 0} of SEAM states (H, S, R) whose exit starts no Hamiltonian path, the S-row-zero gadgets act by 2x2 integer maps [[a, b], [c, d]] on (H, R). The covering density D_k(t) = sum over the distinct nonsingular such maps of order <= k of (a + b t)^(-2) is 0.117, 0.257, 0.486 at t = 2 for k = 4, 5, 6, roughly doubling per order. D_k(t) >= 1 is the necessary density condition for a finite-alphabet self-covering of a 2D cone of states. Pre-registered: D_7(2) lies in [0.75, 1.25], and the order-7 alphabet thickens the realised P-slice"
status: PRE-REGISTERED (predictions P1-P3 committed in a98199e8bb, P4-P5 in 973c1f8312, both under the colliding name HYP-9470; renumbered HYP-9471 because opus-S79 reserved HYP-9470 first, in 395a317b0e; results below)
session: kind-pasteur-2026-11-04-S11
date: 2026-11-04
related: [HYP-9468, HYP-9469, THM-4865, THM-4875, THM-4876]
script: 04-computation/experiments/seam_gadget_density_20261104/
---

# HYP-9471 - when does a finite gadget alphabet fill a plane of seam states?

## Setting (inherited)

- **Seam transfer.** SEAM.md (`05-knowledge/results/seven_return_native_recursions_20261029/SEAM.md`) is the
  universal three-coordinate seam transfer. For a rooted block (X, a), v(X, a) = (H, S, R):
  - H = H(X);
  - S = the number of Hamiltonian paths starting at a;
  - R = sum over U ∋ a of E(X[U], a) H(X \ U).
- **Composition.** For every marked gadget (G, u, z) there is an integer 3x3 matrix with
  v(X ⊙ G, z) = M_G v(X, a), and M_(X⊙Y) = M_Y M_X. The values of H on seam words are therefore an orbit of a
  matrix semigroup.
- **Related work.** opus-S78 (THM-4875, HYP-9468) and Sessions 193-195 use the same calculus or richer ones. The
  Session 195 combined grammar (h, e, R, d), with gadgets of order <= 5 and two order-6 pieces, misses only
  {7, 21} through 20,000.

## Exploratory facts (computed before this file; not predictions)

- **Census of transfer matrices.** The number of distinct transfer matrices of marked gadgets of order exactly k
  (all tournaments, all entry/exit port pairs) is:

  | k | distinct matrices |
  |---|---|
  | 1 | 1 |
  | 2 | 4 |
  | 3 | 12 |
  | 4 | 44 |
  | 5 | 220 |
  | 6 | 1,660 |
  | 7 | 20,500 |

  Order 7 was counted only; nothing below was computed from it.
- **Plane P.** P = {S = 0}. A gadget whose S-row is zero maps every state into P and acts on P by
  (H, R) -> (aH + bR, cH + dR). The singleton seam gives H' = H - S + R, and Redei makes H and H' odd, so R ≡ S
  (mod 2) for every rooted tournament. States in P therefore have H odd and R even.
- **Covering density.** D_k(t) is summed over the distinct nonsingular 2x2 restrictions of order <= k:

  | k | D_k(2) | D_k(3) |
  |---|---|---|
  | 4 | 0.1170 | 0.0596 |
  | 5 | 0.2571 | 0.1315 |
  | 6 | 0.4861 | 0.2498 |

  The realised P-states have R/H >= 2.01 in the runs below, so t = 2 is the dense end of the cone.
- **Seam BFS, order <= 6** (all 1,941 matrices; states with H, R <= X expanded; seeds = small realisable triples):
  - every odd H <= 20,000 except 7 and 21 is attained;
  - distinct states: 240,925, 4,743,006 and 17,081,577 for X = 2,000, 10,000 and 20,000. That is growth about X^1.85,
    a thin orbit in Z^3.
- **Baseline P-slice fill** (order <= 6, X = 5,000). The fill is the fraction of lattice points (H odd, R even) in the
  band that are realised states with S = 0:
  - band H in [801, 999], R in [2H, 2.75H]: 0.2640;
  - band H in [1601, 1799], same slopes: 0.2389;
  - ratio 0.905, so the realised fraction falls with scale.

## The heuristic being tested

Treat the images of the realised set under the maps as a random covering. A point at scale H is then covered with
expected multiplicity about μ f, where μ ≈ D_k(t) and f is the fill at the smaller scale. The self-similar fill
satisfies f = 1 - exp(-μ f).
- This has a positive solution only when μ > 1.
- **Below 1** the orbit is thin: the fill decays with scale and the state count grows like X^δ with δ < 2 in the
  plane.
- **Above 1** the fill tends to a positive constant.

So D_k(t) = 1 is the threshold where a finite alphabet starts to fill a 2D region with positive density. It is a
precondition, not a sufficient condition, for a finite self-covering certificate. Such a certificate would also need
lattice-coset alignment: S78's descent lemma rules out unimodular coverings, and every P-map has even determinant.

## Pre-registered predictions (order <= 7 alphabet)

- **P1.** D_7(2) lies in [0.75, 1.25] and D_7(3) in [0.38, 0.65]. This is the doubling trend continued.
- **P2.** The same seam BFS with all matrices of order <= 7 (X = 5,000, same seeds and bands) has P-slice fill on the
  band H in [801, 999], R in [2H, 2.75H] at least 1.5 x 0.2640 = 0.396.
- **P3.** Under that run, the ratio fill[1601, 1799] / fill[801, 999] exceeds the order-6 value 0.905: the alphabet
  thins out less with scale.

Every outcome will be reported. P2 and P3 are measured with the same program, seeds, bands and X as the baseline.

**Addendum (committed after P1 was read, before any order-8 computation).**
- P1 came out D_7(2) = 0.8479 and D_7(3) = 0.4374, both inside the predicted ranges.
- The growth ratio per order is slowing: x2.20, x1.89, x1.74.
- **P4 (blind).** D_8(2) lies in [1.1, 1.8], so the threshold D = 1 is crossed at order 8, not 7. D_8(3) lies in
  [0.55, 0.95].
- **P5 (blind).** No S-row-zero gadget of order <= 8 has an odd 2x2 determinant. This holds at every order <= 7.

## What this would and would not mean

- **Confirmation would show** that the recursive composition of pieces of order <= 7 is near the critical
  thickness. Above it, a finite gadget alphabet fills a plane of states with positive density, and the
  H-projection is then expected to be complete with overwhelming redundancy.
- **It would not prove** the {7, 21} spectrum. Positive density is not exact filling; a proof would need an exact
  coset-aligned covering of a cone, or an analytic local-global argument.

## Results

(to be filled after the run)
