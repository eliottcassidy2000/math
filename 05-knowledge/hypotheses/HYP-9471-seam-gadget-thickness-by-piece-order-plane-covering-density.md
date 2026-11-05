---
id: HYP-9471
title: "Seam-gadget semigroup thickness by piece order. On the plane P = {S = 0} of SEAM states (H, S, R) whose exit starts no Hamiltonian path, the S-row-zero gadgets act by 2x2 integer maps [[a, b], [c, d]] on (H, R). The covering density D_k(t) = sum over the distinct nonsingular such maps of order <= k of (a + b t)^(-2) is 0.117, 0.257, 0.486 at t = 2 for k = 4, 5, 6, roughly doubling per order. D_k(t0) >= 1 at the dense end t0 of a cone is a weak necessary condition for a finite-alphabet self-covering of that cone; the sharp counting condition is much stronger (audit S11). Pre-registered: D_7(2) lies in [0.75, 1.25], and the order-7 alphabet thickens the realised P-slice"
status: "PARTLY CONFIRMED: P1, P2, P4, P5 CONFIRMED; P3, P6 REFUTED. Pre-registered in a98199e8bb (P1-P3) and 973c1f8312 (P4-P5), both under the colliding name HYP-9470, and in 7297af7bf9 (P6). Renumbered HYP-9471 because opus-S79 reserved HYP-9470 first, in 395a317b0e. The density D_8(2) passes 1 at slope 2, but the realised plane orbit still thins with scale at order 8: the density averaged over its slopes is 0.882. The sharp counting condition for a self-covering cone is <= 0.27 at order 8 (audit S11). Audited PASS WITH CORRECTIONS (MISTAKE-717)"
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

  The realised P-states have R/H >= 2, with equality attained, in the runs below, so t = 2 is the dense end of the cone.
  (Corrected after the audit; first written "R/H >= 2.01".)
- **Seam BFS, order <= 6** (1,941 matrix lines, 1,884 distinct matrices; states with H, R <= X expanded; seeds =
  small realisable triples):
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

**Addendum 2 (committed after P1-P5 were read, before any pure-plane run).**
- **P2, P3 and P5 results.**
  - P2: the order-7 fill is 0.4573 (CONFIRMED).
  - P3: the order-7 ratio is 0.9008 against 0.9051 (REFUTED).
  - P5: no odd determinant among the 22,522 distinct S-row-zero matrices of order <= 8 (CONFIRMED).
- **P4 result.** D_8(2) = 1.4288 and D_8(3) = 0.7395 (CONFIRMED). D_8(2.5) = 0.9999.
- **Why P3 failed.** In the run used for P2 and P3, the plane is fed by the whole 3D orbit through rank-2 maps. Its
  decay with scale then tracks the 3D orbit's exponent: a ratio of 0.90 corresponds to delta of about 1.83, close to
  the measured 1.85. It does not track the 2D density D_k. The test of the 2D heuristic is the pure plane dynamics.
- **P6 (blind).** `plane_bfs.c` runs the BFS of (H, R) under the 2x2 restrictions alone, with the twelve realisable
  S = 0 seeds of order <= 5, X = 5,000, and the same two bands. Predictions:
  - (a) the order-<=8 ratio fill[1601, 1799] / fill[801, 999] exceeds the order-<=7 ratio;
  - (b) the order-<=7 ratio is below 0.95 and the order-<=8 ratio is at least 0.95.
  - Both are read off one run per alphabet (orders <= 6, <= 7 and <= 8).

## Results

| prediction | outcome | data |
|---|---|---|
| P1: D_7(2) in [0.75, 1.25], D_7(3) in [0.38, 0.65] | CONFIRMED | D_7(2) = 0.8479, D_7(3) = 0.4374 |
| P2: order-7 slice fill >= 0.396 | CONFIRMED | 0.4573 (order 6: 0.2640) |
| P3: order-7 decay ratio > 0.9051 | REFUTED | 0.9008 |
| P4: D_8(2) in [1.1, 1.8], D_8(3) in [0.55, 0.95] | CONFIRMED | D_8(2) = 1.4288, D_8(3) = 0.7395, D_8(2.5) = 0.9999 |
| P5: no odd 2x2 determinant to order 8 | CONFIRMED | 0 of 22,522 distinct S-row-zero matrices; M_(z,u)(G) even in all 20,787 order-8 ones |
| P6 (a): pure-plane ratio, order 8 > order 7 | REFUTED | 0.7174 vs 0.7175 |
| P6 (b): order-7 ratio < 0.95 <= order-8 ratio | REFUTED | (b) fails at order 8 |

Order 8 has 6,880 classes and 427,412 distinct transfer matrices. The cumulative distinct counts by order are
1, 5, 16, 56, 264, 1,884, 22,243, 448,915. The S-row-zero ones are 0, 1, 5, 16, 57, 267, 1,919, 22,522, of which
0, 0, 1, 5, 29, 191, 1,632, 21,023 restrict to nonsingular 2x2 maps.

| k | D_k(2) | D_k(2.5) | D_k(3) | D_k(4) |
|---|---|---|---|---|
| 4 | 0.1170 | 0.0812 | 0.0596 | 0.0361 |
| 5 | 0.2571 | 0.1787 | 0.1315 | 0.0798 |
| 6 | 0.4861 | 0.3388 | 0.2498 | 0.1521 |
| 7 | 0.8479 | 0.5923 | 0.4374 | 0.2668 |
| 8 | 1.4288 | 0.9999 | 0.7395 | 0.4519 |

The growth ratio per order at t = 2 is 2.20, 1.89, 1.74, 1.69.

**Pure plane dynamics (P6)** (`plane_bfs.c`, X = 5,000, twelve S = 0 seeds):

| alphabet | maps | states | fill[801, 999] | fill[1601, 1799] | ratio |
|---|---|---|---|---|---|
| order <= 6 | 191 | 20,149 | 0.0256 | 0.0222 | 0.8691 |
| order <= 7 | 1,632 | 46,521 | 0.0615 | 0.0442 | 0.7175 |
| order <= 8 | 21,023 | 124,006 | 0.1506 | 0.1081 | 0.7174 |

The fill rises 2.4-fold per order, but the decay with scale does not soften.

**Why the heuristic failed (exploratory, after the fact).**
- The order-8 plane states with H in [801, 1799] number 54,994.
- Their slopes have median 2.64 and 90th percentile 4.15, with histogram 2: 36,235, 3: 12,363, 4: 5,221,
  5: 1,150, 6: 25.
- D_8 averaged over those realised slopes is 0.8819 < 1. The audit's exact mean is `aud_slope_variants.py`; the first
  version subsampled the slopes and gave 0.8825. The threshold is crossed only at the dense end
  (t <= 2.5), and the maps push states toward larger slopes, where the density is lower. So the realised plane orbit
  is still subcritical at order 8.
- Continuing the per-order ratio (about 1.6-1.7) would put the slope-averaged density above 1 at order 9.
  Order 9 (191,536 classes, 81 port pairs) was not computed.

**Verdict.** The thickness of the seam-gadget semigroup grows steadily with piece order. Measured by D_k it
roughly doubles per order, and it crosses the 2D covering threshold at the dense end of the cone at order 8. No
alphabet of order <= 8 yet makes the realised plane orbit stop thinning. A finite self-covering certificate for the
{7, 21} spectrum on this plane is therefore not available at order 8. Even above the threshold, exact coset alignment
would still be needed: every plane map has even determinant (THM-4879 (E), P5).

## Audit corrections (audit S11, PASS WITH CORRECTIONS; MISTAKE-717)

Audit report: `04-computation/experiments/seam_gadget_thickness_20261104_audit_S11/AUDIT.md`. Every pre-registered
number reproduced exactly with the audit's own code, and git confirms the registration order.

- **P1 blindness.** The order-7 matrix list existed when P1 was registered (its count is quoted above). P1's
  blindness therefore rests on this file's statement that D_7 had not been computed. git cannot prove it.
- **The threshold framing.** D_k(t) >= 1 at one slope is necessary but weak. The sharp necessary condition for
  self-covering the lattice points of a cone t0 <= R/H <= t1 comes from counting points with H <= X:

      (1/(t1 - t0)) Σ_i ∫ (a_i + b_i t)^(-2) dt >= 1,

  where the integral runs over the slopes t in [t0, t1] whose image slope A_i.t also lies in [t0, t1].
  - For the order-<=8 alphabet this is at most 0.2634 on all 40 cones tested with 2 <= t0 <= 4
    (`aud_cone_condition.out`). So order 8 misses it by about a factor of 4.
  - "Crosses 1 at order 8" therefore means only that the weakest condition now holds at the dense end. This agrees
    with P6 and with the session's own slope-preserving sub-alphabet check.
- **Counts.** "1,941 matrices" was a count of matrix lines including cross-order repeats; 1,884 are distinct. The
  results do not depend on the difference.
- **P5 is now a theorem.** THM-4879 (F) proves M_(z,u)(G) even for every gadget that maps the plane into itself, at
  every order. Every plane map therefore has even determinant.
