---
id: THM-4509
title: "Friedman's fences: the Corner Lemma. At every junction, ends + through-fences + sectors >= pi is at least 3, with equality exactly at T, X, Y and L junctions. With the face angle identity, the sum over fields of (convex corners - 3) is at most n - 3. With polygon isoperimetry, A(n) + 2 mu sqrt(pi A(n)) <= 2(mu+rho) n - 6 rho, so lim A(n)/n <= (6 - P_5)/(8 - P_5) = 0.5224525, an LP optimum made of regular pentagons and unit squares. Angle-potential certificates cannot go below (4 - 12^(1/4))/(6 - 12^(1/4)) = 0.5167670; for convex fields they give 0.5168084. Beating the grid's 1/2 needs fields with at most 4 and fields with at least 5 convex corners together"
status: >
  PROVED + INDEPENDENTLY AUDITED (junction inequality, face angle
  identity, Corner Lemma, the finite bound and lambda <= 0.5224525,
  Theorem B1, Proposition C1); PROVED computer-assisted (Theorem B3:
  convex fields, lambda <= 0.5168084); FINITE-EXACT (exact configurations;
  all 48 records n = 3..50 on Friedman's page satisfy the finite bound);
  EMPIRICAL (the non-convex potential LPs return to 0.5225; relaxations
  and Cairo-type searches never beat 1/2); CITED (L'Huilier / Zenodorus
  polygon isoperimetry; Friedman's table). OPEN: whether lambda = 1/2.
  Setting: THM-4505 (6). n unit fences; two fences meet in at most one
  point, which is an end of at least one; every end lies on another
  fence; every field (bounded face) has area <= 1; A(n) is the supremum
  of the total field area; lambda = lim A(n)/n = sup A(n)/n.
  Junction data: e ends, t in {0,1} through-fences, sectors of angles
  summing to 2 pi, r = #{sectors >= pi}, k = #{sectors < pi},
  alpha = 2 - r.
  (1) Junction inequality: 3 alpha - k <= e, i.e. e + t + r >= 3. It is
  tight exactly for T (one stem), X (two stems from opposite sides of a
  through-fence), Y (three ends, all sectors convex) and L (two ends,
  bent).
  (2) Face angle identity. A field with kappa convex corners and h holes
  has sum_convex theta + sum_reflex (theta - pi) = pi (kappa - 2 + 2h).
  Globally, 2C = sum_j (k_j - alpha_j) + 2H + 2c_o.
  (3) Corner Lemma: sum over fields of (kappa_f - 3)
  = sum_j (3 alpha_j - k_j)/2 - kappa_o - 3H - 3c_o <= n - 3.
  Equality holds for the unit square and the n = 5 and n = 6 records.
  (4) Every field satisfies a <= mu P_kappa sqrt(a) + 2 rho (kappa - 3),
  with P_kappa = 2 sqrt(kappa tan(pi/kappa)), mu = 1/(8 - P_5),
  rho = (4 - P_5)/(2(8 - P_5)).
  (5) Hence A + 2 mu sqrt(pi A) <= 2(mu+rho) n - 6 rho for every
  configuration: A(4) <= 1.0768, A(12) <= 4.3661, A(50) <= 22.0163 (5-10%
  below THM-4505's F2), and lambda <= (6 - P_5)/(8 - P_5) = 0.5224525.
  This improves 1/sqrt(pi) = 0.5642 and Hales's honeycomb value
  12^(-1/4) = 0.5373.
  (B1) Every angle-potential certificate (mu, rho, Phi(theta)), which
  generalizes the dual of (3)-(5), has 2(mu+rho) >= (4 - 12^(1/4))/(6 -
  12^(1/4)) = 0.5167670. So this method cannot prove lambda = 1/2. The
  forcing constraints are T(60,120), T(90,90), the unit square, the
  regular hexagon of area 1 and the tiny triangle. The LP attains the
  value with a fractional tiling (hexagons at triangular roundabouts plus
  squares) that fence lengths forbid.
  (B3, computer-assisted) If every field has a convex outer boundary, then
  A + 2 mu sqrt(pi A) <= 0.5168084 n.
  (C1) If every field has at most 4 convex corners, then A < n/2. If every
  field has at least 5, then A <= (n-3)/2. So a configuration beating 1/2
  needs both kinds. Pythagorean fence tilings have density < 1/2, and the
  Cairo tiling has no fence decomposition.
source: orchestrator derivation (junction inequality, Corner Lemma, the 0.5224525 bound as an LP optimum and its dual), 2026-09-26, answering the owner's prompt to relate Friedman's fence problem to the repository; independently audited and extended by the fencelim lane (X junctions, B1, B3, C1-C6) of session collatz-procgen-20260922; promoted by the session orchestrator after a second, exact audit
depends_on:
  - 01-canon/theorems/THM-4505-small-graphs-encoding-arithmetic-square-sum-zigzag-fences-collatz-alphabet.md (setting, Proposition F1, Theorem F2)
related:
  - 05-knowledge/results/procgen_kuratowski_20260925_tait_kempe_triples.md (Euler/Kuratowski counting; the Corner Lemma is its angle form)
  - 01-canon/theorems/THM-4486-min-max-cycle-density-game.md (LP-duality certificates; ANALOGY only)
note: 05-knowledge/results/procgen_fencelim_20260926_fence_density.md
scripts: 04-computation/experiments/procgen_fencelim_20260926_{run,geom,audit,records,potlp,typed}.py
script_audit: 04-computation/experiments/procgen_fencelim_20260926_orchestrator_check.py
output: 05-knowledge/results/procgen_fencelim_20260926.out
output_sha256: 7e0c605b6823eddebf544f41326eaaffe7ae942dd74df8aec5659fbac5d129f9
output_audit: 05-knowledge/results/procgen_fencelim_20260926_orchestrator_check.out
output_audit_sha256: 7bbe2d5833f225ad0da38c448a160ff96f58dbdd0948106c82a68ebc2b0d4815
hash_basis: raw bytes
audit: >
  Two audits.
  First, the fencelim lane wrote complete proofs of (1)-(5), covering
  straight joins, several stems on a fence, bridges, pinches, holes and
  several components. It found that X junctions are tight too (missing
  from the orchestrator's list; the bound is unaffected). It checked the
  identities exactly on 26 configurations plus the reproduced n = 6
  record.
  Second, the orchestrator's independent code
  (procgen_fencelim_20260926_orchestrator_check.py; its own exact
  half-edge face computation in rational coordinates; the lane's geometry
  code was not read) confirms:
  - the tight junction types {T, X, Y, L};
  - the global corner identity, the Corner Lemma and Proposition F1's
    field count on 9 configurations, including a Y, an X, T's, a hole and
    two components; equality holds for the unit square, the square with a
    chord and the n = 5 record;
  - B(4), B(12), B(50) and the 16 records quoted in THM-4505;
  - Theorem B1's five-constraint LP optimum;
  - C1's ingredients;
  - a Monte Carlo spot check of the B3 certificate. The lane's Phi was
    obtained as data from its LP, and this script's own constraint code
    tested it at 200000 random real junction and convex-face
    configurations (min slacks 2.2e-8 and 6.0e-9).
  The lane's runner was re-run (20.5 s, 19 [OK], 0 failures). Its output
  is identical to the committed .out up to the timing line.
  Scope notes:
  - B3's exact reduction to grid points (cells as polytopes with grid
    vertices; tangent-plane and concavity bounds for the faces) was read
    and spot-checked but not independently re-implemented.
  - B3 needs every field to have a convex outer boundary.
  - The non-convex obstruction (a notched regular pentagon, scored by its
    hull) is OPEN.
---

# THM-4509 — Friedman's fences: the Corner Lemma and the density bound

**PROVED + INDEPENDENTLY AUDITED.** Full note: [procgen_fencelim_20260926_fence_density](../../05-knowledge/results/procgen_fencelim_20260926_fence_density.md).

## 1. The discrete half: corners are counted by Euler's formula

In a fence configuration, fields get corners only where fences end.
- A junction where one fence ends on another (T) gives two corners.
- Three ends meeting (Y) give three.
- Every junction satisfies `ends + through + (sectors >= pi) >= 3`, with equality exactly at T, X, Y and L.

Summing the angle identity of every face turns this into the **Corner Lemma**: the fields together have at most `3·#fields + n - 3` convex corners. It is the angle form of THM-4505's Euler field count, `#fields = n + c - V0`.

## 2. The continuous half: isoperimetry prices each corner

A field with `kappa` convex corners and area `a` needs perimeter at least that of the regular `kappa`-gon of area `a`. Pricing perimeter at `mu` and corners at `2 rho` gives the per-field inequality. Summing it with the Corner Lemma gives `lambda <= 0.5224525`.

The optimum of this LP mixes regular pentagons and unit squares:
- pentagons are the cheapest way to beat the square's perimeter-to-area ratio;
- the corner budget allows about 91% of the fields to be pentagons.

Regular pentagons do not tile, so the truth is strictly lower. Any finer pricing that looks only at angles (Theorem B1) still stalls at `0.5167670`. For convex fields it reaches `0.5168084` (B3).

## 3. Where 1/2 stands

A pattern beating the grid must mix fields with at most 4 and at least 5 convex corners (C1). Every construction tried stays below 1/2:
- Pythagorean pinwheels;
- Cairo-type pentagons, which cannot be cut into unit fences.

Deciding `lambda = 1/2` needs accounting at the level of individual fences: each side of each fence is cut into face sides that sum to exactly 1. That is what the potentials cannot see.
