---
id: THM-4839
title: "Bridges into the triangle-group world: for N(P) >= 4, K4, K5, K7 are the only complete unit Cayley graphs of imaginary-quadratic residue fields; Petersen, K5, K3,3 are the smallest non-planar simple skeletons of regular maps of types {5,3}, {4,4}, {6,3} (hemi-dodecahedron; Gaussian {4,4}_(2,1); Eisenstein {6,3}_(1,1)), with K7 = {3,6}_(2,1), Heawood = {6,3}_(2,1) (chiral) and K4,4; Euclid's parametrisation squared is the Klein-four cover X(4) -> X(2), so Pythagorean Frey curves are the X(4)-subfamily of the Z/2 x Z/4 curves, with rational 4-torsion point (a(a+c), ac(a+c))"
status: PROVED (elementary) + CITED (Coxeter-Moser regular maps; Darmon-Granville 1995 section 3 and 7.1 for X(4) -> X(2); Kubert 1976; Halbeisen-Hungerbuehler 2022) + FINITE-EXACT (explicit isomorphisms; group law on 6000 Berggren cases; the audit re-checked all fundamental discriminants |D| <= 2000); independently audited SOUND WITH CORRECTIONS (MISTAKE-693)
session: opus-2026-10-28-S68
date: 2026-10-28
note: 05-knowledge/results/attractors_modular_triangle_world_abc_genealogy_20261028.md
script: 04-computation/experiments/attractor_bridges_regular_maps_belyi_frey_20261028.py (+ .out)
audit: 04-computation/experiments/attractors_census_bridges_20261028_audit_S68/AUDIT.md
related: [HYP-9437, THM-4838, THM-4837, THM-4825, THM-4808, THM-4510, THM-3333, THM-4529, THM-4701, THM-261, THM-4650]
---

# THM-4839 - bridges into the triangle-group world

**(a) Complete unit Cayley graphs (N(P) >= 4).** For an imaginary quadratic field K and a prime P of O_K with
N(P) >= 4, the Cayley graph of O_K/P with connection set the image of O_K^x is complete iff (K, N(P)) is one of
(Q(sqrt-3), 4), (Q(i), 5), (Q(sqrt-3), 7). These are K4, K5 and K7.
- For N(P) in {2, 3} the units +-1 always give K2 or K3.

*Proof.* Surjectivity onto (O_K/P)^x needs N - 1 <= |O_K^x| in {2, 4, 6}, so N in {4, 5, 7}.
- N = 4 needs 2 inert with units of order 3 mod P: only Z[omega].
- N = 5 needs 5 split and four units: Z[i].
- N = 7 needs 7 split and six units: Z[omega].

Checked over d in {1,2,3,7,11,19,43}, p < 60. The audit checked all 611 fundamental discriminants |D| <= 2000, p <= 200.

**(b) Regular maps (Coxeter-Moser's sense, convention verified by the audit).**
- Z[i]/(2 + i) with the unit steps is K5 = {4,4}_(2,1), and Z[i]/(2 + 2i) is K4,4 = {4,4}_(2,2).
- Z[omega]/(3 + omega) with the six unit steps is K7 = {3,6}_(2,1).
- The honeycomb modulo 2 + omega (norm 3) is K3,3 = {6,3}_(1,1). Modulo 3 + omega (norm 7) it is the Heawood graph =
  {6,3}_(2,1).
- {3,6}_(2,0) is K4 with doubled edges.
- The dodecahedron modulo the antipodal map is the Petersen graph, i.e. the hemi-dodecahedron {5,3}_5 on RP^2. Its map
  group [5,3]/{+-1} = A_5 = PSL_2(F_5) acts regularly on its 60 flags.
- So Petersen, K5 and K3,3 are the smallest non-planar simple skeletons of regular maps of types {5,3}, {4,4} and {6,3}.
  The dual types {3,5} and {3,6} give K6 (hemi-icosahedron) and K7.
- The {4,4}_(2,1), {3,6}_(2,1) and {6,3}_(2,1) maps are **chiral**; the others listed are reflexible.
- K5, K3,3, K7, Heawood and K4,4 fail the girth form of Euler's bound E <= g(V - 2)/(g - 2).

**(c) The Klein-four cover (CITED: classical).** lambda(t) = ((1 - t^2)/(1 + t^2))^2 is the Galois Belyi map of degree 4
with deck group {t -> +-t^(+-1)} and passport (2^2, 2^2, 2^2).
- It is the modular cover **X(4) -> X(2)**. With t = theta_4^2/(theta_3^2 + theta_2^2), lambda(t(tau)) = lambda(tau);
  and Gamma(2)-bar, free of rank 2, has a unique normal subgroup with quotient V_4, namely Gamma(4)-bar.
- It is the r = 2 dihedral cover of Darmon-Granville 1995, section 7.1, and p. 516.
- With t = n/m, Euclid's (a, b, c) = (m^2 - n^2, 2mn, m^2 + n^2) gives lambda = (a/c)^2 and 1 - lambda = (b/c)^2.

**(d) Pythagorean Frey curves.** For a^2 + b^2 = c^2, the curve y^2 = x(x - a^2)(x + b^2) contains
Q = (a(a + c), ac(a + c)), with 2Q = (a^2, 0), so Q has order 4.
- Proof: the halving criterion holds because a^2 and a^2 + b^2 = c^2 are squares; directly,
  y^2 = a(a+c) * ac * c(a+c).
- These are exactly the Legendre curves with lambda = k^2 and 1 - lambda = k'^2: the fibres over the rational points of
  X(4). They have Z/2 x Z/4 inside E(Q) and E[4] inside E(Q(i)).
- X(4) is a double cover of Kubert's Z/2 x Z/4 curve X_1(2,4), so these curves are a **thin subfamily** of the
  Z/2 x Z/4 curves. For example, y^2 = x(x+1)(x+4) has torsion exactly Z/2 x Z/4 and is not Pythagorean.
- CITED: Kubert 1976. Related: Halbeisen-Hungerbuehler, J. Number Theory 233 (2022), where every Z/2 x Z/8 curve
  arises from a pythagorean pair. THM-4650 has Z/2 x Z/4 torsion elsewhere in the repository.

**Measured companions (FINITE-EXACT).**
- 60 abc hits among the depth-<=9 Berggren triples (a^2, b^2, c^2); best q = 1.2659 at (36207, 18424, 40625).
- Every Collatz clock K/L (convergent or semiconvergent of log2 3) from 3/2 on is an abc hit.
  - Proof: for L >= 2, |L log2 3 - K| <= 0.2451 < log2(6/5). So |2^K - 3^L| < max/6, hence 6 rad(d) < max. The audit
    checked all 40 clocks with L <= 20000.
  - Baker makes log(max)/log(6|2^K - 3^L|) -> 1. "q -> 1" itself would need abc and is not known.

**Novelty.**
- The repository already had:
  - Petersen = hemi-dodecahedron (`noble_polyhedra_hill_fissary_20261002.md`, `odd_zeta_parallels_petersen_20261001.md`,
    `petersen_product_4polytope_attempt_20261001.md`, `pythagorean_form_conic_graphs_rp_margin_20261018.md`);
  - K7 as the chiral map {3,6}_(2,1) with Heawood dual (`camion_busch_gaps_polyhedra_collatz_20261002.md` section 6.4);
  - the N(P) = 7 case of (a) (`torus_voltage_and_halfplane_forcing_20261020/geometric_forcing_bridges.md`, G2).
- Literature pointers for (a): the Gaussian and Eisenstein-Jacobi networks of Martinez-Beivide et al. (IEEE Trans.
  Comput. 57, 2008; Probl. Inf. Transm. 44, 2008).
- New to the repository:
  - (a) for N(P) in {4, 5} as a classification;
  - the assembly of (b) by map type, with chirality;
  - the X(4) reading of Euclid's parametrisation linking Berggren, the square-sum chains and the abc lineages, typed
    CITED for the map itself;
  - (d) as stated, with its thin-subfamily correction.
