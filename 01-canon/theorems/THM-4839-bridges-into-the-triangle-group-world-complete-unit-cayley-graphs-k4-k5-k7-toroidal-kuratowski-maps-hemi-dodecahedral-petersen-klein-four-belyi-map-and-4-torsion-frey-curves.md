---
id: THM-4839
title: "Bridges into the triangle-group world: K4, K5, K7 are the only complete unit Cayley graphs of imaginary-quadratic residue fields; K5 = {4,4}_(2,1), K3,3 = {6,3}_(1,1), K7 = {3,6}_(2,1), Heawood = {6,3}_(2,1), K4,4 = {4,4}_(2,2) from Gaussian and Eisenstein ideals and Petersen = the hemi-dodecahedron, so the owner's {Petersen, K5, K3,3} are regular-map quotients of the triangle groups (2,3,5), (2,4,4), (2,3,6); Euclid's parametrisation squared is the Klein-four Belyi map onto the Legendre line and every Pythagorean Frey curve has the rational 4-torsion point (a(a+c), ac(a+c))"
status: PROVED (elementary) + CITED (Coxeter-Moser regular maps; Euclid; Frey-Hellegouarch; Kubert 1976) + FINITE-EXACT (explicit isomorphism checks, group-law checks on 6000 Berggren cases); independently audited (see the note's section 7)
session: opus-2026-10-28-S68
date: 2026-10-28
note: 05-knowledge/results/attractors_modular_triangle_world_abc_genealogy_20261028.md
script: 04-computation/experiments/attractor_bridges_regular_maps_belyi_frey_20261028.py (+ .out)
related: [HYP-9437, THM-4838, THM-4837, THM-4825, THM-4808, THM-4510, THM-3333, THM-4529, THM-4701, THM-261]
---

# THM-4839 - bridges into the triangle-group world

**(a) Complete unit Cayley graphs.** Let K be imaginary quadratic and P a prime of O_K with N(P) >= 4. The Cayley graph
of O_K/P with connection set the image of O_K^x is complete iff (K, N(P)) is one of
(Q(sqrt-3), 4), (Q(i), 5), (Q(sqrt-3), 7). These are K4, K5 and K7.

*Proof.* Completeness means the units surject onto (O_K/P)^x, which needs N - 1 <= |O_K^x| in {2, 4, 6}. So
N in {4, 5, 7}.
- N = 4 needs 2 inert with units of order 3 mod P: only Z[omega], since 2 ramifies in Z[i].
- N = 5 needs 5 split and four units: Z[i], where i maps to an element of order 4.
- N = 7 needs 7 split and six units: Z[omega], where -omega has order 6.
- If O_K^x = {+-1}, the image has at most 2 elements.

Checked over d in {1,2,3,7,11,19,43}, p < 60.

**(b) Regular maps and the owner's trio.**
- Gaussian ({4,4}): Z[i]/(2 + i) with the unit steps is K5 = {4,4}_(2,1), and Z[i]/(2 + 2i) is K4,4 = {4,4}_(2,2).
- Eisenstein ({3,6} and its dual {6,3}):
  - Z[omega]/(3 + omega) with the six unit steps is K7 = {3,6}_(2,1);
  - the honeycomb modulo 2 + omega (norm 3) is K3,3 = {6,3}_(1,1);
  - the honeycomb modulo 3 + omega (norm 7) is the Heawood graph = {6,3}_(2,1).
- K4 = {3,6}_(2,0), with every edge doubled because the six units collapse in pairs mod 2.
- The dodecahedron modulo the antipodal map is the Petersen graph, the hemi-dodecahedron on RP^2. Its map group is
  A_5 = PSL_2(F_5), a quotient of the (2,3,5) triangle group.
- Every identification is an explicit graph isomorphism (with Kneser K(5,2) and the Fano incidence graph). K5, K3,3,
  K7, Heawood and K4,4 fail the girth form of Euler's bound E <= g(V - 2)/(g - 2).
- Hence Petersen, K5 and K3,3 are regular-map quotients of the triangle groups (2,3,5) (spherical), (2,4,4) and (2,3,6)
  (Euclidean: the Gaussian and Eisenstein grids of THM-4825).

**(c) The Klein-four Belyi map.** lambda(t) = ((1 - t^2)/(1 + t^2))^2 is a Galois Belyi map of degree 4 with deck group
{t -> +-t^(+-1)} = Z/2 x Z/2 and passport (2^2, 2^2, 2^2) over (0, 1, oo).
- With t = n/m, Euclid's triple (a, b, c) = (m^2 - n^2, 2mn, m^2 + n^2) gives lambda = (a/c)^2 and
  1 - lambda = (b/c)^2.
- So the abc triples (a^2, b^2, c^2) attached to primitive Pythagorean triples are the rational points of the Legendre
  line in the image of this cover.

**(d) 4-torsion.** For a^2 + b^2 = c^2, the Frey curve y^2 = x(x - a^2)(x + b^2) contains Q = (a(a + c), ac(a + c)), with
2Q = (a^2, 0) and hence 4Q = O.
- Proof: the halving criterion for 2-torsion points applies because a^2 - 0 = a^2 and a^2 + b^2 = c^2 are squares.
  Direct check: y^2 = a(a+c) * ac * c(a+c) = (ac(a+c))^2.
- These curves form the genus-zero family with torsion containing Z/2 x Z/4 (Kubert 1976). Berggren's tree therefore
  parametrises a modular family of Frey curves.

**Novelty.**
- The components of (b) are classical (Coxeter-Moser, ch. 8) and (d) is classical.
- (a) as stated is a small observation, not found in the repository.
- New to the repository: the assembly as a bridge between the finite-geometry and modular hubs measured by HYP-9437,
  and the Klein-four reading of Euclid's parametrisation linking Berggren, the square-sum chains and the abc lineages.
- Measured companions (FINITE-EXACT):
  - 60 abc hits among the depth-<=9 Berggren triples (a^2, b^2, c^2), best q = 1.2659 at (36207, 18424, 40625);
  - every Collatz clock K/L from 3/2 to 84/53 is an abc hit, with 27/17 at q = 1.086.
