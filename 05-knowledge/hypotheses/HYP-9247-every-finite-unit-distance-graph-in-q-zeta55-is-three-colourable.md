---
id: HYP-9247
title: "Every finite unit-distance graph in Q(zeta55) is three-colourable (the field graph has chromatic number exactly three): triangle-free field, every explicit eleventh-root context to 38,514 vertices three-colourable, no additive colouring with fewer than 19 colours and no residue colouring at the smallest residue field (109) with fewer than 7; if true, no context in the field changes the six-colour relation at the eight fixed ports (THM-4614) and the six-colour program must leave the field or the ports"
status: >
  UPDATE 2026-10-16 (opus-2026-10-16-S38, THM-4661): the gate is now general. In every odd-degree carrier over a
  conjugation-stable field without sqrt(-3), a three-colourable exterior preserves every old six-colouring (any old set,
  not only the eight ports), and a killing exterior has at least twenty vertices. EI19 is excluded from every odd-degree
  extension of L; the fish needs the same point of C3 over F_4 in every odd-degree conjugation-stable carrier over K.
  Probes of the cubic carrier L(theta) with the width-three unit are three-colourable but Cartesian (no rational relation
  between its ten conjugates and the roots of unity), so they add no evidence beyond L.
  UPDATE 2026-10-13 (opus-2026-10-13-S25, THM-4628): the residue no-go is now PROVED at every conjugation-fixed prime, not only at 109:
  the residue graph is Cay(F_(q^2), N^1) with Kloosterman eigenvalues, so residue colourings need at least 1 + (q+1)/(2 sqrt q) colours
  (7 at 109; 12, 15, 16, 20, 21, 23, 24 at 439, 659, 769, 1319, 1429, 1759, 1979; at least 788, 2266, 27737 at the degree-ten
  residue fields of 19, 29, 79, previously unexamined); every fixed prime has even residue degree and q >= 109. The remaining
  content of this file (three-colourability of finite unit graphs in L; non-additive, non-residue rules) is unchanged.
  OPEN. PROVED: the lower bound 3 (unit five-cycles) and the additive no-go (additive colourings need a
  modulus whose primes are conjugation-fixed in L, the least being 19). NUMERICAL with margin: the residue
  graph at the smallest residue field, 109, is at least 7-chromatic by the Hoffman bound (floating-point
  eigenvalues, bound 6.60); the residue fields of 19, 29, 79 (degree ten) are unexamined. FINITE-EXACT:
  three-colourings of every tested context (8,436; 2,278; 10,598; 846; 1,304; 1,408; 4,974; 6,000; 38,514 vertices).
  Independently audited 2026-10-12 (the five default probes, the Hoffman bound and the prime criterion
  reproduced; the residue scope corrected; MISTAKE-593). VERIFIED (same day, addendum): the smallest known
  triangle-free four-chromatic unit graph EI17 has no placement in L with coordinate height below 1e63
  (LLL at 2900 digits); its exact field is OPEN, the published S20 certificate being refuted. PROVED (same day):
  EI19 has no faithful placement in L (forced angles with sqrt(-2) and sqrt(-35)); the fish graph's
  one-parameter family is OPEN.
  The conjecture is the upper bound 3. Its refutation needs a four-chromatic triangle-free unit-distance
  graph with coordinates in Q(zeta55); the known four-chromatic triangle-free unit graphs (Wormald, Chilakamarri,
  O'Donnell, Hochberg-O'Donnell) have not been placed in this field, and the repo's four-critical diamond chains
  use a 60-degree rotation and triangles, so they do not live in it.
session: opus-2026-10-12-S24
source: 05-knowledge/results/six_colour_core_doubling_20261012.md (section 4)
scripts:
  - 04-computation/experiments/six_colour_core_doubling_20261012.py (probes; + .json)
related:
  - THM-4614 (the doubling extension theorem that makes this conjecture decisive for six-colour forcing at the ports)
  - 05-knowledge/results/cyclotomic_unit_component_20261011.md (triangle-freeness; the unit-path component is the localization at conjugation-exchanged primes)
  - 05-knowledge/results/cyclotomic_palette_escape_20261011.md (C3: 109 is the minimal fixed-prime residue field)
  - 05-knowledge/results/rhombus_overlap_elliptic_20261011.md (G2: chi(Q(zeta5, sqrt(-11))) <= 5 by the residue colouring at 19)
  - 05-knowledge/results/eleven_thirds_boundary_20261012.md (the nonlinear three-colouring of the thirds lattice)
---

# HYP-9247 — every finite unit-distance graph in Q(zeta55) is three-colourable

## Statement

Let `L = Q(zeta55)` with the unit-distance graph `x ~ y` iff `(x-y)(bar x - bar y) = 1`. **Conjecture:** every finite subgraph is three-colourable; equivalently (de Bruijn-Erdos) `chi(L) = 3`.

## Why it matters here

THM-4614 proves that a finite context `G` in `L` containing the eight ports changes the six-colour boundary relation only if `G-B` is not three-colourable. So the conjecture, if true, closes the eleventh-root six-colour program at these ports completely; if false, its witness is exactly the first object the program needs.

## Evidence

* **Lower bound 3 (PROVED):** unit five-cycles (`alpha` powers) and the length-eleven closed walk `5 + 3v + 3 bar v = 0` of the thirds lattice.
* **Triangle-free (inherited, PROVED):** a unit equilateral triangle needs `zeta3`, which is not in `L`.
* **Every tested context is three-colourable (FINITE-EXACT):** the mixed-generator graph (8,436 vertices, 24,238 edges), the thirds context (2,278; 6,405), their union (10,598; 30,506), the apex context (846; 2,187) and four balls over the tenth and eleventh roots of unity, the ten Mersenne directions `(2-gamma^j)/(2-gamma^(-j))` and the thirds unit: radius two (1,304 and 1,408 vertices) and radius three (4,974 and 6,000 vertices; 14,224 and 17,170 complete unit edges), and the radius-four ball over the fifth roots of unity and the Mersenne directions (38,514 vertices, 143,614 complete unit edges; the declared `--large` run of the script). All colourings are verified on every complete unit edge, with the edges found by the two-prime filter plus the exact norm.
* **No cheap upper bound (PROVED for the additive part, NUMERICAL with margin for the Hoffman part):** an additive colouring `S -> Z/m` of the unit-path component `S` (a ring containing `1/p` for every prime `p` with `-1` outside the subgroup generated by `p` modulo 55, hence for 2, 3, 5, 7, 11, 13, 17) needs `m` supported on conjugation-fixed primes, the least being 19; the residue colouring at the conjugation-fixed prime with the smallest residue field, `ell = 109`, is the Cayley graph of `F_(109^2)` on its 110 norm-one residues, whose least eigenvalue `-19.64` gives `chi >= 7` by Hoffman (the residue fields `F_(ell^10)` of 19, 29, 79 are unexamined). So any colouring with at most six colours, and in particular a three-colouring, must be a non-additive rule, as the periodic rule `c(A,B) = c_5(A mod 5) + B mod 3` on the thirds lattice already is.

## Why it might fail

Triangle-free four-chromatic unit-distance graphs exist in the plane (Wormald 1979, 6,448 vertices, girth five; Hochberg-O'Donnell's 23-vertex fish, 1995, girth four; Chilakamarri's 47-vertex moth, 1995; O'Donnell's 40/56-vertex graphs and his arbitrary-girth construction, 2000; the current smallest girth-four examples have 17, 19 and 21 vertices, Exoo-Ismailescu 2016 and Soifer 2024; see MathWorld's tables, read 2026-10-12). Their coordinates involve specific rotation angles; whether some four-chromatic triangle-free unit graph has all coordinates in `Q(zeta55)` is the open question. `L` has infinitely many unit directions (every `u/bar u`), so it is a rich field, and vertex-transitive balls are a weak probe of chromatic number four, which in known examples is carried by delicate odd-cycle arrangements rather than by density.

## Cheapest decisive tests

* Place a known small four-chromatic triangle-free unit graph in `L`. **Done for EI17 (17 vertices, 31 edges, minimally rigid; [note](../results/ei17_placement_q_zeta55_20261012.md)):** in the gauge `v17 = 0`, `v11 = -1` its 31 edges carry 13 directions (nine unit 4-cycles are parallelograms), and integer LLL at 2900 digits on the `zeta55` power basis finds no relation for any non-gauge direction below the random height `1e69`; so no placement with coordinate height below about `1e63` exists (VERIFIED). The exact decision is OPEN: the published field certificate (arXiv:2607.19995, degree-20 polynomial with Galois group `S20`) does not vanish at the published coordinates (nearest root `8.8e-8` away) and cannot be cited. **EI19 (19 vertices) is excluded outright (PROVED, [note](../results/ei19_fish_placement_q_zeta55_20261012.md)):** two support-three relations among its direction classes force the ratios `(-1 +- 2 sqrt(-2))/3` and `(1 +- sqrt(-35))/6` in every faithful realization, and neither square root lies in `Q(zeta55)`. **The 23-vertex fish is a one-parameter family** (42 = 2*23 - 4 edges; 24 classes, 11 relations, no forced angle); whether its curve has an `L`-point is OPEN (no `L`-point at 117 `L`-unit parameter values). **EI21 (21 vertices, 38 = 2*21 - 4 edges; [note](../results/ei21_kappa_q_20261012.md)) is, like the fish, a one-parameter family** (18 classes, 8 relations, no relation of support three, two of support four); 117 `L`-unit parameter values give 524 faithful realizations with no further direction in `L`; OPEN. Next: the genus and `L`-points of the fish and EI21 curves; the exact six-angle elimination for EI17.
* Mycielski-type search: unit five-cycles in `L` with shadow vertices (the pentagram apexes are the shadows of the regular pentagon) and a common unit neighbour of the shadows; SAT over the exact unit graph of such configurations.
* An additive proof of three-colourability is impossible, and a residue-type proof is impossible with the smallest residue field (above; the residue fields of 19, 29, 79 are unexamined); a proof would need a global nonlinear rule, for example a rule on the exact localization `S` built from the periodic thirds rule and the `alpha`-power rule `a+b+c+d mod 3` on `Z[alpha]`.
