---
id: THM-4614
title: "In Q(zeta55) no unit-distance context with a three-colourable exterior changes the six-colour boundary relation at the eight fixed ports: the points at unit distance from two ports are exactly the ports and nine rhombus completions (five pentagram apexes, four pentagon-rhombus apexes), no point of the field outside the ports is adjacent to three ports, and the doubling lemma (three colour classes to three disjoint colour pairs; a blocking vertex kills six of the ninety options, so fewer than fifteen never block) extends every proper six-colouring of the ports to every finite G with G-B three-colourable; a forcing context therefore needs a four-chromatic subgraph of the triangle-free field graph, and none is known"
status: >
  PROVED: statements 2 and 3 (the doubling lemma and the extension theorem), given statement 1.
  PROVED + FINITE-EXACT: statement 1 (the census). Its completeness rests on the squareness criterion
  D(D-4) for the apexes of a pair at squared distance D; eleven pairs receive exactly constructed apexes
  verified at unit distance (nine of them supply the nine new points, two supply only ports), six pairs
  receive a quadratic non-residue witness at a split prime, the rest are unit edges (triangle-free field)
  or beyond distance two.
  FINITE-EXACT: statement 4 (the 570-pattern SAT confirmation and the five three-colourable probe contexts).
  NUMERICAL with margin: the Hoffman bound at the minimal fixed prime (floating-point eigenvalues; the bound
  6.60 exceeds 6 by 0.6). Independently audited 2026-10-12 (adversarial audit with its own exact arithmetic,
  edge finder and SAT: every census entry, witness prime, apex, pattern count and probe reproduced; the port
  quantifier and the residue-colouring scope corrected; CORRECT WITH FIXES, fixes applied, MISTAKE-593).
session: opus-2026-10-12-S24
source: 05-knowledge/results/six_colour_core_doubling_20261012.md
scripts:
  - 04-computation/experiments/six_colour_core_doubling_20261012.py (+ 05-knowledge/results/six_colour_core_doubling_20261012.json)
related:
  - 05-knowledge/results/eleven_thirds_boundary_20261012.md (the 570-pattern bare relation of the old 8,436-vertex context)
  - 05-knowledge/results/eleven_port_layers_20261011.md (P1: one K-port per new point; P2: one-layer contexts)
  - 05-knowledge/results/pentagon_context_number_fields_20261011.md (B1: five-colourable contexts of the twenty-port pentagon gadget)
  - 05-knowledge/results/cyclotomic_unit_component_20261011.md (the field graph is triangle-free)
  - 05-knowledge/results/cyclotomic_palette_escape_20261011.md (C3: 109 is the minimal fixed-prime residue field)
  - HYP-9247 (three-colourability of every finite unit-distance graph in Q(zeta55))
---

# THM-4614 — three-colourable contexts cannot change the six-colour port relation in Q(zeta55)

## Setting

`L = Q(zeta55)`, `alpha = zeta5`, `gamma = zeta11`, `eta = (-3+sqrt(-55))/8`. The unit-distance graph on `L` has an edge between `x` and `y` when `(x-y)(bar x-bar y) = 1`. The eight ports are

    B = {P_0,...,P_4} union {0, -eta, 3/4, -bar(eta)},   P_j = sum_(i<j) alpha^i,

the unit pentagon and the unit rhombus sharing the vertex `0 = P_0`; `G[B]` has nine edges. For a finite unit-distance graph `G` in `L` with `B` a subset of `V(G)`, `R_6(G;B)` is the set of proper six-colourings of `G[B]` that extend to proper six-colourings of `G`. There are 570 proper colourings of `G[B]` up to relabelling (30/180/247/113 using 3/4/5/6 colours) and 327,600 labelled ones.

## Statement

1. **(Census of double-port vertices.)** A point `x` of `L` not in `B` is at unit distance from two ports if and only if it is one of the following nine points, each the fourth vertex of a unit rhombus spanned by two port edges at a common port:
   * the five pentagram apexes `P_j + alpha^(j+1)` (`j = 0,...,4`), i.e. `alpha`, `1+alpha^2`, `1+alpha+alpha^3`, `1+alpha+alpha^2+alpha^4`, `2+alpha+alpha^2+alpha^3`, completing the rhombi `P_j, P_(j+1), P_(j+2)`;
   * the four pentagon-rhombus apexes `1-eta`, `1-bar(eta)`, `P_4-eta`, `P_4-bar(eta)`, completing the rhombi at `0` spanned by `P_1` or `P_4` with `-eta` or `-bar(eta)`.
   Each of the nine has exactly two port neighbours; no point of `L` outside `B` is at unit distance from three ports (the port `0` sees four); the nine are pairwise at non-unit distance. All other port pairs have no apex in `L`: the nine unit edges (an apex would close a unit triangle, impossible in `L`), the pairs `(0,3/4)` and `(-eta,-bar(eta))` whose apexes are ports, the pairs `(P_2,-eta)`, `(P_3,-eta)` at distance above two, and the six pairs `(P_1,3/4)`, `(P_2,3/4)`, `(P_3,3/4)`, `(P_4,3/4)`, `(P_2,-bar(eta))`, `(P_3,-bar(eta))`, for which `D(D-4)` is a quadratic non-residue at a split prime of `L` (witnesses 661, 661, 661, 2531, 991, 331).

2. **(Doubling lemma; any graph.)** Let `H` be a finite graph, `B` a subset of `V(H)`, `c` a proper 3-colouring of `H-B`, and suppose every vertex of `H-B` has at most two neighbours in `B`. Let `t` be a proper 6-colouring of `H[B]`. Call an *option* a perfect matching `M` of the six colours into three pairs together with a bijection `sigma` from the three classes of `c` to the pairs of `M` (90 options). A vertex `x` of `H-B` with two `B`-neighbours of distinct `t`-colours *blocks* the option when `sigma(c(x)) = t(N_B(x))`. If some option is blocked by no vertex, `t` extends to a proper 6-colouring of `H`. Each vertex blocks exactly six options, so fewer than fifteen blocking vertices never exhaust the options.

3. **(Extension theorem.)** Let `G` be a finite unit-distance graph in `L` with `B` a subset of `V(G)` and `G-B` three-colourable. Then every proper six-colouring of `G[B]` extends to `G`:

       R_6(G;B) = R_6(G[B];B).

4. **(Confirmations and the shape of the obstruction; FINITE-EXACT, except the eigenvalue bound, which is NUMERICAL with margin.)**
   * For each of the 570 patterns, no assignment of classes to the nine apexes kills all 90 options (SAT, 570 instances), even with all nine present at once.
   * The apex context (ports, the nine apexes, and one layer of 52 directions: the fifth and eleventh roots of unity, the ten Mersenne directions `(2-gamma^j)/(2-gamma^(-j))`, the thirds unit `(-5+sqrt(-11))/6`, and their negatives; 846 vertices, 2,187 unit edges) is three-colourable, peels to the eight ports at palette six, and has the full relation with counts 30/180/247/113.
   * Four balls of radius two and three over subsets of these directions (1,304 / 1,408 / 4,974 / 6,000 vertices; 2,534 / 2,738 / 14,224 / 17,170 complete unit edges) and the radius-four ball over the fifth roots of unity and the Mersenne directions (38,514 vertices, 143,614 edges; the declared `--large` run) are three-colourable.
   * At the conjugation-fixed prime of `L` with the smallest residue field, `ell = 109` (residue field `F_(109^2)`, 110 norm-one residues, triangle-free), the Cayley graph of the residue field on the norm-one residues has least eigenvalue `-19.64` against degree 110, so the Hoffman bound gives chromatic number at least 7: the residue colouring with the smallest residue field cannot give six colours, let alone show port freedom. The other conjugation-fixed primes (19, 29, 79, residue fields `F_(ell^10)`) were not examined.

## Proof

**1.** Let `p != q` be ports with `D = (q-p)(bar q-bar p)` (an element of the real subfield). The apexes of the pair, the points `x` of the plane with `|x-p| = |x-q| = 1`, exist iff `D <= 4` and are

    x = p + (q-p) (1/2 +- s/(2D)),      s^2 = D(D-4),

because `(x-p)/(q-p)` is the unit-modulus number `1/2 +- i sqrt(4-D)/(2 sqrt D)` and `-(4-D)/D = D(D-4)/D^2`. Since `q-p` and `D` lie in `L`, `x` lies in `L` iff `s` lies in `L`. The ports lie in `K = Q(alpha, sqrt(-11))`, `D` lies in the real subfield of `K`, and a square root of such an element that lies in `L` lies in `K`, because `L/K` has odd degree five while `K(s)/K` has degree at most two. So squareness is decided in `K`: a square root is proposed by an integer relation on the eight-term basis of `K` (PSLQ at 120 digits) and *verified exactly* by squaring in the 40-dimensional power basis; a non-square is certified by a split prime `ell = 1 mod 55` with a ring map `Z[zeta55][1/den] -> F_ell` at which the image of `D(D-4)` is a non-residue (if `D(D-4) = s^2` in `L`, then `s` is integral at that prime and its image squares to the image). The five pentagon pairs at squared distance `phi^2 = (3+sqrt5)/2` (where `D(D-4) = -(5+sqrt5)/2 = (alpha-alpha^(-1))^2`) and the four pentagon-rhombus pairs `(P_1,-eta)`, `(P_1,-bar eta)` (squared distance `5/4`, `D(D-4) = -55/16`) and `(P_4,-eta)`, `(P_4,-bar eta)` (squared distances `2+2Re(P_4 bar eta) = 3.995` and `2+2Re(P_4 eta) = 0.468`) receive both apexes; one of the two is a port and the other is the rhombus completion listed. The pairs `(0,3/4)` and `(-eta,-bar eta)` also receive both apexes, which are ports. The six excluded pairs receive the witnesses listed. The pairs at unit distance have `D(D-4) = -3`, not a square in `L` since `sqrt(-3)` is not in `L` (adjoining `zeta3` doubles the degree 40). The remaining two pairs have `D > 4`. The three further properties (two port neighbours each, none shared, no unit distance among the nine) are read off the complete unit graph on the seventeen points, computed exactly. ∎

**2.** Give each `x` in `H-B` a colour in `sigma(c(x)) \ t(N_B(x))`. This set is nonempty: `|t(N_B(x))| <= 2`, with equality only when `x` has two `B`-neighbours of distinct colours, and then it is empty only if `sigma(c(x)) = t(N_B(x))`, which is the blocked case. Two adjacent vertices of `H-B` lie in different classes of `c`, hence receive colours from disjoint pairs; an edge from `x` to a port `p` is proper because `t(p)` was removed from the list of `x`; edges inside `B` are proper because `t` is. For the count: a given pair `{a,b}` lies in exactly three of the fifteen matchings (the other four colours pair in three ways), and a given class is sent to a given pair of `M` by exactly two of the six bijections. ∎

**3.** By statement 1, every vertex of `G-B` has at most two port neighbours, and at most nine vertices of `G-B` have two. By statement 2 with a proper 3-colouring of `G-B`, at most `9*6 = 54 < 90` options are blocked, so some option extends `t`. ∎

**4.** Computation: `apex_census`, `doubling_blockable` (one SAT instance per pattern over the nine apexes, with the ninety option clauses), `fast_graph` (bucketed version of the engine's two-prime filter followed by the exact norm on every survivor; it is checked against the engine's pairwise finder on the seventeen marked points), `colourable` (SAT, verified on every edge), `PE.boundary_packet` (inherited relative-six-core peeling and the 570 lifts), and the eigenvalue control (character sums over the 110 norm-one residues of `F_(109^2) = F_109(sqrt d)`). ∎

## Why this is the obstruction, and what it reroutes

Every eleventh-root context filed so far has an explicit three-colouring: the 8,436-vertex mixed-generator graph, the 2,278-vertex thirds context, their 10,598-vertex union with 88 seam edges, the apex context and the four balls above. Statement 3 says their full boundary relations are not accidents of their generators or sizes: with the ports fixed, no three-colourable exterior can remove a pattern. The only lever a three-colourable exterior has is a double-port vertex, there are nine of them in the whole field, and fifteen would be needed.

A context that changes `R_6(G;B)` must therefore contain a four-chromatic finite unit-distance graph of `L`. The field graph is triangle-free, additive colourings of its unit-path component need at least 19 colours (every prime below 19, and 2, 3, 5, 7, 11, 13, 17 in particular, divides some unit-direction denominator), and the residue colouring with the smallest residue field is at least 7-chromatic; so any colouring of `L` with at most six colours is a non-additive rule, and HYP-9247 conjectures that no four-chromatic subgraph exists. The general form of statement 2 (palette `m`, exterior chromatic number `k`, at most `r` port neighbours per exterior vertex: unconditional extension when `m >= k(r+1)`, and the blocking analysis at `m = 2k`, `r <= 2`) applies to any port set once its double-port census is known.
