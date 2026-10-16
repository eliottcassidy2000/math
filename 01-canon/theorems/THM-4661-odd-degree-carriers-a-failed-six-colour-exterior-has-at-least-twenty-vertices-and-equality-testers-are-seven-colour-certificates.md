---
id: THM-4661
title: "Odd-degree carriers and the six-colour boundary: a vertex-minimal exterior that deletes an old six-colouring is L-critical with lists of size five or six, its low vertices span a forest of edges and odd cycles of length at least five, it is not three-colourable and not five-choosable, and it has at least twenty vertices (nineteen exterior vertices always preserve every old six-colouring; the earlier bound was seventeen; an exhaustive SAT search over the necessary conditions raises this to twenty-two, and every twenty-two-vertex model it finds is five-choosable); a six-colour equality tester is exactly six-chromatic, cannot live in a five-colourable field plane, and in any field decides the plane at seven colours; EI19, the fish and EI21 gain nothing from odd-degree extensions of Q(zeta5, sqrt(-11)) or Q(zeta55), in particular from the cubic carrier Q(zeta55, theta), and EI21 carries the fish's genus-one curve, so one question C3(F_4) decides the necessary condition for both"
status: >
  PROVED: statements 1-4 (elementary, given the CITED classical inputs: the Borodin / Erdos-Rubin-Taylor
  degree-choosability theorem, Riemann-Roch on a genus-one curve, the de Bruijn-Erdos compactness theorem and the
  Moser spindle principle; and the INHERITED repo inputs: Lemma 1 and Corollary 2a of hn_elliptic_bridge_20261016c,
  the forced EI19 angles, equations (1) and the genus-one curve C3 of the fish note).
  FINITE-EXACT: statement 5 (CaDiCaL verdicts over all attachment patterns, not proof-checked; Alon-Tarsi
  coefficients computed exactly). The counting of statement 2(d) is also re-derived by the n = 19 run of the same search.
  Self-audited; independent audit pending (S38).
session: opus-2026-10-16-S38
source: 05-knowledge/results/six_colour_odd_carrier_twenty_vertices_20261016.md
related:
  - 05-knowledge/results/hn_elliptic_bridge_20261016c.md (clock_duality_oct15: Lemma 1, Theorem 2 = the seventeen-vertex bound, Corollary 2a, inequality (2); statement 1 improves Theorem 2)
  - 05-knowledge/results/elliptic_crossedges_20261016c.md (clock_duality_oct15: G7, the cubic carrier K' = K(theta), L' = L(theta))
  - 05-knowledge/results/cubic_crossedge_packets_20261016d.md (clock_duality_oct15: C1 chi(K') <= 5, C2 the 64-term cubic circuit; received during S38)
  - THM-4614 (the doubling lemma behind the three-colour gate)
  - HYP-9247 (every finite unit graph in Q(zeta55) is three-colourable)
  - 05-knowledge/results/existence_count_amplifiers_20261011e.md (the equality-tester target e > 0, d = 0)
  - 05-knowledge/results/constructive_interfaces_and_bridges_20261011.md (the six/seven table)
  - 05-knowledge/results/ei19_fish_placement_q_zeta55_20261012.md (forced EI19 angles)
  - 05-knowledge/results/ei21_kappa_q_20261012.md (EI21: one-parameter family, left OPEN in S24; statement 4(c) reduces it to C3(F_4))
  - THM-4662 (C3(F_4) is empty: the fish and EI21 are excluded from Q(zeta55) and from every odd carrier over K)
  - 05-knowledge/results/fish_marked_elliptic_cover_20261013.md (the fish's square-class equations (1), the curve C3, its degree-eight divisor)
---

# THM-4661 — odd-degree carriers: twenty exterior vertices, the three-colour gate, and equality testers

## Setting

`F` is a subfield of `C` with `bar F = F` and `sqrt(-3)` not in `F`; `M` is a finite extension of `F` of **odd**
degree (not necessarily conjugation-stable). Points are elements of `M`, two points are adjacent when they are at
distance one. An **old** set `A ⊂ F` (possibly infinite) carries a proper `q`-colouring `phi` of its unit graph; an
**exterior** set `B ⊂ M \ F` is finite. `phi` **extends to** `B` if some proper `q`-colouring of the unit graph on
`A ∪ B` (all old-old, old-new and new-new unit pairs) agrees with `phi` on `A`.

Inherited facts (hn_elliptic_bridge_20261016c, Lemma 1 and Corollary 2a; clock_duality_oct15):

- **(I1)** a point of `M \ F` has at most one unit neighbour in `F` (a common unit neighbour of two points of `F` is
  quadratic over `F`, and `[M:F]` is odd);
- **(I2)** the unit graph on `M` is triangle-free (`sqrt(-3)` would lie in an odd-degree extension of `F`);
- **(I3)** two distinct points have at most two common unit neighbours (two unit circles).

For a finite `H ⊂ B` put `t(b) = |N_A(b)|` (0 or 1 by I1) and `L(b) = [q] \ phi(N_A(b))`, so `|L(b)| = q - t(b)`.
Then `phi` extends to `H` exactly when the unit graph on `H` is `L`-colourable. A vertex `b` of `H` is **low** if
`deg_H(b) = |L(b)|` and **high** if `deg_H(b) > |L(b)|`.

## Statement

1. **(Exterior bound.)** If `q = 6` and `|B| <= 19`, then every proper old six-colouring extends to `B`.
   For every `q >= 4`, `|B| <= q + ceil(q(q-2)/2)` suffices (for `q = 6` this simple form gives 18).
2. **(Minimal obstructions.)** Let `H ⊂ B` be vertex-minimal such that `phi` does not extend to `H`. Then
   (a) `H` is connected, not `L`-colourable, every proper induced subgraph is `L`-colourable, and
       `deg_H(b) >= |L(b)| >= q - 1` for every `b`;
   (b) every component of the subgraph spanned by the low vertices is a tree of blocks each of which is an edge or
       an odd cycle of length at least five; in particular `H` has a high vertex;
   (c) `H` is not `floor(q/2)`-colourable (for `q = 6`: not three-colourable) and not `(q-1)`-choosable;
   (d) for `q = 6`: `|H| >= 20`, `H` has minimum degree at least five, every vertex without an old neighbour has
       degree at least six, and some vertex has degree at least six.
3. **(Equality testers.)** Let `G` be a finite unit-distance graph (edges are unit pairs; further unit pairs may be
   present or absent) with non-adjacent vertices `x, y` such that `G` is six-colourable and every proper
   six-colouring gives `x` and `y` the same colour (an **equality tester**; existence count `e > 0`, separating count
   `d = 0` in the notation of existence_count_amplifiers_20261011e).
   (a) `chi(G) = 6`.
   (b) A finite unit-distance graph that is not six-colourable exists; hence `chi(R^2) = 7`. Conversely, if
       `chi(R^2) = 7` then an equality tester exists (at distance one).
   (c) No equality tester has all its vertices in a field plane whose unit graph is five-colourable; in particular
       none lies in `K = Q(zeta5, sqrt(-11))` (rhombus_overlap_elliptic_20261011.md G2, residue colouring at 19) or in
       the cubic base `K' = K(theta)` (clock_duality_oct15's C1 in cubic_crossedge_packets_20261016d: the same
       residue colouring extends because the cubic has a simple root modulo 19).
   (d) In the setting above with `q = 6`: if `x, y ∈ F` and the unit graph on `A = V(G) ∩ F` has a proper six-colouring
       that separates `x` and `y` (for instance if it is five-colourable), then `V(G) \ F` contains a vertex-minimal
       failed exterior, so `|V(G) \ F| >= 20` and statement 2 applies to it. For the cubic carrier, `F = K'` is
       five-colourable (C1), so every equality tester whose old vertices lie in `K'` has such an exterior in `L' \ K'`.
4. **(Odd-degree descent of three known exteriors.)**
   (a) The 19-vertex Exoo-Ismailescu graph EI19 has no faithful placement in any finite odd-degree extension of
       `L = Q(zeta55)`, in particular not in the cubic carrier `L' = L(theta)`, `2 theta^3 + 2 theta^2 - 2 theta - 1 = 0`.
   (b) Let `M ⊇ K` be conjugation-stable with `[M:K]` odd (for example `K' = K(theta)`, `L` or `L'`). If the
       Hochberg-O'Donnell fish has an injective unit-edge placement in `M`, then the genus-one curve `C3` of
       fish_marked_elliptic_cover_20261013 has a point over the real quartic field `F_4 = Q(tau)`,
       `tau^4 - 55 tau^2 + 605 = 0`. This is the necessary condition already known for `L` (that note's (5)); odd-degree
       carriers over `K` add nothing to it.
   (c) **EI21 carries the same curve.** In every faithful realization of the 21-vertex Exoo-Ismailescu graph
       (House of Graphs 51378, labels of the S24 script), with the hub `v_15` at the origin, `a, b, c, d = v_13, v_14,
       v_19, v_20` and `P = a - d`, `Q = b - c`: `|P|^2 + |Q|^2 = 1` and `Re(P bar Q) = 0`. Consequently, with `X = |P|^2`,
       any placement in a conjugation-stable `M ∋ sqrt(-11)` gives `z_0, z_1, z_2 ∈ M^+` with `11 z_0^2 = X(1-X)`,
       `11 z_1^2 = X(4-X)`, `11 z_2^2 = (1-X)(X+3)`, `0 < X < 1`: a point of the fish's `C3` over `M^+`. So (b) holds
       verbatim for EI21. In particular EI21 has a placement in `L` only if `C3(F_4)` is nonempty, and it has none in
       `Q(zeta11, sqrt5)` or in any conjugation-stable `M ∋ sqrt(-11)` with `M^+` of odd degree over `Q(sqrt5)` (the fish
       note's (4): `C3` has no point over any odd-degree extension of `Q(sqrt5)`).
5. **(Exhaustive search, FINITE-EXACT.)** No graph on at most 21 vertices satisfies the necessary conditions
   (N1) triangle-free with at most two common neighbours per pair, (N2) degrees and list sizes as in 2(a),
   (N3) no 4-cycle among low vertices and `4|E(low)| <= 5(#low - 1)`, (N4) a high vertex. Hence every failed six-colour
   exterior has **at least 22 vertices**, and 21 exterior points always preserve every old six-colouring.
   On 22 vertices the conditions, together with (N5) non-three-colourability, are satisfiable (36 of the 137
   attachment patterns rooted at a high vertex of maximum degree; none of the two with root degree seven). Every model found is **five-choosable**, by an
   exact Alon–Tarsi coefficient, so by 2(c) none of them is a failed exterior.

## Proof

**Statement 2(a).** Minimality gives that `H - b` is `L`-colourable for every `b`; if `deg_H(b) < |L(b)|`, colour
`H - b` and then `b`, a contradiction. A non-`L`-colourable component would contradict minimality, so `H` is
connected. `|L(b)| >= q - 1` is (I1).

**Statement 2(b).** We use the theorem of Borodin (1977) and Erdős, Rubin and Taylor (1979): *if `C` is connected,
`L'` a list assignment with `|L'(v)| >= deg_C(v)` for all `v`, and `C` is not `L'`-colourable, then `C` is a Gallai
tree: every block is a complete graph or an odd cycle.* Let `C` be a component of the low subgraph `Λ`. By
minimality `H - V(C)` has an `L`-colouring `psi` (empty if `H = C`). For `v ∈ C` let `L'(v) = L(v)` minus the colours
`psi` gives to neighbours of `v` outside `C`; those neighbours are high, so
`|L'(v)| >= |L(v)| - (deg_H(v) - deg_C(v)) = deg_C(v)`. An `L'`-colouring of `C` would combine with `psi` into an
`L`-colouring of `H`; so `C` is a Gallai tree. By (I2) its complete blocks are single edges, and its cycle blocks have
length at least five. If every vertex were low, `H = Λ` would be such a tree, whose leaf block (or single block) has a
vertex of degree at most two `< q - 1`; so `H` has a high vertex. A tree of such blocks on `c` vertices has at most
`(5/4)(c - 1)` edges (an edge block adds one vertex and one edge, an `m`-cycle adds `m - 1` vertices and `m <= (5/4)(m-1)`
edges for `m >= 5`), hence

    |E(Λ)| <= (5/4)(ℓ - 1),     ℓ = number of low vertices >= 1.                                  (1)

**Counting lemma.** For `v ∈ H` of degree `d`: its neighbours are pairwise non-adjacent (I2), every neighbour `u` has
`deg(u) - 1` further neighbours, all at distance exactly two from `v`, and each vertex at distance two is adjacent to
at most two neighbours of `v` (I3). Hence

    |H| >= 1 + d + |N_2(v)|,     |N_2(v)| >= (1/2) sum_(u ∈ N(v)) (deg(u) - 1) >= d(q - 2)/2.         (2)

**Statement 1, simple form, and 2(d) without the refinement.** A high vertex has degree at least `|L| + 1 >= q`, so
(2) gives `|H| >= 1 + q + ceil(q(q-2)/2)` (19 for `q = 6`). The degree claims of 2(d) are 2(a) with `|L(b)| = 6 - t(b)`.

**Statement 2(d): nineteen is impossible.** Let `q = 6` and suppose `|H| = 19`.

- (i) Every vertex has degree 5 or 6: degree `d >= 7` gives `|H| >= 1 + 7 + 14 = 22` by (2).
- (ii) Every vertex `v` of degree 6 is **tight**: (2) reads `19 >= 1 + 6 + |N_2(v)|` with `|N_2(v)| >= 12`, so every
  neighbour of `v` has degree exactly 5, every vertex of `N_2(v)` has exactly two neighbours in `N(v)`, and
  `V(H) = {v} ∪ N(v) ∪ N_2(v)`.
- (iii) A vertex of degree 5 has `|L| <= deg = 5`, hence `|L| = 5` and it is low. A high vertex has degree 6 and
  `|L| = 5`. By (ii) all neighbours of a high vertex have degree 5, hence are low; high vertices are pairwise
  non-adjacent and exactly `6h` edges join high to low vertices (`h` = number of high vertices, `ℓ = 19 - h`).
- (iv) Summing low degrees and using (1): `5ℓ <= 2|E(Λ)| + 6h <= (5/2)(ℓ - 1) + 6h`, i.e. `50 <= 8.5 h`, so `h >= 6`.
- (v) Fix a high vertex `z` and put `W = N_2(z)`, `|W| = 12`. By (ii) at `z` and (iii), every other high vertex lies
  in `W`. Let `W_6, W_5` be the vertices of `W` of degree 6 and 5; `|W_6| >= h - 1 >= 5`. Each `w ∈ W` is non-adjacent
  to `z` and has exactly two neighbours in `N(z)`, so it has `deg(w) - 2` neighbours in `W`; by (ii) at `w ∈ W_6` its
  four `W`-neighbours lie in `W_5`. Counting the `W`-edges at `W_5`:

      3|W_5| = 4|W_6| + 2 e(W_5),   i.e.   2 e(W_5) = 36 - 7|W_6| <= 1,

  so `e(W_5) = 0` and `7|W_6| = 36`, impossible.

Hence `|H| >= 20`, which is statement 2(d) and the refined statement 1.

**Statement 2(c).** If `H` has a proper colouring with classes `X_1, ..., X_m`, `m <= floor(q/2)`, give `X_i` its own
pair `P_i` of colours (pairwise disjoint). Each `b ∈ X_i` has at most one old neighbour (I1), so `P_i \ phi(N_A(b))`
is nonempty; choosing from it colours `H` properly (classes are independent, palettes disjoint), a contradiction.
This is the doubling lemma of THM-4614 with (I1) removing every blocking vertex. If `H` were `(q-1)`-choosable it
would be `L`-colourable since `|L| >= q - 1`.

**Statement 3.** (a) If `G` had a proper five-colouring, recolouring `y` with the sixth colour (its neighbours carry
at most the other five) gives a proper six-colouring separating `x` and `y`. (b) Let `d = |x - y|`. If `d >= 1/2`,
rotate a copy of `G` about `x` until the image `y'` of `y` satisfies `|y - y'| = 1`; in every six-colouring of the
union, `y` and `y'` both receive the colour of `x` but are adjacent. If `d < 1/2`, choose points
`p_0, ..., p_m` with `|p_i - p_(i+1)| = d`, `|p_0 - p_m| = 1` (`m >= 1/d`) and place a congruent copy of `G` with ports
`p_i, p_(i+1)`; all `p_i` share a colour while `p_0 p_m` is an edge. Extra unit pairs among the copies only add
constraints. Conversely, if `chi(R^2) = 7`, de Bruijn–Erdős gives a finite non-six-colourable unit-distance graph;
an edge-critical one minus an edge `xy` is an equality tester. (This is the Moser spindle principle, there with three
colours and distance `sqrt 3`.) (c) Restrict a five-colouring of the field plane to `V(G)` and apply (a)'s
recolouring. (d) Restrict an equality tester's colourings: a separating six-colouring of `A` does not extend to
`V(G) \ F`; a vertex-minimal non-extendable subset is a failed exterior.

**Statement 4.** (a) The two support-three relations of EI19 force `(-1 ± 2 sqrt(-2))/3` and `(1 ± sqrt(-35))/6` among
the direction ratios of every faithful realization (ei19_fish_placement_q_zeta55_20261012), so a placement field
contains `sqrt(-2)` and `sqrt(-35)`. A quadratic subfield `Q(sqrt a)` of an odd-degree extension `M` of `L` satisfies
`[L(sqrt a):L] | [M:L]`, so `sqrt a ∈ L`; the quadratic subfields of `L` are `Q(sqrt 5)`, `Q(sqrt(-11))`,
`Q(sqrt(-55))`. (b) `M` contains `sqrt(-11)`, so `M = M^+(sqrt(-11))` with `M^+ = M ∩ R` and every
conjugation-anti-invariant element of `M` is `sqrt(-11)` times an element of `M^+`. The derivation of equations (1)
of fish_marked_elliptic_cover_20261013 uses only this, so a placement gives `D ∈ M^+` and `z_0, z_1, z_2 ∈ M^+` with
`11 z_0^2 = D(1-D)`, `11 z_1^2 = D(4-D)`, `11 z_2^2 = (1-D)(D+3)`, `0 < D < 1`: a point `P ∈ C3(M^+)`. Here
`[M^+ : F_4] = [M:K]` is odd, so the closed point of `C3` over `F_4` under `P` has odd degree `e`. The fibre of the
degree-eight map `C3 -> P^1_D` over `D = 2` is an `F_4`-rational effective divisor of degree 8. Since
`gcd(8, e) = 1`, `C3` carries an `F_4`-rational divisor of degree one, and Riemann–Roch on the genus-one curve `C3`
turns it into an `F_4`-rational point. (c) In an injective placement every unit 4-cycle is a parallelogram. The
four 4-cycles `(10,13,15,14)`, `(13,15,19,17)`, `(14,15,20,16)`, `(15,19,18,20)` give `v_10 = a + b`, `v_17 = a + c`,
`v_16 = b + d`, `v_18 = c + d` (hub at 0), and the edges `10-18`, `16-17` give `|P + Q| = |P - Q| = 1`, i.e.
`|P|^2 + |Q|^2 = 1`, `Re(P bar Q) = 0`. As in the fish note, `P bar Q` is anti-invariant, so `P bar Q - bar P Q =
2 sqrt(-11) z_0` with `z_0 ∈ M^+` and `11 z_0^2 = |P|^2 |Q|^2 = X(1-X)`; the unit pair `a, d` at squared distance `X`
gives `(a bar d - bar a d)^2 = -X(4-X)`, and the unit pair `b, c` at squared distance `1 - X` gives
`(b bar c - bar b c)^2 = -(1-X)(X+3)`; dividing these last two anti-invariant differences by `sqrt(-11)` yields
`z_1, z_2 ∈ M^+` with the displayed equations (the factor 2 belongs only to the earlier definition of `z_0`). `X = 0` or `1` would identify
`v_13 = v_20` or `v_14 = v_19`. The argument of (b) then applies to this point of `C3`. ∎

**Statement 5.** six_colour_exterior_twenty_sat_20261016.py roots the search at a high vertex `z` of maximum degree
`d` among high vertices; all degrees are then at most `d` (low vertices have degree `|L| <= 6`), and by (2) at most
`floor((n-1)/3)`. The neighbours of `z` are independent; every second neighbour is adjacent to one or two of them,
never two second neighbours to the same pair (I3). The attachment patterns are enumerated up to isomorphism (pair
graphs from the networkx atlas, single attachments up to the pair graph's automorphisms), and the remaining edges
are SAT variables with (N1)-(N4) and the lazy cut for (N5). For `n = 19, 20, 21` every pattern is UNSAT (1, 7 and
38 patterns; the `n = 19` run reproduces the hand proof of 2(d) without needing (N5)). For `n = 22` the 36 models are
checked by alon_tarsi_five_choosability_20261016.py: for an orientation with all outdegrees at most four, the
coefficient of `prod x_v^(outdeg v)` in `prod_(uv in E) (x_u - x_v)` is nonzero, so by the Combinatorial
Nullstellensatz every list assignment with lists of size five extends.

## Remarks

- Statement 1 improves clock_duality_oct15's Theorem 2 (seventeen) by two; its special exclusions of five-cores of
  sizes 16 and 17 are not needed and five-cores of sizes 18 and 19 are not excluded (only as failed exteriors).
  The palette version `q + ceil(q(q-2)/2)` improves their `q(q-1)/2` for `q >= 4`; the parity step (v) also works for
  other even `q` but is stated only for six.
- Statement 3 makes the target precise: an equality tester in the cubic carrier (of any size) would prove
  `chi(R^2) = 7`, and it would in particular be a six-chromatic unit-distance graph; no six-chromatic unit-distance
  graph is known. Its exterior must already pass statements 2 and 5: at least twenty-two vertices, minimum degree five, a
  degree-six vertex, and **chromatic number at least four**, i.e. a four-chromatic triangle-free unit-distance graph
  inside `L' \ K'`. None is known in `L` (HYP-9247 conjectures there is none). Among the four smallest known
  triangle-free four-chromatic unit graphs, statement 4 shows that EI19 cannot supply one in any odd-degree extension
  of `L`, and that the fish and EI21 can only if `C3(F_4)` is nonempty. THM-4662 proves `C3(F_4)` empty by a two-descent,
  so neither can; only EI17 (rigid, field open) is untouched.
- Statements 4(b) and 4(c) are existence reductions. `C3(F_4)` is decided (empty) in THM-4662, which therefore
  excludes the fish and EI21 from `Q(zeta55)` and from every conjugation-stable odd-degree carrier over `K`.
