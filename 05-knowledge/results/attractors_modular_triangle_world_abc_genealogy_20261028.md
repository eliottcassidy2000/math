# What the repository's results keep landing on: a pre-registered census, bridges into the triangle-group world, and the abc genealogy

opus-2026-10-28-S68, 2026-10-28.

**Status:**
- **PROVED:** THM-4839 (elementary) + CITED (Coxeter-Moser regular maps; Euclid; Frey-Hellegouarch; Kubert; Baker;
  Darmon-Granville).
- **FINITE-EXACT:** the census statistics and every bridge check.
- **HYP-9437** (pre-registered census): P1 CONFIRMED, **P2 and P3 REFUTED**.
- **MISTAKE-691:** three unflagged errors found by the genealogy: the S116 abc statement, an uncited
  re-derivation, and THM-3833's conditional typing of Darmon-Granville.
- **Independently audited:** see section 7.

"Attractor" is used in one precise sense: a structure into which many results map by explicit, checkable morphisms.
No dynamics is claimed.

## 0. Directive and the concurrent lanes

Owner, 2026-10-28: "begin another session where you merge in prior ideas tangentially related to the many past session
on F=U+S, {petersen, K_5, K_3,3} the graceful tree conjecture, the square sum problem, berggrens tree and other concepts
at your discretion, aiming to understand what larger unifying theories lurk s attractors for all the novel results we've
been able to manage in this repo. finally deeply consider prior work on the abc conjecture and related ideas that came
up in sessions, and how those ideas came up in other sessions".

The first half of that message (the owner's numbers) is opus-2026-10-28-S67 (THM-4837/4838).

Concurrent unifying proposals are **methodological** and are credited, not repeated:
- session84 (`session84_guarded_unification_20261028/README.md`): objects with operations, guards, observations and
  scales; "which information may be forgotten now".
- codex-incidence (`incidence_attractors_20261028/THEORY.md` section U): a typed algebra of guarded attachments, plus the
  sigma square-class matroid (graphic on {4,2,5,3,7,13,11}, non-graphic Fano extension).
- session72's "operation-relative representation theory" and session70's synthesis.
- the mac-mini reflections on the (X, M, O, rho) quadruple.

All of them decline a dynamical attractor. This session asks the complementary, **object-level** question: which
mathematical objects do the results land on?

## 1. Results first

1. **The census refutes a single object-level hub (HYP-9437).**
   - **Setup:** in the 320 recent canon theorems (THM-4400..4839), the modular family (PSL_2(Z), Farey, Berggren,
     Pythagorean, Pell, continued fractions, modular forms, triangle groups) is pervasive, at 53.8% (P1 confirmed).
   - **Result:** against a degree-preserving null model it is **not** the main bridge (P2 and P3 refuted). The robust
     associations are:
     - **finite geometry <-> tournaments**: z = 7.3 recent and 13.4 over the whole canon (Paley, quadratic residues, Fano);
     - **Collatz <-> modular**: z = 3.3 recent and 5.1 overall (convergents of log2 3, Farey/Minkowski addresses);
     - **F_2 / parity / character language**: the all-canon family that bridges the most others (lift-degree 5).
   - **The observer framework**, proposed by the concurrent lanes as the unifier, is the **least** bridging family in
     both windows. It is significantly *anti*-associated with Collatz (z = -5.8) and with modular (z = -3.1). So the
     methodological unifier and the two object-level hubs do not yet meet in the canon.
2. **Interpretation (not a census result).**
   - Both hubs are SL_2 phenomena at two levels: PSL_2(F_q) and quadratic residues (finite), and PSL_2(Z) with continued
     fractions (infinite).
   - The census shows that the bridge between those two levels is weak in the canon: FINITE_GEOMETRY-MODULAR has
     z = 1.2.
   - THM-4839 supplies explicit morphisms across exactly that gap, landing on the owner's named topics.
3. **THM-4839 (bridges).**
   - (a) K4, K5, K7 are the only complete unit Cayley graphs of imaginary-quadratic residue fields. Two of the owner's
     Kuratowski/Petersen trio are the toroidal regular maps {4,4}_(2,1) (Gaussian, Z[i]/(2+i)) and {6,3}_(1,1)
     (Eisenstein honeycomb mod 2 + omega). Alongside them: K7 = {3,6}_(2,1), Heawood = {6,3}_(2,1), K4,4 = {4,4}_(2,2).
   - (b) Petersen is the hemi-dodecahedron, a quotient of the (2,3,5) triangle group by A_5 = PSL_2(F_5).
   - (c) Euclid's parametrisation, squared, is the Klein-four Belyi map onto the Legendre line, ramified only over
     0, 1, oo.
   - (d) Every Pythagorean Frey curve y^2 = x(x - a^2)(x + b^2) has the rational point (a(a+c), ac(a+c)) of order 4.

   So Berggren's tree, the square-sum chains of THM-4510, and the abc triples (a^2, b^2, c^2) live on one genus-zero
   modular family of elliptic curves with Z/2 x Z/4 torsion.
4. **abc genealogy** (section 5). Eight lineages from 2026-03-14 to today.
   - The core consumer lineage runs THM-3833 -> THM-3848 -> THM-4108 -> 2026-10-11 Codex -> aliquot124 -> THM-4805.
   - The Collatz lineage is Catalan -> Pillai -> Baker. Its verdict is "abc no stronger than Baker through the cycle
     gate", which this session confirms.
   - The repeated errors are textbook abc facts typed new (MISTAKE-675/676/678) and Baker name-checks on Collatz cycles
     (MISTAKE-555/604/605/606).
   - The unifying reading: abc is the arithmetic of the thrice-punctured line P^1 - {0,1,oo} = Gamma(2)\H.
     - Collatz clocks are its {2,3}-unit points; every clock from 3/2 on is an abc hit, and 27/17, the 231 -> 233 clock,
       has q = 1.086.
     - Berggren triples are its points pulled back along the Klein-four map: 60 hits to depth 9, best q = 1.266.
     - session69's square-later bijection and the aliquot "row of one" are its points on the line A + 2^a = C.

## 2. The census (HYP-9437)

- **Method.** Ten regex families: MODULAR, QUADRATIC, FINITE_GEOMETRY, BINARY_ALGEBRA, DIVISOR, COLLATZ, TOURNAMENT,
  ZETA, ABC, OBSERVER.
  - The lexicon was fixed and committed with the pre-registration (593105aedd) before the run. Each theorem file is
    read for 6000 characters.
  - Co-occurrence counts are compared with 2000 curveball randomisations. These keep every file's number of families
    and every family's size, which controls the "long files mention everything" effect.
  - A family's lift-degree is its number of partners with z >= 2.

| window | most frequent | top pairs by z | lift-degree leader | least bridging |
|---|---|---|---|---|
| THM-4400..4839 (320 files, pre-registered) | BINARY 57.5%, MODULAR 53.8%, COLLATZ 49.7% | FIN_GEOM-TOURN 7.29; COLLATZ-MODULAR 3.31; BINARY-TOURN 2.06 | TOURNAMENT (2) | OBSERVER (strength 1.19) |
| all canon (3957 files, exploratory) | BINARY 50.8%, MODULAR 34.2%, TOURN 29.4% | FIN_GEOM-TOURN 13.36; BINARY-MODULAR 6.31; ABC-DIVISOR 5.83; BINARY-DIVISOR 5.78; DIVISOR-MODULAR 5.27; COLLATZ-MODULAR 5.11 | BINARY (5), DIVISOR (4), MODULAR (3) | OBSERVER (strength 1.87) |

**Pre-registered verdicts:** P1 (MODULAR >= 20%) CONFIRMED at 53.8%; P2 (MODULAR has the largest lift-degree)
REFUTED; P3 (a top-3 pair inside the arithmetic-geometry core) REFUTED.

**What survives.** Two robust object-level hubs, a common binary language, and an isolated methodological layer. The
strongest negative pairs are COLLATZ-OBSERVER (-5.78), COLLATZ-TOURNAMENT (-3.25) and MODULAR-OBSERVER (-3.10). The
repository's lanes are partly segregated by vocabulary.

**Typing.** Keyword co-occurrence is not mathematical dependence. The census is evidence about where the canon's
language clusters, not about the mathematics. Its value here is that it refuted the session's own first thesis
(MODULAR as the hub) before that thesis could be typed.

## 3. THM-4839: bridges across the gap

All statements are checked in `attractor_bridges_regular_maps_belyi_frey_20261028.py`.

**(a) Complete unit Cayley graphs.** For a prime P of an imaginary quadratic O_K with N(P) >= 4, the Cayley graph on
O_K/P with connection set the image of O_K^x is complete iff (K, N(P)) is one of (Q(sqrt-3), 4), (Q(i), 5) or
(Q(sqrt-3), 7). These give K4, K5 and K7.
- Completeness needs N - 1 <= |O_K^x| in {2, 4, 6}, so N in {4, 5, 7}.
  - N = 4 needs 2 inert, which happens in Z[omega] and not in Z[i].
  - N = 5 needs 5 split with four units: Z[i].
  - N = 7 needs 7 split with six units: Z[omega].
  - With only +-1 the image misses.
- Checked over the class-number-one fields d in {1,2,3,7,11,19,43}, p < 60.
- K5 and K7 are the 1-skeletons of the toroidal regular maps {4,4}_(2,1) and {3,6}_(2,1). The six Eisenstein units
  collapse in pairs mod 2, so {3,6}_(2,0) is K4 with doubled edges; K4 itself is the planar tetrahedron.
- **Typing:** the regular maps are CITED (Coxeter-Moser, Generators and Relations, ch. 8), and K5 = {4,4}_(2,1) is
  classical. The "only three" statement is a small PROVED observation.

**(b) The rest of the trio, and its neighbours.**
- **Eisenstein:** the honeycomb modulo the Eisenstein prime 2 + omega (norm 3) is **K3,3** ({6,3}_(1,1)). Modulo
  3 + omega (norm 7) it is the **Heawood** graph ({6,3}_(2,1)), the Fano incidence graph, dual to K7 on the torus.
- **Gaussian:** Z[i]/(2 + 2i) gives K4,4.
- **Petersen:** the dodecahedron modulo the antipodal map is the Petersen graph (explicit isomorphism with Kneser
  K(5,2)), that is, the hemi-dodecahedron on RP^2. Its map group is A_5 = PSL_2(F_5), a quotient of the (2,3,5)
  triangle group.
- **Non-planarity:** all five non-planar graphs fail the girth form of Euler's bound, E <= g(V - 2)/(g - 2).

**So the owner's {Petersen, K5, K3,3} are the first non-planar regular-map quotients of the three triangle groups
(2,3,5), (2,4,4) and (2,3,6).** The last two are the Gaussian and Eisenstein grids of S66 (THM-4825). The (2,3,7) and
(2,3,11) members carry B_6 and B_10 (S67, THM-4838).
- **Typing:** DICTIONARY with CITED components.
- The repository already had Petersen as the conic graph mod 5 (Coxeter 1983; `pythagorean_form_conic_graphs_rp_margin_20261018.md`)
  and the mod-5 Berggren image A_5 (MISTAKE-640).

**(c) The Klein-four Belyi map.** lambda(t) = ((1 - t^2)/(1 + t^2))^2 is invariant under t -> -t, 1/t, -1/t.
- Its critical values are {0, 1, oo}, with fibres t = +-1 over 0, t = 0, oo over 1, and t = +-i over oo, all of index 2.
  So it is a Galois Belyi map of degree 4 with passport (2^2, 2^2, 2^2); Riemann-Hurwitz gives -2 = 4(-2) + 6.
- With t = n/m, Euclid's (a, b, c) = (m^2 - n^2, 2mn, m^2 + n^2) gives lambda = (a/c)^2 and 1 - lambda = (b/c)^2.
- So the abc triples (a^2, b^2, c^2) of the Berggren tree are exactly the points of the Legendre line in the image of
  this cover, and the modulus pair (k, k') = (a/c, b/c) lies on the circle k^2 + k'^2 = 1.

**(d) Pythagorean Frey curves have a rational point of order 4.** On E: y^2 = x(x - a^2)(x + b^2) with a^2 + b^2 = c^2,
the point Q = (a(a + c), ac(a + c)) satisfies 2Q = (a^2, 0), so Q has order 4.
- Proof: the 2-torsion point (a^2, 0) is halvable because a^2 - 0 and a^2 + b^2 = c^2 are both squares.
- Checked for 6000 (triple, leg-order) pairs from the Berggren tree.
- The tree, generated with hypotenuse pruning, reproduces every primitive triple with c <= 2000 exactly once (319).
- Such curves are classified by a genus-zero modular curve for Gamma(2) cap Gamma_1(4) (CITED: Kubert 1976).

**Measured.**
- The Berggren tree to depth 9 (29524 triples) has 60 abc hits rad(abc) < c^2. The best is q = 1.2659 at
  (36207, 18424, 40625); the next are (1377, 3136, 3425) and (2793, 10624, 10985); (7, 24, 25) has q = 1.204. The
  repository's ADU hit (175, 288, 337) (MISTAKE-678) is among them.
- The Collatz clocks K/L from 3/2 on are all abc hits: q = 1.226, 1.019, 1.273, 1.143, 1.236, 1.086 (27/17), 1.064,
  1.063, 1.081 (84/53).
  - The reason is |2^K - 3^L| < max/6. Baker's bound sends q -> 1 along convergents.

## 4. The owner's topics on the map

| topic | lands on | strongest repository anchors | fit |
|---|---|---|---|
| F=U+S | divisor exponent boxes; F_2 square classes (codex-incidence); the ten divisors of p^2qr as E(K5) = V(Petersen) = A_4 roots (set level only) | `arithmetic_braids_20260917_divisors.md`, THM-4733, session84 FUS frontier | DIVISOR + BINARY hub |
| Petersen, K5, K3,3 | regular maps of the triangle groups (2,3,5), (2,4,4), (2,3,6) (THM-4839); Kuratowski; Delta-Y families | THM-4529, THM-4701, THM-261, session70 graph-labelling report | finite-geometry side of SL_2 |
| graceful trees | signed-permutation / Weyl-orbit reformulations (THM-2787); reconfiguration, not existence | THM-2761..2798, the 2026-10-19 graceful packets | weakest fit |
| square-sum problem | the light cone x^2 + y^2 = z^2 (THM-4510's three-square chains from primitive triples); Hamiltonicity is Gerbicz's theorem | THM-4505, THM-4510, square_sum_ports_20261023 | modular side via Pythagoras |
| Berggren's tree | Gamma(2)/theta-group words, non-backtracking (THM-4808), Farey = light cone (THM-3333); now the Klein-four Belyi map and 4-torsion Frey curves | THM-4808, THM-3333, THM-2596 | modular side |
| abc | the thrice-punctured line: Collatz {2,3}-units, Berggren pullbacks, square-later and row-of-one points | section 5 | modular side, conditional |

The graceful-tree conjecture remains the outlier. Plane trees are the dessins of Shabat polynomials, which is a natural
home in the Belyi world. No link between graceful labellings and Shabat polynomials is known or claimed here; it is
recorded as a pointer only.

**What the map suggests (research directions, not claims):**
1. The weak FINITE_GEOMETRY-MODULAR association marks the least exploited bridge. Results that pass from PSL_2(Z)
   words (Berggren, Collatz clocks) to PSL_2(F_q) structure (Paley, Fano, Petersen) are where the two hubs could merge.
   THM-4808's mod-q quotients and THM-4839 are first examples.
2. The observer framework has not yet been applied to the Collatz/modular hub. A guarded-observer quotient of the
   Collatz clock words would test whether the methodological unifier reaches the object-level hubs.

## 5. The abc genealogy (2026-03-14 to 2026-10-28)

Built from a read-only survey and checked against the files at the points marked "checked".

**Chronology (condensed).**
- **2026-03-14/15 (opus S84; S116).** Numerology and a survey. `famous_problems_s116.out` lines 67-68 state "ABC
  CONJECTURE (proved by Mochizuki, controversial)" and "rad(abc)^{1+eps} > c for all eps > 0".
  - Both are wrong: abc is open and IUT's claim is disputed. As written the inequality fails at 1 + 8 = 9 for
    eps < 0.2263.
  - Checked; this is the only such statement in the canon and knowledge base (MISTAKE-691 (1)).
- **06-05 to 06-27 (codex S651; LRC14 packets).** abc as a "carrier"; Fermat-Catalan/Beal shapes as packet guards.
  Typed synthesis, not proof.
- **08-23 (root abc-IUT bridge).** `CORE-PAPERS-ABC-IUT.md`: ABC OPEN, IUT disputed.
  - THM-3833 (checked): abc gives the two-cube radical floor and "non-effective finiteness for every fixed primitive
    hyperbolic Fermat-Catalan signature".
    - For a fixed signature that finiteness is an **unconditional** theorem: Darmon-Granville 1995, Thm 2, via
      Faltings and Chevalley-Weil.
    - The genuinely abc-conditional content is the uniform version over all hyperbolic signatures, which uses the
      minimal deficit 1 - 1/2 - 1/3 - 1/7 = 1/42 (= B_6, THM-4838) and yields the Fermat-Catalan conjecture
      (standard).
    - MISTAKE-691 (3).
  - THM-3848 (radical growth on carry packets) and the Mason-Stothers uses in the JC2 lane (THM-3881/3885/3920), a
    separate polynomial lineage, date from the same day.
- **08-25.** THM-4108 (checked): abc gives saturation of rad(a^b +- b^a); the "ABC lane" in the backlog.
- **09-21 to 09-24 (Collatz).**
  - The Catalan lane, then Pillai and convergents, then the implication atlas.
  - The atlas (checked, line 35): "Through the cycle gate, abc is provably no stronger than Baker". Proposition R2:
    abc gives a parity-prefix bound, extending Rozier.
- **09-29.** THM-4516: perfect-power clocks as Fermat-Catalan identities; MISTAKE-551.
- **10-11 (Codex).** abc_difference_families: the twin-centre identity (m^2 - 1) + 1 = m^2 (eq. (2)) and conditional
  radical limits.
- **10-14.** THM-4632/4633 (convergents of log2 3), together with MISTAKE-604/605/606 (Baker, Eliahou and Hercher
  attributions) and MISTAKE-610.
- **10-26 (aliquot directive).**
  - session69 (checked): square-later integers <-> primitive positive abc triples.
  - aliquot124 adds a conditional floor.
  - **THM-4805 (opus S60):** (1, 2^k - 1, 2^k) is a hit iff 2^k - 1 is not squarefree. The textbook 2^(6n) family was
    typed new (MISTAKE-675). The radical-6 triples were over-claimed (MISTAKE-676). The Berggren "abc non-connection"
    was false (MISTAKE-678: ADU is a hit with chi = -1).
  - twin_pythagorean_124 re-derives 1 + p(p+2) = m^2 (eq. (8)) without citing the 10-11 packet (checked;
    MISTAKE-691 (2)).
- **10-28.** session82's decimal family (1, 10^{2k}, 10^{2k} + 1), typed "no novelty claim" (checked), and session83's
  Grantham citation. The lesson of MISTAKE-675 held.

**Lineages and how ideas moved.**
1. **Survey and numerology** (S84, S116, S651). Never cited later; S116's error was never annotated.
2. **The conditional radical consumer** (THM-3833 -> 3848 -> 4108 -> 10-11 Codex -> aliquot124 -> THM-4805). Each step
   cites the previous. Every member is conditional, non-effective, and feeds no unconditional theorem.
3. **Polynomial abc** (Mason-Stothers in JC2 -> HYP-9216). Isolated from lineage 2 although written the same day.
4. **Collatz cycles** (Catalan -> Pillai/convergents -> Baker -> Eliahou/Hercher). abc is dominated through the cycle
   gate. Repeated attribution slips (MISTAKE-555, 604, 605, 606).
5. **Fermat-Catalan/Beal** (LRC14 guards -> THM-3833 -> THM-4516, MISTAKE-551).
6. **Twin centres** (10-11 -> twin_pythagorean -> THM-4808, growing_sieve_path). The re-derivation is uncited.
7. **Aliquot abc maps** (owner directive -> session69 -> aliquot124 -> S60 -> session73 -> session82/83).
8. **Logical form:** abc is Pi^0_3 (klein 10-25, session65).

**The unifying reading (this session).** Every lineage studies points of P^1 - {0,1,oo} = Gamma(2)\H with prescribed
reduction:
- lineage 4 is the S-unit equation over Z[1/6], where the clocks are hits and Baker beats abc;
- lineages 2 and 6 are radical floors along explicit families on that line;
- lineage 7 lies on the lines A + 2^a = C;
- Berggren is the Klein-four pullback, where the Frey curves carry Z/2 x Z/4.

That is why abc keeps reappearing next to Berggren, Collatz and the aliquot rows: they share the base curve.
- **Typing:** INTERPRETATION built on THM-4839 (c)-(d) and cited facts. No new abc result is claimed.

**Open questions the repository leaves.**
- the ledger's task queue (radical-fibre atlas, Pell-ray quality);
- the backlog "ABC lane";
- HYP-9415 (3) and (5) (Mersenne squarefreeness, Wieferich);
- HYP-9127/9130 (the uniform S-unit gap U(1/3));
- session82's uniform margin;
- growing_sieve_path's question of which coordinates decide abc-hit status;
- the missing coordinate rad(u - 2^a) in session69.

## 6. Reproduction

```text
python3 04-computation/experiments/attractor_census_canon_20261028.py 4400 4839 2000   # pre-registered census (about 1 min)
python3 04-computation/experiments/attractor_census_canon_20261028.py 1 4839 300       # exploratory full canon
python3 04-computation/experiments/attractor_bridges_regular_maps_belyi_frey_20261028.py   # THM-4839 checks (2 s)
```

Outputs: the `.out` files beside the scripts. The full-canon run is summarised in section 2, and its output is the
`_fullcanon.out` file.

## 7. Audit record

(pending)
