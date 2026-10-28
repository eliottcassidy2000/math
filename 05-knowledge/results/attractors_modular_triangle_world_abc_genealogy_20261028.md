# What the repository's results keep landing on: a pre-registered census, bridges into the triangle-group world, and the abc genealogy

opus-2026-10-28-S68, 2026-10-28.

**Status.**
- PROVED: THM-4839 (elementary).
- CITED: Coxeter-Moser; Darmon-Granville 1995; Kubert; Halbeisen-Hungerbuehler; Elkies; Baker.
- FINITE-EXACT: every bridge check.
- **HYP-9437** (pre-registered census): P1 CONFIRMED, **P2 and P3 REFUTED**. The registered lexicon has a bug, disclosed in
  section 2; the verdicts stand, but the readings were corrected after audit.
- MISTAKE-691 (genealogy findings, amended after audit) and MISTAKE-693 (this session's own errors).
- Independently audited, **SOUND WITH CORRECTIONS**, all applied (section 7).

"Attractor" is used in one precise sense: a structure into which many results map by explicit, checkable morphisms.
No dynamics is claimed.

## 0. Directive and the concurrent lanes

Owner, 2026-10-28: "begin another session where you merge in prior ideas tangentially related to the many past session
on F=U+S, {petersen, K_5, K_3,3} the graceful tree conjecture, the square sum problem, berggrens tree and other concepts
at your discretion, aiming to understand what larger unifying theories lurk s attractors for all the novel results we've
been able to manage in this repo. finally deeply consider prior work on the abc conjecture and related ideas that came
up in sessions, and how those ideas came up in other sessions".

The first half of that message is opus-2026-10-28-S67 (THM-4837/4838).

Concurrent unifying proposals are **methodological** and are credited, not repeated:
- session84's guarded observers (`session84_guarded_unification_20261028/README.md`);
- codex-incidence's typed algebra of guarded attachments and its sigma square-class matroids
  (`incidence_attractors_20261028/THEORY.md` section U);
- session72's operation-relative representation theory, and session70;
- the mac-mini (X, M, O, rho) reflections.

All decline a dynamical attractor. This session asks the complementary, object-level question.

## 1. Results first

1. **The census refutes the modular family as the bridging hub** (HYP-9437: P2 and P3 REFUTED as registered; P1 passes).
   - **The bug.** The registered MODULAR share, 53.8%, is inflated by a lexicon bug: the pattern `Hecke`
     (case-insensitive, no word boundary) matches "checked". With a word boundary the share is 23.1%.
   - **What is robust.** The one association that survives seeds, reading length, front matter and lexicon corrections
     is **finite geometry <-> tournaments**: z = 7.3 recent and 13.4 over the whole canon, the strongest pair in all 26
     audit variants (Paley, quadratic residues, Fano).
   - **What is an artifact.** The registered Collatz <-> modular pair (z = 3.3 recent, 5.1 overall) comes from the same
     bug: 1.4 and 0.8 after correction. The binary family's all-canon lead also does not survive the corrected lexicon.
   - **The OBSERVER proxy** is mostly the word "carrier". It is anti-associated with Collatz in every variant; its other
     properties depend on that word.
2. **The bridge that is missing.** Over the whole canon, finite geometry and modular vocabulary are *anti*-associated
   (z = -4.3). The repository's PSL_2(F_q) world (Paley, Fano, Petersen) and its PSL_2(Z) world (Farey, Berggren,
   continued fractions) rarely meet in one theorem. THM-4839 supplies explicit morphisms between them that land on the
   owner's named topics.
3. **THM-4839.**
   - (a) For N(P) >= 4, K4, K5 and K7 are the only complete unit Cayley graphs of imaginary-quadratic residue fields.
   - (b) Petersen, K5 and K3,3 are the smallest non-planar simple skeletons of regular maps of types {5,3}, {4,4} and
     {6,3}: the hemi-dodecahedron (map group A_5 = PSL_2(F_5)), the Gaussian {4,4}_(2,1) and the Eisenstein
     {6,3}_(1,1). These are the grids of THM-4825, and the {5,3} type belongs to the spherical (2,3,5) family of
     THM-4838.
   - (c) Euclid's parametrisation, squared, is the Klein-four cover X(4) -> X(2) (classical; Darmon-Granville).
   - (d) The Pythagorean Frey curves are the X(4)-subfamily of the curves with Z/2 x Z/4 torsion, with rational 4-torsion
     point (a(a+c), ac(a+c)).

   So Berggren's tree, the square-sum chains of THM-4510 and the abc triples (a^2, b^2, c^2) live on the genus-zero
   modular curve X(4).
4. **abc genealogy** (section 5).
   - Eight lineages from 2026-03-14; each member of the conditional-radical lineage cites THM-3833, though the chain is
     not linear.
   - The Collatz lineage is Catalan -> Pillai -> Baker. The atlas's verdict is "abc provably no stronger than Baker
     through the cycle gate".
   - The repeated errors are textbook abc facts typed new (MISTAKE-675), over-reached claims (MISTAKE-676/678) and
     Baker/Eliahou attributions on Collatz cycles (MISTAKE-555/605/606).
   - The unifying reading, Elkies' and Vojta's viewpoint, is that abc is the arithmetic of P^1 - {0,1,oo} = X(2):
     - Collatz clocks are Pillai-type points 2^K - 3^L = d on it; every clock from 3/2 on is an abc hit, and 27/17 has
       q = 1.086;
     - Berggren triples are its X(4)-points, with 60 hits to depth 9 and best q = 1.266;
     - session69's square-later bijection lies on the lines A + 2^a = C.

## 2. The census (HYP-9437)

**Method (registered).**
- Ten regex families over the first 6000 characters of each canon file.
- 2000 curveball randomisations, which preserve both degree sequences.
- A family's lift-degree is its number of partners with z >= 2.
- The lexicon was fixed in 593105aedd before the run.

| window | registered top pairs (z) | registered lift-degree leader | after the audit's lexicon corrections |
|---|---|---|---|
| THM-4400..4839 (320 files) | FIN_GEOM-TOURN 7.29; COLLATZ-MODULAR 3.31; BINARY-TOURN 2.06 | TOURNAMENT (2) | FIN_GEOM-TOURN stays strongest; COLLATZ-MODULAR 1.4; MODULAR share 23.1% (20.6% fully corrected) |
| all canon (3957 files, exploratory) | FIN_GEOM-TOURN 13.36; BINARY-MODULAR 6.31; ABC-DIVISOR 5.83 | BINARY (5) | FIN_GEOM-TOURN 13.8-14.9; COLLATZ-MODULAR 0.75; BINARY lift-degree 2; FIN_GEOM-MODULAR -4.3 |

- **Pre-registered verdicts:** P1 CONFIRMED (53.8% registered; still 23.1% / 20.6% after correction), P2 REFUTED,
  P3 REFUTED. These verdicts keep their registered values under the corrections.
- **Robustness (audit, exploratory).**
  - FINITE_GEOMETRY-TOURNAMENT is the strongest pair in all 26 variants. COLLATZ-OBSERVER is negative in all variants.
  - COLLATZ-MODULAR, the binary lead and MODULAR-OBSERVER do not survive.
  - P2 would flip to CONFIRMED at BODY_CHARS = 12000.
  - Several associations live in the front matter: stripping it sends COLLATZ-MODULAR to 0.11.
  - A random spot check of 10 files found 15 of 31 family assignments false. MODULAR was wrong in 6 of 7; FINITE_GEOMETRY
    and TOURNAMENT were clean.
- **Reproducibility.**
  - The curveball shuffles sets of strings, so results depend on PYTHONHASHSEED: hash seeds 0, 1, 2 give
    FG-T z = 7.38, 7.13, 7.17, with unchanged verdicts.
  - Reproduce the registered corpus at 593105aedd, or pass ID_MAX = 4838 on later checkouts.
- **What survives:** one robust association (finite geometry with tournaments, the repository's founding Paley
  programme) and the separation of Collatz notes from "carrier"/"observer" vocabulary.
- **Typing.** Keyword co-occurrence is lexical evidence only. The census refuted this session's first thesis (modular as
  hub) before it could be typed. The audit then showed that a further apparent hub was a lexicon artifact.

## 3. THM-4839: bridges across the gap

Checked in `attractor_bridges_regular_maps_belyi_frey_20261028.py` and independently in the audit folder.

- **(a)** For N(P) >= 4, the unit Cayley graph of O_K/P is complete iff (K, N) is (Q(sqrt-3), 4), (Q(i), 5) or
  (Q(sqrt-3), 7).
  - The units must surject, so N - 1 <= 6. The session checked d in {1,2,3,7,11,19,43}, p < 60; the audit checked all
    fundamental discriminants |D| <= 2000, p <= 200.
  - For N(P) in {2, 3}, K2 and K3 always occur.
- **(b)** K5 = {4,4}_(2,1), K4,4 = {4,4}_(2,2), K7 = {3,6}_(2,1), K3,3 = {6,3}_(1,1) and Heawood = {6,3}_(2,1);
  {3,6}_(2,0) is K4 with doubled edges.
  - The dodecahedron modulo the antipodal map is Petersen, the hemi-dodecahedron {5,3}_5, whose map group A_5 acts
    regularly on 60 flags.
  - Petersen, K5 and K3,3 are the smallest non-planar simple skeletons of regular maps of types {5,3}, {4,4}, {6,3}. The
    dual types give K6 and K7.
  - {4,4}_(2,1), {3,6}_(2,1) and {6,3}_(2,1) are chiral.
  - Repository prior work: Petersen = hemi-dodecahedron in four notes of 2026-10-01..18; K7 = {3,6}_(2,1) in
    camion_busch section 6.4; the N = 7 case in torus_voltage G2.
- **(c)** lambda(t) = ((1 - t^2)/(1 + t^2))^2 is the classical Klein-four cover X(4) -> X(2) (Darmon-Granville 1995,
  section 7.1, r = 2), with passport (2^2, 2^2, 2^2). Euclid's t = n/m gives lambda = (a/c)^2 and 1 - lambda = (b/c)^2.
- **(d)** Q = (a(a+c), ac(a+c)) has order 4 on y^2 = x(x - a^2)(x + b^2).
  - The Pythagorean Frey curves are the fibres over X(4): a thin subfamily of Kubert's Z/2 x Z/4 curves.
    y^2 = x(x+1)(x+4) is a Z/2 x Z/4 curve that is not Pythagorean. Halbeisen-Hungerbuehler 2022 is related.
  - The Berggren tree with hypotenuse pruning reproduces all 319 primitive triples with c <= 2000.

**Measured.**
- 60 abc hits among the depth-<=9 Berggren triples; best q = 1.2659 at (36207, 18424, 40625).
- Every clock K/L from 3/2 on is a hit: q = 1.226, 1.019, 1.273, 1.143, 1.236, 1.086 (27/17), 1.064, 1.063, 1.081
  (84/53).
  - Proof: |L log2 3 - K| <= 0.2451 < log2(6/5) for L >= 2; the audit checked all 40 clocks with L <= 20000.
  - Baker makes log(max)/log(6|2^K - 3^L|) -> 1; q -> 1 would need abc.

## 4. The owner's topics on the map

| topic | lands on | strongest repository anchors | fit |
|---|---|---|---|
| F=U+S | divisor exponent boxes; F_2 square classes (codex-incidence); the ten proper nontrivial divisors of p^2qr as E(K5) = V(Petersen) = the 10 positive roots of A_4 (set level only) | `arithmetic_braids_20260917_divisors.md`, THM-4733, session84 FUS frontier | divisor / binary |
| Petersen, K5, K3,3 | smallest non-planar regular-map skeletons of types {5,3}, {4,4}, {6,3} (THM-4839); Kuratowski; Delta-Y families | THM-4529, THM-4701, THM-261, session70 report, noble_polyhedra (hemi-dodecahedron) | the finite-geometry side |
| graceful trees | signed-permutation / Weyl-orbit reformulations (THM-2787); reconfiguration, not existence | THM-2761..2798, the 2026-10-19 graceful packets | weakest fit |
| square-sum problem | the light cone x^2 + y^2 = z^2 (THM-4510's three-square chains); Hamiltonicity is Gerbicz's theorem | THM-4505, THM-4510, square_sum_ports_20261023 | modular side via Pythagoras |
| Berggren's tree | Gamma(2)/theta-group words and non-backtracking (THM-4808); Farey = light cone (THM-3333); now X(4) and its Frey curves | THM-4808, THM-3333, THM-2596 | modular side |
| abc | the thrice-punctured line X(2): Pillai-type clock points, X(4)-points, square-later lines | section 5 | modular side, conditional |

The graceful-tree conjecture remains the outlier. Plane trees are the dessins of Shabat polynomials, a natural Belyi
home, but no link to graceful labellings was found by the session or the audit. It is recorded as a pointer only.

**Directions (not claims):**
1. The anti-association of finite-geometry and modular vocabulary marks the least exploited bridge: results passing
   from PSL_2(Z) words (Berggren, Collatz clocks) to PSL_2(F_q) structure (Paley, Fano, Petersen). THM-4808's mod-q
   quotients and THM-4839 are first examples.
2. The observer framework could be tested on Collatz clock words, the family its vocabulary avoids most.

## 5. The abc genealogy (2026-03-14 to 2026-10-28)

Built from a read-only survey, checked against the files, and corrected after audit.

**Chronology (condensed).**
- **2026-03-14/15 (opus S84; S116).** `famous_problems_s116.out` lines 67-68, from the generator
  `04-computation/famous_problems_s116.py` line 76, state "ABC CONJECTURE (proved by Mochizuki, controversial)" and
  "rad(abc)^{1+eps} > c for all eps > 0".
  - abc is open, and as written the inequality fails at 1 + 8 = 9 for eps < 0.2263.
  - This is the only such statement in the canon and knowledge directories (MISTAKE-691 (1)).
- **06-05 to 06-27 (codex S651; LRC14 packets).** abc as a "carrier"; Fermat-Catalan/Beal shapes as packet guards.
- **08-23 (root abc-IUT bridge).** `CORE-PAPERS-ABC-IUT.md`: ABC OPEN, IUT disputed.
  - **THM-3833**: abc gives the two-cube radical floor and, in Theorem 4.1, finiteness with a height bound (18) for each
    fixed hyperbolic Fermat-Catalan signature. Its (24) shows that 1/42 is the least hyperbolic curvature.
  - Its implication is true, but it does not record that fixed-signature finiteness is unconditional (Darmon-Granville
    1995, Thm 2). The uniform version from (18) and (24) is DG 1995 section 5.2 (MISTAKE-691 (3)).
  - THM-4516 (09-29) and the 10-11 packet do cite Darmon-Granville.
  - Same day: THM-3848, and the separate Mason-Stothers lineage in JC2 (THM-3881/3885/3920).
- **08-25.** THM-4108 (rad(a^b +- b^a) saturation under abc); the "ABC lane" in the backlog.
- **09-21 to 09-24 (Collatz).** The Catalan lane, then Pillai and convergents, then the implication atlas. Atlas line 35:
  "Through the cycle gate, abc is provably no stronger than Baker". Proposition R2 is a parity-prefix bound extending
  Rozier.
- **09-29.** THM-4516 (perfect-power clocks as Fermat-Catalan identities; cites Darmon-Granville); MISTAKE-551.
- **10-11 (Codex).** abc_difference_families: the twin-centre identity (m^2 - 1) + 1 = m^2 (eq. (2)), conditional
  radical limits, and Darmon-Granville section 5.2 (line 34).
- **10-14.** THM-4632/4633 (convergents of log2 3), together with:
  - MISTAKE-605/606 (Eliahou credited, and the cycle bound attributed to Baker instead of Eliahou/Barina);
  - MISTAKE-604 (Gersonides and the trivial cycle);
  - MISTAKE-610.
- **10-26 (aliquot directive).**
  - session69: square-later integers <-> primitive positive abc triples **with B = 2^a and C odd**.
  - aliquot124 adds a conditional floor, citing THM-3833.
  - **THM-4805 (opus S60):** (1, 2^k - 1, 2^k) is a hit iff 2^k - 1 is not squarefree. The textbook 2^(6n) family was
    typed new (MISTAKE-675).
  - Separately, in mac-mini lanes:
    - MISTAKE-676: novelty over-claimed for repository-internal radical-6 facts.
    - MISTAKE-678: a Berggren "abc non-connection" asserted beyond its evidence; ADU is a hit with chi = -1.
  - twin_pythagorean_124 duplicates the 10-11 identity as its eq. (8) without citing it. No novelty is claimed there
    (MISTAKE-691 (2)).
- **10-28.** session82's decimal family, typed "no novelty claim", and session83's Grantham citation.

**Lineages.**
1. Survey and numerology (S84, S116, S651). Never cited later; S116's error was never annotated.
2. The conditional radical consumer (THM-3833, 3848, 4108, the 10-11 packet, aliquot124, THM-4805). Each member cites
   THM-3833, but the chain is not linear: aliquot124 cites THM-3833 and not the 10-11 packet. Every member is
   conditional and non-effective.
3. Polynomial abc (Mason-Stothers in JC2 -> HYP-9216). Isolated from lineage 2.
4. Collatz cycles (Catalan -> Pillai/convergents -> Baker -> Eliahou/Hercher), with repeated attribution slips
   (MISTAKE-555, 605, 606).
5. Fermat-Catalan/Beal (LRC14 guards -> THM-3833 -> THM-4516, MISTAKE-551).
6. Twin centres (10-11 -> twin_pythagorean -> THM-4808, growing_sieve_path), with an uncited duplication.
7. Aliquot abc maps (owner directive -> session69 -> aliquot124 -> S60 -> session73 -> session82/83).
8. Logical form: abc is Pi^0_3 (klein 10-25, session65).

**The unifying reading** is INTERPRETATION, from the classical Elkies (1991) / Vojta viewpoint. Each lineage studies
rational points of P^1 - {0,1,oo} = X(2) with prescribed reduction:
- lineage 4 studies Pillai-type points 2^K - 3^L = d, where only 3/2 is a {2,3}-unit solution;
- lineages 2 and 6 study radical floors on explicit families;
- lineage 7 lies on the lines A + 2^a = C;
- Berggren is the X(4)-pullback, whose Frey curves carry Z/2 x Z/4.

That shared base curve is why abc keeps reappearing next to Berggren, Collatz and the aliquot rows. No new abc result is
claimed.

**Open questions left by the repository:**
- the ledger's task queue (radical-fibre atlas, Pell-ray quality);
- the backlog "ABC lane";
- HYP-9415 (3) and (5);
- HYP-9127/9130 (the uniform S-unit gap U(1/3));
- session82's uniform margin;
- growing_sieve_path's question of which coordinates decide hit status;
- session69's missing coordinate rad(u - 2^a).

## 6. Reproduction

```text
python3 04-computation/experiments/attractor_census_canon_20261028.py 4400 4838 2000   # registered corpus (320 files); PYTHONHASHSEED=0 for stable output
python3 04-computation/experiments/attractor_census_canon_20261028.py 1 4838 300       # exploratory full canon
python3 04-computation/experiments/attractor_bridges_regular_maps_belyi_frey_20261028.py
```

The audit's corrected-lexicon, robustness and independent THM-4839 scripts are in
`04-computation/experiments/attractors_census_bridges_20261028_audit_S68/`.

## 7. Audit record

Independent audit `04-computation/experiments/attractors_census_bridges_20261028_audit_S68/AUDIT.md`:
**SOUND WITH CORRECTIONS**, all applied here, in THM-4839, HYP-9437, MISTAKE-691 and MISTAKE-693.

- **Pre-registration integrity confirmed.**
- **MAJOR.**
  - (1) The census lexicon bug (`Hecke` matches "checked"). Registered verdicts stand; readings corrected; one robust
    association remains.
  - (2) THM-4839 (d) named the wrong modular family. The Pythagorean curves are the X(4)-subfamily of Kubert's
    Z/2 x Z/4 curves.
- **MINOR.**
  - Baker does not give q -> 1.
  - "First non-planar" holds only by map type, and three of the maps are chiral.
  - The N(P) >= 4 hypothesis was missing from the title.
  - The field list was incomplete.
  - The A_5 phrasing was garbled.
  - Internal prior art (hemi-dodecahedron, the chiral K7, torus_voltage G2) and Darmon-Granville section 7.1 for (c)
    were uncredited.
  - Census fragility and hash-seed dependence; "carrier" dominates OBSERVER.
  - MISTAKE-691 (3) was unfair to THM-3833, and (2) is a credit omission.
  - Genealogy paraphrases: session69's scope, clocks not being unit-equation points, "Baker beats abc", and the
    MISTAKE-604/676/678 descriptions.
- **NOTES.**
  - A_4 has 10 positive roots.
  - Gaussian/EJ networks literature.
  - A proof that every clock is a hit.
  - Elkies/Vojta for the thrice-punctured-line view.
  - Delta(2,4,4) is not a quotient of PSL_2(Z).
  - The S116 generator script.
  - The THM-4839 filename exceeded Windows MAX_PATH and has been shortened.
