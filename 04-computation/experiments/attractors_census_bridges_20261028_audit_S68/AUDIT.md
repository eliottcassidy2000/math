# Independent audit: S68 attractor census, THM-4839 bridges, abc genealogy, MISTAKE-691

Referee: independent adversarial audit (Claude, 2026-10-28), worktree `math-wt-s68-attractors` at 4e32454073.
Nothing outside this folder was edited; nothing was committed.

**VERDICT: SOUND WITH CORRECTIONS.** The proved mathematics of THM-4839 (a)-(d) is correct, and the pre-registration of
HYP-9437 is honest; but the census lexicon contains a substring bug ("Hecke" matches "checked") that manufactures the
note's headline numbers (MODULAR "53.8%", the "Collatz <-> modular" hub), THM-4839 (d) names the wrong modular curve
for the Pythagorean Frey curves, and several sentences of the note and of MISTAKE-691 need precision fixes.

Files under audit:
- N = `05-knowledge/results/attractors_modular_triangle_world_abc_genealogy_20261028.md` (the note)
- T = `01-canon/theorems/THM-4839-bridges-into-the-triangle-group-world-...-4-torsion-frey-curves.md`
- H = `05-knowledge/hypotheses/HYP-9437-the-modular-triangle-group-family-...-pre-registered-census.md`
- M = `01-canon/MISTAKES.md`, MISTAKE-691 (lines 21087-21118)
- the two scripts and three `.out` files in `04-computation/experiments/attractor_*_20261028*`

What was redone independently (own code, in this folder):
- `audit_thm4839_independent.py` (+ `.out`): every THM-4839 claim rebuilt from scratch (fields of every class number,
  genuine maps with faces and chirality, theta-function identification of the Belyi map, torsion, abc counts, clocks).
- `audit_census_independent.py`: a re-implementation of the census (only the registered regex dictionary is read, via
  `ast`, from the committed script), a robustness grid, per-term match counts and a random spot check.
- `rerun_session_census.py`: the session's census script run unchanged under three hash seeds.
- Darmon-Granville 1995 read in full at the cited theorem pages (McGill scan, rendered locally); web prior-art search.

---

## Findings

### MAJOR

**MAJOR-1. The census lexicon matches "checked": the 53.8% MODULAR share, the Collatz-modular "hub" and the binary
"bridging leader" are artifacts. (N lines 41-50, 92-93, 98-100; H lines 64, 79-85; results INDEX line 8.)**
- Evidence. The registered MODULAR pattern contains the alternative `Hecke` with no word boundary, compiled with
  `re.I` (`attractor_census_canon_20261028.py` line 32, line 52). It matches "checked", "checker", "cross-checked",
  "checkerboard". In the 320-file window, 114 files contain "checked" and only 2 contain the word Hecke
  (`audit_census_terms.out`: `Hecke` hits 127 files, 98 of them matched by no other MODULAR term).
  - With word boundaries on `Hecke` alone, MODULAR falls from 172 files (53.8%) to 74 (23.1%). 58 of the 98 spurious
    files are COLLATZ files (Collatz notes record verifications as "checked"), so the files in both COLLATZ and
    MODULAR fall from 102 to 44.
  - COLLATZ-MODULAR z: registered 3.37 (own implementation; session 3.31) -> **1.44** with only `Hecke` fixed ->
    **1.41** with the full audit lexicon. Full canon: 4.92 -> **0.75**.
  - BINARY_ALGEBRA also matches "characteristic/characterisation" (`character`, 45 files matched only by it) and file
    names (`parity`); with the corrected lexicon its full-canon lift-degree is 2, not 5, and the full-canon leaders are
    DIVISOR, MODULAR and ZETA (3 each).
  - Random spot check of 10 window files (seed 68, `audit_census_spotcheck.out`): 31 family assignments, 15 false
    positives by content (48%). MODULAR was assigned to 7 of the 10 files; 5 of those 7 came from "checked"-type words,
    and 6 of the 7 are false (the 7th, THM-4686 on CM supercongruences, is genuinely modular but was caught only through
    "checked"). FINITE_GEOMETRY (2/2) and TOURNAMENT (5/5) were clean.
- What survives (robustness grid below): FINITE_GEOMETRY-TOURNAMENT is the strongest pair in **all 26 variants**
  (z = 5.7-7.4 recent, 13.8-14.9 full canon). COLLATZ-OBSERVER is negative in all variants (-3.0 to -6.0 in the window,
  -2.4 and -2.6 over the full canon). P1-P3 keep their registered verdicts under the `Hecke` fix and under the full
  corrected lexicon (P1 at 23.1% / 20.6%).
- Correction. The registered verdicts stand as registered; the readings must change.
  - N lines 41-50, replace item 1 by: "**The census refutes the modular family as the bridging hub (HYP-9437: P2 and P3
    REFUTED as registered; P1 passes).** The registered MODULAR share, 53.8%, is inflated by a lexicon bug: the pattern
    `Hecke` (case-insensitive, no word boundary) matches 'checked'; with a word boundary the share is 23.1%. The one
    association robust to seeds, reading length, front matter and lexicon corrections is finite geometry <->
    tournaments (z = 7.3 recent, 13.4 overall: Paley, quadratic residues, Fano). The registered Collatz <-> modular
    pair (z = 3.3 recent, 5.1 overall) is an artifact of the same bug (1.4 and 0.8 after correction), and the binary
    family's all-canon lead does not survive the corrected lexicon. The OBSERVER family, mostly the word 'carrier' in
    this window, is anti-associated with Collatz in every variant; its other properties depend on 'carrier'."
  - N lines 51-56: delete "Both hubs are SL_2 phenomena at two levels"; keep the FINITE_GEOMETRY-MODULAR weakness as
    the reason for the THM-4839 bridges (the full canon gives FG-MODULAR z = -4.3, an anti-association, which the note
    should also report).
  - N lines 98-100 ("What survives. Two robust object-level hubs ..."): "What survives: one robust association (finite
    geometry with tournaments) and the isolation of 'carrier'/'observer' vocabulary from the Collatz notes."
  - H line 64: append "(registered lexicon; `Hecke` matches 'checked' in 98 files; 23.1% with a word boundary, verdict
    unchanged)". H lines 79-85 ("Reading"): same replacement as N lines 98-100. H line 4 (status): replace "(the modular
    family is pervasive but is not the bridging hub)" with "(not the bridging hub; the registered 53.8% share is
    inflated by a lexicon bug, 23.1% corrected)".
  - Results INDEX line 8 and hypotheses INDEX line 2: replace "the hubs are finite geometry-tournaments and
    Collatz-modular" and "COLLATZ-MODULAR z = 3.3 (5.1)" accordingly.
  - Leave the registered script unchanged; cite this folder's corrected-lexicon runs as exploratory.

**MAJOR-2. THM-4839 (d) names the wrong modular family: the Pythagorean Frey curves are parametrised by X(4), a thin
subfamily of Kubert's Z/2 x Z/4 family X_1(2,4). (T lines 52-53; N lines 66-67 and 152.)**
- Claims. T: "These curves form the genus-zero family with torsion containing Z/2 x Z/4 (Kubert 1976)". N line 152:
  "Such curves are classified by a genus-zero modular curve for Gamma(2) cap Gamma_1(4)".
- Evidence.
  - The order-4 point is correct (verified exactly for 6560 (triple, leg-order) pairs; torsion is exactly
    Z/2 x Z/4 for the first six curves, with gcd #E(F_p) = 8; no Z/2 x Z/8 occurs among the 6560).
  - The Z/2 x Z/4 family is strictly larger. y^2 = x(x+1)(x+4) has (0,0) halvable (2(2,6) = (0,0)), so its torsion
    contains Z/2 x Z/4, and it is exactly Z/2 x Z/4. But none of its six lambda-values {4, 1/4, -3, -1/3, 4/3, 3/4} has
    both lambda and 1 - lambda squares, so it is not a Pythagorean Frey curve, nor a twist of one.
  - The Pythagorean Frey curves are the Legendre curves with lambda = k^2 and 1 - lambda = k'^2 (k = a/c, k' = b/c),
    i.e. exactly the image of the Klein-four map of (c). That map **is the modular cover X(4) -> X(2)**:
    - with t(tau) = theta_4^2/(theta_3^2 + theta_2^2) one has lambda(t(tau)) = lambda(tau) to 2e-31;
    - t is invariant under six elements of Gamma(4) and transforms by t -> 1/t, -t, 1/t, -1/t under four elements of
      Gamma(2) \ +-Gamma(4) (`audit_thm4839_independent.out`, section D);
    - this is forced in any case: Gamma(2)-bar is free of rank 2, so it has a unique normal subgroup with quotient V_4,
      namely Gamma(4)-bar.
  - X(4) is a double cover of X_1(2,4) = X(Gamma(2) cap Gamma_1(4)), so the Pythagorean curves are the image of a
    degree-2 map of genus-0 curves inside Kubert's family: a thin subset. Over Q(i) they have full 4-torsion, since all
    three 2-torsion points are halvable there (-a^2, -b^2, -c^2 become squares).
- Correction.
  - T lines 52-53: "These curves are the Legendre curves with lambda = k^2 and 1 - lambda = k'^2, i.e. the fibres over
    the rational points of the Klein-four cover of (c), which is the modular cover X(4) -> X(2) (t is a Hauptmodul for
    Gamma(4)). They have Z/2 x Z/4 inside E(Q) and E[4] inside E(Q(i)), and are parametrised by X(4) (genus 0), a double
    cover of Kubert's Z/2 x Z/4 curve X_1(2,4) = X(Gamma(2) cap Gamma_1(4)). They are therefore a thin subfamily of the
    Z/2 x Z/4 curves: y^2 = x(x+1)(x+4) is one that is not Pythagorean. CITED: Kubert 1976 (the Z/2 x Z/4
    parametrisation); Darmon-Granville 1995, p. 516 and section 7.1 (the (2,2,2) cover X(4) -> X(2))."
  - N line 152: the same in one sentence. N line 67: "... live on the genus-zero modular curve X(4): elliptic curves with
    torsion containing Z/2 x Z/4."
  - Related literature to cite: Halbeisen-Hungerbuehler, "Pairing pythagorean pairs", J. Number Theory 233 (2022)
    467-480. To a pythagorean pair (a, b) they attach Gamma_{a,b}: y^2 = x(x + a^2)(x + b^2), with torsion
    Z/2 x Z/4. Their Proposition 1, via Kubert p. 217, shows that **every** curve with torsion Z/2 x Z/8 is
    Gamma_{a^2,b^2} for a pythagorean pair. Pythagorean parametrisations of Kubert families are therefore in print,
    although they are not this family.

### MINOR

**m1. "Baker's bound sends q -> 1 along convergents" is not a theorem (N line 160; bridges script line 504/out line 59).**
Baker gives |2^K - 3^L| >= max(2^K, 3^L) K^(-C). That makes log(max)/log(6d) -> 1, which is a *lower* bound for
q = log(max)/log(6 rad d). q -> 1 needs rad(2^K - 3^L) >= max^(1-o(1)), which is abc-type and open (even the squarefree
part of 2^n - 1 is not controlled). Replace with: "Baker's bound makes log(max)/log(6|2^K - 3^L|) -> 1, so the clocks are
hits of quality close to 1 unless |2^K - 3^L| has a large square part; q -> 1 itself would follow from abc and is not
known."

**m2. "the first non-planar regular-map quotients of the three triangle groups" (N lines 133-134; T lines 38-39, title).**
- For (2,3,5), the hemi-icosahedron {3,5}_5 has skeleton K6 (6 vertices, non-planar), smaller than Petersen.
- For (2,3,6), the first non-planar {3,6} map is K7 and the first non-planar {6,3} map is K3,3.
- The small-map scan in section B of `audit_thm4839_independent.out` confirms the statement only by type.
- K5 = {4,4}_(2,1), K7 = {3,6}_(2,1) and Heawood = {6,3}_(2,1) are **chiral** (orientably regular, not reflexible;
  checked by the monodromy group acting regularly on darts and by the absence of a mirror automorphism), while
  {4,4}_(2,2), {6,3}_(1,1) and {3,6}_(2,0) are reflexible. "Regular" is Coxeter-Moser's sense.
- Correction: "Petersen, K5 and K3,3 are the smallest non-planar simple skeletons of regular maps of types {5,3}, {4,4}
  and {6,3} (the dual types {3,5} and {3,6} give K6 and K7); the {4,4}_(2,1), {3,6}_(2,1) and {6,3}_(2,1) maps are chiral."

**m3. "K4, K5, K7 are the only complete unit Cayley graphs of imaginary-quadratic residue fields" omits N(P) >= 4.**
(T title line 3; N line 58; results INDEX line 8; hypotheses INDEX line 2.) Every field with a prime of norm 2 or 3 gives
K2 or K3, because the units +-1 already fill F_2^x and F_3^x. The body of T (line 14) has the hypothesis. Add "with
N(P) >= 4" everywhere.

**m4. "Checked over the class-number-one fields d in {1,2,3,7,11,19,43}" (N line 118).** The class-number-one list also
contains 67 and 163. The proof needs no class-number restriction. Audit check: all 611 fundamental discriminants with
|D| <= 2000 (class numbers 1 to 56), all p <= 200, explicit Cayley graphs for N(P) <= 60. Exactly (Q(sqrt-3), 4),
(Q(i), 5) (both primes over 5) and (Q(sqrt-3), 7) (both primes over 7). Write "checked over all fundamental
discriminants |D| <= 2000, p <= 200 (audit)", or name the seven fields without "the class-number-one fields".

**m5. "Petersen is the hemi-dodecahedron, a quotient of the (2,3,5) triangle group by A_5 = PSL_2(F_5)" (N line 61) is
garbled.** The (2,3,5) rotation group *is* A_5. The correct statement is the one the audit verified: the antipodal
quotient of the dodecahedral map has 6 pentagons, 10 vertices and 15 edges on RP^2. Exactly 60 of the 120
automorphisms of the Petersen graph preserve its face set; they act regularly on the 60 flags, with element orders
{1:1, 2:15, 3:20, 5:24} (A_5). Replace with: "Petersen is the hemi-dodecahedron {5,3}_5, whose map group
[5,3]/{+-1} = A_5 = PSL_2(F_5) acts regularly on its 60 flags." The same applies to T lines 34-35 ("a quotient of the
(2,3,5) triangle group" is vacuous).

**m6. Internal prior art for (a)/(b) is not credited (N lines 137-138, T "Novelty").**
- Petersen = hemi-dodecahedron is already in the repository:
  - `noble_polyhedra_hill_fissary_20261002.md` line 182 ("The antipodal map turns the dodecahedron into the
    hemi-dodecahedron, whose graph is the Petersen graph K(5,2)");
  - `odd_zeta_parallels_petersen_20261001.md` line 42;
  - `petersen_product_4polytope_attempt_20261001.md` line 80;
  - `pythagorean_form_conic_graphs_rp_margin_20261018.md` line 155 (the hemi-dodecahedron's automorphism group A_5 of
    order 60).
- K7 on the torus as "the chiral regular map {3,6}_(2,1)", with dual the Heawood graph:
  `camion_busch_gaps_polyhedra_collatz_20261002.md` section 6.4, lines 285-299.
- The N(P) = 7 case of (a): `torus_voltage_and_halfplane_forcing_20261020/geometric_forcing_bridges.md` section 3
  (G2), lines 32-48. There the triangular lattice modulo an index-7 sublattice is K7, and "the six steps have the six
  nonzero images modulo seven".
- Related: THM-4650 (Z/2 x Z/4 torsion at the rational-diagonal rectangles).
- Add these to the "repository already had" sentence. THM-4839's novelty line "(a) ... not found in the repository"
  should read "(a) for N(P) in {4, 5}; the N(P) = 7 case is torus_voltage G2".

**m7. Literature prior art for (c) is not cited.** In Darmon-Granville 1995 (the paper MISTAKE-691 cites):
- Proposition 3.1 (p. 525): Galois covers of signature (p,q,r) have degree 2/chi; for (2,2,2) this is degree 4,
  genus 0.
- p. 516: "the covering X(2p) -> X(2) which is of signature (p,p,p)" (p = 2: X(4) -> X(2)).
- Section 7.1 (p. 536): dihedral covers pi_{2,2,r} with group {t -> zeta t, t -> 1/t}. For r = 2 the invariant
  u = (t^2 + t^-2)/2 satisfies lambda(t) = (u - 1)/(u + 1) exactly.

So "Euclid's parametrisation squared is the Klein-four Belyi map" is the r = 2 dihedral case of DG section 7.1 and the
cover X(4) -> X(2). T's "New to the repository: ... the Klein-four reading of Euclid's parametrisation" may stand, but
type the map CITED (classical; DG 1995) and say it is X(4) -> X(2).

**m8. Census reproducibility and fragility (N section 6, H "Results").**
- (i) The `.out` is not bit-reproducible: the curveball shuffles lists built from Python sets of strings, whose order
  depends on PYTHONHASHSEED. The session script, run unchanged on the registered corpus with hash seeds 0, 1 and 2, gives:
  - FG-T z = 7.38, 7.13, 7.17 (`.out`: 7.29);
  - at hash seed 0, ABC-ZETA (2.15) displaces BINARY-TOURN from the top 3.
  - The verdicts are unchanged.
- (ii) Running the documented command at HEAD reads 321 files, including THM-4839 itself; the registered run had 320
  (the corpus of 593105aedd). On Windows with this checkout path the script crashes on THM-4839's 205-character file
  name (MAX_PATH). Document "run at 593105aedd, or exclude THM-4839".
- (iii) P2's refutation depends on the arbitrary BODY_CHARS. At BODY = 12000, MODULAR ties FINITE_GEOMETRY at
  lift-degree 3, so P2 would be CONFIRMED (ties count). Several associations live in the front matter (titles,
  status, `related:` file names): with the front matter stripped, COLLATZ-MODULAR = 0.11 and BINARY-TOURN = 0.05.
- (iv) "5N swaps" (H line 40) are 5N *attempted* trades; trades with a = b or nothing to exchange are no-ops.

State these caveats next to the verdicts.

**m9. The OBSERVER family is mostly the word "carrier" (N lines 48-50).**
- In the window it matches `carrier` (42 files, 37 only), `observer` (14), `guarded` (8), and nothing else: Myhill,
  Nerode, "minimal automaton" and "deterministic quotient" have 0 hits. "Carrier" is general repository vocabulary,
  not the concurrent lanes' observer framework.
- Without `carrier` (strict lexicon), OBSERVER is no longer the least bridging family (QUADRATIC is), OBSERVER-TOURNAMENT
  = +2.14 and MODULAR-OBSERVER = -1.09. Only COLLATZ-OBSERVER stays significantly negative (-3.0).
- Correction: "The OBSERVER proxy (mostly 'carrier') is anti-associated with Collatz in every variant; its isolation
  from modular vocabulary and its last place are not robust to dropping 'carrier'."

**m10. MISTAKE-691 (3) mischaracterises THM-3833 (M lines 21087, 21099-21105, 21110-21111).** THM-3833 was read in
full; it never mentions Darmon, Granville, Faltings, Beukers or Chevalley, nor do its script/output,
CORE-PAPERS-ABC-IUT.md or the 08-23 reflection.
- Its Theorem 4.1 is a *true* abc-conditional implication.
- It is more than finiteness: (18) gives the height bound (z^r)^(kappa/2) <= K_eps, which DG's Faltings-based
  theorem does not give.
- THM-3833 section 5, (24), already proves that 1/42 is the least hyperbolic curvature. The uniform statement is one
  line from (18) + (24); it is DG 1995 section 5.2 (p. 533: "1/p+1/q+1/r <= 41/42 ... z^r <= kappa_{1/83}^{83} ...
  This implies the 'Fermat-Catalan' Conjecture"). DG themselves note (p. 527) "Theorem 2 may be deduced directly from
  the abc-conjecture".
- The real defect is the omission of the unconditional theorem in a file titled "what ABC would actually buy". Its own
  section 1 does flag the analogous S-unit case.
- The repository already cited DG's fixed-signature finiteness in THM-4516 (2026-09-29, line 14: "one of finitely many
  (Darmon-Granville)") and DG section 5.2 in `abc_difference_families_20261011c.md` line 34.

Correction:
- Header: "... and an unconditional finiteness theorem omitted from an abc-conditional statement (THM-3833)".
- Item (3): "THM-3833 states abc => fixed-signature finiteness (true) without recording that this finiteness is
  unconditional (DG 1995, Thm 2). The abc-conditional content of its Theorem 4.1 is the height bound (18). The uniform
  version, which follows from (18) and THM-3833's own (24), is DG 1995 section 5.2."
- "How it was caught": credit THM-4516 and the 10-11 packet, which already cite DG.

**m11. MISTAKE-691 (2) is a credit omission, not an "error" (M lines 21096-21098).**
- Verified: no file of the `twin_pythagorean_124_20261026` packet, `SOURCES.json` included, cites
  `abc_difference_families_20261011c`. Its eq. (8)-(9) duplicate that packet's eq. (2)-(3) and its "quality < 1"
  remark.
- But the identity is a one-line difference of squares, and the packet's README line 5 disclaims novelty ("no
  historical novelty claim").
- Retype as "uncited duplication (no novelty claimed)" and drop "error" from the header for this item.

**m12. Genealogy: the consumer lineage is not a citation chain (N lines 69, 230-231, "Each step cites the
previous").** `aliquot_power_124_20261026` does not cite the 10-11 Codex packet; it cites THM-3833 directly
(`chains/RESULT.md` line 22). The other links hold:
- THM-3848 cites THM-3833;
- THM-4108 depends on THM-3833 and relates THM-3848;
- the 10-11 packet cites THM-3833 and THM-4108;
- THM-4805 cites aliquot124 and THM-3833.

Write "each member cites THM-3833; the chain is not linear".

**m13. Paraphrase errors in the genealogy.**
- (a) N line 218: session69's bijection is between square-later integers and primitive positive abc triples **with
  B = 2^a >= 2 and C odd** (`session69-power/REPORT.md` lines 262-263), not all primitive triples. Line 78 has it right.
- (b) N line 75 ("Collatz clocks are its {2,3}-unit points") and line 242 ("lineage 4 is the S-unit equation over
  Z[1/6]") are imprecise. A clock triple (2^K, |2^K - 3^L|, 3^L) has a non-unit middle term (13, 139, 7153, ...). Only
  3/2 (1 + 8 = 9) is a {2,3}-unit-equation solution; the others are Pillai-type 2^K - 3^L = d.
- (c) "Baker beats abc" (line 242) and "which this session confirms" (line 71) overstate the atlas's proved "abc is
  provably no stronger than Baker" (line 35). The session computed qualities and proved nothing about the gate.
- (d) N line 72 / lines 219-222: MISTAKE-678's abc item was a non-connection asserted beyond its evidence, not a
  textbook fact typed new. MISTAKE-676's abc item concerned repository-internal facts (THM-3833, CORE-PAPERS,
  THM-4640). MISTAKE-604 concerns Gersonides / the trivial cycle, not Baker. MISTAKE-676 and 678 are mac-mini lanes,
  not THM-4805 (S60); move them out of the THM-4805 bullet.

### NOTE

- n1. Section 4 table (N line 166): "the ten divisors of p^2qr" are the ten *proper nontrivial* divisors
  (d(p^2qr) = 12), and "A_4 roots" are the 10 positive roots (root lines) of A_4, which has 20 roots.
- n2. T line 4 and N line 12 say "independently audited" while N section 7 is "(pending)". Set these after this audit's
  corrections are applied.
- n3. (a) literature: the graphs are the Gaussian networks of Martinez, Beivide et al., IEEE Trans. Comput. 57 (2008),
  and the Eisenstein-Jacobi networks of Martinez, Stafford, Beivide and Gabidulin, Probl. Inf. Transm. 44 (2008).
  Their diameter-1 members are K5 = G(2+i) (the "dense Gaussian network" alpha_k = k + (k+1)i at k = 1) and K7. No
  literature novelty is claimed for (a), but add the pointer.
- n4. Convention check passed. The session's (b,c) is Coxeter-Moser's (vertex counts b^2+c^2 and b^2+bc+c^2). In
  Wedd's/Conder's torus tables the same maps appear as {3,6}(1,3) (K7), {6,3}(1,3) (Heawood) and {6,3}(0,2) (K3,3),
  with coordinates (|b-c|, b+c); {4,4}(2,1) = K5 agrees.
- n5. "Every clock from 3/2 on is an abc hit" (N lines 75 and 158; T line 62) is true for **all** convergents and
  semiconvergents, not only up to 84/53, but the note gives no proof.
  - Proof: for clocks with L >= 2 the error |L log2 3 - K| is at most 0.2451 (attained at 5/3, which is
    |3 log2 3 - 5|), and 0.2451 < log2(6/5) = 0.2630. Hence |2^K - 3^L| < max/6, and so 6 rad(d) < max.
  - The audit verified all 40 clocks with L <= 20000. Add the proof, or state the checked range.
- n6. DG's proof of Thm 2 (pp. 526-527) uses Beckmann's ramification lemma, Hermite-Minkowski and Faltings. "Via
  Faltings and Chevalley-Weil" (M line 21102) names the principle and is acceptable.
- n7. "abc is the arithmetic of the thrice-punctured line" (N lines 74 and 240) is the classical Elkies (1991) / Vojta
  viewpoint. Cite it; the INTERPRETATION typing (line 248) is correct.
- n8. "27/17, the 231 -> 233 clock" (N line 75): 27/17 has 12 recorded windows
  (`collatz_procgen_20260922_order_laws.md` line 104), of which 231 -> 233 is one.
- n9. H lines 23-25 put Delta(2,4,4) among the modular group's "triangle-group quotients". Delta(2,4,4) is not a
  quotient of PSL_2(Z) = C2 * C3 (its abelianisation Z/2 x Z/4 is not a quotient of Z/6); it is a quotient of
  Delta(2,4,oo), which contains Gamma_0(2). The (2,3,6) and (2,3,p) cases are fine.
- n10. The S116 statement also sits in the generator `04-computation/famous_problems_s116.py` line 76, outside the four
  directories MISTAKE-691 (1) scopes. Its "only" claim is otherwise verified (grep of `mochizuki` and of abc/IUT
  "proved" phrasings in 01-canon, 05-knowledge, 07-reflections and 00-navigation).
- n11. "No link between graceful labellings and Shabat polynomials is known" (N line 174). None was found by the audit
  either (searches including arXiv 1811.07614 and 2510.10192). Write "none found".
- n12. The lexicon also has false negatives: THM-4686, a CM-supercongruence/2-isogeny theorem, matched no genuine
  MODULAR term in its first 6000 characters. The `Baker` hits behind ABC-QUADRATIC (z = 1.87) are four Heegner-Baker-Stark
  class-number citations, not abc.

---

## Recomputed numbers

| quantity | session | audit (own code) | status |
|---|---|---|---|
| complete unit Cayley graphs, N(P) >= 4 | (Q(sqrt-3),4), (Q(i),5), (Q(sqrt-3),7); 7 fields, p < 60 | same; all 611 fundamental D, \|D\| <= 2000 (h <= 56), p <= 200 | correct |
| {4,4}_(2,1), {4,4}_(2,2) | K5, K4,4 | K5 (V,E,F = 5,10,5; chiral); K4,4 (8,16,8; reflexible) | correct |
| {3,6}_(2,0), {3,6}_(2,1) | doubled K4, K7 | 4,12,8 with multiplicity 2; K7 (7,21,14; chiral) | correct |
| {6,3}_(1,1), {6,3}_(2,1) | K3,3, Heawood | K3,3 (6,9,3; reflexible); Heawood (14,21,7; chiral) | correct |
| girth-Euler bound | 5 graphs violate | K5, K3,3, K7, Heawood, K4,4, Petersen, K6 all violate | correct |
| hemi-dodecahedron | Petersen, A_5 | Petersen; 60 of 120 automorphisms keep the 6 faces; flag-regular; orders {1:1,2:15,3:20,5:24} | correct |
| Klein-four map | passport (2^2,2^2,2^2), R-H | same; numerical monodromy V_4; equals X(4) -> X(2) | correct (+) |
| Frey point (a(a+c), ac(a+c)) | order 4, 6000 pairs | order 4, 6560 pairs; torsion exactly Z/2 x Z/4 (first six); no Z/2 x Z/8 | correct |
| Berggren abc hits, depth <= 9 | 60; best q = 1.2659 (36207, 18424, 40625) | 60; q = 1.265900 | correct |
| next hits | (1377,3136,3425), (2793,10624,10985) | q = 1.2425, 1.2182 | correct |
| (7,24,25); ADU (175,288,337) | q = 1.204; among hits | 1.2040; present, q = 1.0424 | correct |
| primitive triples c <= 2000 | 319 | 319 | correct |
| Collatz clocks 3/2 ... 84/53 | 1.226, 1.019, 1.273, 1.143, 1.236, 1.086, 1.064, 1.063, 1.081 | 1.2263, 1.0190, 1.2728, 1.1433, 1.2359, 1.0861, 1.0637, 1.0633, 1.0814 | correct |
| clocks with 2 <= L <= 20000 | "from 3/2 on" | all 40 have \|2^K-3^L\| < max/6 (hits) | correct; proof in n5 |
| max hyperbolic 1/p+1/q+1/r | 41/42 only at (2,3,7) | 41/42, unique (2,3,7) | correct |
| 1+8=9 failure of rad^{1+eps} > c | eps < 0.2263 | eps <= 0.22629439 | correct |
| census corpus | 320 files | 320 (HEAD minus THM-4839; 321 at HEAD) | correct (see m8) |
| census observed n_ij, shares | `.out` | identical | correct |
| z (registered): FG-T, C-M, BA-T | 7.29, 3.31, 2.06 | 7.12, 3.37, 2.13 (10 seeds: 6.93-7.37, 3.30-3.43, 2.07-2.20) | within Monte Carlo noise |
| MODULAR share | 53.8% | 53.8% registered; **23.1%** with `\bHecke\b`; **20.6%** corrected lexicon | artifact (MAJOR-1) |

## Census implementation check

- Curveball: correct (Strona et al.). Each trade keeps both row sums and column sums; asserted in the audit
  implementation. Samples are successive states of the chain with 5N attempted trades between them; the chain is
  uniform on the fixed-margin bipartite graphs (Carstens 2015).
- z = (n_obs - null mean)/null sd (ddof = 1), lift-degree = #{partners with z >= 2} and strength = sum of positive z:
  all as described.
- P1 = MODULAR share >= 0.20, P2 = MODULAR's lift-degree equals the maximum, P3 = some top-3 pair lies inside the core:
  all as registered.
- Verdicts are read correctly: P1 CONFIRMED, P2 REFUTED (TOURNAMENT 2 vs MODULAR 1), P3 REFUTED.
- Pre-registration integrity. `git show --stat 593105aedd` contains the script and the HYP. `git diff 593105aedd
  4e32454073` shows the script unchanged, and in the HYP only the status line changed and the Results section was
  appended; the Protocol and Predictions text is untouched. The `.out` reports 320 files, the corpus at 593105aedd
  (THM-4839 was still a RESERVED stub), so the run used the registered corpus. The commit order (10:28:44 pre-register,
  10:40:23 results) is consistent with the claim. Git cannot prove the run happened after the commit, but nothing
  contradicts it.

## Robustness grid (`audit_census_robustness.out`; exploratory, the registered verdict stands)

Strict lexicon = registered with `\bHecke\b`, `\bparity\b`, `\bcharacters?\b`, Markov only for
number/triple/equation/tree/spectrum/level, Steiner only for system/triple/quadruple/rule, Baker only for
theorem/bound/theory/method/estimate, and Minkowski not followed by "sum" or "bound".

| run | FG-T | C-M | BA-T | BA-FG | C-O | M-O | FG-M | max lift-degree | least | P1 P2 P3 |
|---|---|---|---|---|---|---|---|---|---|---|
| registered (320) | 7.12 | 3.37 | 2.13 | 1.87 | -5.86 | -3.17 | 1.24 | 2: TOURNAMENT | OBSERVER | T F F |
| registered, seeds 1-9 | 6.93-7.37 | 3.30-3.43 | 2.07-2.20 | 1.84-1.98 | -5.51 to -5.92 | -3.11 to -3.30 | 1.13-1.20 | 2: TOURNAMENT | OBSERVER | T F F |
| THM-4839 included (321) | 7.09 | 3.40 | 2.18 | 1.72 | -5.80 | -3.20 | 1.27 | 2: TOURNAMENT | OBSERVER | T F F |
| BODY 3000 | 6.36 | 2.59 | 0.89 | 2.62 | -4.62 | -3.14 | -0.33 | 3: FIN_GEOM | OBSERVER | T F F |
| BODY 12000 | 6.67 | 2.78 | 1.20 | 1.86 | -5.80 | -4.05 | 2.08 | 3: MODULAR, FIN_GEOM | OBSERVER | T **T** F |
| BODY whole file | 6.19 | 2.50 | 1.13 | 2.17 | -5.52 | -4.05 | 2.25 | 4: FIN_GEOM | OBSERVER | T F F |
| front matter stripped | 5.74 | 0.11 | 0.05 | 2.05 | -4.06 | -1.43 | 0.98 | 3: FIN_GEOM | OBSERVER | T F F |
| only `Hecke` fixed | 7.39 | **1.44** | 2.24 | 1.93 | -5.75 | -2.30 | 1.17 | 2: TOURNAMENT | OBSERVER | T F F |
| strict (seeds 20261028, 1, 2, 3) | 6.88-7.10 | **1.38-1.46** | 1.79-1.87 | 2.50-2.59 | -5.65 to -5.97 | -2.18 to -2.28 | 0.99-1.07 | 2: FIN_GEOM, BINARY | OBSERVER | T F F |
| strict, no 'carrier' | 6.72 | 1.00 | 1.28 | 2.11 | -3.00 | -1.09 | 0.75 | 2: FIN_GEOM, BINARY, TOURN | QUADRATIC | T F F |
| strict, front matter stripped | 5.66 | -0.05 | -0.10 | 2.09 | -4.06 | -0.44 | 0.96 | 3: FIN_GEOM | OBSERVER | T F F |
| strict, BODY 3000 | 6.22 | 0.83 | 0.85 | 3.28 | -4.84 | -2.69 | 0.43 | 3: FIN_GEOM | OBSERVER | **F** F F |
| strict, BODY 12000 | 6.53 | 0.54 | 1.03 | 2.72 | -5.23 | -1.43 | 0.41 | 3: FIN_GEOM | OBSERVER | T F F |
| full canon, registered (300 null) | 13.80 | 4.92 | -3.46 | -1.69 | -2.41 | -2.36 | -4.12 | 5: BINARY | OBSERVER | T F F |
| full canon, strict (300 null) | 14.89 | **0.75** | 1.33 | 0.53 | -2.62 | -1.99 | -3.63 | 3: MODULAR, DIVISOR, ZETA | OBSERVER | **F T** F |

Qualitative findings asked about:
- FINITE_GEOMETRY-TOURNAMENT strongest: **survives everywhere.**
- COLLATZ-MODULAR positive: **does not survive** the `Hecke` fix, front-matter stripping or the corrected lexicon.
- OBSERVER least bridging: survives unless 'carrier' is dropped.
- OBSERVER anti-associated with COLLATZ: **survives everywhere** (-2.4 to -6.0).
- OBSERVER anti-associated with MODULAR: not robust (-0.4 to -4.1).

## abc genealogy and MISTAKE-691: the "checked" claims

| claim (N section 5 / M) | file evidence | verdict |
|---|---|---|
| S116 `.out` lines 67-68 say abc "proved by Mochizuki" and "rad(abc)^{1+eps} > c for all eps" | `famous_problems_s116.out` lines 67-68 verbatim; git: added 2026-03-15 (kind-pasteur S116) | correct |
| the inequality fails at 1+8=9 for eps < 0.2263 | 6^(1+eps) > 9 iff eps > 0.226294 | correct |
| only such statement in 01-canon / 05-knowledge / 07-reflections / 00-navigation | grep: only the S116 `.out` (+ quotes in MISTAKES, the note, INDEX); the generator `.py` line 76 is outside scope | correct (n10) |
| THM-3833 types fixed-signature finiteness as abc-conditional; no DG/Faltings acknowledgement | status lines 8-10, Theorem 4.1; no Darmon/Granville/Faltings/Beukers/Chevalley in THM-3833, its script/out, CORE-PAPERS-ABC-IUT, 08-23 reflection | wording correct; framing unfair (m10) |
| DG 1995 Thm 2 is unconditional fixed-signature finiteness | DG p. 515, Theorem 2 (read in the McGill scan) | correct |
| abc => Fermat-Catalan (uniform) is standard; 41/42 | DG section 5.2, p. 533 (41/42, kappa_{1/83}^{83}); also DG p. 527 | correct; cite section 5.2 |
| 41/42 maximal, equality only at (2,3,7) | exhaustive check; also THM-3833 (24) | correct |
| twin_pythagorean eq. (8) re-derives abc_difference_families eq. (2), uncited | `arithmetic/README.md` line 164 vs `abc_difference_families_20261011c.md` line 87; no citation anywhere in the packet | correct; severity (m11) |
| Collatz atlas line 35 | `procgen_atlas_20260924_collatz_implication_atlas.md` line 35 verbatim; Prop. R2 line 34 | correct |
| session69 square-later <-> primitive abc triples | REPORT lines 250-283: bijection onto triples with B = 2^a, C odd; missing coordinate rad(u - 2^a) | restriction dropped (m13a) |
| THM-4805: (1, 2^k-1, 2^k) hit iff 2^k-1 not squarefree; MISTAKE-675/676/678 | statement correct (2 rad(2^k-1) <= 2(2^k-1)/3 < 2^k when non-squarefree); 675 = 2^(6n) family typed new; 676 = radical-6 novelty (mac-mini, HYP-9416); 678 = ADU hit with chi = -1 (mac-mini, THM-4808) | attribution blurred (m13d) |
| THM-3848 / Mason-Stothers JC2 lane same day | THM-3848, 3881, 3885, 3920 all sourced 2026-08-23, Mason-Stothers present | correct |
| session82 "no novelty claim"; session83 Grantham | DECIMAL_TRIPLE.md line 108; session83 README line 94 | correct |
| "Each step cites the previous" | aliquot124 cites THM-3833, not the 10-11 packet | wrong (m12) |

## Prior-art search log

- **Darmon-Granville**, Bull. LMS 27 (1995) 513-543, `math.mcgill.ca/darmon/pub/Articles/Research/12.Granville/pub12.pdf`
  (scanned PDF, rendered locally). Relevant passages:
  - Thm 2, p. 515;
  - p. 516, X(2p) -> X(2) of signature (p,p,p);
  - Prop. 3.1, p. 525;
  - proof of Thm 2, pp. 526-527 (Beckmann, Minkowski, Faltings; "may be deduced directly from the abc-conjecture");
  - section 5.2, p. 533 (abc => Fermat-Catalan with 41/42);
  - section 7.1, p. 536 (dihedral (2,2,r) covers; r = 2 is the Klein-four map).
- **Halbeisen-Hungerbuehler**, "Pairing pythagorean pairs", J. Number Theory 233 (2022) 467-480: Gamma_{a,b} with
  Z/2 x Z/4; Prop. 1, every Z/2 x Z/8 curve is Gamma_{a^2,b^2} (Kubert p. 217). Also Halbeisen-Hungerbuehler-Zargar,
  arXiv:2405.12989 (Gamma_{a^h,b^h}, torsion Z/2 x Z/4).
- **Gaussian / Eisenstein-Jacobi networks**: Martinez, Beivide, Stafford, Moreto, Gabidulin (IEEE TC 2008); Martinez,
  Stafford, Beivide, Gabidulin (PIT 2008); dense Gaussian networks alpha_k = k + (k+1)i. Unit Cayley graphs of
  Z[i]/(alpha) and Z[omega]/(alpha); K5 at alpha = 1 + 2i (associate of 2 + i).
- **Regular toroidal maps**: Coxeter-Moser ch. 8; Wedd's/Conder's census (weddslist.com/rmdb): K5 = {4,4}(2,1) chiral,
  K7 = {3,6}(1,3) chiral with C7:C6, Heawood = {6,3}(1,3), K3,3 = {6,3}(0,2) (coordinate change, n4). Hemi-dodecahedron
  = Petersen is classical (McMullen-Schulte, *Abstract Regular Polytopes*).
- **abc hits from 2^K vs 3^L**: no source found stating the clock family as such; related are de Weger arXiv:2602.08051
  (approximation gain, convergents of surds) and Sankaran arXiv:2606.08416 (powers-of-2-and-3 family). The session
  claims no novelty here.
- **Graceful trees vs Shabat polynomials**: no link found (arXiv:1811.07614, arXiv:2510.10192).
- **Repository greps** (01-canon, 05-knowledge, 07-reflections, 00-navigation):
  - "hemi-dodecahedr": 4 earlier notes;
  - `{3,6}_`: `camion_busch` (K7 chiral map);
  - "K7 ... torus": `torus_voltage` G2 and others;
  - "Belyi" together with Pythag/Euclid/Legendre/X(4): only S68 files;
  - "Frey curve": only S68 files;
  - "Kubert": the supersingular note;
  - "Z/2 x Z/4": THM-4650;
  - `1 + p(p+2)`: twin_pythagorean, THM-4808 (and the 10-11 packet in the m^2 form);
  - "Darmon/Granville": THM-4516, the 10-11 packet, the Collatz necklace notes.
- **Novelty typing.**
  - (a) "small PROVED observation, not found in the repository" must cite torus_voltage G2 for N = 7 and the
    Gaussian/EJ network literature.
  - (b) is classical and partly in the repository (m6).
  - (c) must be typed CITED/classical (DG 1995 section 7.1; X(4) -> X(2)).
  - (d)'s order-4 point is classical (the halving criterion); its modular-family sentence is wrong (MAJOR-2).

## How to rerun (from the repository root)

```text
python 04-computation/experiments/attractors_census_bridges_20261028_audit_S68/audit_thm4839_independent.py      # ~5 s, all asserts
PYTHONIOENCODING=utf-8 python 04-computation/experiments/attractors_census_bridges_20261028_audit_S68/audit_census_independent.py run     # registered, ~10 s
PYTHONIOENCODING=utf-8 python .../audit_census_independent.py terms                 # per-alternative counts
PYTHONIOENCODING=utf-8 python .../audit_census_independent.py spot --k 10 --seed 68 # spot check with contexts
PYTHONIOENCODING=utf-8 python .../audit_census_independent.py robust                # grid, ~5 min -> audit_census_robustness.out
PYTHONHASHSEED=0 python .../rerun_session_census.py 4400 4839 2000                  # the session's script, unchanged, minus THM-4839
```

Outputs in this folder:
- `audit_thm4839_independent.out`;
- `audit_census_R0_registered.out`, `audit_census_robustness.out`, `audit_census_terms.out`,
  `audit_census_spotcheck.out`;
- `rerun_session_census_hashseed{0,1,2}.out`.

The audit census code is deterministic for a given seed, independent of PYTHONHASHSEED. On Windows it opens long paths
with the `\\?\` prefix.
