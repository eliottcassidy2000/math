# Independent audit of opus-2026-11-06-S84 (THM-4890, THM-4891, HYP-9482, MISTAKE-728)

Auditor: independent Claude instance, 2026-11-06. All code in this folder is the auditor's own; no session module is
imported. Audited state: checkpoint 4a4d5cbfad, the addendum pre-registration 25fde96277 (HYP-9482 P9-P11), and the
uncommitted working-tree text of HYP-9482 and the README as of about 20:00 (the P9-P11 outcomes).

## Verdict: PASS WITH CORRECTIONS

The proved statements are correct, and every number checked was reproduced. There are 2 majors, 13 minors and 4
upgrades.
- **The majors are about how results are read, not about proofs.**
  - THM-4890 (D) is correct as a statement about a random model. The session reads it as "decorrelation". Exact
    computation shows the true micro values are positively correlated with the macro strips, by a factor of about 3.7.
  - THM-4891's title and index lines say the gaps "need zone-wide coherence". What is proved is weaker: neither face
    type alone excludes 7 or 21.
- **Main upgrade (U1).** The balance lemma's sharp form is now proved.
  - The sharp form is |A - S| <= 1 for every collision pair, with no exceptional region.
  - The windows always have the form 1...1 against 0...0, and X Xi < 0.
  - HYP-9482 P5 becomes a theorem for all n, and the session's next obligation 2 is closed.
  - Checked on 258,617,162 ordered collision pairs (n <= 28). A - S takes only the values -1, 0, 1.

## 1. THM-4890

**(a) Inherited identities (A), (B).** These are Session 193's laws (1) and (2); the substitution was checked against
PROOF.md lines 35-43. `aud_blocks.out`:
- g is computed by the MSB-first state machine and checked against the recursion.
- (A): exhaustive for m <= 8 and h < 2^9 (260,610 cases); 20,000 random big-integer cases with lengths up to 150.
- (B) with its parity term eps: exhaustive over every word of length 2..13 and every split with a, w >= 1, s >= 0
  (548,862 cases); 20,000 random cases with a, w, s <= 100.
- Also checked: 1 <= g(x) <= x + 1, with equality exactly at x = 0 and x = 2^k, k >= 1; C_w is an involution;
  g(C_w x) = g(x) for w <= 14 (leading zeros allowed); v -> v* maps Y_n onto Y_n.
- 0 violations throughout.

**(b) Balance lemma (C).** The proof is correct line by line.
- Ranges: a, s >= 1. The prefix h lies in [2^(a-1), 2^a) and tau = C_s(sigma) in [2^(s-1), 2^s).
- X and Xi are odd: the windows differ in bit 0 and bit w-1, and bit 0 of C_w(x) is 1 - x_(w-1).
- 2^(A+w-2) <= h|Xi| < 2^(A+w) and 2^(S+w-2) <= tau|X| < 2^(S+w).
- |hXi + tauX| <= 2^w + 1.
- Both sign cases are correct, and so is the corollary |a - s| <= |A - S| + |dXi - dX|.

Exhaustive check with the auditor's C code (`aud_pairs.c`, `aud_pairs_4_24.out`, `aud_pairs_25_28.out`):
- every ordered collision pair, converse pairs included, for n = 4..28;
- every intermediate inequality of the proof, the statement, the session's code test, P5 and (C');
- 0 violations.

The stronger facts the data show:
- A - S takes only the values -1, 0 and +1 at every n <= 28 (n = 28: 7,716,910 / 113,737,020 / 7,716,910 pairs).
- So |A - S| = 2 and the "max(A, S) <= 3" exception never occur.
- The same-sign case X Xi > 0 never occurs either.
- Both facts are proved in U1 below.

**(c) Converse symmetry (C').** Correct. Checked for every pair to n = 28: the image is a collision pair in Y_n, and
(a, s, dX, dXi) maps to (s, a, dXi, dX).

**(d) Randomised micro self-energies (D).**
- **Hypotheses.** All hold.
  - n even gives B = P/C = 8 and n odd gives B = 16. C = 2^m >= 16 for n >= 12.
  - THM-4887's strip lemma needs h' > 2W with W = C; here h' >= P >= 8C.
  - The per-pair bound (C/beta + 1)(2 beta C/h' + 3) <= 2C^2/P + 3C/beta + 5 was checked in exact arithmetic for
    beta <= B (`aud_D_constants.out`).
- **Arithmetic.**
  - The sum of 1/beta is THM-4887 step 4 verbatim. The l = l' pairs are at most 2P^2 = 2B PC.
  - P(F(l') - F(l) = t) = (C - |t|)/C^2 <= 1/C.
  - The bracket is at most 75.72, attained at n = 13. The total is at most 108.72 <= 109 for n = 12..80.
- **Exact expectation, not sampled** (`aud_randmodel.out`, n = 12..22, plus Monte Carlo at n <= 20).
  - The random-model energy is 1.64-1.99 |X_n|, so the constant 109 is loose by a factor of 55-66.
  - Type 3 is 0.04-0.55 PC, against the bound of 16-32 PC.
  - The deterministic split E_n = PC + T3 + DC matches a sorted energy exactly at every n.
  - The equivalence (DC) <=> linear energy is correct. It is also tautological: DC is E_n minus at most 33 |X_n|.
- **Deterministic DC against the random model.** This is the content of MAJOR-1.
  - DC is 3.31, 3.67 and 3.69 times the random-model expectation at n = 12, 20 and 22.
  - Against the model with the true marginal law of g(l) - 1, it is 1.64, 1.83 and 1.84 times. This reproduces the
    session's uncommitted `dc_ratio_*.out` digit for digit.
  - At n = 20, 0.965 |X_n| of DC is converse-partner pairs (y 1^k, y* 1^k), forced by g(v*) = g(v)
    (`aud_dc_split.out`). That alone exceeds the random model's entire type-4 expectation of 0.592 |X_n|.

**(e) Reduction (E).** Correct. Checked for n = 4..28 directly on all of X_n:
- the layer ranges g(Y_j) in [2^(j-2) + 1, 2^j - 1], with 0 violations;
- at most two layers per value;
- E_n <= 1 + 2 sum E^Y_j, with E_n/(1 + 2 sum) at most 0.558.

## 2. THM-4891

**(a) Proofs.** (i), (ii), (iii) and the recursion f(n) <= C(n,2) f(n-2) are correct.
- In a hexagon, "opposite edges" are E_k and E_(k+3) of the 6-cycle. They swap the same pair, with the third element
  on opposite sides, so "same way" means equal direction vectors x(head) - x(tail).
- Independent check (`aud_perm.py` part A, `aud_perm_AB.out`):
  - The Cayley graph of S_n was built directly, and 2-faces were generated from commuting and braid pairs of adjacent
    transpositions.
  - Union-find with parity was run on geometric opposite-edge relations, for n = 3..7.
  - The class counts are exactly C(n,2) 2^(n-2) for squares, C(n,2)(n-2)! for hexagons and C(n,2) for both.
  - Each class is exactly (pair, set before), (pair, relative order of the others) or the zone. There were no
    contradictions and no non-parallel opposite pairs.

**(b) Sink spectra at n = 3, 4.** Exhaustive enumeration of all class assignments, with acyclicity by Kahn's algorithm
(`aud_perm_enum.c`). Identical to the session's table, multiplicities included:
- n = 4, squares: {0..6}, acyclic {1..6}.
- n = 4, hexagons: {0..9, 12}, acyclic {1..9, 12}. 7 sinks occur for exactly 48 of 4,096 orientations, all acyclic.
- n = 3: squares {0..3}, acyclic {1..3}; hexagons {1, 3}.

**(c) n = 5.** A different encoding (edge variables with coherence equalities), a different solver (Glucose4) and a
different cardinality encoding (sequential counter) were used (`aud_perm_C.out`).
- Squares {0..25}, 26..60 UNSAT. Hexagons {0..54, 56, 60}; 55, 57, 58, 59 UNSAT.
- The decisive cases were rechecked with Minisat22 and the kmtotalizer encoding: squares 7, 21, 25 SAT and 26 UNSAT;
  hexagons 7, 21, 54, 56, 60 SAT and 55, 57, 58, 59 UNSAT.
- The auditor's own CEGAR found acyclic witnesses with 7 and 21 sinks for both relaxations. Each was verified
  independently, including "coherent on the stated face type and not on the other".
- The session's witnesses:
  - its two printed tables (n4_L2_7 and n5_L1_7) were rebuilt and verified;
  - for the three witnesses printed only as sink lists, each sink set is exactly realised by an acyclic orientation
    coherent on the stated face type.

**(d) Bound and values.**
- f(3) = 3, f(4) = 6 and f(5) = 25 are confirmed; 26 is UNSAT with two solvers.
- f(6) >= 90 holds by a Glucose4 / sequential-counter witness, independently checked (`aud_perm_D.out`). The witness
  can even be taken acyclic. With the bound, f(6) = 90.

## 3. HYP-9482

**(a) Pre-registration.** Consistent.
- 21b02648ec (19:01:43) holds the predictions and the in-sample n = 14..20 output only.
- The P1-P8 text is unchanged since then: the diff touches only the status line and appended outcomes.
- `s84_energy.py` and `s84_model.py` are byte-identical between 21b02648ec and 4a4d5cbfad.
- `outofsample_21_26.out` first appears in 4a4d5cbfad (19:17:14).
- The addendum 25fde96277 (19:19:29) precedes `outofsample_28.out` (mtime 19:22:11).
- Git cannot exclude uncommitted earlier runs. The auditor found no Y_n statistics for n >= 21 in the repository from
  before the registration (grep for E^Y and the model constants). The X_n energies for n = 22..28 do exist (see
  MINOR-4).

**(b) Independent recomputation** (`aud_pairs.c`). Every registered quantity matches `outofsample_21_26.out` and
`outofsample_28.out` to the printed precision.

| n | E^Y/\|Y\| | nonconv | values | q=1 | q=2 | q=3 | q=4 | K=8 | K=16 | max | cross | E/\|X\| |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 22 | 2.942013 | 0.942989 | 0.401259 | 0.317588 | 0.068725 | 0.011994 | 0.001733 | 0.902970 | 0.902656 | 14 | 0.324018 | 3.260717 |
| 24 | 2.934855 | 0.935343 | 0.401433 | 0.318090 | 0.069175 | 0.011823 | 0.001635 | 0.907489 | 0.906218 | 16 | 0.322637 | 3.257433 |
| 26 | 2.929600 | 0.929844 | 0.401942 | 0.319300 | 0.068806 | 0.011727 | 0.001649 | 0.908719 | 0.907978 | 18 | 0.321487 | 3.252571 |
| 28 | 2.924796 | 0.924918 | 0.402267 | 0.319909 | 0.068779 | 0.011634 | 0.001614 | 0.909357 | 0.908753 | 16 | 0.320490 | 3.247338 |

- The level law and the sharp centred-window check (0 violations) also match. The odd n = 27 gives 2.926717.
- The E/|X| column also equals HYP-9479's 3.2607, 3.2574, 3.2526 and 3.2473.

**(c) Model constants** (`aud_model.out`, mpmath, 30 digits).
- The density of XY, by convolution, is 4 ln(4t) on [1/4, 1/2] and 4 ln(1/t) on [1/2, 1].
- The integral of f^2 is 24 - 32 ln 2.
- E^Y/|Y| -> 14 - 16 ln 2 = 2.909645; the non-converse part and the smoothed slope are both 12 - 16 ln 2.
- values/|Y| -> 3/4 - (ln 2)/2 = 0.403426.
- Cross layers: 12 ln 2 - 8 = 0.317766 = 4 - 4 ln 2 - (12 - 16 ln 2). E/|X|: 6 - 4 ln 2 = 3.227411.
- The multiplicity formula agrees with direct integration of the Poisson mixture for q = 1..12.
- Tail at |Y| = 2^24: 282.0, 24.15, 1.843 and 0.1268 values of multiplicity 12, 14, 16 and 18.
  - So P(max >= 18) = 0.126. A 100-replicate Monte Carlo of the model gives 0.12.
  - MISTAKE-728's 0.13, 1.84, 24.2 and 282 are right.
- P11's 0.54 values of multiplicity >= 18 and 0.034 of >= 20 at 2^26 are right (0.5408 and 0.0337).

**(d) Verdicts.** P1, P2, P3, P4, P5, P6(a) and P8 CONFIRMED, P6(b) and P7 REFUTED, and P9-P11 CONFIRMED. All are
stated correctly. Only third-decimal ratio rounding is off (MINOR-5), and P5 is now a theorem (U1).

## 4. Prior art

- **The n!/2^floor(n/2) bound.** This is Session 188's elementary tournament bound with the same orbit proof
  (`session188_structural_bridges_20261104/equal_nine_ambient/PROOF.md` section 1, eq. (1)). It is not credited
  (MINOR-2).
- **The XY limit law and its density.** These are proved in Session 16 (`BINARY_PRODUCT_DENSITY.md` section 3). Not
  credited (MINOR-1).
- **No prior art found** for:
  - the balance lemma (repo grep for collision windows, common prefix/suffix and differing-bit statements);
  - the constants 14 - 16 ln 2, 12 - 16 ln 2, 6 - 4 ln 2, 12 ln 2 - 8 and 0.40343 (only numerical coincidences in
    unrelated outputs);
  - square-only or hexagon-only coherence, context tournaments, chosen sides of Boolean squares, and relaxed sink
    spectra (repo and web, delegated search).
- **Related but distinct.**
  - Session 210's relaxations (convex-linear and acyclic only) are credited.
  - Session 211 and Session 212, and mac-mini's THM-4900 and HYP-9490, have no square or hexagon relaxation. The
    "permutohedral face formula" there is the coherent gate fibration.
  - Higher Bruhat orders B(n, 2) (Manin-Schechtman, Ziegler) use the same square/hexagon split, on paths rather than
    orientations.
  - Joswig-Kaibel-Koerner, arXiv math/0012204, has "unique sink on every 2-face".
  - Plott path independence.
  - None of these states (i)-(iii) or the bound.
- **Name clash.** "Context tournament" already means a tournament on the context or exterior vertices in five
  repository notes (MINOR-3).

## 5. Findings

### MAJOR

**MAJOR-1. THM-4890 (D) is read as "decorrelation", but the true micro values are correlated with the strips (factor
about 3.7).**
- **What (D) proves.** It is correct for i.i.d. uniform micro values F(l), and only for them.
- **What the session concludes from it.**
  - THM-4890 title and (D): "So the macro lattice ... creates no excess collisions", and "decorrelation criterion
    (DC)".
  - README section 0 (1): "The macro lattice ... is harmless. The whole difficulty is a digit-reversal-versus-
    multiplication decorrelation".
  - README section 4: "(DC) asks that ... occurs with frequency O(1/C), as it does for random values".
  - README section 9, item 1: the micro self-energies "do not correlate with the macro strips".
  - Both INDEX lines: "linear energy ... is one decorrelation criterion".
- **What exact computation shows** (`aud_randmodel.out`, `aud_dc_split.out`).
  - The true type-4 count DC is 3.31-3.69 times the model's expectation for n = 12..22: 2.222 against 0.602 |X_n| at
    n = 22.
  - With the true marginal law of g(l) - 1 it is 1.64-1.84 times. That ratio is still creeping up: 1.825, 1.832 and
    1.844 at n = 18, 20 and 22.
  - The converse involution alone forces 0.965 |X_n| type-4 collisions at n = 20. That is more than the model's
    whole type-4 expectation.
  - For n >= 14 the model's energy, 1.64-1.85 |X_n|, even lies below the deterministic floor
    2|X_n| - 2^(n/2) - 1 (THM-4886).
- **What follows.** Linear energy needs this ratio to stay bounded. It does not need decorrelation, and the data refute
  decorrelation.
- **Fixes.**
  - **THM-4890 title.** Replace "So linear energy, which stays OPEN, is equivalent up to constants to a decorrelation
    between micro self-energies and macro strips." with:
    "Linear energy (OPEN) is equivalent, up to constants, to criterion (DC): the collisions with h != h' and l != l'
    number O(|X_n|). In the random model their expectation is at most 76 |X_n|. The true count is 3.3-3.7 times the
    model value at n = 12..22, with about |X_n| of it forced by the converse involution. So (DC) is a bounded-
    correlation statement, not a decorrelation statement."
  - **(D).** Replace "So the macro lattice {g(h) + h c} creates no excess collisions." with:
    "So, in the random model, the macro lattice {g(h) + h c} creates no excess collisions. For the true micro values
    the type-4 count is about 3.7 times the model value (audit, n = 20, 22)."
  - **(D).** Change "decorrelation criterion (DC)" to "criterion (DC)".
  - **README section 0 (1).** Replace the two "What this means" sentences with:
    "In the random model the macro lattice creates no excess collisions. The true micro values give about 3.7 times
    the random-model count of collisions with h != h' and l != l', about |X_n| of them forced by the converse
    involution. Linear energy needs this ratio to stay bounded (criterion (DC), section 4)."
  - **README section 4 ("The criterion") and section 9, item 1.** Replace "do not correlate with" and "decorrelation"
    by:
    "correlate with the macro strips only boundedly: the type-4 count is at most K times its random-model expectation,
    for an absolute K. The observed ratio is 3.3-3.7 against uniform micro values and 1.6-1.8 against their empirical
    law, still slowly rising at n = 22."
  - **Both INDEX lines.** Change "decorrelation criterion" to "bounded-correlation criterion (DC); the true count is
    about 3.7 times the random model".

**MAJOR-2. THM-4891 claims the gaps "need zone-wide coherence"; only "neither face type alone" is proved.**
- **Where it appears.** THM-4891's title, README section 0 (3), both INDEX lines and the commit message.
- **What is proved.** Neither square coherence alone nor hexagon coherence alone excludes 7 or 21. Both together are
  zone-wide coherence, by (i).
- **What is not proved.** That every condition excluding 7 and 21 implies zone-wide coherence. Intermediate
  relaxations are untested: all hexagons plus some squares, or coherence on contexts of bounded size. README section 9,
  item 4 lists them as open.
- **Fixes.**
  - **Title.** Replace "The exclusion of 7 and 21 needs both face types, i.e. zone-wide coherence" with:
    "Neither face type alone excludes 7 or 21. Both together are zone-wide coherence (tournaments), which does. Whether
    an intermediate relaxation still excludes them is OPEN."
  - **(v), third bullet.** Replace "So the exclusion of 7 and 21 (THM-1370) is a property of the global zone structure.
    A coherence condition on 2-faces of only one type cannot carry it." with:
    "So no coherence condition on 2-faces of a single type carries the exclusion of 7 and 21 (THM-1370). It is not
    shown that full zone-wide coherence is necessary."
  - **Remarks.** Replace "so the gaps are not local" with "so neither face type alone carries the gaps". Mark the
    THM-4890 comparison as "(analogy)".
  - **README section 0 (3) and both INDEX lines.** Replace "the gaps need zone-wide coherence" and "a global
    (zone-wide) phenomenon" with "neither face type alone keeps 7 and 21 out".

### MINOR

1. **The XY law is proved, and uncredited (HYP-9482 Model; README sections 0 (2) and 5).** g(v)/2^n -> XY on reduced
   words, with density 4 ln(4t) and 4 ln(1/t), is proved in Session 16
   (`feedback_word_hamilton_spectrum_20261020/BINARY_PRODUCT_DENSITY.md` section 3, (7)-(8), uniform CDF error
   O(2^(-n/2))). Add after "The limit law of g(v)/2^n on Y_n is XY":
   "(PROVED: Session 16, BINARY_PRODUCT_DENSITY.md section 3). Only the Poisson clustering of converse classes at the
   unit scale is heuristic."
   Add Session 16's product law to README section 1.

2. **The n!/2^floor(n/2) bound is Session 188's (THM-4891 (ii)).** Add to the THM-4891 status, INHERITED:
   "the bound n!/2^floor(n/2) and its orbit proof are Session 188's tournament bound
   (session188_structural_bridges_20261104/equal_nine_ambient/PROOF.md (1)). The orbit argument applies verbatim to
   context tournaments, since an orbit fixes the pair at positions (2i-1, 2i) and the set before it. New here: equality
   at n = 4 and 6 (6 and 90 sinks), where tournaments reach only 5 and 45."

3. **Name clash.** "Context tournament" already means a tournament on the context or exterior vertices:
   - spectrum_basins_20261020e/five_atoms.md:115;
   - session33 BOUNDARY_OBSERVABILITY.md:143;
   - session73 REPORT.md:165;
   - square_sum_spectrum_overnight_20261031f/TOURNAMENT-PRICES.md:165;
   - HYP-2731:130.

   Rename to "set-context tournament", or add "(not the earlier 'context tournament' on exterior vertices)".

4. **Pre-registration disclosure (HYP-9482 outcomes; README section 0 table).** E_n/|X_n| for n = 22..28 was published
   in HYP-9479 (S83) before 21b02648ec: n = 26 and 28 in 8acbc16ebc (12:31), the full sequence 3.2607 ... 3.2473 in
   a2de7874e0 (13:41). Since E_n = sum_j E^Y_j + cross + 1, the E/|X| column is a
   reproduction, and P1 and P8 were jointly constrained by known data. Add:
   "The E/|X| column reproduces HYP-9479 (published before 21b02648ec) and is not out of sample; P1 and P8 are not
   independent of it."

5. **Ratios computed from rounded levels.** Exact counts give:
   - n = 24, a = 8 and 9: 0.464 and 0.375 (not 0.465 and 0.376);
   - n = 26, a = 8: 0.498 (not 0.497); a = 10: 0.373 (not 0.370);
   - n = 20, a = 7: 1052/2822 = 0.373 (not 0.372).

   Correct the HYP-9482 status ("ratio 0.373-0.380"), the P6 outcome list and MISTAKE-728 ("0.380, 0.375, 0.373"). In
   the THM-4890 remark, L_2 at n = 26 is 0.259 (0.259498), not 0.260. The verdicts are unchanged. MISTAKE-728 may add
   that the in-sample boundary ratio was below 0.45 at every even n: 0.370, 0.383 and 0.373 at n = 16, 18 and 20.

6. **README section 7.** "Three of the six pairs flip" should read "Four of the six pairs (01, 02, 12, 23) flip". All
   48 seven-sink hexagon-coherent orientations at n = 4 have exactly four flipping pairs (`aud_perm_flips.out`).

7. **P5's "proved form" disagrees with THM-4890 (C).** P5 quotes it as "a <= 4 + dXi and s <= 4 + dX", and
   `s84_energy.py` tests that weaker form, while THM-4890 (C) has 3 + dXi and 3 + dX. Replace the P5 outcome with:
   "CONFIRMED, and PROVED for all n by the audit's sharp lemma (U1): |a - s| <= 1 + |dX - dXi|."

8. **(D) novelty typing.** The type-4 bound is (1/C) times THM-4887's strip total (steps 3-5) at W = C. Type it as
   "PROVED (new statement; method and constants from THM-4887 steps 3-5)". Add that the exact random-model energy is
   1.64-1.99 |X_n| for n = 12..22, so 109 is not sharp.

9. **Witness reproducibility.** `witnesses.out` prints full tables for only 2 of the 5 witnesses, and the session's
   "independent check" reuses the session's class convention. Print all five tables. Record that the audit verified
   all five independently (`aud_perm_C.out`).

10. **Unarchived data claim.** README section 4, "sup_D N_t(c, c', D)(c + c' + 1)/2^t <= 4.7 (t = 10, 12)", has no
    committed script or output. Archive one, or mark it "(exploratory, not archived)".

11. **"Second moments are Poisson-pair" is too strong (README section 0 (2); HYP-9482 Reading).** The model is
    heuristic. The non-converse energy is 2.2 % above it at n = 26 and 1.7 % at n = 28. The upper tail is heavier:
    - at multiplicity 12: 1.26 and 1.16 times the model (n = 26, 28);
    - at 14: 1.24 and 1.35 times;
    - at 16: 1.63 and 1.36 times (3 and 10 values).

    Use: "Second moments match the Poisson-pair model (heuristic) to about 2 %. The upper tail is 15-35 % heavier at
    multiplicities 12-14."

12. **HYP-9482 P1 wording.**
    - "|E^Y_(n+2) - E^Y_n|/|Y|" is ill-posed, since the two |Y| differ by a factor of 4. Say it was evaluated as
      |E^Y_(n+2)/|Y_(n+2)| - E^Y_n/|Y_n||.
    - "has decreased since: 2.9328, 2.9349, 2.9303, 2.9296" is not monotone. Say "decreasing with an odd/even
      oscillation".

13. **Stale lines after U1.**
    - THM-4890 status: "this is not proved".
    - THM-4890 (C) corollary: "within 1 ... unless a, s <= 3".
    - README section 3 ("It is not proved") and section 9, item 2.
    - Results INDEX: "sharper |A - S| <= 1 unproved".

    Update all of them to U1. With no carries, the window midpoint lies within 1/2 of the word's midpoint, with no
    exception.

### UPGRADES

**U1. Sharp balance lemma (PROVED here).**
- **Statement.** Let v != v' in Y_n with g(v) = g(v'). Then:
  - X Xi < 0 and |A - S| <= 1;
  - if w >= 2, the word with the larger window value has window 1...1 and the other has 0...0;
  - if w = 1, ceil(h/2) = ceil(tau/2);
  - hence |a - s| <= 1 + |dX - dXi|.

  Also a sign rule. Take any v = h|x|sigma and v' = h|x'|sigma in Y_n with x = 1...1, x' = 0...0 and w >= 2,
  colliding or not. Then A >= S + 2 implies g(v) < g(v'), and S >= A + 2 implies g(v) > g(v').
- **(L0) Exact local formula.** For every w-bit string z (leading zeros allowed), with b_1...b_w its bits MSB first:

      g(z) = z C_w(z)/2^w + 1 + L(z),  L = #10(z) - R(z),  R = sum over i < j with b_i = 0, b_j = 1 of 2^(i-j-1).

  - Read z with Session 193's machine. Each 0 at position i adds A_i + b_(i-1). The sum of A_i over the zeros is the
    j < i part of z C_w(z)/2^w = sum over b_j = 1, b_i = 0 of 2^(i-j-1). The sum of b_(i-1) over the zeros is #10(z).
  - Group the zeros into maximal runs. A run contributes [preceded by a 1] - R_run, with 0 <= R_run < 1.
  - Hence L(z) <= #10(z), and L(z) > -1.
- **Step 1 (signs).** Swap so that X > 0, so x_(w-1) = 1 and x'_(w-1) = 0. For w = 1, Xi = -X automatically. For
  w >= 2, suppose Xi > 0.
  - Then x is even and at least 2^(w-1), and x' is odd and below 2^(w-1).
  - So g(x) >= 1 + 2^(w-2) and g(x') = g((x' - 1)/2) <= 2^(w-2).
  - (B) then reads h Xi + tau X = g(x') - g(x) + [tau odd] + [h odd] <= 1, while the left side is at least h + tau
    >= 2. Contradiction.
- **Step 2 (w = 1).**
  - (B) with g(1) = g(0) = 1 gives tau + [tau odd] = h + [h odd].
  - So h and tau lie in {2k-1, 2k}, and their bit lengths a and s differ by at most 1.
  - Since dX = dXi = 0, A = a and S = s.
- **Step 3 (w >= 2: the windows are 1...1 against 0...0).**
  - X >= 3 and |Xi| >= 3, so |Xi| >= 2^(w-1-dXi) + 1 and X <= 2^k - 1, where k = w - dX.
  - y = C_w(x) < 2^(w-1) and x' < 2^(w-1).
  - (L0) in (B) gives an identity valid for all such configurations:

        g(v) - g(v') = T X - H |Xi| - delta,
        T = tau + y/2^w,  H = h + x'/2^w,  delta = L(x') - L(x) + [h odd] - [tau odd].

  - If A >= S + 2, then H|Xi| >= 2^(A+w-2) + 2^(a-1) >= 2^(S+w) + 2^(a-1), and T X < (2^s - 1/2)(2^k - 1).
    So H|Xi| - T X > 2^(a-1) + 2^s + 2^(k-1) - 1/2 >= 2^(k-1) + 5/2.
  - Also -delta < #10(x) + 2 <= 3 + (k-1)/2. If dX >= 1, x = 1 0^(dX-1) z with z a k-bit word ending in 1.
  - Hence g(v) - g(v') < 1/2 + (k-1)/2 - 2^(k-1) < 0.
- **Step 4 (S >= A + 2).** Apply Step 3 to (v*, v'*), which (C') makes a collision pair with A and S exchanged. QED.
- **Checks** (`aud_sharp.out`):
  - (L0) exhaustively for w <= 16, and its bounds for w <= 14;
  - the Step 1 bounds for w <= 22;
  - the identity and the sign rule on 300,000 random configurations (125,989 with A >= S + 2 and 126,121 with
    S >= A + 2), biased toward carries;
  - every step on all 1,007,162 collision pairs for n <= 20;
  - the statement on 258.6 M pairs for n <= 28 (`aud_pairs_*.out`).

**U2. Extended checks.**
- (C), (C') and (E) now hold exhaustively for n = 4..28 (previously 10..20).
- The HYP-9482 statistics are reproduced independently for n = 21..28, including the P9-P11 addendum.

**U3. Exact expectation of the (D) random model.**
- For n = 12..22 the expected energy is E_F = 1.64-1.99 |X_n|, and the type-4 part is 0.44-0.60 |X_n|.
- With the empirical law of g(l) - 1 the type-4 part is 0.88-1.20 |X_n|.
- These numbers make (DC) quantitative (MAJOR-1).

**U4. THM-4891 verified by a fully independent construction.**
- The geometric face relations (i) and the dictionaries hold for n <= 7.
- The n = 5 spectra were rederived with two other solvers and encodings.
- The f(6) = 90 witness can be taken acyclic.

## 6. Scripts and outputs

All files are in this folder.
- **Build** with `gcc -O2 -o aud_pairs aud_pairs.c`, likewise for `aud_randmodel.c` (add `-lm`), and
  `gcc -O2 -o aud_perm_enum.exe aud_perm_enum.c`, which `aud_perm.py` part B needs.
- **Requirements:** Python 3 with numpy, mpmath and PySAT.

| script | output | content |
|---|---|---|
| `aud_blocks.py` | `aud_blocks.out` | (A), (B), the range of g, C_w, v* |
| `aud_pairs.c` | `aud_pairs_4_24.out`, `aud_pairs_25_28.out` | (C), sharp form, P5, (C'), (E), all HYP-9482 statistics, n = 4..28 |
| `aud_sharp.py` | `aud_sharp.out` | proof of U1 checked step by step |
| `aud_randmodel.c` | `aud_randmodel.out` | (D): exact random-model expectation (uniform and empirical law), Monte Carlo, DC, E_n split |
| `aud_dc_split.py` | `aud_dc_split.out` | converse-forced share of DC |
| `aud_D_constants.py` | `aud_D_constants.out` | (D) arithmetic (76, 109) and the per-pair strip bound |
| `aud_model.py` | `aud_model.out` | HYP-9482 constants, multiplicity law, tails, model Monte Carlo |
| `aud_perm.py`, `aud_perm_enum.c` | `aud_perm_AB.out`, `aud_perm_C.out`, `aud_perm_D.out` | THM-4891: classes n <= 7, spectra n = 3, 4, SAT n = 5, witnesses, f(6) >= 90 |
| `aud_perm_flips.py` | `aud_perm_flips.out` | flipping pairs of the 48 seven-sink orientations (MINOR-6) |
