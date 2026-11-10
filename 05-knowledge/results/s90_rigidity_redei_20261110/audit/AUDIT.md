# Independent audit of the S90 packet (THM-5124, THM-5125, HYP-9657, README)

Auditor: independent adversarial sub-agent, 2026-11-10. Own code throughout (nothing imported from the packet).
State audited: THM-5124 and THM-5125 as committed in `b52e250ee9` (unchanged during the audit); HYP-9657 and the
README as in the working tree at 12:27 local (sha1 `419b69f0c1`, `ec60832ec9`), i.e. including the items added
after the pre-audit checkpoint ((C) 3, (B'), (D); README sections 3.3 scope check, 4 last bullet, 5.5, 6).

## Verdict: PASS WITH CORRECTIONS

The mathematics of THM-5124 (statements 1-4, both proofs of 3, Remark 4) and of THM-5125 (statements 1-4) is
correct; every proof step was checked and every number in the two theorems was reproduced by independent code, with
exhaustive verification one order beyond the packet for statements 2 and 3. The corrections concern (i) novelty and
attribution (THM-5125 (1) and (4a) are known in the repository and in the literature; THM-5122 is misdescribed and
the arbitrary-dimension cubes are THM-5117's), (ii) HYP-9657 (A), whose universal sentences are false on the
packet's own windows, and (iii) scope and wording.

## MAJOR findings

**M1. THM-5125 statement 1 (a)-(b) is prior art, in the repository and in the literature; status and credit are
wrong.** Files: THM-5125 title, status ("PROVED (opus-2026-11-10-S90): statement 1 ... the finite realisation
inside the class 3 mod 8 is KNOWN ... the countable statement is the tournament analogue of ... the Rado graph"),
statements 1 (a), (b), proof 1; README section 2.
- `05-knowledge/results/analogy_reduction_atlas_20261025/QR_TOURNAMENT.md` (klein, 2026-10-25) defines the same
  object with the same orientation ("the Legendre tournament `L` on `P_3`", `p -> q` iff `(p/q) = +1`).
- `05-knowledge/results/root27_machine_bridges_20261026/anchor/LEGENDRE_SIDECAR.md` (2026-10-26) gives the full
  extension property by CRT + Dirichlet + reciprocity, the back-and-forth, and the conclusion "the prime Legendre
  tournament is the countable random tournament", typed "CITED known result", citing Cameron, *The Random Graph
  Revisited* (2001), p. 3, and G. A. Jones, *Infinite Paley graphs*, arXiv:1912.02709 (sources inspected there; not
  re-inspected by this audit).
- HYP-9375 (PROVED 2026-10-29, `tournament_open_problem_routes_20261029/PROOFS.md` P1) proves quasirandomness of the
  prime prefixes, which contains finite universality with frequencies.
- `pythagorean_hourglass_census_tower_20261025.md` item 5 already phrases the Selmer dimension as
  `dim ker(A^2 + A + I)` "of the Legendre tournament".
Correction: retype 1 (a), (b) KNOWN / CITED with these four sources (plus the bsd note for the class 3 mod 8); use or
cross-reference the existing name "Legendre tournament"; keep as this session's content only the reformulation (c)
("residue chains"), the Heegner data and the least-realiser table.

**M2. THM-5124 misdescribes THM-5122 and asserts an existence statement that is false at order 6.** Files: THM-5124
status ("THM-5122 (Session246) constructs additive cubes of arbitrary dimension up to the negative inertia in the
EXCHANGE (ear) cube of one fixed core"), Remark 2 ("additive cubes exist up to dimension `n_-(W)`"); README
section 1 ("THM-5122 shows that in the exchange cube of one core (reversing which side of a cut a vertex is on)
additive cubes exist"), section 6.
- THM-5122 / Session246 `switches/PROOF.md` proves the affine criterion, the UPPER bound
  `k <= min(n_-(W), floor(n/2))` and exhibits ONE two-dimensional cube (345, 347, 347, 349); it says explicitly that
  the inertia bound "neither produces such directions".
- Affine cubes of every dimension `k` are THM-5117's (Session245): `k` disjoint pair exchanges in the ordinal
  composition of one h13 and `k - 1` h37 modules, order `7k + 1`, spacing `50 * 37^(k-2)`. Recounted here for `k = 2`
  (increments -74, 124, mixed difference 0; `audit_thm5117_cube.out`). So `k` is about `n/7`, far from the bound.
- "Exist up to `n_-(W)`" is false: exhaustively, no core of order 6 carries three pairwise orthogonal disjoint
  exchanges with all three increments nonzero (720 affine (core, matching) pairs, each with a zero increment), while
  `floor(n/2) = 3` and `n_-(W) = 5` for the THM-5122 witness core. At order 7 they exist, e.g. core
  `[68, 49, 18, 103, 105, 69, 6]`, cut `0010111`, exchanges `0->3, 1->6, 2->5`, counts
  `873, 1267, 801, 1195, 1001, 1395, 929, 1323` (`audit_exchange.out`, `audit_misc.out`).
- "EXCHANGE (ear) cube" and "reversing which side of a cut a vertex is on" conflate two cubes. On the ear cube
  proper (single flips) every 2-face has strictly negative interaction `-W_ij <= -2 H(C - i - j) < 0`
  (THM-4114 (33), Session124, Session192 `cut_sumsets/PROOF.md`: "cannot produce an additive cube of dimension two");
  additivity needs exchange directions.
Correction: "THM-5122: criterion, capacity bound and one 2-cube; THM-5117: affine exchange cubes of every dimension
(order `7k + 1`); the capacity bound is not attained at order 6 (audit)"; state the trichotomy presence (`>= 0`, zero
only in 1, 3, 3, 5) / ear (`< 0`) / exchange (signed).

**M3. HYP-9657 (A): the universal sentences are false on the packet's own windows.** Files: HYP-9657 title, (A),
"Reading"; README section 3 item 3 (iii), section 4 third bullet. Independent recomputation reproduces every number
printed by `f1_mass.py` (both windows), and then shows:
- "For every odd prime `p` ... largest for `c = 0` and `c = -1`": false at the first untested prime. On both
  `[2^18, 2^21)` and the pre-registered `[2^21, 2^23)` the classes 0 and -1 are NOT the two largest for
  `p = 29, 31, 37, 53, 61, 73, 89, 127` (e.g. `p = 29`: 10.995, 10.959 against other classes up to 11.092; out of
  sample 12.185, 12.206 against 12.231); at `p = 73, 127` they are below average. Among the primes up to 61 and
  73, 89, 127, only `p <= 23` and `41, 43, 47` show the pattern in both windows (`audit_f1_classes_*.out`).
- "the remaining `p - 2` classes are equal up to sampling error": false at `p = 7` on the packet's own output
  (10.48-10.91 in sample, 11.47-12.19 out of sample; 27 and 80 standard errors of a class mean, with the
  independence yardstick stated in `audit_f1_classes.py`), and significantly heterogeneous at
  `p = 11, 13, 17, 19, 31`. The mod-7 bias of the one-apex image is known (THM-4878 (5), HYP-9480 (A), mass note
  section 7).
- "no Kronecker-symbol dependence": at `p = 7` the split by `(-N/7)` is -0.25 (17 s.e.) in sample and -0.44 (53 s.e.)
  out of sample (`N = 1, 2, 4` heavier than `N = 3, 5`); since the squares mod 7 are the powers of 2 this cannot be
  told apart from a base-2 effect. At the other primes the split is not robust (`p = 11`: -4.6 s.e. in sample,
  +1.1 out of sample).
- "about `1/p` more solutions": the qualitative reading is confirmed (the class-0 excess comes entirely from the
  sheet `s >= 1` and the class-(-1) excess entirely from `s = 0`, `audit_f1_split_23_21.out`), but the size is not
  `1/p`: relative excess `0.56, 0.35, 0.23, 0.18, 0.036, 0.26, 0.045, 0.043, -0.002, 0.005, 0.008, 0.030, 0.098` at
  `p = 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43` against the generic `0.375, 0.208, 0.146, 0.092, 0.077,
  0.059, 0.053, 0.044, 0.034, 0.032, 0.027, 0.024, 0.023`.
The pre-registered tests Q1-Q3 are unaffected (they name `p = 3, 5, 7, 11, 17`) and are CONFIRMED independently.
Correction: restrict (A) to the primes tested, record the absent primes, delete "equal up to sampling error" and "no
Kronecker dependence" (or restrict them to `p = 5`), retype the `1/p` reading REFUTED as a size law and keep the
sheet decomposition; adjust the title and README 3.3 / 4.

**M4. HYP-9657 (A)/(B) and the principle "mass, not congruence": the mechanism is largely S82's and is not
credited.** `one_apex_join_representation_mass_pi2_over_12_20261106.md` section 5 already records "a bounded local
factor at every small odd prime, the `N = 0` class being over-represented", `r1` by `N mod 3 = (18.1, 12.5, 14.1)`
(classes 0 and -1 above class 1 for F1), mean 0.796 on primes and 0.813 on `N` coprime to `3...17` (F2), and the
rotating mod-7 bias of `g`; the S78 note records the residue profile of the F1 misses. New here: the class `-1` at
`p >= 5`, its attribution to the sheet `s = 0`, the Rabinowitsch stratum. Correction: cite S82 section 5 and the S78
profile in (A)/(B) and in README 3.3 and 6; replace the citation "(THM-4875 (D), THM-4888)" for the mod-3 bias (THM-4888
is about O/Z banks).

**M5. THM-5125 statement 4 (a) and the term "cycle gas" are known in the repository (attribution only; the
mathematics is correct).** `per(I + A)`, `det(I + A)` as sums over cycle collections: THM-467 (2026-06-10,
refutation table and closing remark), THM-1950 / Session137 `GLOBAL_PROOF.md`; `H = I(Omega, -1) mod 3`: THM-4754,
THM-4756; `I(Omega, -1)` as an Euler characteristic: THM-467; "cycle gas": Session37 (2026-10-23), Session159. README
3.5 cites "THM-5125 (4a)" for the mod-3 fact. Statements 4 (b), (c) (the `F_2`/`F_4` balance criteria) were not found
earlier in the repository; their first sentence (coefficients of `det(xI + A)` as sums over cycle covers of
`k`-subsets) is the classical coefficient theorem for digraphs (Harary, Sachs) reduced mod 2. Correction: type 4 (a)
KNOWN with these sources, cite the coefficient theorem in 4 (b), keep the balance criteria as the session's content.

## MINOR findings

1. THM-5124 statement 4: `H(S(A)) = H(A) + D(A)` fails for empty `A` (gives 2, true value 1); the cube formula
   `S -> c (H(A[S]) + D(A[S]))` gives `2c` at `S = empty`. Write `D(A) + [A nonempty] H(A)`.
2. THM-5124 statement 4 is the full-cut case of the inherited closure identity: `S(A)` is the protected closure of
   `A` with every vertex selected, and `F(1) = h + 2 a.1 - sum_(i<j) W_ij = H(A) + D(A)` (verified on every `A` of
   order `<= 5`); `D` is Session197's "labelled two-bin covers". Credit.
3. THM-5124 proof of statement 2, last step: `f(W1 - z, empty) = 3c` needs the line "`a` in `W1 - z` and
   `f({a}, empty) = 3c`".
4. THM-5124 status: "explains the smallest gap mechanism" - it shows that three parallel switches give 9c or 15c,
   not why 7 is not a value.
5. THM-5124 Remark 4: "(THM-5122 (6))", "(3)-(4)" are equation numbers of Session246 `switches/PROOF.md`. The
   remark itself is correctly derived and scoped (constants 1.1696, 0.8270, exponent 0.6460 checked). Addendum: a
   unit-Lipschitz cube has at most `k + 1` values.
6. THM-5125 statement 3: "Its sub-tournament counts" refers to the seven-prime tournament, not to the six primes.
7. THM-5125 statement 3 / proof 3: `a(35) = 1901942427` IS globally least (the only candidate outside the box,
   `3*7*11*19*23*31*607 = 1899519699`, has 151 chains; all 668 candidates with at least 7 prime factors were
   examined). The certification argument silently uses that 8 primes give more than `3 * 10^9`.
8. THM-5125 Remark 1, README 3.2: "the action of an element of order 3 on the 2-torsion of the curve" is not proved
   or cited; `E_n` has no automorphism of order 3. What holds at matrix level: for all `p_i = 3 mod 8` the kernel of
   the Monsky matrix is stable under `(x, y) -> (y, x + y)`, which acts as `A + I` (`audit_monsky_tau.out`; fails for
   general `n`). State that, or type ANALOGY.
9. THM-5125 Remark 2, README 3.3 (ii): "(q + 1)/4 prime" needs `q >= 7` (`q = 3` gives 1). HYP-9657 has it right.
10. THM-5125 statement 2 is the classical inertness of the primes below `(q + 1)/4`; type it as a restatement.
    (Audit observation: the domination property fails for all 156 non-Heegner primes `q = 3 mod 4`, `23 <= q < 2100`.)
11. `redei_model.py` docstring: "`[p_j/p_i] = 1`" means Legendre symbol -1 in the cited note; the code uses +1
    (harmless by the orientation independence, verified on every tournament of orders 3-6).
12. HYP-9657 setting: "THM-4885 gives the mean of `r`" - no: `r` counts triples `(h1, h2, s)`, THM-4885 counts
    pairs `u <= v`. Means 12.19 against 13.01 on `[2^21, 2^23)` (THM-4885 predicts 13.01); `r(11) = 2, r(19) = 4`
    are 1 and 3 as pairs. Zero sets agree. With the pair count the top-two pattern already fails at `p = 13` on the
    pre-registered window (class -1: 13.048, another class: 13.051).
13. HYP-9657 (B): "the Heegner numbers are rough": 11 and 19 have `3 | 12`, `5 | 20`; only 43, 67, 163 (and 7) are.
14. HYP-9657 (C) 3 / README 4: the equality cases of `2 mu >= max(d)` through order 7 are 14 profiles, not two
    (e.g. `(5, 4, 4, 2)`, `(11, 8, 8, 4)`); the order-7 line of `interaction_gap_refuted.py` is a printed constant
    (confirmed here, every pair); docstring typo `(15, 14, 10, 6)`; the order-10 sampled minimum `14/34` is beaten by
    `18/50` (`audit_gap_witness.out`); the minima found are `(2^j + 1)/(3 * 2^j + 1)`, suggesting the limit 1/3.
15. HYP-9657 Non-consequences: "that list was not re-analysed here" contradicts the added (B').
16. README 1: the "no free-staircase family" sentence is in section 1 (route 2) of the mass note, not section 2.
17. README 3.1: the corrected "odd" map is already THM-467 ("even lengths cancel, odd double") and LEM-004;
    "`H = U_T(1, 0, 0, ...)`" must mean the variables `x_i` (every `p_k -> 1`), not `p_1 = 1, p_3 = 0`.
18. README 3.3: "This is why 7, 43, 67, 163 appear among the misses ... and of the O/Z/J readout" - a statistical
    tendency, exploratory for O/Z/J.
19. README 3.4: 679 and 2719 are unconditional exceptions (Jones-Pall; Gupta 1941); Ono-Soundararajan (GRH) is
    completeness only. "Misses 29 odd numbers": FINITE-EXACT through 50,001 (THM-4884).
20. README 5.1 "Dead": HYP-9480 (C) is OPEN; the window part of THM-4876 (5) is FINITE-EXACT through order 8.
21. README 5.3: Bourgain-Kontorovich's exceptional set is not power-saving (that is Magee-Oh-Winter with the
    Bourgain-Kontorovich-Magee appendix); "best: THM-4887" is improved asymptotically by Session216
    `density/RESULT.md`, `A(X) >= X^(2/3) / (60 (ln X)^(2/3))`.
22. README 6: "`log H` is ... multiplicative on order joins" - additive; the bold slogan of section 1 omits
    "beyond 1, 3, 3, 5".

## Verified independently (method, universe)

- THM-5124 (`audit_rigidity.py`; path tables by depth-first enumeration of all directed paths, cross-checked
  against permutations on all 1,099 labelled tournaments of orders 1-5; structure (ii) through strong components):
  pair lemma on every labelled tournament and pair, orders 3-7 (24, 240, 2160, 23520); statement 2 on every
  ordered partition, orders 3-7 (48, 1440, 38880, 1095360 instances; values, unique active vertices, hubs in `B`,
  singleton components; sub-selection additivity to order 6); statement 3 on every module triple, orders 3-7
  (3.57 * 10^9 at order 7; none); three pairwise additive vertices, orders 5-7 (9c transitive, 15c cyclic, never
  7c); statement 4 on every `A` of order 0-5; random tournaments of orders 8, 9.
- Proofs: every step of THM-5124 (pair lemma both directions, Steps 1-2, all arcs of `Y => Q => X`; statement 2
  including beta = c, uniqueness and the 5-cycle for every orientation; both proofs of statement 3, Session214
  eq. (8) does cover induced extensions of any size; statement 4; Remark 4 with Busch's floor) and of THM-5125
  (extension property with `r = 3 mod 8`, norm form and signs, identities 4 (a)-(c) rederived, orientation
  independence) was checked by hand.
- THM-5125 (`audit_redei.py`, Jacobi symbols by reciprocity; `audit_cyclegas.py`, packings as permutations, Ryser,
  exact determinants, characteristic polynomials mod 2 coefficientwise, elimination over `F_4`): arcs, scores, 8
  triangles, `H = 71`, prefixes, 29, value set, Monsky nullity 4 in three conventions; domination; every entry
  `a(N)`, `N <= 45`, with no box; identities and both bridges on all labelled tournaments of orders 3-6 (first
  failures: order 5, 120 labelled tournaments, all with `H = 13`, the first being `[18, 8, 3, 5, 14]`, the masks of
  Session245's h13 core; pure adjacency: order 6, `H = 17, 29`).
- HYP-9657 (`audit_f1.py`, from the canon form `F1(u, v)`, checked against 256 explicit joins and 64 apex
  triples): the psi-coordinates, both windows, Q1-Q3, the Heegner values, the eight zeros among 30; (C) 3 and (B')
  (`audit_new_items.py`). The pre-registration commit `cf65422aaa` (11:39:04) is on origin/main and precedes the
  output file (11:40).
- Nine of the ten packet scripts were re-run (all but `redei_model.py`; `f1_mass.py` on the in-sample window);
  every `.out` reproduces (`rerun_*.out`; `rigidity_check.out` up to its timing field).

## Upgrades the packet may take from this audit

1. Statements 2 and 3 of THM-5124 are FINITE-EXACT through order 7 (packet: order 6).
2. `a(35) = 1901942427` is globally least.
3. No core of order 6 carries an active additive exchange 3-cube; order 7 does (explicit example in M2).
4. The class-0 excess of F1 lives on the sheet `s >= 1` and the class-(-1) excess on `s = 0`; sizes by prime in M3.
5. Monsky kernel, primes `3 mod 8`: stable under `(x, y) -> (y, x + y)`, acting as `A + I` (two-line proof from
   eq. (5) of the cited note; checked on 300 random `n`).
6. Second-order gap: `2 mu >= max(d)` holds for every pair at orders 5-7; ratios `3/7, 5/13, 9/25` at orders 8-10.

## Files

`audit_rigidity.py`, `audit_exchange.py`, `audit_thm5117_cube.py`, `audit_redei.py`, `audit_cyclegas.py`,
`audit_monsky_tau.py`, `audit_f1.py`, `audit_f1_classes.py`, `audit_f1_split.py`, `audit_new_items.py`,
`audit_gap_witness.py`, `audit_misc.py`, each with its `.out`; `rerun_*.out` are the packet's scripts re-run.

## Post-audit note (12:33 local)

While this report was being finalised the author began applying it. The working-tree THM-5124 at 12:33 (status
field, statement 4, last step of statement 2, Remarks 1, 2, 4) was re-read: the applied corrections for M2 and for
minor items 1-5 are accurate. One residual precision: the strict inequality `W_ij >= 2 H(C - i - j) > 0` on the ear
cube is THM-5122 (5) / Session124 (THM-4114 (33) gives only `<= 0`). THM-5125, HYP-9657 and the README had not been
corrected at that time; M1, M3, M4, M5 and minor items 6-22 refer to them as described under "State audited".

## Not checked

Monsky's formula itself and the Selmer interpretation; the cited pages of Cameron and Jones; `redei_model.py` was
not re-run (its output is reproduced by other code); the order-8 census for the second-order gap; the heuristic
program (D) beyond its elementary bound `k! (k+1)^n` and the differencing identity; the 2^35 coverage it cites.
