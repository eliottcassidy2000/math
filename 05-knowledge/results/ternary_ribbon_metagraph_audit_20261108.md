# Audit of THM-4952, THM-4953, HYP-9541 and the note ternary_ribbon_metagraph_lines_20261108.md

Independent, adversarial audit (2026-11-08, same machine, separate scratch directory
`/private/tmp/claude-501/-Users-e-Documents-GitHub-math/99e0ab32-35be-49e2-9a09-5ca4d960241b/scratchpad/audit4952/`).
No repository file other than this report was created or modified. The literal placeholders PENDING-N9,
PENDING-S2-1E5, PENDING-S3 and PENDING-AUDIT were ignored, as instructed; the three censuses behind them were left
running untouched.

## 1. Scope

Audited: `01-canon/theorems/THM-4952-ribbon-odd-cycle-packing-polynomial-...md` (statements (A)-(F), status);
`01-canon/theorems/THM-4953-every-metagraph-line-...md` ((A)-(D), census paragraph, status);
`05-knowledge/hypotheses/HYP-9541-...md` (typing against the theorems); sections 2-6 of
`05-knowledge/results/ternary_ribbon_metagraph_lines_20261108.md`; the scripts and receipts in
`04-computation/experiments/ternary_ribbon_metagraph_20261108/`. Background canon (THM-4880, THM-4940, THM-4744,
THM-643/644/4729, THM-4930, THM-4865 (c), THM-4922, HYP-9453) was read for definitions and cited, not re-audited.

## 2. Method

1. Reran the author's receipts: `ribbon_packing.py selftest`, `rankdrop.py`, `spec3_small.py 7`,
   `readiness.py 7 99` and `readiness.py 6 45`, `lines_formula.py 7`, `./lines_census 7 H`,
   `./ribbon_levels_mod16 ... 2 3000 20 aud_pk3 3000000 0 3` (outputs written to the scratch directory only), and
   diffed the committed `formula_n{5..8}.txt` against the `#cuts` sections of `lines_n{5..8}_p0.txt`.
2. Wrote independent code for every proved statement that is finitely checkable: a cycle-enumeration brute force
   for I(Omega(T), x) (explicit simple directed cycles found by DFS from their minimum vertex, then an
   independent-set DP over vertex subsets) and a backtracking Hamiltonian-path counter (`aud_packing.py`); a
   flint-based rank/determinant study of the span-2 functionals (`aud_rank2.py`, `aud_rank3.py`, two inline scripts);
   a complete enumeration of the closed tilings at n = 5, 6, 7 with own strong-component, cut-set, grid-symmetry
   and anti-automorphism code plus an own implementation of the Moebius formula and a pynauty node census at
   n = 5, 6 (`aud_lines.py`); an own packing census of every tournament of order <= 8 in bases 2..7 with own
   strong-connectivity test (`aud_spec_small.py`); an own near-join check and own root-response sets at orders <= 6
   (`aud_nearjoin.py`); a copy of `ribbon_levels_mod.c` instrumented to compare consecutive pruned state SETS, not
   cardinalities (`fixpoint_check.c`).
3. Read every proof argument in the note and the theorem files line by line; checked OEIS identifications via the
   JSON API; checked every census number quoted in THM-4953 against `lines_analysis.txt`, `readiness_n8.txt`,
   `pk3_s2_1e4_levels.txt`, `spec3_orders_le_8.txt` and the `pk*_intrep.txt` headers.

Computation time used: about 40 minutes of wall clock on a machine already saturated by the author's eight
`lines_census 9` processes.

## 3. Findings

Severity: MAJOR = a false or unproved claim typed as proved; MINOR = wording, attribution, receipt.

**1. MAJOR (proof gap in THM-4952 (B), repaired by the audit; the statement is true).**
Sentence at fault (THM-4952 status): "the generic rank 14 certified at eighteen specialisations (a nonzero 15 x 15
minor would have at most 15 roots)"; note section 2: "the rank is 14 at eighteen distinct specialisations (a nonzero
15 x 15 minor has degree <= 15, so it would have at most 15 roots)".
Why it is wrong: the rows of the functional matrix are phi.M_v for |v| <= 4; only the letter matrices are linear in
x, the rows are products and their entries have degree up to 2 (depth 4) and 3 (depth 5), so a 15 x 15 minor can
have degree 30 (45). Eighteen specialisations therefore do not show that every 15 x 15 minor vanishes
identically, i.e. the upper bound "generic rank <= 14" was unproved. Moreover the committed `rankdrop.py` tests
only thirteen generic points (7/3 and -3, -1, 1/2, 1, 3/2, 5/2, 3, 4, 5, 6, 7, 10); the values 11, 12, 13, 17, 19,
23 listed in (B) are not in the receipt.
Audit repair: rank 14 at the 65 integer points 3..67 for the depth-<= 4 matrix (341 rows, entry degree <= 2, minor
degree <= 30) and at the 61 points 3..63 for the depth-<= 5 matrix (1,365 rows, entry degree <= 3, minor degree
<= 45). Hence all 15 x 15 minors of both matrices vanish identically, the depth-4 span equals the depth-5 span over
Q(x), so it is letter-invariant and is the whole functional space: the minimal dimension over Q(x) is exactly 14.
(Depths <= 2 and <= 3 give only 10 and 13; depth 4 is genuinely needed.) An own nonsingular 14 x 14 minor (rows
0, 2, 3, 4, 10, 11, 14, 15, 16, 19, 47, 48, 63, 255 of the depth-4 matrix, columns 0-7, 9-14) has determinant
-x^13 (x - 2)(x + 2), computed by interpolation at 31 points. The exact ranks of the full letter-invariant
Z-modules, computed with the author's `packing_intrep.build` stabilisation at x = -2, 2, 0, -1, 1, are 13, 13, 1,
14, 14. Corrected sentence: "the generic rank is 14: the depth-4 and depth-5 functional matrices have rank 14 at
65 resp. 61 integer specialisations, more than the degree bounds 30 resp. 45 of their 15 x 15 minors, so the
depth-4 span is letter-invariant; the full module has rank 13 at x = +-2 and 1 at x = 0 (exact Z-module
stabilisation), and a nonsingular 14 x 14 minor has determinant -x^13 (x - 2)(x + 2)". Receipt: add the
65/61-point check to `rankdrop.py` or cite this audit.

**2. MINOR (receipt does not establish the stated premise; conclusion independently true).**
Sentence (note section 4, mirrored in THM-4952 (D)): "The pruned state sets are a fixed point from length 19
(1,197,559 states, unchanged through length 24); since the level n + 1 is the image of level n ... the state set of
every length is this fixed point". `ribbon_levels_mod.c` records only cardinalities; equal counts at consecutive
levels do not imply equal sets, which is what the fixed-point argument needs. The audit instrumented a copy of the
program to test set equality: S_19 != S_18 (counts differ), S_20 = S_19, S_21 = S_20, S_22 = S_21 as sets (zero
missing states), so the premise holds and the "for all lengths" conclusion stands. Note also that the conclusion is
secured without any fixed point: 10 = 1 + 3*3 and 31 = 1 + 3*10 are excluded for every tournament by THM-4940 (3)
and (4), and the strong miss 64 by THM-4952 (E). Correction: state that set equality was verified (cite this audit's
`fixpoint.log`), or add the set test to the program; mention the THM-4940 fallback.

**3. MINOR (over-claim relative to the evidence; true, but the justification must be supplied or the clause
restricted).** Sentence (THM-4952 title (E) and note section 4): "the first composite attained value without a
strong realiser in any base" / "base 3 is the first base in which a composite attained value has no strong
realiser". The committed evidence covers bases 2-5 at orders <= 8 only. The audit's own census of all 12,357
tournaments of orders <= 8 gives, for bases 6 and 7, the attained values without a strong realiser at those orders
as {1, 175} and {1, 120}; the squares 49 = H_6(C3 => C3) and 64 = H_7(C3 => C3) ARE strongly attained; 175 and 120
exceed f(9) = 75 so they are undecided but irrelevant to "first". For b >= 8 every value below 64 is 1 + ab with
a <= 7, attained iff a lies in the intercyclic spectrum A (THM-4940 (2)); A ∩ [1, 7] = {1, 2, 4, 5, 6, 7} and each
of these a is the odd-cycle count of a strong intercyclic tournament of order <= 5 (C3, the strong 4-tournament,
and the strong 5-tournaments with H = 9, 11, 13, 15), so each such value has a strong realiser. Hence the claim is
correct, but the file must either add this argument (bases 6-7 census + block-lemma argument for b >= 8) or say
"in bases 2 to 5".

**4. MINOR (transcription error in THM-4953's census paragraph).** Sentence: "First creators (THM-4922; nodes 6,
33, 328, 5,457): lines FC-FC 24, 371, 13,723, 932,047". At n = 5 the FC-FC count is 18 (6 blue + 12 black,
`lines_analysis.txt`; the note's table correctly says "18 (6)"); 18 + 14 + 0 = 32 lines, whereas 24 + 14 + 0 = 38.
Correction: "FC-FC 18, 371, 13,723, 932,047".

**5. MINOR (attribution).** THM-4953 (B): "L_k = 0 for k >= 2 on both ends (the closed end is strong), so lines
between two non-strong nodes do not exist", and the k = 1 column 18, 391, 14,389, 983,928. Both are already canon:
THM-4729 (3) (non-strong nodes are pairwise non-adjacent and loopless; every line at a non-strong node ends at a
strong node) and THM-4729 (4) (strong-strong lines = 2^(m-1) - NS(n), NS = 14, 121, 1995, 64648), whose numbers
the k = 1 column reproduces exactly. The new content of (B) is the full cut-set distribution N(S) and its blue
refinement. Correction: cite THM-4729 (3)-(4) as KNOWN/CITED for these two facts in the status and in (B).

**6. MINOR (inconsistency between hypothesis file and note).** HYP-9541 (2): "grows without bound ... but
subexponentially"; note section 7 (2): "grows without bound but at most linearly in n". Choose one formulation
(the data 2, 14, 16, 26 support neither specifically); the hypothesis file is the one of record.

**7. MINOR (wording, THM-4952 (B)).** "The 13 at x = 2 is the irreducible dimension of THM-4880 (G) (the path
series' rank 14 counted the empty word)". The identification is justified, but by two facts that should be
stated: OCF makes the packing and path series agree on every nonempty word, so the Hankel rank on nonempty
prefixes is THM-4880 (G)'s 13; the computed rank 13 of the full packing module at x = 2 adds that the empty-prefix
row (packing value 1, path value 0) is dependent on the others.

**8. MINOR (wording, THM-4952 (E)).** "In base 2 the strong spectrum omits only 1 (THM-4880, orders <= 8 here)".
THM-4880 (C) concerns the ribbon family (values = multiplicative closure of strong values); for all tournaments the
statement is finite data at orders <= 8 (audit census: the only base-2 value at orders <= 8 without a strong
realiser is 1). Cite it as FINITE data, not THM-4880.

**9. MINOR (receipt caution for the pending items, not a claim audited here).** `pk3_s3_5e3.log` records "state cap
reached at n=11" and `pk3_s3_5e3_levels.txt` shows exactly 5,000,000 states at length 11: that run was truncated
by its state cap, so PENDING-S3 must not be filled from this file. The replacement process now running uses the
default 20,000,000-state cap with 70 int32 coordinates, which exceeds the memory of this 8 GB machine.
`pk3_s2_1e5` stands at length 16 with 12,872,609 states (no fixed point yet).

No other faults were found. Typing is otherwise consistent: PROVED parts of THM-4952 are (A) (the frontier method
correctly credited as KNOWN, the tournament automaton new), (B) (after finding 1), (C), (D) (after finding 2), (E);
(F) is FINITE-EXACT; THM-4953 (A)-(D) are elementary and verified; HYP-9541 is OPEN with NUMERICAL sub-clauses,
and its PROVED references point to the right theorem parts. No statement overreaches from n <= 8 data beyond the
conjectures explicitly typed NUMERICAL in HYP-9541 (whose own caveat "the conjecture rests on HYP-9520's thinning,
not on monotone data" is appropriate).

## 4. Independent computations (all in the scratch directory; files named)

- **(A) automaton.** `aud_packing.py`: own cycle-enumeration brute force vs `poly_with_matrices` vs the author's
  set-partition brute force on 100 random span-2 and 100 random span-3 words, n <= 8: 200/200 agree; x = 2 equals
  the own backtracking Hamiltonian-path count in every case. Author's selftest rerun: 4 / 20 / 112 states, 360/360
  agree. Reading of `step()`: heads must lie in the window (age < w, killed at age w after the step that could
  still feed them with span w), tails become free at age w, closing requires fragment parity 0 (odd cycle), a
  singleton can never close (its tail condition and head condition contradict), merging uses head of the tail-
  fragment and tail of the head-fragment, each packing corresponds to exactly one accepted path. Correct and
  complete.
- **(B) dimension.** See finding 1: ranks at 65 + 61 points, own minor determinant -x^13 (x - 2)(x + 2), exact
  module ranks 13 / 13 / 1 / 14 / 14 at x = -2 / 2 / 0 / -1 / 1. `rankdrop.py` rerun reproduces the author's output
  (determinant, roots, ranks 1, 13, 13, one relation at x = 2). Base-3 letter matrices of `pk3_span2_intrep.txt`:
  characteristic polynomials x^13 (x - 1), x^11 (x - 2)(x^2 + x + 1), x^3 (x - 1)(x + 1) q_9, x (x - 1)(x + 1) q_11
  with q_9, q_11 exactly as in the note; spectral radii 1, 2, 2.446327, 2.718009; the only ordered two-letter blocks
  with characteristic polynomial x^13 (x - c) are --.--, --.1-, 1-.-- with c = 1. The author's rank-14 base-3
  representation agrees with the own brute force at x = 3 on 60 random words (n <= 9). Entry bounds read from the
  files: 2, 5, 6, 21 for x = 2, 3, 4, 5 (span 2), 4 for span 3 at x = 3, ranks 13, 14, 14, 14, 14, 14, 14 (x = 2, 3,
  4, 5, 1, 6, 7) and 70 / 61 (span 3, x = 3 / 2), as stated.
- **(C) tower.** `aud_packing.py`: for n <= 12 every directed cycle of 1-^n is an interval traversed
  i -> j -> j-1 -> ... -> i (checked cycle by cycle), the polynomial satisfies f(n) = f(n-1) + f(n-2) + (x-1) f(n-3),
  H_3 = (3*2^n + c_n)/7 with c_n = 4, 1, -5, and H_2 equals the own path count 1, 1, 3, 5, 9, 17, 31, 57, 105, 193,
  355, 653. OEIS (JSON API): A122552 = 1, 1, 1, 4, 7, 13, 28, 55, 109, 220, 439, 877, 1756, ... with
  a(n) = a(n-1) + a(n-2) + 2a(n-3); A000213 = tribonacci 1, 1, 1, 3, 5, 9, 17, 31, 57, 105, ... Both identifications
  correct.
- **(D) census.** `ribbon_levels_mod16` rerun with N = 3000: misses 10, 31 (all) and 10, 31, 64 (strong), states
  stabilising at 199,726 from length 17. `fixpoint_check.c` (instrumented copy) with N = 10^4 to length 22: counts
  identical to the author's table through length 19 and set equality S_20 = S_19 = S_21 = S_22 verified (finding
  2). Soundness of the pruning: appending letters realises T_u as an induced subtournament of T_uv, so THM-4940 (1)
  (monotone under vertex deletion) is exactly the right statement; deduplication on (reduced state, condensation
  state) loses nothing because the reduced state determines H_3 of every continuation and the condensation state
  determines strongness of every continuation; the letter order of `span2_condensation.txt` (`letters(w)` of
  ribbon.py) matches the intrep file.
- **(E) 64.** `aud_spec_small.py`: own packing polynomials of all 12,357 tournaments of orders <= 8 (gentourng),
  own strongness test. Base 3: 382 values, misses 10, 31 then only values >= 1069; attained but by no strong
  tournament: {1, 64}; the realisers of 64 at orders <= 8 all have alpha = (1, 6, 5) and are non-strong; 16 is
  realised by (1, 2, 1) non-strong, (1, 5) non-strong and (1, 5) strong. Base 2: {1}; base 4: misses 13, 57, 109
  then >= 1669, not-strong {1}; base 5: misses 16, 51, 136, 371 then >= 2286, not-strong {1}; bases 6, 7 as in
  finding 3. H <= H_3 termwise (alpha_k 2^k <= alpha_k 3^k) and f(9) = 3*5^2 = 75 (THM-4744's Busch floor) check.
  `spec3_small.py 7` rerun agrees with the committed order-7 numbers.
- **(F) readiness.** `aud_nearjoin.py`: on 20 random rooted tournaments (orders 3-7) the near-join
  t -> y -> z -> t, all other vertices beating y and z, has exactly h + 2 R_t Hamiltonian paths with R_t computed
  by the formula sum_X h(T[X]) s_t(T - X); own R-sets over all tournaments of orders <= 6 (h <= 45) agree row by row
  with `readiness.py 6 45`; `readiness.py 7 99` rerun is consistent with `readiness_n8.txt` (readiness only gained
  at order 8). `readiness_n8.txt` has 97 rows with h >= 3, of which 50 admit R = h + 1 and 55 admit R = h + 2; the
  26-element "neither" list of the theorem equals the file's.
- **THM-4953 (A), (B), (C), (D).** `aud_lines.py`, own enumeration of all closed tilings at n = 5, 6, 7 (32, 512,
  16,384 lines): zero violations of: components of both ends are intervals; k = |cuts| + 1; closed end strong; cut
  sets disjoint; k + k' <= n + 1 with equality exactly once (the transitive line); no end with n - 1 components;
  non-cut runs of length >= 2; every blue line has a grid-symmetric open end and both ends admit the
  anti-automorphism rho(i) = n + 1 - i. The own Moebius formula (all and blue) equals the own enumeration cut set by
  cut set at n = 5, 6, 7; k-tables 18/10/3/0/1, 391/101/15/4/0/1, 14389/1816/153/20/5/0/1 with blue
  6/0/1/0/1, 27/1/3/0/0/1, 239/0/15/0/1/0/1 as in the theorem; N({1}) = N({n-1}) = 5, 50, 903 against
  2^C(n-2,2) = 8, 64, 1024; non-strong shares 0.4375, 0.2363, 0.1218. Fibre 3 carries 3, 4, 5 = n - 2 lines, fibre 5
  carries 10, 15, 20 = 5(n - 3), fibre 9: 8, 37, 57, fibre 15: 1, 34, 94; the open end has H = 1 on exactly one line;
  the complement 1^m of the transitive tiling obeys i -> j iff j = i - 1 or j >= i + 2 and has H = 9, 17, 31 at
  n = 5, 6, 7 (own path counts; A000213). Author's `lines_formula.py 7` and `lines_census 7 H` reruns are identical
  to the committed files, and the committed formula files equal the committed census `#cuts` sections at
  n = 5, 6, 7, 8.
- **Census paragraph (h).** `aud_lines.py` pynauty census at n = 5, 6 over all 2^m tilings: nodes 12 (8 SC), 56
  (12 SC); fibre size = H/|Aut| at every node; blue endpoints odd and black even at SC nodes, no blue and black odd
  at NS nodes; S-classes SC-SC blue 8 / black 12, SC-NS 12 (n = 5) and 32 / 10, SC-NS 156, NS-NS 314 (n = 6);
  value colours (nodes 12/0 and 48/8; lines blue-blue 32 and 316, mixed 0 and 176, black-black 0 and 20); maximum
  multiplicity 2, 6; self-lines 4, 8; same-value lines 7, 53. The BLACK list in `lines_analysis.py` has 83 members,
  the union of THM-4930's three exact lists (29 + 33 + 21) below 1131, and no node at n <= 8 has H >= 1131. The
  parity law of (C) follows from THM-643's per-node parities by counting black endpoints at NS nodes
  (2 x NS-NS + SC-NS = sum of odd numbers). All n = 7, 8 numbers in THM-4953's census paragraph were checked against
  `lines_analysis.txt` and match except finding 4.

Skipped (time budget): an own n = 8 node census (1,048,576 lines with pynauty; the author's `lines_nodes_n8.json`
derived numbers were checked for internal consistency and against THM-4729's NS(8) = 64648 only); own readiness at
order 8; own recomputation of the span-3 ranks 70 and 61 (read from file headers and build logs); the 10^5 and
span-3 level censuses (pending); THM-4880's condensation automaton and THM-4940's monotonicity (cited canon).

## 5. Verdict

**PASS WITH CORRECTIONS.** One major correction (finding 1: the generic-rank certificate in THM-4952 (B) rests on a
false degree bound and on more test points than the receipt contains; the statement is true and the audit supplies
a valid certificate) and eight minor ones (findings 2-9: a fixed-point receipt that records counts, not sets; an
"in any base" clause whose justification must be added; one transcription error (24 -> 18); attribution of the
k = 1 / non-strong facts to THM-4729 (3)-(4); one hypothesis/note inconsistency; two wording items; one caution on
the truncated span-3 run). Every proved statement of THM-4952 (A)-(E) and THM-4953 (A)-(D) was reproduced with
independent code or argument; the FINITE-EXACT and NUMERICAL claims checked agree with the committed receipts.
