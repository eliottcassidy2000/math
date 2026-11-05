# Audit of THM-4880 (ribbon tournaments) and its source note — 2026-11-04

**Auditor:** independent adversarial audit (Claude Code subagent, 2026-11-04 evening), run inside the worktree
`math-wt-recursive-20261104`; no repository file was modified (all reruns were made on a copy of the scripts
directory in the auditor's scratchpad, listed in "Independent computations").
**Verdict (short):** PASS WITH CORRECTIONS — one major proof gap (H2, repaired by the audit with an exact
computation; the conclusion stands), the rest minor (wording, a mislabelled table row, script labels, receipts,
attribution). No false theorem statement, no even `H`, no arithmetic slip in the quoted data.

## Scope

1. `01-canon/theorems/THM-4880-ribbon-tournaments-bounded-span-back-arcs-give-a-finite-state-transfer-for-hamiltonian-paths-strong-connectivity-is-regular-the-span-2-series-is-an-irreducible-14-dimensional-integer-representation-47-needs-span-3.md`
   — statements (A)–(H), the status block and the FINITE-EXACT paragraph.
2. `05-knowledge/results/ribbon_tournaments_bounded_span_transfer_20261104.md` — sections 3–8 (proofs and data),
   the width table, the level table, the junction ranks. The placeholders LEVELS-FINAL, THREE-LETTER-RESULT,
   SPAN3-RESULT, DFS-BOUND, SPAN3-P2, THREE-LETTER-P3, AUDIT-SECTION were ignored as instructed.
3. Scripts and receipts in `04-computation/experiments/ribbon_spectrum_20261104/` (ribbon.py, condensation.py,
   intrep.py, exact_algebra.py, algebra.py, minrec.py, recurrences.py, witness47.py, width.py, probe47.py,
   digitscan2.py, junctions.py, ribbon_enum.c, ribbon_levels.c, and the data files span2_matrices.txt,
   span2_intrep.txt, span2_condensation.txt, test_s2_misses.txt, lv2_5e4_levels.txt, lv2_1e5_levels.txt,
   intrep_w3.log, lll_span3.log).

## Method

* Every proof in sections 3–5 of the note was re-derived by hand (Theorem A profile argument, Theorem B merge
  lemma line by line, Corollary C, Lemma D's algebra, Proposition E's Busch bound, Proposition F's
  commuting-operator argument, Theorem G's conclusion from dimension 169, Theorem H's two inferences).
* Every script named in the task was rerun from an isolated copy; all outputs were compared with the numbers
  quoted in the theorem and the note.
* Independent code (written from the definitions, not from ribbon.py) was used for: a Hamiltonian-path counter
  (memoised subset recursion, cross-checked by a permutation count for `n <= 7`), a strong-connectivity test,
  a condensation routine, a re-implementation of the merge-lemma automaton from the note's text, the integer
  representation read directly from `span2_intrep.txt`, exact (python-flint / sympy) ranks and characteristic
  polynomials, the 47 census, the end-twin test, the run sequences, the full digit-block distribution, the
  junction ranks, and an independent DFS census to 50,000.
* Time used: about 14 minutes of CPU in total (the 50,000 DFS took 643 s; width.py, recurrences.py and
  digitscan2.py a few minutes together). Nothing was skipped.

## Findings

Severity: **major** = a false or unproved claim typed as proved; **minor** = wording, attribution, receipt.

1. **MAJOR (proof gap; conclusion true; repaired by the audit).** Theorem (H): "Every block of length `<= 4` has
   `rank(M_B - I) >= 25` on the 26-dimensional reachable space. Hence no bounded-carry digit system and no
   arithmetic-progression gadget exist among short blocks"; note §8 (H2): "For every block of length `<= 4`,
   `rank(M_g - I)` on the reachable space is `>= 25` of 26 (recurrences.py); no block inserts an arithmetic
   progression of values." The rank datum is correct (reproduced: minimum 25, attained by `--`, `--.--`, …),
   but it does not prove the heading "No identity-plus-nilpotent block": a single nilpotent Jordan block of size
   26 has `rank(M - I) = 25`, so rank `>= 25` is exactly compatible with `M = I + N`, `N` nilpotent. What does
   prove it (audit, exact, python-flint): on the 14-dimensional quotient, `charpoly(A_B) != (x - 1)^14` for all
   340 blocks; more precisely the eigenvalue 1 has algebraic multiplicity 0 for 289 blocks and 1 for the
   remaining 51, never `>= 2`. Since the Hankel rank is 14 the reachable space maps onto `Q^14` under `Psi`, so
   `M_B - I` nilpotent on the reachable space would force `charpoly(A_B) = (x - 1)^14`; hence no block of length
   `<= 4` is identity-plus-nilpotent. Multiplicity `<= 1` also gives the arithmetic-progression statement
   rigorously: a non-constant progression `H(u B^k) = H(u) + k c` needs a Jordan block of size `>= 2` at the
   eigenvalue 1. **Correction:** replace the rank sentence as the proof by "the eigenvalue 1 of every block
   matrix `A_B` (length `<= 4`) is at most simple (exact charpoly), so no block is unipotent on the reachable
   space and no block inserts a non-constant arithmetic progression"; keep rank `>= 25` only as a datum.
2. **MINOR (wording that changes the content).** Theorem (H): "the residual `H(uB) - alpha H(u)` over 8,040 prefix
   states is bounded only when it vanishes"; note §8: "takes a bounded set of values only when it is identically
   0". Over a finite set of prefixes every residual set is bounded. The computed fact is: exactly 22 of the
   4,664 (block, alpha) pairs have a single residual value (then 0), namely the ordinal-sum blocks with
   `H(B) = alpha` (2 of length 3 with `alpha = 3`; of length 4, 8 with `alpha = 3` and 12 with `alpha = 5`; audit
   list below), and every other pair takes at least 227 distinct residual values (over the 7,873 distinct
   reduced prefix states, i.e. the 8,040 distinct profile vectors). Say "constant only when it vanishes; every
   other pair has `>= 227` distinct residuals". Also the status line "340/342 blocks": digitscan2.py scans 583
   blocks (363 three-letter blocks of length `<= 5` plus 220 blocks of length `<= 4` containing an empty letter)
   times 8 multipliers; 340 is the number of junctions. "342" corresponds to nothing in the scripts.
3. **MINOR (file name / index wording, false as written).** The theorem's file name says "the-span-2-series-is-an-
   irreducible-14-dimensional-integer-representation", and the INDEX entries say "the span-2 representation is
   irreducible (full matrix algebra M_13(Q))". The 14-dimensional representation is reducible: rows 1 and 2 of
   every letter matrix in `span2_intrep.txt` coincide, so the hyperplane `y_1 = y_2` (the span of the
   nonempty-word states) is invariant. What is irreducible is the 13-dimensional subrepresentation on that
   hyperplane, which the YAML title's (G) states correctly. Rename the file (e.g. "…a-rank-14-integer-
   representation-irreducible-on-nonempty-words…") and fix both INDEX lines before the commit.
4. **MINOR (attribution: KNOWN as a method).** Theorem A's profile transfer is the standard frontier-based /
   bounded-bandwidth transfer-matrix dynamic programme for Hamiltonian paths (the "connectivity state" DP used for
   Hamiltonian paths and cycles in grid graphs and in graphs of bounded pathwidth; Knuth's simpath, TAOCP 4A
   §7.1.4; frontier-based search). The novelty is the tournament specialisation (tails outside the window are
   free, at most one unfed head, 7/27/119 states) and everything built on it (ranks, integer quotient, G, H, the
   census). Likewise Theorem B's regularity is known in principle: strong connectivity is MSO-definable and the
   tournament is interpreted from the word, so the language is regular (Büchi–Elgot–Trakhtenbrot); the explicit
   merge-lemma automaton (4/17/121 states) is the new content. Add a "KNOWN (method)" line to the status for A and
   B; nothing else is silently re-derived (one-arc join, near-join, tower, Busch floor, insertion lemma,
   multiplicativity are all cited).
5. **MINOR (Corollary C, construction stated for span 2 only).** Note §5: "the ordinal sum `u => v` of two ribbon
   words is the ribbon word `u . -- . (second letter of v intersected with {1}) . v_3 ...` (the first letter of
   `v` becomes empty and its second letter loses the span-2 bit…)". For general `w` (the corollary is stated for
   span `w`) the `k`-th letter of `v` must be intersected with `{1, …, k - 1}` for every `k <= w` (for `w = 3`
   the third letter loses its span-3 bit, which would point into `u`). The claim is true with that correction.
6. **MINOR (Theorem B wording).** Note §4: "`C_{j_1 + 1}, ..., C_m` (all inside `A`, hence singletons or a pair of
   window vertices)". A strong component is never a pair. For `w = 2` these components are singletons; in general
   they are components contained in `A` (`|A| <= w`; for `w = 3` a 3-cycle on the window is possible). The
   merge lemma itself is correct (checked line by line; initial-stretch monotonicity, order preservation, the
   `j_1 < i_0` case where `v` becomes its own component before `C_{i_0}`, and the dead rule (iii) all hold).
7. **MINOR (table row wrong at two entries).** Note §6, row "strong, width 2 | 1 (C3 has width 1) | 1 | 5 | 28 |
   123 | 418": the entries at orders 3 and 4 must be 0 — width.py itself prints `strong: {1: 1}` at orders 3 and
   4 (`C3 = --.1-.1-` and the strong 4-tournament `--.1-.1-.1-` have width 1; the parenthetical already says so).
   All other rows and the strong-by-width lists (1/28/6; 1/123/226/3; 1/418/4483/1106) reproduce exactly.
8. **MINOR (script label off by one; the note is right).** `intrep.py` prints "words of length <= 2: module rank
   14" and "<= 3: 14", but from `k = 2` on the printed length is one less than the true suffix length (the frontier
   is multiplied twice between prints). Exact ranks (audit): suffix functionals of length `<= 0, 1, 2, 3, 4` have
   Q-rank 1, 5, 11, 14, 14. So the note's "generated already by the suffixes of length `<= 3` (length `<= 2` gives
   rank 11)" and the theorem's "generated by suffixes of length `<= 3`" are correct, and the 14 listed semantic
   suffixes do have rank 14. Consequence for `intrep_w3.log`: "length <= 3: rank 62" means length `<= 4`
   (and "<= 2: 51" means `<= 3`). Fix the print label (`k + 1`) or the log's reading.
9. **MINOR (receipts of the level enumeration).** (i) `ribbon_levels.c` writes a header with 11 column names but 10
   numbers per row: "strong_values_n" is never printed (`nstrong_vals` is unused). (ii) The strong bitmap is
   pre-seeded with 1, 3, 5, so at `n <= 3` the strong column credits values not yet attained at that length
   (5 at `n = 3`); harmless from `n = 4` on, where `cum_all = cum_strong` at every level is what justifies "strong
   values coincide with all values at every level". (iii) At the start of the audit `lv2_5e4_levels.txt` was being
   regenerated (rows `<= 14` only, no misses line); by the end of the audit the run had finished and the file
   contains rows 15–17 (5,046,447 / 13,566,206 / 27,797,197 states; solid [49, 45457], [49, 49999]; cumulative
   misses 12, 4, 4) and the lines "odd misses <= N (all words): 7 21 47 50001" and the same for strong words —
   all table values for `n = 8..16` match. Commit the completed file. (iv) `lv2_1e5_levels.txt` has no misses
   lines, but its row `n = 16` (cumulative misses 4 for all and for strong words below 100,000) already proves
   coverage to 100,000: the four misses can only be 7, 21, 47 and the boundary artefact 100001. A free
   strengthening of the FINITE-EXACT bound if the misses lines are regenerated.
10. **MINOR (precision in (E)).** Title: "every realiser of 47 has width 3" is verified for orders 7–8 only, but it
   holds for every order and the argument should be written: the strong realisers of 47 are exactly the two
   order-7 classes (Busch: order `<= 8`; width.py shows the four order-8 realisers are non-strong, i.e. ordinal
   sums of an order-7 realiser with one vertex), and ordinal sums with transitive pieces keep width exactly 3
   (`>= 3` by the induced-subtournament bound of (C), `<= 3` by the closure). Also say explicitly in §6 that the
   four order-8 classes with `H = 47` are not strong.
11. **MINOR (receipt).** Note §5: "the profile states differ, and later letters separate them: twins_state.py".
   The script prints nothing for this pair by design (it only reports state-equal or fully H-equivalent pairs),
   so the receipt is the absence of a line. The count equality itself is Lemma D (proved; audit algebra agrees)
   plus the random tests.
12. **MINOR (interpretation of a receipt).** Note §7: "first-attainment lengths 1:1 3:1 4:1 … 14:787 15:370 16:32"
   (from `test_s2_misses.txt`). These are the lengths at which the depth-first search first met each value, not
   the minimal realising lengths (the DFS takes `--` to full depth before backtracking). The minimal-length
   counts are the level differences of the level table (e.g. 109 - 42 = 67 new values at `n = 8`, not 25).
   Rephrase as "DFS first encounters".
13. **Remark (no fault).** Among the 96 rank-1 junctions, 94 contain an uncrossed cut of the vertex order; the
   other two, `--.12.-2.--` and `--.12.-2.1-`, are ordinal sums along a non-interval cut (`u + y | x z t v`,
   where `y` is the junction's second vertex) — the MISTAKE-717 phenomenon again, worth a footnote in §8.

Checked and found correct (no correction needed): Theorem A's proof (heads of non-initial fragments lie in the
window, free tails, at most one unfed fragment, the four attachment rules, ageing, the kill rule; the sorted
tuple loses nothing because no two fragments share a (head, tail) pair); the merge lemma; Corollary C's
multiplicative-closure and prime statements; Lemma D's algebra (the two "c a"/"c b" sums share
`sum_{S'} h(X[S']) s_a(X - b - S')` and the difference is `e_b(X) - (h(X - b) - e_a(X - b)) = 0`); Busch's floor
values `f(3q) = 3·5^(q-1)`, `f(3q+1) = 5^q`, `f(3q+2) = 9·5^(q-1)`, `f(9) = 75`, hence strong realisers of 47 have
order `<= 8`; the Krylov + Berlekamp–Massey argument (shift operators commute; `q(E)H` satisfies the proved
order-16/21 recurrence and vanishes on 112 consecutive terms, hence identically); the recurrence of the
`{2}`-run on all printed terms (and at `n = 8` with `H(0) := 1`); the growth constants 2.288291, 2.562972, and the
span-3 growths 2.68 and 3.32 (dominant roots 2.68076 and 3.32260 of the Krylov polynomials; the term ratios
oscillate around them); the conclusion of Theorem G from dimension 169 (`= 13^2`, so the generated algebra is
all of `M_13(Q)`, hence no common invariant subspace and no smaller quotient) and the 13-vs-14 explanation
(`y_1(u) = H(u --) = H(u) = y_2(u)` for nonempty `u`, `y(empty) = e_1`); the exact characteristic polynomials of
the four letters on the 13-dimensional space; the counts 996, 500, 86,778 (= 254 + 21,844 + 64,680), 300,
8,040, 340, 19,764,437, 6,133,659; the solid-prefix factors (1.56–2.87, quoted as 1.6–2.9); both run sequences
absent from the OEIS (API queries returned no match; a control query for the tribonacci numbers returned
A000213).

## Independent computations

Scratchpad directory of the audit (scripts in `audit_scripts/`, reruns in `ribbon_copy/`, outputs in `out/`):
`/private/tmp/claude-501/-Users-e-Documents-GitHub-math/99e0ab32-35be-49e2-9a09-5ca4d960241b/scratchpad/`.

| item | what was run | result |
|---|---|---|
| (a) | `python3 ribbon.py` (rerun) | 7 / 27 / 119 states, 4 / 8 / 14 one-fragment states; 332 random words per `w` verified; ALL CHECKS PASSED |
| (a) | `check_A_transfer.py`: own subset-recursion counter (own tournament builder) vs `ribbon.count_with_matrices` | 200 random words (100 of span 2, 100 of span 3, `n <= 9`): all equal, all odd; own permutation count agreed on the 152 words with `n <= 7` |
| (b) | `python3 intrep.py 2` (rerun on a copy) | rank 14, largest entry 5, 500 random words verified; the produced `span2_intrep.txt` is byte-identical to the repository file |
| (b) | `check_B1_rep.py`: `span2_intrep.txt` read independently, 14-dim integer recursion vs own counter | 100 random words (`n <= 10`) agree; entries in `[-4, 5]`; rows 1 and 2 of every letter matrix identical; `y_1 = y_2` on all 100 |
| (b) | `check_B2_ranks.py`: exact ranks (flint) of suffix functionals from `span2_matrices.txt` | ranks 1, 5, 11, 14, 14 for `|v| <= 0..4`; reachable rank 26; Hankel block rank 14; the 14 listed semantic suffixes have rank 14 |
| (b) | `python3 minrep.py 1 2` (rerun) | Hankel rank 4 (`w = 1`), 14 (`w = 2`); forward 26 / backward 14 |
| (c) | `python3 condensation.py` (rerun) | 4 / 17 / 121 states (2, 2, 3 accepting); 0 mismatches on 254 / 21,844 / 64,680 words |
| (c) | `check_C_strong.py`: own strong test and own condensation; own re-implementation of the merge-lemma automaton from the note's text | counterexample `--.1-.-2.1-`: back arcs 2->1, 3->1, 4->3, not strong, condensation {2} < {1,3,4}; 0 mismatches vs condensation.py on all 87,380 span-2 words of length `<= 8` and 6,000 random span-3 words (`n <= 10`); my automaton has 4 / 17 / 121 states and 0 mismatches vs my brute force; 4,408 of the 21,840 span-2 words of length 2–7 have a non-interval strong component |
| (e) | `check_D_twins.py`: own counter, 300 random `X` of orders 3–8 (+ `c`) | 0 failures; control without the sole-exit hypothesis: 273 of 273 differ |
| (f) | `check_E_47.py`: Busch values; own exhaustive census of normalised span-2 words of length `<= 8` with own counter and own strong test; the two span-3 witnesses | `f(3..13) = 3, 5, 9, 15, 25, 45, 75, 125, 225, 375, 625`; largest order with `f <= 47` is 8; no strong span-2 word of order `<= 8` has `H = 47`; odd numbers `<= 75` not attained strongly at orders `<= 8`: 7, 21, 47; both witnesses: `H = 47`, strong |
| (f) | `python3 probe47.py`, `python3 witness47.py`, `python3 width.py` (reruns) | no strong word with `H` in {7, 21, 47} for `n <= 9` (11: 9 words, 13: 9, 67: 7); three-letter solid intervals [49,165], [49,257], [129,639], [249,1833], [397,3033] with overlaps; witnesses `H = 47` by transfer and brute force; width table orders 3–8 reproduced exactly, order 7 = 1/20/194/238/3, strong 1/123/226/3, `H = 47`: two strong width-3 classes at order 7, four non-strong width-3 classes at order 8 |
| (g) | `python3 minrec.py`, `python3 recurrences.py` (reruns); `check_F_runs.py` (own counts) | minimal recurrences of orders 8 (`n >= 9`) and 12 (`n >= 13`) with the stated polynomials, checked to `n = 119`; Krylov orders 7 / 16 / 21 (span 2) and 8, 17, 55, 22, 72, 89, 86 (span 3); own counts of `{2}^n`, `{1,2}^n`, `{1}^n`, `{3}^n`, `{1,2,3}^n` for `n <= 10` equal the quoted terms; recurrence verified by hand on the printed terms; dominant roots 2.28829113889 and 2.56297156736, both polynomials irreducible over Q; span-3 dominant roots 2.68076 and 3.32260 |
| (h) | `python3 exact_algebra.py`, `python3 algebra.py` (reruns); `check_G_H.py` (exact charpolys on `y_1 = y_2`) | nonempty-word span 13; algebra dimension 169 of 169 with four letters and with the three nonempty letters; exact charpolys `x^12 (x-1)`, `x^10 (x^3-x^2-x-1)`, `x^5 (deg 8)`, `x (deg 12)` as quoted |
| (i) | `python3 digitscan2.py` (rerun); `check_G_H.py` (full distribution and exact unipotency test); eigenvalue-1 multiplicities | 8,040 prefix vectors; the 15 smallest residual sets are all `{0}`; full distribution: 22 pairs with one residual (all ordinal-sum blocks: `--.--.-2`, `--.1-.1-` (alpha 3); `--.--.--.-2`, `--.--.1-.1-`, `--.--.-2.--`, `--.--.12.12`, `--.1-.--.12`, `--.1-.1-.--`, `--.1-.-2.1-`, `--.1-.12.-2` (alpha 3); `--.--.-2.1-`, `--.--.-2.-2`, `--.--.-2.12`, `--.--.12.1-`, `--.1-.--.-2`, `--.1-.1-.1-`, `--.1-.1-.-2`, `--.1-.1-.12`, `--.1-.-2.-2`, `--.1-.-2.12`, `--.1-.12.1-`, `--.1-.12.12` (alpha 5)), next smallest size 227; 583 blocks x 8 multipliers; `rank(M_g - I) >= 25` reproduced; no block of length `<= 4` has `charpoly(A_B) = (x-1)^14`; eigenvalue 1 has multiplicity 0 (289 blocks) or 1 (51 blocks) |
| (j) | `./ribbon_enum span2_matrices.txt span2_condensation.txt 2 5000 18 audit_s2` (rerun; output `ribbon_copy/audit_s2_misses.txt`) | nodes 19,764,437, strong leaves 6,133,659, attained 2497 of 2501, misses 7, 21, 47, 5001 — identical to `test_s2_misses.txt` except the timing |
| (j) | **own census**: `./ribbon_enum … 2 50000 21 out/s2_5e4` (maxdepth 21 = largest order with Busch `f <= 50,000`) | 1,637,211,036 nodes, 513,524,921 strong leaves, attained 24,997 of 25,001, misses exactly 7, 21, 47, 50001; 643 s; an undeduplicated 27-coordinate confirmation of the level enumeration's 50,000 claim |
| (j) | level files | `lv2_5e4_levels.txt` rows `n = 8..16` (states, values, solid prefix, cumulative misses) all equal the note's table; completed file ends with misses 7 21 47 50001 for all and for strong words; `cum_all = cum_strong` at every level `n >= 4` |
| (k) | `python3 junctions.py` (rerun); `check_K_junctions.py` (exact flint ranks) | 1: 96, 2: 27, 3: 153, 4: 34, 8: 9, 9: 6, 10: 11, 12: 4 (sum 340), both floating and exact |
| (l) | OEIS API (curl, JSON) for `1,1,3,5,11,29,65,149,343,777,1781,4081` and `1,1,1,5,15,41,105,261,659,1689,4347,11165` (and the versions without the leading 1s) | no match; control query for the tribonacci numbers returned A000213 |

## Verdict

**PASS WITH CORRECTIONS.** Statements (A)–(G) and the FINITE-EXACT census are correct and reproduced (the
50,000 bound additionally confirmed by an independent undeduplicated DFS). Statement (H) is correct but its second
half (no identity-plus-nilpotent block / no arithmetic progression) was not proved by the evidence given; the
audit's exact eigenvalue computation (multiplicity of 1 at most simple for all 340 blocks) proves it, and the
theorem should cite that instead of `rank(M_B - I) >= 25` (finding 1). Apply findings 2–12 (wording, the file
name's "irreducible 14-dimensional", KNOWN-method attribution for the profile DP and for regularity, the general-
`w` ordinal-sum rule, the width-table entries at orders 3–4, the intrep.py labels, the level-file header, and the
precise form of the digit-block statement) before the final commit.
