# Independent audit of THM-4873 / HYP-9466 (opus-2026-11-04-S77, C2 digit heredity)

**Scope.** Independent adversarial audit, 2026-11-04, worktree `C:/tmp/wt-s77`. I read in full: the note
`05-knowledge/results/ternary_parent_c2_digit_heredity_20261104.md`; THM-4873; HYP-9466 (at d96e81539f and at HEAD); the
three S77 programs and `c2_digit_heredity_20261104.out`; the S77 index lines; and the prior art in section 5. Every number
below comes from code written for this audit. The S77 programs were only compiled and run, in a scratchpad, to test their
reproduction commands.

## Verdict: PASS WITH CORRECTIONS

- Every computed claim reproduces exactly, and the pattern holds out of sample at |B| = 8 (550,400 grandchildren).
- The end-count lemma is correct; its proof needs wording fixes.
- Four MAJOR corrections are needed:
  - M1: provenance of the heredity no-go;
  - M2: a false bold sentence;
  - M3: the preparation route is mis-typed;
  - M4: the claim that a ternary parent "must come from fresh bases" is false.
- Section 4 upgrades the corrected (iii) to PROVED at every order.

## 1. Reproduced (independent code)

| claim | S77 | audit |
|---|---|---|
| classes, orders 1..8 | 7,412 | 1, 1, 2, 4, 12, 56, 456, 6,880 (191,536 at order 9); canonical form = min code over score-sorted relabelings |
| DP cross-check | - | (h, e_a, B_a, strong) multisets identical to klein's nauty-based `triples_4..8.txt` (12 / 42 / 255 / 2,739 / 51,327 rows) |
| (B, a) pairs and grandchildren, orders 1-7 | 3,613; 31,973 | 3,613 (1, 2, 6, 16, 60, 336, 3,192); 31,973 (3, ..., 28,728) |
| children with a digit-1 / digit-2 mark, \|B\| = 2..7 | 1..6 / 1, 4, 4, 4, 4, 4 | same. Digit 1 sits only at y (B transitive, a not the sink). Digit 2 sits at a for (TT_n, a penultimate), and at a and y for (TT_j => C3, a in C3) |
| grandchildren with digit 1 / digit 2 | 0 / 3 | 0 / 3, but each of the 3 has two digit-2 marks: y2 and the old mark u (M2) |
| smallest positive digit, \|B\| = 7 | >20 in 28,463 | histogram 3 / 16 / 246 / 28,463 identical |
| C2 law, sigma_x = 0, sigma_y = h - e_a | 0 failures | 0 failures at both levels for \|B\| <= 8; permutation brute force on all 285 labelled (B, a) with \|B\| <= 4 |
| out of sample, \|B\| = 8 | - | 55,040 pairs: 7 children with digit 1, 4 with digit 2. 550,400 grandchildren: 0 with digit 1, 3 with digit 2 |
| end-count lemma (ends and starts) | PROVED | 0 violations over all 2,131,019 labelled tournaments of order <= 7 and all classes of order <= 9 |
| sweep misses, odd v <= 1,500 | {1, 7, 11, 13, 21, 47, 67} | same; every odd v in [69, 3,141] is attained (= klein's `coverage_8`) |
| realised odd h; with digit 1; with digit 2 | 320 / 65 / 79 | 320 (largest 661) / 65 / 79 |
| explicit digit-1 family, "about 35%" | - | 5,410 / 15,000 odd values <= 30,000 (36.1%); complete dyadic windows 0.36-0.39; 116 of 200 odd values <= 400 missed. Closed form = DP (t <= 11) = permutations (t <= 5); sigma_a = 1 for all t <= 11 |
| pre-registration | - | d96e81539f committed 16:20:21; Session191 rebased onto it 16:21:06; `.out` committed 16:23:10; the program is unchanged. OK |

S77's own reproduction commands work and match: `./any 8 6 1500` (0.3 s) and `./fam 6 14 30000` (2 min 15 s).

## 2. MAJOR findings

**M1 (provenance): the heredity no-go for strong inputs is inherited prior art and is not cited.**
- `two_entry_funnels_20261103/READY.md` section 4 (2026-11-03, audited) proves the following for the reset, which is the
  converse of C2. The reset of any strong input of order >= 4 has no 1- or 2-ready vertex under arbitrary re-marking.
  The proof contains the native bound h - s >= 3, and section 4 records the C3 boundary: an H = 9 output with two
  2-ready roots.
- `multiport_defect_grammar_20261103/construction/RESULT.md` is the very file S77 cites for law (4). It says "Our C2
  operation is exactly its converse ... its failure to remain 1/2-ready are inherited". Its section 4 dualises the
  no-go: "C2 cannot produce either next nonzero digit at any auxiliary".
  - Eq. (12) gives S_y' = H(C2) + h - e_a (that is, sigma_y = h - e_a) and S_a' = 4S - H(B - a).
  - Eqs. (8) and (10) give sigma_x = 0.
- `nonzero_preparation_20261103/algebra/ONE_EXIT_AND_PERSISTENCE.md` section 1 says the same.
- **Fix** (note l. 17-18 and 63-66; THM status l. 4-12 and `related` l. 15-22):
  - cite READY section 4, RESULT section 4 with eqs. (8), (10) and (12), and ONE_EXIT section 1;
  - narrow the claimed novelty to the end-count lemma, non-strong bases and the grandchild level.

**M2 (correctness): note l. 91-92 is false.** The note says "So **no old mark of any grandchild ever has digit 1 or 2**".
- In each of the 3 digit-2 grandchildren C2(TT_{n-1} => C3, u), with u in the 3-cycle, the old mark u also has digit 2,
  beside y2.
- This holds for every \|B\| = 1..8 (`aud_heredity_*.out`) and under permutation brute force (`aud_bruteforce.out`). It is
  READY section 4's C3 boundary (gap list (2, 14, 12, 0, 2) in the funnels AUDIT).
- `c2_digit_heredity_20261104.c` l. 287-290 never record which mark carries digit 2, so the sentence was inferred from a
  count.
- **Fix:** "each of the three has exactly two digit-2 marks, the native y and the second mark u". Make the same change in
  THM (iii) and in the HYP-9466 results.

**M3 (typing): "All three known self-reproducing routes stop at generation two".** Affected: note l. 96-103 and typing row
l. 137; THM Reading l. 69-73; `results/INDEX.md` l. 149.
- H-neutral preparation is an iff criterion (PREPARATION_REDUCTION section 3). Its closure is "open either way"
  (HYP-9454 (2)).
- klein's data show the route alive past generation two (THM-4862 (f), section 7). All four level-2 leaves of the 1429 tree
  are digit-1 ready; both +4 branches close weak-allowed; 7 of the 8 level-3 values are ready at order <= 15.
- RESULT section 3's zero channel (digit 0) and READY section 5's R-stable families reproduce at every depth.
- Only two routes are shown to stop at generation two: the C2 self-marks, and Session191's fixed-core same-clock fans.
- **Fix:** keep those two routes and retype preparation as OPEN.

**M4 (logic): "a ternary parent for full spectrum coverage must come from fresh bases" is false.** Affected: THM title l. 3
and l. 74; note l. 96-99 and 121-122.
- Many Hamiltonian path counts have no ternary parent at all.
- Two facts locate small digits:
  - a mark outside the last strong component has digit >= H, since B_u(M => R) = B_u(M) K(R) with K(R) >= 2h(R);
  - in the last component S, sigma_u(T) = h(T - S) sigma_u(S).
- So a digit-1 or digit-2 mark of a tournament with H >= 3 lives in a strong S with H(S) = H. This is check (L) in
  `aud_proofcheck.out`.
- By Busch's bound, every strong S with H <= 42 has order <= 8. So all classes of order <= 8 decide every predecessor
  h <= 42.
- Result (`aud_sweep_8.out`, last line): 20 odd N <= 130 have no ternary parent. 14 of them are still attained as h + 2B_a
  with a digit >= 3: 23, 25, 49, 55, 63, 65, 77, 83, 85, 103, 107, 109, 125 and 127.
- Example: 77 = 3·25 + 2, and no tournament with H = 25 carries a digit-1 mark.
- This is an existence obstruction, independent of heredity. For digit 1 it already follows from klein's complete R_1 gaps
  (25, 27, 35, 41, 43), because digit-1 marks are sole exits by (D).
- klein's numerations do not use ternary parents: they allow any digit on 3h + 2 sigma = N.
- "Any proof ... must use fresh bases with growing digits" (THM l. 74) is unproved. Near-join, tower and two-base joins
  are other routes.
- **Fix:** "full coverage cannot come from ternary parents; it needs digits >= 3 on 3h + 2 sigma = N (HYP-9454 (6))".

## 3. MINOR findings

- **m1** (note l. 52-53). The m = 1 parenthetical is muddled. The end counts are then H and 0, which excludes H - 1 and
  H - 2 for H >= 3; "H is odd" is not the reason.
- **m2** (note l. 54-55; THM l. 36). "Open a Hamiltonian cycle before that vertex" yields a path that *starts* there; cut
  after the vertex instead. Cite Camion.
- **m3** (note l. 56; THM l. 37-41).
  - The phrase "leaves at least (m - 1) M <= d" is garbled.
  - The case a outside C (e_a = 0, so d = H) is missing.
  - THM (i) omits the case m = 1, d = 2.
- **m4** (note l. 38). "not found elsewhere in the repository" holds for the statement, but credit:
  - READY section 4 (the strong case, gap >= 3);
  - the inequality h - s >= q (`orbit_spectrum_session_20261103/coverage/NOTE.md` section 1;
    `source_selection_closure_20261103/README.md` section 3);
  - the component product in PREPARATION_REDUCTION section 2, FUNNELS F1, and Session191 UNIT_START_BOUNDARY sections 1
    and 3 (16:18, before the S77 checkpoint).
- **m5** (note l. 81-88; THM l. 60-65). The tables omit the \|B\| = 1 row (1 pair; 3 grandchildren with digit 2), although
  the totals include it.
- **m6** (note l. 90; THM (iii)).
  - Counts are (class, labelled mark) pairs. The 3 TT => C3 children are one rooted class, and so are the 3 digit-2
    grandchildren.
  - The 4 digit-2 children carry digit 2 at the old mark a (and at y in three of them). Say so.
- **m7** (note l. 96-97). "Apart from the transitive seeds, no C2 descendant carries digit 1" has no stated range; the
  S77 data reach only order 11. Cite section 4 or qualify the claim.
- **m8** (note l. 108-112).
  - The program includes non-strong members (b_t = 1), which FUNNELS F1 excludes; they add only H = 1.
  - "no single explicit digit-1 family is a ternary parent" over-generalises.
  - The comment at `ternary_parent_digit_families_20261104.c` l. 98 claims a check of H = h(C) + F(C, r) - 1 that the
    code never performs (l. 100). The identity is true: the audit's closed form agrees.
- **m9** (note l. 115-117).
  - "h <= 1,500" is vacuous (the maximum is 661).
  - The any-mark digit-1 set is {1} together with klein's strong sole-exit R_1 (64 values), as READY section 1 predicts.
  - Digit 2 adds 13 values, all from READY section 2's exceptional families: 5, 11, 13, 19, 25, 31, 39, 47, 49, 55, 95,
    97, 115.
- **m10** (note l. 150-152). There are no committed outputs for `anymark_c2_digits_20261104` or
  `ternary_parent_digit_families_20261104`. Commit them.
- **m11** (`c2_digit_heredity_20261104.c` and `anymark_c2_digits_20261104.c`).
  - Lines 1-24 are the header of an unrelated "fewcycle" program (`Usage: fewcycle NMAX K`).
  - Lines 123-198 are dead code (`all_pairs_sub`, `ham_cycles`, `is_strong_mask`, `fib`, `Rec`).
- **m12** (note l. 124-125; THM l. 74).
  - The any-vertex line statement from 69 is HYP-9454 (6), not (1).
  - Its FINITE-EXACT range is [69, 86,229]. The figure 8,079,323 is the union of klein's three numerations.
- **m13** (note l. 29; THM l. 71). "densities about 0.75" overstates R_2, whose complete windows are 0.52-0.79.
- **m14.** THM-4862 (f) already uses "heredity" for one-vertex extensions of ready classes. Say "C2 digit heredity".

## 4. Upgrade: the corrected (iii) is PROVED at every order (identities checked for \|B\| <= 8)

**Statement.** Let T1 = C2(B, a).
- T1 has a digit-1 mark iff B is transitive and a is not its sink (the mark is y).
- T1 has a digit-2 mark iff B = TT_n with a penultimate (mark a), or B = TT_j => C3 with a in the C3 (marks a and y).
- Hence no grandchild has a digit-1 mark.
- The digit-2 grandchildren are exactly those with B transitive, a its sink and u in {a, x, y}. Each has digit 2 at y2
  and at u.

**(D) Out-degree bound.** sigma_v >= d(d - 1) for every vertex v of out-degree d. This is READY (7), dualised; a direct
proof follows.
- sigma_v counts the covers (P1, P2) in which P2 starts at v and P1 is nonempty and ends at an out-neighbour of v. The
  remaining covers correspond to the Hamiltonian paths.
- Fix out-neighbours w != w' and a Hamiltonian path q_1..q_m of T - {v, w, w'}. Among the splittings
  (q_1..q_k w, v w' q_(k+1)..) and (q_1..q_k w', v w q_(k+1)..), the valid ones are the transitions of the bit pair
  ([q_k -> w], [q_k -> w']) along a walk from (1, 1) to (0, 0).
- Leaving (1, 1) contributes at least 1, and entering (0, 0) contributes at least 1 (2 if both happen at once). So each
  pair {w, w'} contributes at least 2 covers.

**Proof of the statement.**
1. In T1, x and y have out-degree 1, a has d_B(a) + 1, and every other vertex v has d_B(v) + 2. By (D), a small digit can
   only sit at x, at y, at the sink z of B (z != a), or at a with d_B(a) <= 1.
2. sigma_x = 0. sigma_y = h - e_a, which the end-count lemma settles.
3. sigma_z >= h(B - z) + B_a(B - z) + 2 >= 4: the first two terms are the paths of C2(B - z, a) that end at x or at y.
4. The mark a.
   - If d_B(a) = 0, then sigma_a = 0, because x is the sink of T1 - a.
   - If a -> w is its only out-arc, then sigma_a = e_w(U) + B_w(U) with U = B - a (RESULT (12) plus the one-exit law).
   - B_w(U) >= h(U), with B_w(U) = 1 iff U is transitive with sink w, and then e_w = 1. Hence sigma_a != 1, and
     sigma_a = 2 iff U is transitive with w its sink or its penultimate vertex.
5. For grandchildren, apply 1-4 with base T1. T1 is never transitive (x -> a -> y -> x) and has no sink, so
   T1 = TT => C3 iff h + 2B_a = 3, i.e. iff B is transitive with a its sink.

**Checks.** `aud_proofcheck.out` covers all 58,652 pairs (B, a) with 2 <= \|B\| <= 8: RESULT (12), the sink bound and the
a -> w identity give 0 failures, and small digits occur only at y or a. READY section 4 (dualised) also checks out: no
child of a strong B with 4 <= \|B\| <= 8 has a small digit (50,779 cases).

## 5. Prior art and same-day overlap

- **Heredity no-go for strong inputs (all orders, all marks): already PROVED.**
  - READY section 4: commit 4ed40bd587, 2026-11-03 10:46.
  - RESULT section 4: commit 8923c5fa5e, 2026-11-03 21:03.
  - ONE_EXIT section 1: 2026-11-04 00:44.
- **End-count lemma.** No earlier statement found. I searched for H - 1 / H - 2 end counts, "all but one", "h - s >= q"
  and "last strong component". The ingredients are listed in m4. It is an elementary corollary of the component product;
  I know of no literature source.
- **Small digits are rare.** klein THM-4862 (b) (digit = cycle count) and (f) (complete R_1 gaps, which already imply M4
  for digit 1).
- **Same-day overlap.** I read all 52 origin/main commits since 15:30 (through 38ee304be2, 16:25:44).
  - Only Session191 is adjacent.
  - SECOND_GENERATION ((61, 63)) is credited correctly.
  - UNIT_START_BOUNDARY uses the same technique and is not credited (m4).
  - CYCLE_CROSSING is a generation-1 fan whose continuation is open.
  - two_apex and deletion_growth do not overlap.
  - Session190 CP3 section 5.1 is cited correctly.
  - CP4, kind-pasteur S10, mac-mini THM-4871 and opus S76 do not overlap.
  - Nothing after 16:23 states the lemma or the no-go.

## 6. Files

- `aud_c2.c` (`gcc -O3 -march=native -o aud_c2 aud_c2.c`; modes `classes`, `heredity 1 7`, `heredity 8 8`,
  `endcount 7 9`, `sweep 8 1500`, `proofcheck 8`, `triples N`, `family 11 0`); `aud_family.py` (closed form vs DP vs
  brute force); `aud_bruteforce.py` (permutation-only checks, including the M2 cases).
- Outputs: `aud_heredity_1_7.out`, `aud_heredity_8.out`, `aud_endcount.out`, `aud_sweep_8.out`, `aud_proofcheck.out`,
  `aud_family.out`, `aud_family_dp_t11.txt`, `aud_bruteforce.out`.
