# The ternary parent and C2 digit heredity: the end-count lemma, the complete small-digit classification of C2 children and grandchildren, and why full coverage cannot come from ternary parents

**Session:** opus-2026-11-04-S77. Worktree `C:/tmp/wt-s77`. Canon: THM-4873. Pre-registration: HYP-9466. Corrected after
the independent audit (section 7; MISTAKE-713).

**Owner prompt (verbatim).** "start another session on the ternary parent for full spectrum coverage"

## 0. The problem and what is already known

**The target.** Chan-Pak Conjecture 6.2 (TCS 1015 (2024) 114776, arXiv:2308.10214): every odd number except 7 and 21 is
the Hamiltonian path count of some tournament. OPEN.

**The C2 law.** For a tournament B and a mark a, add x with x -> a only (every other vertex beats x), then y beaten by
all of B with y -> x. The child has

    H = h(B) + 2 B_a(B) = 3 h(B) + 2 sigma_a,    sigma_a = B_a - h,  B_a = sum_{S subset V - a} h(B[S]) s_a(B[V - S]).

This is `multiport_defect_grammar_20261103/construction/RESULT.md` law (4). klein's THM-4865 states it at any vertex.

**The ternary parent.** Every odd N has a unique odd predecessor under 3h, 3h + 2, 3h + 4. A *ternary parent* of N is a
realiser of that predecessor carrying a mark of the needed digit 0, 1 or 2.

**Prior art on C2 digit heredity** (S77 audit, M1; not cited in the first version):
- `two_entry_funnels_20261103/READY.md` section 4 (2026-11-03, audited) proves that the reset, the converse of C2, of any
  strong input of order >= 4 has no 1- or 2-ready vertex under any re-marking. Its proof contains the native bound
  h - s >= 3, and section 4 records the C3 boundary.
- RESULT.md section 4 dualises this: "C2 cannot produce either next nonzero digit at any auxiliary". Its eqs. (8), (10)
  and (12) give sigma_x = 0 and sigma_y = h - e_a.
- `nonzero_preparation_20261103/algebra/ONE_EXIT_AND_PERSISTENCE.md` section 1 says the same.

**New here:**
- the end-count lemma;
- the extension to non-strong bases;
- the complete classification of small digits in C2 children and grandchildren at every order (proof completed by the
  audit, section 2).

**Related computations, credited.**
- klein's root-response numerations attain every odd number in [69, 8,079,323]. That figure is the union of three
  numerations (THM-4865).
- The any-vertex line statement HYP-9454 (6) is FINITE-EXACT on [69, 86,229].
- The full-digit sweep covers [307, 26551] (THM-4862 (d)).
- The densities of the ready spectra are 0.65-0.77 for R_1 and 0.52-0.79 for R_2 (THM-4862 (f)).
- Session191 builds two-port ternary fans by fixed-core wrappers; its second generation fails at (61, 63).

## 1. The end-count lemma (PROVED, elementary)

**Lemma.** Let T be a tournament with H = H(T) Hamiltonian paths.
- (i) If H >= 3, no vertex is the end of exactly H - 1 Hamiltonian paths.
- (ii) A vertex is the end of exactly H - 2 Hamiltonian paths iff T = TT_j => C3 and the vertex lies in the 3-cycle.

The same holds for starts, via the converse.

*Proof.*
- Let C be the last strong component of T, of order m, and M = prod H(earlier components) >= 1. Every Hamiltonian path
  ends in C. So e_v = 0 for v outside C, and e_v(T) = M e_v(C) for v in C; also H = M H(C).
- **A vertex outside C** has e = 0, which equals H - d only if d = H >= 3. So it is excluded for d = 1, 2 when H >= 3.
- **m = 1.** The end counts are H (at the single vertex of C) and 0. Neither equals H - 1 or H - 2 when H >= 3.
- **m >= 3.** C is strong, hence Hamiltonian (Camion 1959). Each of its vertices ends at least one Hamiltonian path of
  C: cut a Hamiltonian cycle just after that vertex. So every vertex of C ends at least M paths of T.
- If a vertex a ends H - d paths, the other m - 1 vertices of C together end at most d, so (m - 1) M <= d.
  - For d = 1 this forces m <= 2, impossible.
  - For d = 2 it forces m = 3 and M = 1: C = C3, the earlier components are single vertices, and H = 3. Conversely,
    each vertex of the 3-cycle of TT_j => C3 ends exactly 1 = H - 2 paths.

QED.

**Credit for the technique.** The component-product and cycle-cutting technique also appears in:
- READY section 4 (the strong case, with gap h - s >= 3);
- the inequality h - s >= q (`orbit_spectrum_session_20261103/coverage/NOTE.md` section 1;
  `source_selection_closure_20261103/README.md` section 3);
- PREPARATION_REDUCTION section 2 and FUNNELS F1;
- Session191's UNIT_START_BOUNDARY sections 1 and 3 (16:18 the same day).

## 2. Small digits of C2 children and grandchildren (PROVED at every order)

**Native digits** (RESULT eqs. (8), (10), (12); klein's tower lemmas): sigma_x = 0 and sigma_y = h(B) - e_a(B).

**Theorem (C2 digit heredity, complete).** Let T1 = C2(B, a).
- T1 has a digit-1 mark iff B is transitive and a is not its sink. The mark is y.
- T1 has a digit-2 mark iff B = TT_n with a its penultimate vertex (digit 2 at the old mark a), or B = TT_j => C3 with a
  in the C3 (digit 2 at a and at y).
- Consequently, no C2 grandchild has a digit-1 mark. The digit-2 grandchildren are exactly
  C2(C2(TT_n, sink), u) with u in {a, x, y}, the 3-cycle of TT => C3. Each has exactly two digit-2 marks: the native y2
  and the old mark u.

The proof is the S77 audit's (AUDIT.md section 4), checked here and recorded with credit.

- **(D) Out-degree bound** (READY (7), dualised): sigma_v >= d(d - 1) for every vertex v of out-degree d.
  - sigma_v counts covers (P1, P2) with P2 starting at v and P1 nonempty ending at an out-neighbour of v.
  - For two out-neighbours w != w' and a Hamiltonian path q_1..q_m of T - {v, w, w'}, consider the splittings
    (q_1..q_k w, v w' ...) and (q_1..q_k w', v w ...). Their validity follows the bit pair ([q_k -> w], [q_k -> w'])
    from (1, 1) to (0, 0).
  - Leaving (1, 1) and entering (0, 0) give at least 2 covers per pair {w, w'}.
- **Where small digits can sit.** In T1, x and y have out-degree 1, a has d_B(a) + 1, and every other vertex v has
  d_B(v) + 2. By (D), a digit <= 2 can sit only at x, at y, at the sink z of B (z != a), or at a with d_B(a) <= 1.
- **sigma_x and sigma_y.** sigma_x = 0, and sigma_y = h - e_a is settled by the end-count lemma.
- **The sink.** sigma_z >= h(B - z) + B_a(B - z) + 2 >= 4. The first two terms count the paths of C2(B - z, a) that end
  at x or at y.
- **The mark a.**
  - If d_B(a) = 0, then sigma_a = 0, since x is the sink of T1 - a.
  - If a -> w is its only out-arc, then sigma_a = e_w(U) + B_w(U) with U = B - a.
  - B_w(U) >= h(U), with B_w(U) = 1 iff U is transitive with sink w, in which case e_w = 1. So sigma_a != 1, and
    sigma_a = 2 iff U is transitive with w its sink or its penultimate vertex.
- **Grandchildren.** Apply the above with base T1. T1 is never transitive (x -> a -> y -> x) and has no sink. So
  T1 = TT => C3 iff h + 2B_a = 3, i.e. B is transitive with its sink as mark.

**FINITE-EXACT checks.** All classes of order 1-7 here, and of order 8 by the audit. Counts are (class, labelled mark)
pairs:

| \|B\| | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|---|
| children | 1 | 2 | 6 | 16 | 60 | 336 | 3,192 |
| with a digit-1 mark | 0 | 1 | 2 | 3 | 4 | 5 | 6 |
| with a digit-2 mark | 0 | 1 | 4 | 4 | 4 | 4 | 4 |
| grandchildren | 3 | 8 | 30 | 96 | 420 | 2,688 | 28,728 |
| with a digit-1 mark | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| with a digit-2 mark | 3 | 3 | 3 | 3 | 3 | 3 | 3 |

- The 4 digit-2 children carry digit 2 at the old mark a, and in three of them also at y.
- The 3 digit-2 grandchildren form one rooted class.
- At |B| = 8 the audit found 7 and 4 small-digit children, and 0 and 3 grandchildren among 550,400. Its checks
  (`aud_proofcheck.out`) cover all 58,652 pairs with |B| <= 8.

## 3. What this means for the ternary parent

- **C2 self-marks do not reproduce small digits.** This settles the C2-self-mark route at every order. So does
  Session191's fixed-core same-clock fan, at its second generation.
- **H-neutral preparation is still OPEN.** It is an iff criterion, membership of h in R_s (PREPARATION_REDUCTION
  section 3), and klein's data show it alive past generation two (THM-4862 (f), section 7). The zero-digit channel
  (RESULT section 3) and READY section 5's R-stable families reproduce at every depth.
- **Full coverage cannot come from ternary parents at all** (S77 audit, M4).
  - A digit-1 or digit-2 mark of a tournament with H >= 3 lives in a strong component S with H(S) = H. A mark outside
    the last strong component has digit >= H, and in the last component sigma_u(T) = h(T - S) sigma_u(S).
  - By Busch's bound, every strong S with H <= 42 has order <= 8. So the classes of order <= 8 decide every
    predecessor h <= 42.
  - Result: 20 odd N <= 130 have no ternary parent. 14 of them are nevertheless attained as h + 2B_a with a digit >= 3:
    23, 25, 49, 55, 63, 65, 77, 83, 85, 103, 107, 109, 125 and 127.
  - For example 77 = 3 * 25 + 2, and no tournament with H = 25 carries a digit-1 mark.
- **So coverage needs digits >= 3** on the line 3h + 2 sigma = N (HYP-9454 (6)), or other routes: near-join, tower,
  two-base joins. Chan-Pak completeness stays OPEN.

## 4. Supporting data

**Explicit digit-1 family** (`ternary_parent_digit_families_20261104`; output `..._20261104.out`).
- Construction: unit-start rooted tournaments with a single-vertex core and a free transitive tail, conversed and given
  a sole exit.
- The program also includes non-strong members (b_t = 1), which FUNNELS F1 excludes; they add only H = 1.
- The values cover 36.1% of the odd numbers up to 30,000 (audit).
- The identity H = h(C) + F(C, r) - 1 is true (the audit's closed form agrees), although the program's comment claims
  a check it does not perform.

**Any-mark one-step sweep** (`anymark_c2_digits_20261104`; output committed).
- From all 7,412 classes of order <= 8 and every mark, the odd numbers missed below the maximum 661 (and in the stated
  window up to 1,500) are {1, 7, 11, 13, 21, 47, 67}. This re-derives klein's height-two numeration.
- The any-mark digit-1 set is {1} together with klein's strong sole-exit R_1 (64 values), as READY section 1 predicts.
- Digit 2 adds 13 values, all from READY section 2's exceptional families.

## 5. Typing

| claim | type |
|---|---|
| End-count lemma | PROVED (elementary; technique credited) |
| Native digits sigma_x = 0, sigma_y = h - e_a | inherited (RESULT (8), (10), (12)); re-verified |
| Complete small-digit classification of C2 children; no digit-1 grandchild; the digit-2 grandchildren | PROVED at every order (S77 audit section 4, with READY (7)); FINITE-EXACT orders 1-8 |
| Heredity no-go for strong inputs | inherited (READY section 4; RESULT section 4) |
| No ternary parent for 20 odd N <= 130 | FINITE-EXACT (audit; decided by orders <= 8 via Busch) |
| Preparation closure | OPEN (not a no-go) |

## 6. Pre-registration record

HYP-9466 (commit d96e81539f) pre-registered P1-P3 at |B| = 7 after the exploratory run at |B| <= 6. All hold. The
first version's reading of P3 ("no old mark has digit 2") was wrong: the old mark u also has digit 2 in each of the
three grandchildren (M2).

## 7. Audit

Independent adversarial audit:
[`AUDIT.md`](../../04-computation/experiments/c2_digit_heredity_20261104_audit_S77/AUDIT.md). Verdict:
**PASS WITH CORRECTIONS.**

**Reproduced from scratch.**
- Everything in the session, with the audit's own code and klein's nauty-based rows.
- An out-of-sample check at |B| = 8.
- The end-count lemma on 2,131,019 labelled tournaments and on all classes to order 9.

**Majors, applied.**
- **M1.** The prior art: READY section 4, RESULT section 4, ONE_EXIT section 1.
- **M2.** A false sentence about the old marks.
- **M3.** Preparation is not shown to stop.
- **M4.** No ternary parent exists for 20 odd N <= 130.

**Minors, applied.**
- The lemma proof wording (Camion; the outside and m = 1 cases).
- The |B| = 1 row.
- Count conventions.
- Credits.
- HYP-9454 numbering and range.
- The R_2 densities.
- Committed outputs.
- Program headers.

**Upgrade.** The audit proved the corrected classification at every order (section 2).

**Logged as MISTAKE-713.**

## 8. Reproduction

    gcc -O3 -march=native -o her c2_digit_heredity_20261104.c && ./her 7
    gcc -O3 -march=native -o any anymark_c2_digits_20261104.c && ./any 8 6 1500
    gcc -O3 -march=native -fopenmp -o fam ternary_parent_digit_families_20261104.c && ./fam 6 14 30000
