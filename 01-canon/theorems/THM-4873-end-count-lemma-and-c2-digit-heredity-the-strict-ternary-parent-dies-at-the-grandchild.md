---
id: THM-4873
title: "The end-count lemma and complete C2 digit heredity. (i) In a tournament with H >= 3 Hamiltonian paths no vertex ends exactly H - 1 of them, and a vertex ends exactly H - 2 iff the tournament is TT_j => C3 and the vertex lies in the 3-cycle (dually for starts). (ii) The C2 child of (B, a) has a digit-1 mark iff B is transitive and a is not its sink (at y), and a digit-2 mark iff B = TT_n with a penultimate (at a) or B = TT_j => C3 with a in the C3 (at a and y). Hence no C2 grandchild has a digit-1 mark, and the digit-2 grandchildren are exactly C2(C2(TT_n, sink), u) with u in the 3-cycle (digit 2 at y2 and at u). (iii) Full coverage cannot come from ternary parents: 20 odd N <= 130 have none"
status: >
  PROVED: (i) (opus-2026-11-04-S77, elementary). (ii) at every order: the proof is the S77 independent audit's (AUDIT.md
  section 4), using the out-degree bound sigma_v >= d(d - 1) (READY (7), dualised), the end-count lemma and
  RESULT eqs. (8), (10), (12). (ii) is FINITE-EXACT at orders 1-8; order 7 was pre-registered as HYP-9466 and CONFIRMED.
  FINITE-EXACT: (iii), from the S77 audit, using Busch's bound to decide every predecessor h <= 42 at orders <= 8.
  INHERITED: the heredity no-go for strong inputs (two_entry_funnels_20261103/READY.md section 4;
  multiport_defect_grammar_20261103/construction/RESULT.md section 4; ONE_EXIT_AND_PERSISTENCE section 1); the C2 law at
  any mark (RESULT law (4); THM-4865).
  NOT CLAIMED: Chan-Pak completeness, which is OPEN. H-neutral preparation is also OPEN, not a no-go.
  AUDIT: PASS WITH CORRECTIONS (4 majors, 14 minors, all applied; MISTAKE-713);
  04-computation/experiments/c2_digit_heredity_20261104_audit_S77/AUDIT.md.
session: opus-2026-11-04-S77
source: 05-knowledge/results/ternary_parent_c2_digit_heredity_20261104.md
related:
  - THM-4862
  - THM-4865
  - HYP-9454
  - HYP-9466
  - THM-4872
  - 05-knowledge/results/two_entry_funnels_20261103/READY.md
  - 05-knowledge/results/multiport_defect_grammar_20261103/construction/RESULT.md
  - 05-knowledge/results/session191_tournament_growth_20261104/ternary_fan/PROOF.md
  - 05-knowledge/results/nonzero_preparation_20261103/README.md
verification:
  - 04-computation/experiments/c2_digit_heredity_20261104.c
  - 04-computation/experiments/c2_digit_heredity_20261104.out
  - 04-computation/experiments/anymark_c2_digits_20261104.c
  - 04-computation/experiments/anymark_c2_digits_20261104.out
  - 04-computation/experiments/ternary_parent_digit_families_20261104.c
  - 04-computation/experiments/ternary_parent_digit_families_20261104.out
  - 04-computation/experiments/c2_digit_heredity_20261104_audit_S77/AUDIT.md
---

# THM-4873 - the end-count lemma; complete C2 digit heredity; no ternary parents for full coverage

## (i) End-count lemma

Let C be the last strong component of T, of order m, and M = prod H(earlier components). All paths end in C, and
e_v = M e_v(C) there.
- A vertex outside C ends 0 paths.
- If m = 1, the end counts are H and 0.
- If m >= 3, C is Hamiltonian (Camion). Cutting a Hamiltonian cycle after any vertex shows that each vertex of C ends at
  least M paths.

So a vertex ending H - d paths forces (m - 1) M <= d:
- d = 1 is impossible for H >= 3;
- d = 2 forces T = TT_j => C3 with the vertex in the 3-cycle (H = 3).

## (ii) Complete small-digit classification under C2

The S77 audit's proof (AUDIT.md section 4):
- **(D) Out-degree bound.** sigma_v >= d(d - 1) for out-degree d (READY (7), dualised). So in T1 = C2(B, a) a small
  digit can sit only at x, at y, at the sink z of B, or at a with d_B(a) <= 1.
- **x and y.** sigma_x = 0 and sigma_y = h - e_a. The latter is decided by (i).
- **The sink.** sigma_z >= h(B - z) + B_a(B - z) + 2 >= 4.
- **The mark a.**
  - If d_B(a) = 0, then sigma_a = 0.
  - If a -> w only, then sigma_a = e_w(B - a) + B_w(B - a). This is never 1, and it equals 2 iff B - a is transitive
    with w its sink or penultimate vertex.
- **Grandchildren.** T1 is never transitive and has no sink. So T1 = TT => C3 iff B is transitive with its sink as
  mark.

## (iii) No ternary parents for full coverage (FINITE-EXACT)

- A digit-1 or digit-2 mark of a tournament with H >= 3 lies in a strong component S with H(S) = H.
- By Busch's bound such S with H <= 42 have order <= 8.
- So 20 odd N <= 130 have no ternary parent from any base. 14 of them are attained with a digit >= 3: 23, 25, 49, 55,
  63, 65, 77, 83, 85, 103, 107, 109, 125, 127.
- Coverage therefore needs digits >= 3 on 3h + 2 sigma = N (HYP-9454 (6), FINITE-EXACT on [69, 86,229]; klein's three
  numerations together reach 8,079,323), or other routes.
