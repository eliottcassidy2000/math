---
id: THM-4875
title: "Sharp old-vertex digits for C2 children, closed one-arc gadget entries, and one-apex joins to 2^32. (A) In C2(B, a), the sink z != a of B has digit exactly 3h(B-z) + 2R_a(B-z) - e_a(B-z) (>= 4; the bound is THM-4873 (ii)), equal to 4 iff (TT_n, a penultimate); every other old vertex other than the mark has digit >= 9, equal to 9 iff (TT_n, a sink, w penultimate). (B) Closed subset-sum entries for the one-arc gadget transfer matrix M(G, u, z); the transfer law and the adjoint M(G^op, z, u) = J M^T J are KNOWN (SEAM.md, 2026-10-29). (C) A descent lemma against unimodular cone-covering, with its corrected and limited consequence. (D) Every odd number up to 2^32 except 7 and 21 is the Hamiltonian path count of a one-arc join of a one-apex tournament with a one-apex tournament or a rooted tournament of order <= 7 (FINITE-EXACT; Chan-Pak OPEN)"
status: >
  PROVED:
  (A2) the exact sink digit and its equality case;
  (A3) the non-sink bound >= 9 with its equality case (proved by the independent audit, U1; HYP-9468 P3's all-order
  conjecture);
  (B) the closed entries (new presentation);
  (C) the descent lemma and its corrected consequence (audit U2).
  KNOWN:
  (A1) and the bounds >= 4 / >= 6 and (A5) (THM-4873 (ii), S77 audit section 4; (A5) also RESULT (12));
  (B) the general transfer and the Lorentzian adjoint (seven_return_native_recursions_20261029/SEAM.md eqs. (1)-(6);
  golden_hamilton_transfer_20261103/RECURRENCE.md (G1)-(G4); session29 ADDENDUM section 3 for 18 tiles).
  FINITE-EXACT:
  (A) on all 7,411 classes of order 2-8 and their 58,652 marks (order 8 pre-registered as HYP-9468 P3), and on every
  labelled B of order 2-6;
  (B) closed entries against direct counts for every labelled gadget of order <= 5, adjoint for all 1,206,353 marked
  gadgets of order <= 6;
  (D) pre-registered as HYP-9468 P1 and reproduced independently by the audit.
  INHERITED: the one-arc join law (session29 ADDENDUM sections 2 and 5; THM-4866 (B)); the one-apex triple (session25;
  session179 spectrum/REPORT.md (1); g = session124's f; Session191 one_apex/PROOF.md; THM-4876 (6)); C2's native digits
  (RESULT eqs. (8), (10), (12)); the impossibility of 7 and 21 (THM-1370).
  NOT CLAIMED: Chan-Pak completeness, which is OPEN.
  AUDIT: PASS WITH CORRECTIONS (3 majors, 16 minors, 4 upgrades; all applied; MISTAKE-716);
  04-computation/experiments/old_mark_one_apex_join_20261104_audit_S78/AUDIT.md.
session: opus-2026-11-04-S78
source: 05-knowledge/results/old_mark_bound_lorentzian_gadgets_one_apex_join_coverage_20261104.md
related:
  - THM-4873
  - THM-4876
  - THM-4866
  - THM-4865
  - THM-4862
  - THM-1370
  - HYP-9468
  - HYP-9454
  - HYP-9453
  - 05-knowledge/results/seven_return_native_recursions_20261029/SEAM.md
  - 05-knowledge/results/golden_hamilton_transfer_20261103/RECURRENCE.md
  - 05-knowledge/results/session29_framed_profiles_20261022/session29-profiles/ADDENDUM.md
  - 05-knowledge/results/session191_tournament_growth_20261104/one_apex/PROOF.md
  - 05-knowledge/results/session191_tournament_growth_20261104/two_apex/TWO_APEX.md
  - 05-knowledge/results/square_sum_spectrum_overnight_20261030/README.md
  - 05-knowledge/results/session193_triangulated_proof_20261105/one_arc_joins/PROOF.md
  - 05-knowledge/results/two_entry_funnels_20261103/READY.md
  - 05-knowledge/results/multiport_defect_grammar_20261103/construction/RESULT.md
verification:
  - 04-computation/experiments/old_mark_digit_bound_20261104.c
  - 04-computation/experiments/old_mark_digit_bound_20261104.out
  - 04-computation/experiments/nonsink_digit_nine_check_20261104.py
  - 04-computation/experiments/nonsink_digit_nine_check_20261104.out
  - 04-computation/experiments/gadget_adjoint_20261104.c
  - 04-computation/experiments/gadget_adjoint_20261104.out
  - 04-computation/experiments/rooted_triple_bank_20261104.c
  - 04-computation/experiments/rooted_triple_bank_order7_20261104.txt
  - 04-computation/experiments/chanpak_one_apex_join_cover_20261104.c
  - 04-computation/experiments/chanpak_one_apex_join_cover_20261104.out
  - 04-computation/experiments/join_witness_check_20261104.c
  - 04-computation/experiments/join_witness_check_20261104.out
  - 04-computation/experiments/one_apex_triple_check_20261104.py
  - 04-computation/experiments/old_mark_one_apex_join_20261104_audit_S78/AUDIT.md
---

# THM-4875 - sharp old-vertex digits; closed one-arc gadget entries; one-apex joins to 2^32

Notation. The root response is R_t = sum_{X subset V - t} h(T[X]) s_t(T[V - X]) (klein's B_t), the digit is
sigma_t = R_t - h, and the rooted triple is (h, e_t, R_t). C2(B, a) adds x -> a (beaten by the rest of B), then y,
beaten by all of B, with y -> x.

## (A) Old-vertex digits of a C2 child

Let T1 = C2(B, a), with |B| >= 2.

- **(A1)** (KNOWN, THM-4873 (ii) (D)) For every vertex v,
  sigma_v >= E_v + 2 sum over pairs {p, q} subset N+(v) of h(T - {v, p, q}) >= d(d - 1),
  where E_v = h(T - v) - e_v(T).
- **(A2)** (PROVED) If z != a is the sink of B, then sigma_z(T1) = 3h(B - z) + 2R_a(B - z) - e_a(B - z).
  - This is >= 4; the bound itself is KNOWN from THM-4873 (ii).
  - Equality holds iff B is transitive and a is its penultimate vertex.
- **(A3)** (PROVED; audit U1) Every non-sink old vertex w != a has sigma_w(T1) >= 9. Equality holds iff B is transitive,
  a is its sink and w is its penultimate vertex.
  - For B = TT_m + w with w -> u_j only (j < m) and a = u_m, sigma_w = 4 + 2[j = m - 1] + 3 * 2^(m-j).
- **(A4)** Hence the old vertices other than the mark have digit >= 4. The value 4 occurs only at the sink of (TT_n, a
  penultimate), and the next possible minimum, 9, only at the penultimate vertex of (TT_n, a sink).
- **(A5)** (KNOWN, THM-4873 (ii); RESULT (12): sigma_a = 2R_a(B) - h(B) - h(B - a)) The mark's digit is:
  - 0 if a is the sink of B;
  - e_w(B - a) + R_w(B - a), never 1, if a -> w is its only out-arc;
  - at least (d + 1)d if d_B(a) = d >= 2.

**Proofs.** Source note, section 1.
- (A2): case split on the start of P2. In T1, x's only out-neighbour is a and y's is x.
- (A3): (A1) for d_B(w) >= 2. For d_B(w) = 1, (A1) at the three out-neighbours {c, x, y}, plus one extra cover in the
  transitive case.

## (B) Closed entries of the one-arc gadget transfer

The general law is KNOWN (SEAM.md eqs. (1)-(6), 2026-10-29). For a gadget G with entry u and exit z, and a base (T, v),
let N = T => G with u -> v reversed. Then (h, e, R)(N, z) = M(G, u, z)(h, e, R)(T, v) and
M(G^op, z, u) = J M(G, u, z)^T J with J = diag(1, -1, 1).

New here are the closed entries:

    M = [[h_G, -s_u, E_u], [e_z, -P_(u->z), E_(u,z)], [R_z, -J_1, J_2 + K]],

with E_u, E_(u,z), J_1, J_2 and K the subset sums defined in the source note, section 2.

## (C) Descent lemma (PROVED, elementary)

**Lemma.** Integer matrices of determinant ±1 that strictly increase a linear functional with nonnegative integer values
on a closed cone C cannot cover C by their images.

**Corrected consequence (audit U2).**
- Hypotheses: C has nonempty interior inside {0 <= e <= h, R >= (1 + c)h}, M_i are unimodular gadget matrices, and Λ is
  a full-rank lattice coset.
- Conclusion: no finite F gives C ∩ Λ subset F ∪ union M_i(C ∩ Λ).

**Scope.** No cone-intersect-coset of triples is realisable, so this excludes no proof of Chan-Pak. It only restates
the heuristic of THM-4865 (f) and HYP-9453 section 2.

**Equality.** On realisable triples a gadget keeps h fixed iff the base port is a sink and G is transitive, entered at
its source.

## (D) Coverage (FINITE-EXACT, pre-registered HYP-9468 P1, reproduced by the audit)

Every odd N in [1, 2^32] except 7 and 21 is the Hamiltonian path count of one of:
- a one-arc join of two one-apex tournaments, H = (1 + u)(1 + v) + g(u)g(v) - e(u)e(v);
- a one-arc join of a one-apex tournament with a rooted tournament of order <= 7, H = g(v)y_h - e(v)y_e + (1 + v)y_R.

Here g(0) = 1, g(2v + 1) = g(v), g(2v) = g(v) + v + (v mod 2), and e(v) = [v = 0 or v even].

- The plain counts (F0) of the pre-registered run are redundant: g(2t) = F1(0, t), and the audit's run without F0
  covers [1, 2^28].
- The previous verified range was 99,999,999.
- The bilinear family alone has miss density 0.00087 in [2^31, 2^32).
- It lies inside the two-apex family, whose completeness is the OPEN conjecture of
  square_sum_spectrum_overnight_20261030.

**Addendum (2026-11-06, after close-out; recorded by opus-2026-11-06-S83).** HYP-9478 P5-P6 extend (D) to 2^35. Every odd N <= 2^35 except 7 and 21
is covered, already by the second family alone (rooted tournaments of order <= 7). Source:
`chanpak_f1f2_cover35_bucketed_20261106.out`. This is one run, not independently recomputed at 2^35.
