---
id: THM-4875
title: "The old-mark digit lower bound for C2 children (every old vertex other than the mark has digit >= 4, with equality only for (TT_n, penultimate mark, sink); the sink's digit is exactly 3h(B-z) + 2R_a(B-z) - e_a(B-z)); explicit one-arc gadget transfers with the Lorentzian adjoint M(G^op, z, u) = J M(G, u, z)^T J; a descent lemma against unimodular cone-covering induction; and one-apex one-arc joins realise every odd number up to 2^32 except 7 and 21 (FINITE-EXACT; Chan-Pak OPEN)"
status: >
  PROVED: (A) at every order. (A1) is READY (7) dualised, with the walk proof written out. (A2) and (A4) are new
  (exact sink digit; sharp minimum 4 and its equality case). (A5) is the S77 audit's mark formula, re-derived. (B)'s
  explicit coefficients and the general adjoint identity; the mechanism is KNOWN (session29 ADDENDUM section 3, 18 spine
  tiles), and so is the tile case of the adjoint (session29 section 5). (C) elementary.
  FINITE-EXACT: (A) on all 7,412 classes of order 1-8 and all 58,652 marks, with order 8 pre-registered as HYP-9468 P3;
  non-sink old digits >= 9 at orders 2-8 (proved bound 6); (B) formula against direct counts for every labelled gadget
  of order <= 5 (133,525 checks), adjoint for all 1,206,353 marked gadgets of order <= 6; (D) pre-registered as
  HYP-9468 P1.
  INHERITED: the one-arc join law (session29 ADDENDUM sections 2 and 5; THM-4866 (B)); the one-apex state recurrence
  (session25; Session191 one_apex/PROOF.md); C2's native digits (RESULT eqs. (8), (10), (12)); THM-4873; the
  impossibility of 7 and 21 (THM-1370).
  NOT CLAIMED: Chan-Pak completeness, which is OPEN.
  AUDIT: pending.
session: opus-2026-11-04-S78
source: 05-knowledge/results/old_mark_bound_lorentzian_gadgets_one_apex_join_coverage_20261104.md
related:
  - THM-4873
  - THM-4866
  - THM-4865
  - THM-4862
  - THM-1370
  - HYP-9468
  - HYP-9454
  - 05-knowledge/results/session29_framed_profiles_20261022/session29-profiles/ADDENDUM.md
  - 05-knowledge/results/session191_tournament_growth_20261104/one_apex/PROOF.md
  - 05-knowledge/results/two_entry_funnels_20261103/READY.md
verification:
  - 04-computation/experiments/old_mark_digit_bound_20261104.c
  - 04-computation/experiments/old_mark_digit_bound_20261104.out
  - 04-computation/experiments/gadget_adjoint_20261104.c
  - 04-computation/experiments/gadget_adjoint_20261104.out
  - 04-computation/experiments/rooted_triple_bank_20261104.c
  - 04-computation/experiments/rooted_triple_bank_order7_20261104.txt
  - 04-computation/experiments/chanpak_one_apex_join_cover_20261104.c
  - 04-computation/experiments/chanpak_one_apex_join_cover_20261104.out
  - 04-computation/experiments/join_witness_check_20261104.c
  - 04-computation/experiments/join_witness_check_20261104.out
  - 04-computation/experiments/one_apex_triple_check_20261104.py
---

# THM-4875 - the old-mark digit lower bound; Lorentzian gadgets; one-apex joins to 2^32

Notation. The root response is R_t = sum_{X subset V - t} h(T[X]) s_t(T[V - X]) (klein's B_t), the digit is
sigma_t = R_t - h, and the rooted triple is (h, e_t, R_t). C2(B, a) adds x -> a (beaten by the rest of B), then y,
beaten by all of B, with y -> x.

## (A) The old-mark digit lower bound (PROVED)

Let T1 = C2(B, a), with |B| >= 2.

- **(A1)** In every tournament, a vertex of out-degree d has sigma >= d(d - 1).
- **(A2)** If z != a is the sink of B, then sigma_z(T1) = 3h(B - z) + 2R_a(B - z) - e_a(B - z) >= 4. Equality holds iff
  B is transitive and a is its penultimate vertex.
- **(A3)** Every other old vertex w != a has sigma_w(T1) >= (d_B(w) + 2)(d_B(w) + 1) >= 6. FINITE-EXACT: >= 9 at
  orders <= 8.
- **(A4)** Hence every old vertex w != a has sigma_w(T1) >= 4, with equality iff (B, a, w) = (TT_n, penultimate, sink).
- **(A5)** The mark: sigma_a(T1) is 0 if a is the sink of B; e_w(B - a) + R_w(B - a) if a -> w is its only out-arc
  (>= 2, never 1); and >= (d + 1)d if d_B(a) = d >= 2.

**Proofs.** Source note, section 1.
- (A1): the bit-pair walk over a Hamiltonian path of T - {v, w, w'}.
- (A2) and (A5): case splits on the start of P2. In T1, x's only out-neighbour is a and y's is x.

## (B) One-arc gadget transfers and the Lorentzian adjoint (PROVED)

For a gadget G with entry u and exit z, and a base (T, v), let N = T => G with u -> v reversed. Then
(h, e, R)(N, z) = M(G, u, z)(h, e, R)(T, v), with

    M = [[h_G, -s_u, E_u], [e_z, -P_(u->z), E_(u,z)], [R_z, -J_1, J_2 + K]]

(coefficients defined in the source note, section 2), and

    M(G^op, z, u) = J M(G, u, z)^T J,   J = diag(1, -1, 1).

The adjoint follows from the two readings of T1 => G => T3 and its converse as double joins.

## (C) Descent lemma (PROVED)

Integer matrices of determinant ±1 that strictly increase a positive integral functional on a closed cone C with
integer points cannot cover C by their images. Gadget transfers strictly increase h on cones of triples with
R >= (1 + c)h. So a finite-base unimodular cone-covering induction for Chan-Pak is impossible.

## (D) Coverage (FINITE-EXACT, pre-registered HYP-9468 P1)

Every odd N in [1, 2^32] except 7 and 21 is the Hamiltonian path count of one of:
- a one-arc join of two one-apex tournaments, H = (1 + u)(1 + v) + g(u)g(v) - e(u)e(v);
- a one-arc join of a one-apex tournament with a rooted tournament of order <= 7, H = g(v)y_h - e(v)y_e + (1 + v)y_R;
- a tournament of order <= 8 or a one-apex tournament.

Here g(0) = 1, g(2v + 1) = g(v), g(2v) = g(v) + v + (v mod 2), and e(v) = [v = 0 or v even]. The previous verified
range was 99,999,999. (RESULT TO BE CONFIRMED.)
