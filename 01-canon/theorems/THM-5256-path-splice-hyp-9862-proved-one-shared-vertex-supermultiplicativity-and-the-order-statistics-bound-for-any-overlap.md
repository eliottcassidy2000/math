---
id: THM-5256
title: "The path splice. For induced subtournaments A, B of a tournament U with W = A n B and h_A(sigma), h_B(sigma) the numbers of Hamiltonian paths of A, B inducing the linear order sigma on W: the number of Hamiltonian paths of U restricting to Hamiltonian paths of both A and B is at least sum over sigma of h_A(sigma) h_B(sigma); hence H(U) >= sum_sigma h_A(sigma) h_B(sigma), and for |W| <= 1, H(U) >= H(A) H(B) (HYP-9862 PROVED, every order). Corollaries: THM-5312 (A) (pq >= p + q + 7); for two shared vertices u, v, H(U) >= h_A(uv) h_B(uv) + h_A(vu) h_B(vu), tight (5) at two triangles on an arc; the excess H(U) - sum is the number of extra shuffles plus the number of Hamiltonian paths of U not restricting to both; in the odd-cycle formula, H(U) = H(A)H(B-v) + H(A-v)H(B) - H(A-v)H(B-v) + Mix(U), so HYP-9862 is the inequality Mix(U) >= (H(A) - H(A-v))(H(B) - H(B-v)) on mixed packings."
status: >
  PROVED (klein-2026-11-16-splice; elementary — the merge lemma for two vertex-disjoint paths, applied segment by segment
  between the shared vertices; the image is decoded by restriction; not independently audited). FINITE-EXACT: all
  tournaments of orders 4-7, all subset pairs with |W| <= 2 (4,143,216 instances at order 7), 0 violations, the shuffle
  count never below the bound (scripts/shuffle_check.py); the join session's order-9 census of HYP-9862 (218,925,648
  instances) is now a corollary. Pre-registered P-1, P-2 confirmed (packet PREREGISTERED.md).
  CITED: THM-002 (odd-cycle formula), THM-5312, HYP-9862 (klein-2026-11-16-join), THM-4114 (Mobius positivity).
session: klein-2026-11-16-splice
date: 2026-11-16
source: 05-knowledge/results/path_splice_rows_continuation_20261116/SPLICE.md
scripts:
  - 05-knowledge/results/path_splice_rows_continuation_20261116/scripts/shuffle_check.py (+ receipts/shuffle_check_4..7.out)
related:
  - 05-knowledge/hypotheses/HYP-9862-one-vertex-supermultiplicativity-induced-subtournaments-sharing-at-most-one-vertex-have-h-of-the-union-at-least-the-product.md
  - 01-canon/theorems/THM-5312-gluing-targets-below-125-one-shared-vertex-forces-p-plus-q-plus-7-corrected-overlap-bound-and-the-conflict-graphs-of-7-and-21.md
  - 01-canon/theorems/THM-002-ocf.md
  - 05-knowledge/results/chanpak_gluing_targets_one_vertex_20261116/ONE_VERTEX.md
---
# THM-5256 — the path splice

**Merge lemma.** Two vertex-disjoint directed paths `X`, `Y` of a tournament have a common shuffle that is a directed
Hamiltonian path of the tournament on `V(X) ∪ V(Y)`. (Induction on `|X|`: insert the last vertex of `X` into the
suffix of `Y` that follows its predecessor, at the last position whose predecessor beats it.)

**Proof of the theorem.** For Hamiltonian paths `P_A = α_0 w_1 α_1 … w_r α_r` of `A` and `P_B = β_0 w_1 β_1 … w_r β_r`
of `B` inducing the same order `(w_1, …, w_r)` on `W`, set `Z = μ_0 w_1 μ_1 … w_r μ_r` with `μ_i` a merge of `α_i` and
`β_i`. The first and last vertices of each nonempty `μ_i` are first and last vertices of `α_i` or `β_i`, so `w_i → μ_i → w_(i+1)`
are arcs of `P_A` or `P_B`; empty `μ_i` means `w_i → w_(i+1)` is an arc of `P_A`. Thus `Z` is a Hamiltonian path of `U`
with `Z|_A = P_A`, `Z|_B = P_B`, and `(P_A, P_B) ↦ Z` is injective. ∎

**Corollaries and the odd-cycle form.** Packet SPLICE.md §3. The identity `H(U) = H(A)H(B − v) + H(A − v)H(B) − H(A − v)H(B − v) + Mix(U)`
follows from THM-002 by sorting the packings of `U` according to whether some cycle meets both `A − v` and `B − v`.

**What is new.** The theorem, its proof, the two-vertex order-statistics bound, the excess formula and the identity.
Known and used: the Hamiltonian-path insertion argument (Rédei), THM-002.
