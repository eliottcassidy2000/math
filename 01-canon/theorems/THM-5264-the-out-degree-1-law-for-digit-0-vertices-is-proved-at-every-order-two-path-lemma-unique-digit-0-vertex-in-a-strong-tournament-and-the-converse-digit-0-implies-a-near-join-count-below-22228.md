---
id: THM-5264
title: "The out-degree-1 law for digit-0 vertices, proved at every order. (L) Two-path lemma: in a strong tournament C with distinct vertices x, x', the vertex set is partitioned by two directed paths, one starting at one of x, x' and the other ending at the other (the ordered form fails). (A) A vertex u with out-degree >= 2 has digit sigma_u = B_u - H >= 1; hence sigma_u = 0 iff u is a sink, or u has out-degree 1 and its out-neighbour t_0 lies outside the terminal strong component of S - u (e_{t_0}(S - u) = 0). In a strong tournament every digit-0 vertex has out-degree exactly 1, its port lies in the initial component of S - z, and H(S) = h(S - z) + B_{t_0}(S - z): every strong tournament with a digit-0 vertex is a dead-end append. (U) A strong tournament of order >= 4 has at most one digit-0 vertex (C_3 has three). (D) Dichotomy for strong S with digit-0 vertex z, S - z = A -> D with D the terminal component: H(S) = h(A)h(D) + B_{t_0}(A) P_D(2), P_D(2) = sum over W of h(D[W]) h(D[W^c]); if D is a single vertex S is literally a near-join T -> M_k -> z with k >= 1, otherwise D is strong (a composite right source). (C) Converse: digit 0 implies a near-join count (k >= 1) for every strong tournament with H <= 22,227, and the converse at all orders is exactly {strong digit-0 counts} = {near-join values}; all 62 strong-D cases of orders <= 8 are near-join values from bases <= 7. (E) {counts} . {near-join values} is contained in {near-join values} (left ordinal multiplication), and (2^k + 1) v is a near-join value for every count v (sink trick). (F) Sterile strong counts (no digit-0 vertex) of orders <= 9 are every odd number in [11, 199] except 17, 21, 27; below 125 the counts all of whose strong realisers are fertile are exactly 3, 5, 9, 17, 27 (permanent by Busch's floors), and those all of whose strong realisers are sterile are 11, 13, 47, 67."
status: >
  PROVED (klein-2026-11-16-outdeg): (L), (A), (U), (D), (E) by elementary arguments on the condensation of S - u (Rédei,
  Camion, THM-5263 (I1)-(I3), THM-5260 (B)); each replayed exhaustively (characterisation on 58,650 vertices of orders
  3-8, 0 mismatches; lemma on 352,462 pairs; uniqueness to order 9; closure 400 cases). (C): PROVED for H <= 22,227 via
  THM-5263 (B) and HYP-9453's census (bases <= 8); the structural statement above 22,227 in the strong-D case is OPEN.
  (F) FINITE-EXACT (orders <= 9; the permanence claims use Busch's strong floors 3, 5, 9, 15, 25, 45, 75, 125). Not
  independently audited. CITED: THM-4865, THM-5260, THM-5263, HYP-9453, THM-4094 (reduction to primes), MISTAKE-776.
session: klein-2026-11-16-outdeg
date: 2026-11-16
source: 05-knowledge/results/outdeg_converse_descent_fifteen_20261116/OUTDEG.md
scripts:
  - 05-knowledge/results/outdeg_converse_descent_fifteen_20261116/scripts/twopath.py (+ receipts/twopath.out)
  - 05-knowledge/results/outdeg_converse_descent_fifteen_20261116/scripts/converse_descent.py (+ receipts/converse_descent.out)
  - 05-knowledge/results/outdeg_converse_descent_fifteen_20261116/scripts/converse_census.py (+ receipts/converse_census.out)
  - 05-knowledge/results/outdeg_converse_descent_fifteen_20261116/scripts/closure_check.py (+ receipts/closure_check.out)
  - 05-knowledge/results/outdeg_converse_descent_fifteen_20261116/scripts/sterile9.c (+ receipts/sterile9.out)
related:
  - 01-canon/theorems/THM-5263-digit-0-vertices-backward-junctions-sole-exits-and-the-near-join-the-primes-11-13-47-67-are-permanent-near-join-exceptions.md
  - 01-canon/theorems/THM-4865-the-sole-exit-tower-is-tribonacci.md
  - 01-canon/theorems/THM-5260-forcing-the-midpoint-but-not-the-plus-2-child-history-blind-gadgets-are-ordinal-sums-the-ordinal-sum-laws-for-root-response-end-count-and-charge-and-the-tower-criterion-for-3h-plus-2.md
  - 05-knowledge/hypotheses/HYP-9453-the-near-join-sweep-toward-chan-pak-completeness-density-of-b-sets-two-strong-realisers-for-every-attained-odd-count-at-least-15.md
  - 05-knowledge/results/outdeg_converse_descent_fifteen_20261116/CONVERSE.md
  - 05-knowledge/results/outdeg_converse_descent_fifteen_20261116/DESCENT.md
---
# THM-5264 — the out-degree-1 law, the unique digit-0 vertex, and the converse

Notation as in THM-4865 and THM-5263: `B_u(S) = Σ_{A ⊆ V∖u} h(S[A]) s_u(S[V∖A])`, `σ_u = B_u − H(S)`; (I1) `σ_u` counts
the two-path covers `(P_1, P_2)`, `P_2` from `u`, `P_1` nonempty, with the backward junction `u → last(P_1)`; (I2) for
`d⁺(u) = 1` with out-neighbour `t_0`, `σ_u = e_{t_0}(S − u)`.

**(L)** Let `E_1 → … → E_q` be the strong components of `C − x′`. If `x ∈ E_1`: `P_a` = Hamiltonian path of `E_1` from `x`
then `Ham(E_2), …, Ham(E_q)`, `P_b = (x′)`. If `x ∈ E_b`, `b ≥ 2`: strongness gives `x′ → e ∈ E_1`;
`P_a = (x′, Ham(E_1) from e, Ham(E_2), …, Ham(E_{b−1}), Ham(E_{b+1}), …, Ham(E_q))`, `P_b` = Hamiltonian path of `E_b` ending at
`x`. ∎ (The ordered form fails for 726 of 352,462 pairs at orders 3–8.)

**(A)** Let `C_1 → … → C_m` be the components of `S − u`, `j` the largest index met by `N⁺(u)`, `x′ ∈ N⁺(u) ∩ C_j`. If some
out-neighbour `x` lies in `C_i`, `i < j`: `P_2 = (u, x, Ham(C_{j+1}), …, Ham(C_m))`,
`P_1 = (Ham(C_1), …, Ham(C_i − x), …, Ham(C_j ending at x′))`. Otherwise all out-neighbours lie in `C_j`, strong of order ≥ 3;
(L) on two of them gives `Q_a` (from one) and `Q_b` (ending at the other, `y`); `P_2 = (u, Q_a, Ham(C_{j+1}), …)`,
`P_1 = (Ham(C_1), …, Ham(C_{j−1}), Q_b)`, `u → y`. Either way (I1) gives `σ_u ≥ 1`. The characterisation follows with (I2);
in a strong `S` there is no sink, and `t_0 ∈ C_1` because every vertex other than `t_0` beats `z`. ∎

**(U)** Digit-0 vertices `z_1 → z_2` with `z_2 → t_2`: `C_m(S − z_1) ∌ z_2` forces `C_m(S − z_1) = {t_2}` (any other member
would beat `z_2`), so `t_2` is a sink of `S − z_1`, `t_2 → z_1` by strongness, and `z_1 z_2 t_2` is a 3-cycle of out-degree-1
vertices beaten by every other vertex: `S = R → C_3`, not strong unless `R = ∅`. ∎

**(D)** `H(S) = h(S − z) + B_{t_0}(S − z)` by (I2) with `e_{t_0} = 0`; `S − z = A → D` with `t_0 ∈ A`, and THM-5260 (B) gives
`h(A)h(D)` and `B_{t_0}(A) P_D(2)`; for `D = M_1`, `P = 2`, and the trailing singletons form the transitive suffix `M_k`. ∎
Replay: 341 pairs at orders 4–8, 0 failures; 279 with a sink in the base, 62 with a strong `D`.

**(C)** By (A) the strong digit-0 counts are dead-end-append counts; by THM-5263 (I3) every near-join value is one. HYP-9453's
census makes every odd `≤ 22,227` other than `1, 7, 11, 13, 21, 47, 67` a near-join value, and THM-5263 (B) shows the seven
have no strong digit-0 realiser; `H` is odd. ∎ The 62 strong-`D` cases have least near-join base orders 3–7
(`converse_census.out`). OPEN: a near-join rewriting of `h(A)h(D) + B_{t_0}(A)P_D(2)` in general.

**(E)** `R → (T → M_k → z)` is the near-join with base `R → T`, and `B_{t_0}(R → T) = h(R) B_{t_0}(T)` (THM-5260 (B)); the padded
sink `s` of `T → s` has `B_s = H(T)` (MISTAKE-776), so `(T → s) → M_k → z` with `z → s` has count `(2^k + 1) H(T)`. ∎

**(F)** `sterile9.c` over all strong classes of orders 3–9; Busch's floors bound the orders at which a count has strong
realisers. The three kinds of counts and the descent system are in DESCENT.md; the typed "fifteen / twenty-five" material
(RT_5 orbits, supersingular primes, the first 25 square-sum nodes) is in FIFTEEN.md and carries no theorem.
