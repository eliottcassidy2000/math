---
id: THM-5263
title: "Digit-0 vertices and the near-join. For a vertex u of a tournament S with root response B_u (THM-4865) and digit sigma_u = B_u - H(S): (I1) sigma_u is the number of ordered two-path covers (P_1, P_2) with P_2 starting at u and u beating the last vertex of P_1 (backward junctions); (I2) if u = z has out-degree 1 with out-neighbour t_0 then sigma_z = e_{t_0}(S - z), so z is a digit-0 vertex iff no Hamiltonian path of S - z ends at t_0; (I3) the near-join T' = T -> M_k -> z with z -> t_0 is the rooted append of z at t_0 to T -> M_k, so H(T') = h(T) + 2^k B_{t_0}(T) and z is a digit-0 vertex of T' (strong when T is strong). Consequences: every prime near-join value has a strong realiser with a digit-0 vertex; the primes 11, 13, 47, 67, whose strong realisers (all of order <= 9) have no digit-0 vertex, are never near-join values from any base — HYP-9453's seven exceptions {1, 7, 11, 13, 21, 47, 67} are permanent; a padded sink is always a digit-0 vertex (MISTAKE-776). Every digit-0 vertex of a strong tournament has out-degree 1 (FINITE-EXACT here at orders 3-8, 50,782 rooted strong tournaments; PROVED at every order by THM-5264), so the counts with a strong digit-0 vertex are the dead-end appends h(B') + B_{t_0}(B') with e_{t_0}(B') = 0, which from bases of order <= 7 realise every odd number below 125 except exactly the seven."
status: >
  PROVED (klein-2026-11-16-nearjoin; (I1)-(I3) elementary from THM-4865, THM-4862 (a) and THM-5260 (B), each verified
  exhaustively on small orders; (B) from (I3) and the exhaustive strong census of orders <= 9 (strongB_*.out of the
  previous packet; the strong floor f(10) = 125 makes it exact below 125); not independently audited). FINITE-EXACT: the
  out-degree-1 statement (orders 3-8) and the dead-end-append coverage below 125 (bases <= 7). CITED: THM-4865, THM-4862,
  THM-5260, THM-5126, spectrum_analogies_nearjoin_20261102 (the near-join and its check B(T', z) = H(T')), HYP-9453.
session: klein-2026-11-16-nearjoin
date: 2026-11-16
source: 05-knowledge/results/nearjoin_digit0_readers_20261116/NEARJOIN.md
scripts:
  - 05-knowledge/results/nearjoin_digit0_readers_20261116/scripts/nearjoin_check.py (+ receipts/nearjoin_check.out)
  - 05-knowledge/results/nearjoin_digit0_readers_20261116/scripts/digit0_structure.py (+ receipts/digit0_structure.out)
related:
  - 01-canon/theorems/THM-5264-the-out-degree-1-law-for-digit-0-vertices-is-proved-at-every-order-two-path-lemma-unique-digit-0-vertex-in-a-strong-tournament-and-the-converse-digit-0-implies-a-near-join-count-below-22228.md
  - 01-canon/theorems/THM-4865-the-sole-exit-tower-is-tribonacci.md
  - 05-knowledge/hypotheses/HYP-9453-near-join-family.md
  - 05-knowledge/results/spectrum_analogies_nearjoin_20261102/README.md
  - 01-canon/theorems/THM-5260-forcing-the-midpoint-but-not-the-plus-2-child-history-blind-gadgets-are-ordinal-sums-the-ordinal-sum-laws-for-root-response-end-count-and-charge-and-the-tower-criterion-for-3h-plus-2.md
---
# THM-5263 — digit-0 vertices, dead-end appends, and the permanent near-join exceptions

**(I1)** `B_u = Σ_A h(S[A]) s_u(S[V∖A])` counts pairs (Hamiltonian path of `S[A]`, Hamiltonian path of `S[V∖A]` from
`u`), `A = ∅` allowed. Cutting a Hamiltonian path of `S` just before `u` gives such a pair with `P_1` empty or
`last(P_1) → u`, bijectively; the remaining pairs have `u → last(P_1)`. ∎ **(I2)** Paths from `z` are `z` alone or
`z → t_0 → …`, so `B_z = h(S − z) + B_{t_0}(S − z)`; THM-4862 (a) gives `H(S) = h(S − z) − e_{t_0}(S − z) + B_{t_0}(S − z)`.
∎ **(I3)** In `T → M_k`, every Hamiltonian path ends in `M_k`, so `e_{t_0} = 0`, and `B_{t_0}(T → M_k) = 2^k B_{t_0}(T)`
by THM-5260 (B); `z` beats only `t_0`. ∎ **(B)** A near-join value `v` prime: `T′` has one nontrivial strong component,
which contains `z` (`z → t_0`, everything `→ z`) and has count `v`, so it is a strong realiser of `v` with a digit-0
vertex; the census has none for 11, 13, 47, 67. ∎
