---
id: THM-5269
title: "The q-slot cover count is the odd-cycle packing polynomial at (2q, q). (A) For a tournament T let P_q(T) be the number of ordered q-tuples of disjoint, possibly empty vertex sets covering V(T), weighted by the product of their Hamiltonian counts (P_1 = h, P_2 the two-cover count of THM-5264/THM-5265, P_3 the three-cover count), and G_T(x, y) the generating polynomial of the spanning packings of T by directed odd cycles of length >= 3 (x) and singletons (y), so that h = G_T(2, 1) is THM-002. Then P_q(T) = G_T(2q, q) for every q >= 1. Consequences: P_q is a polynomial in q of degree n whose coefficients a_m = sum_{c+s=m} 2^c N_{c,s} vanish unless m = n (mod 2); P_2 is even; the signed two-cover sum over W of (-1)^|W| h(W) h(W^c) is 0 for every nonempty T; P_q is multiplicative over ordinal sums; two tournaments of order n with equal h and P_2 can differ in P_3 only by an element of 120 Z (n = 7), 360 Z (n = 8), 120 Z (n = 9). (B) Census: among the 7,412 classes of orders <= 8 there are 978 two-cover counts, 602 shared by classes of different count, and exactly one (h, P_2) collision split by P_3, namely (345, 5640) with P_3 in {44505, 45225}; at order 9, 48,870 collision groups of which 338 split, with differences 240, 2160, 2400, 2640; P(C_3) = (3, 12, 33), P(strong-4) = (5, 32, 117), P(RT_5) = (15, 120, 525), P(P_7) = (189, 2016, 12159). (C) Flexible cores under the one-port alphabet S (sink) and N (near-sink at the root), whose state is (H, g, e_r, B_r): over the 4,378 nontrivial count-preserving flips of orders <= 8 and the 18,880 roots off the flipped arc, the flip is invisible at 1,040, stores a delayed choice in g (read only by the second N) at 216, a first-wrap choice (Delta e = ±2) at 6, and changes B at every endpoint of the flipped arc; permanent +2 pairs of kind K1 (same (g, e, B), counts H, H + 2) exist at 201 even centres below 1,000 and of kind K2 (same (H, g, B), end-counts e, e + 2; permanent after one N) at 95, together at 250 of 499; the number of rooted states a fibre offers has Spearman correlation 0.673 with its number of nontrivial neutral edges; iterated N at the newest vertex is the tribonacci tower of THM-4865."
status: >
  PROVED (klein-2026-11-16-compose): (A) by expanding each slot by THM-002 and assigning the components of a spanning
  packing to slots (replayed on all 530 classes of orders 3-7 with the packings enumerated directly, q = 1, 2, 3, and the
  alternating identity, 0 failures); the lattices by solving the two linear constraints on the coefficient differences.
  FINITE-EXACT: (B) (twocover.py on profiles8.npz; covers9.c + collisions9.py on all order-9 classes), (C) (compose.py
  with nauty labelg). Not independently audited. CITED: THM-002, THM-4865, THM-5260 (B - e = 2E), THM-5264, THM-5265,
  THM-5267, THM-5268, and the concurrent one_port_neutral_continuations_20261117 (the alphabet S, N and its transfers)
  and two_cover_feedback_20261117 (the dyadic tile bank, the frozen-fibre selection theorem, the pair (345, 5640)).
session: klein-2026-11-16-compose
date: 2026-11-16
source: 05-knowledge/results/compose_alphabet_twocover_20261116/TWOCOVER.md
scripts:
  - 05-knowledge/results/compose_alphabet_twocover_20261116/scripts/packing.py (+ receipts/packing.out)
  - 05-knowledge/results/compose_alphabet_twocover_20261116/scripts/twocover.py (+ receipts/twocover.out)
  - 05-knowledge/results/compose_alphabet_twocover_20261116/scripts/covers9.c
  - 05-knowledge/results/compose_alphabet_twocover_20261116/scripts/collisions9.py (+ receipts/collisions9.out)
  - 05-knowledge/results/compose_alphabet_twocover_20261116/scripts/compose.py (+ receipts/compose.out)
related:
  - 01-canon/theorems/THM-002-ocf.md
  - 01-canon/theorems/THM-5268-padding-is-transparent-to-balanced-arcs-the-twenty-one-fibres-below-125-with-no-count-preserving-flip-at-any-order-and-the-holes-as-k-1-near-joins-with-borrowed-responses.md
  - 01-canon/theorems/THM-5267-seams-at-a-prescribed-centre-by-one-port-contact-the-arc-flip-lemma-deletion-monotonicity-the-order-10-extension-graph-and-the-heaviest-parent-tree-of-first-creators.md
  - 01-canon/theorems/THM-4865-the-sole-exit-tower-is-tribonacci.md
  - 05-knowledge/results/compose_alphabet_twocover_20261116/COMPOSE.md
  - 05-knowledge/results/two_cover_feedback_20261117/README.md
  - 05-knowledge/results/one_port_neutral_continuations_20261117/README.md
---
# THM-5269 — the cover counts are the packing polynomial

**(A)** `P_q(T) = Σ_{(W_1, …, W_q)} Π_i h(W_i)`; expand each `h(W_i) = Σ_{packings of W_i} 2^{cycles}` (THM-002). A choice
of a cover and of a packing of each slot is the same as a spanning packing of `T` with each component assigned to a
slot: a cycle contributes `2·q`, a singleton `1·q`. Hence `P_q = Σ_π (2q)^{c(π)} q^{s(π)} = G_T(2q, q)`. ∎ The coefficient
of `q^m` is `Σ_{c+s=m} 2^c N_{c,s}`, and `m = n − Σ (ℓ_i − 1) ≡ n (mod 2)`. For the signed two-cover sum assign each
component a sign: a slot-1 component of odd size contributes `−1`, a slot-2 component `+1`, so each component's factor is
`0`. For the lattices write `Δa` for the difference of coefficient vectors of two tournaments with equal `P_1, P_2`:
`Σ Δa_m = 0` and `Σ Δa_m 2^m = 0`; at order 8 the solutions on `(a_6, a_4, a_2)` are `t(1, −5, 4)`, giving `ΔP_3 = 360t`;
at order 7 `t(1, −5, 4)` on `(a_5, a_3, a_1)` gives `120t`; at order 9 `Δa_3 = −21Δa_7 − 5Δa_5`, `Δa_1 = 20Δa_7 + 4Δa_5`,
giving `ΔP_3 = 1680Δa_7 + 120Δa_5`. ∎

**(B), (C)** are the censuses of TWOCOVER.md §2 and COMPOSE.md §1–3. The selection rule of the concurrent packets
(freeze `P_2`, vary `h`; hold both, store `P_3`) is the statement that `G(4, 2)` and `G(2, 1)` are two evaluations of one
polynomial whose remaining freedom is the lattice above. Chan–Pak OPEN.
