---
id: HYP-9423
title: "After THM-4815-4817: (1) every 2-connected non-transitive tournament has an L-descending native move (a Hamiltonian path joining a pair of D_L); (2) the black distance to the regular class at odd orders is R except at self-converse colour obstructions, finitely many per order; (3) the even-order regular target; (4) cut-or-colour for every endpoint potential; (5) non-rigid prime tournaments are asymptotically (1/12) 2^(C(n-6,2)+2(n-6)+5)/(n-6)!; (6) pre-registered order-14 values: sum over non-rigid primes of (1 - 1/|Aut|) = 3,574,801,896/7, primes with an automorphism of order 7 exist, and about 7.66 * 10^8 non-rigid primes"
status: >
  OPEN (mac-mini-2026-10-27-dipole). (1) FINITE-EXACT through order 10 (all 9,733,056 classes: every strong L-local
  minimum, 0, 1, 3, 11, 49, 242, 1,302 at orders 4-10, has a cut vertex) and random to order 16; the stronger
  'min-score to max-score Hamiltonian path' fails for 2-connected tournaments (1, 26, 590 classes at orders 7, 8, 9),
  so the full descending set is needed. (2) census to order 8 (colour obstructions 0, 2, 1 at orders 5-7 for the
  regular target; order 5 is the pure-blue R5); THM-4815 B6 proves they are self-converse. (3) data: 306 failures of
  ceil(R'/2) at order 8. (4) holds for lambda, delta0, L, R through order 10. (5) heuristic plus census at orders
  10-13. (6) prediction only: no order-14 census exists.
session: mac-mini-2026-10-27-dipole
source: 05-knowledge/results/dipole_principle_and_prime_cycle_index_tower_20261027.md
scripts:
  - 04-computation/experiments/dipole_exactness_20261027/endpoint_checks.c (endpoint_checks_n4_9.out, endpoint_checks_n10.out)
  - 04-computation/experiments/dipole_exactness_20261027/endpoint_random.c (endpoint_random_n11_16.out)
  - 04-computation/experiments/dipole_exactness_20261027/dipole_exactness.cpp (exact_n*.out)
  - 04-computation/experiments/prime_cycle_index_tower_20261027/tower.py, nonrigid_primes.c
related:
  - 01-canon/theorems/THM-4815-dipole-principle-native-move-distance-bounds-endpoint-relation-cut-or-colour-exactness-and-the-p7-parity-correction.md
  - 01-canon/theorems/THM-4816-prime-cycle-index-tower-order-eleven-twelve-thirteen-types-from-non-rigid-primes-and-the-ledger.md
  - 01-canon/theorems/THM-4817-order-twelve-thirteen-fourteen-ledgers-prime-self-converse-intransitivity-and-labelled-counts.md
  - 05-knowledge/hypotheses/HYP-9421-black-descent-exceptions-tc-spectrum-and-the-native-forcing-tail.md
  - 05-knowledge/hypotheses/HYP-9422-regular-tournaments-exactly-m-black-moves-from-a-singleton-plus-strong-core-nuclear-norm-maximizers-off-the-hadamard-orders-and-the-order-11-types.md
---

# HYP-9423 — Endpoint conjectures and the prime-tower predictions

1. **2-connected tournaments descend toward the transitive class.** If `T` is not transitive and `T - v` is strong
   for every vertex `v`, then some Hamiltonian path of `T` runs from `a` to `b` with
   `L(n - 1 - s + e_a - e_b) = L(s) - 1`; equivalently, `T` has a native move lowering
   `L = (1/2) sum |s_(i) - i|`. By THM-4815 B5 the non-strong local minima are exactly the tournaments with a
   source or sink. The strong local minima found so far all have a cut vertex: 0, 1, 3, 11, 49, 242, 1,302 at
   orders 4-10, exhaustively. Random samples at orders 11-16, uniform, near-regular and with planted tied extremes,
   found none either. The natural special case, a Hamiltonian path from a minimum-score to a maximum-score vertex,
   is *not* enough: 1, 26 and 590 two-connected classes at orders 7, 8, 9 have no such path but descend through
   another pair of `D_L`. A proof would need Hamiltonian paths against the score direction, in the spirit of
   Thomassen's 1980 results on prescribed endpoints. Test next: the strong order-11 classes with a cut vertex,
   by the endpoint DP.

2. **Black moves to the regular class.** For odd `n` the black distance to the regular class is `R(T)` except at
   colour obstructions (THM-4815 B6: self-converse, every descending path palindromic), and these are rare: one at
   order 7 (scores `1, 3^5, 5`, `|Aut| = 5`). At order 5 the target `R5` itself is pure blue, so every other class is
   cut off. Test: the class-graph census at order 9 (running at filing; see the note).

3. **The almost-regular target at even orders.** At even `n` a move changes `R' = sum_v (|2 s_v - n + 1| - 1)/2` by at
   most 2, and `ceil(R'/2)` is not attained for 306 classes at order 8. Determine the exact distance to the
   almost-regular class at even orders. Find whether a refined potential, for example `R'/2` corrected by the number
   of vertices at the two middle scores, is attained except at cuts.

4. **Cut or colour for every endpoint potential.** Let a functional of the scores or of the cut function have a unit
   tariff. Conjecture: all of its unrestricted local minima have vertex connectivity at most 1, and all of its black
   local minima that are not unrestricted ones are self-converse. The second half is THM-4815 B6. The first half
   holds for `lambda` (none), `delta0` (non-strong), `L` (through order 10) and `R` at odd orders (none, B8).

5. **Asymptotics of the non-rigid primes.** The number of non-rigid prime tournaments of order `n` is
   `N_nr(n) ~ 2^(C(n-6,2) + 2(n-6) + 5) / (12 (n-6)!)`. The leading term counts the classes with an automorphism
   of type `1^(n-6) 3^2`, which the cycle lemma makes the smallest admissible type. The ratio `N_nr(n)/lead(n)` is
   0.49, 0.52, 0.575, 0.619 at orders 10-13. It should tend to 1, as the prime fraction of the invariant
   tournaments tends to 1. Consequently `N_nr(n)/A(n) ~ (n^6/12) 2^(14 - 4n)`.

6. **Pre-registered order-14 values** (from THM-4817, before any order-14 census):
   `sum over non-rigid order-14 primes of (1 - 1/|Aut|) = 3,574,801,896/7 = 510,685,985.142857...`. Some order-14
   prime has an automorphism of order 7, of type `7^2`. If all such groups have order 7, their number is
   `≡ 6 (mod 7)`, and the contributions of the 3-, 5- and 9-parts sum to an integer. The total non-rigid count lies
   in `[510,685,986, 766,028,977]`, and near the top: about 766.0 million, as nearly all have `|Aut| = 3`. The tower
   also predicts `P(14) = 27,687,212,017,612,715` and `SCP(14) = 767,626,075`. Those follow from `A000568(14)` and
   `A002785(14)` and are not independent predictions.
