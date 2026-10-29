---
id: HYP-9445
title: "After THM-4845/4846/4847: (1) a strong tournament with at least two Hamiltonian cycles has torsion in the cokernel of its two- and three-run observers (FINITE-EXACT: k = 2 through order 8, k = 3 through order 7; false for the one-arc observer k = n-1); (2) which unique-Hamiltonian-cycle (Douglas) tournaments are torsion-free (1, 1, 3, 8, 16, 36 of 1, 1, 3, 8, 21, 55 for orders 3..8, k = 2); (3) the non-star kernel of the second permanental compound from order 6; (4) the mod-3 ghost formula; (5) T5 as the only boundary-killing class; (6) the irreducible blocks of the one-arc determinant; (7) consequences of the gcd Gaussian field for A000568; (8) Fermat: is the strong semicomplete test a primality criterion for Fermat numbers, and is there a tournament-native invariant that sees their compositeness"
status: >
  OPEN (mac-mini-2026-10-29-growth, 2026-10-29; pre-audit). Every evidence statement below is FINITE-EXACT, over
  all isomorphism classes from gentourng with exact Smith forms (PARI matsnf) unless marked otherwise.
session: mac-mini-2026-10-29-growth
source: 05-knowledge/results/tournament_growth_gaussian_gf_runs_fermat_20261029.md
related:
  - 01-canon/theorems/THM-4845-a000568-as-a-gaussian-average-gcd-covariance-cycle-index-exponential-and-the-meromorphic-normalized-series-with-poles-at-powers-of-four.md
  - 01-canon/theorems/THM-4846-run-observers-of-tournaments-the-nauru-orientation-the-permanental-compound-t5-stars-and-the-ghost-lemma-behind-its-four-spurious-z6-factors.md
  - 01-canon/theorems/THM-4847-fermat-primes-barriers-to-proving-exactly-five-the-function-field-analogue-and-why-tournament-counts-are-blind-while-three-state-counts-are-pepin.md
---

# HYP-9445 — run-observer torsion and Hamiltonian cycles, the star law, Fermat routes, Gaussian-GF questions

`M_k(T)` is the ordered `k`-run observer (THM-4846 §0): spanning covers of `T` by `k` directed paths with prescribed
starts and ends, which are the base path of a metagraph tiling cut into `k` runs. "Torsion-free" means every nonzero
invariant factor of `M_k(T)` equals 1.

1. **Hamiltonian cycles force torsion (conjecture).** If `T` is strong and has at least two Hamiltonian cycles, then
   `coker M_k(T)` has nontrivial torsion, for `k = 2` and `k = 3`. Consequently a unimodular `M_2` or `M_3` forces a
   unique Hamiltonian cycle.
   - **Evidence for `k = 2`.** All strong classes with at least two Hamiltonian cycles:
     3 (order 5), 27 (order 6), 332 (order 7) and 5953 (order 8), 6315 in all.
   - **Evidence for `k = 3`.** All such classes through order 7: 3, 27 and 332, 362 in all.
   - **Further evidence.**
     - It holds for `k = 4` at order 6 (27 classes).
     - It holds for `k = 3` at order 8 on the delicate classes: all 60 with exactly 2 Hamiltonian cycles and all 124
       with exactly 3 (`out/hc_small_k3_n8.out`).
     - The unique-cycle classes there split 36 torsion-free and 19 with torsion, as at `k = 2`.
   - **False for the one-arc observer `k = n − 1`.** One class at order 5 and two at order 6 have singular,
     torsion-free `M_{n−1}`. Each has exactly two Hamiltonian cycles.
   - **Beyond the census.** Random labelled tournaments and low-cycle tournaments (2 or 3 Hamiltonian cycles,
     found by a random walk of arc flips) were tested with `k = 2`. There were 116 + 584 strong samples at order 9
     (360 of them low-cycle) and 344 at order 10 (200 low-cycle). Every one has torsion (`out/random_torsion.out`;
     sampling evidence only).
   - **The sign sector already carries it.** The two-run sign sector `S_2` (signed covers; THM-4846 (F), a Jacobi
     ratio of cycle-cover minors) has torsion for every strong class with at least two Hamiltonian cycles: 3, 27,
     332 and 5953 classes at orders 5 to 8.
     - For unique-cycle classes `S_2` is torsion-free through order 7 (32 of 32).
     - At order 8 it is torsion-free for 54 of 55. The exception is `1010111111111111111101111101`, with Smith form
       `1^21 2^4 0^3`.
     - The plus sector `P_2` (unordered covers) has torsion for both kinds of class (`p2_vs_hc.py`,
       `out/s2_vs_hc_n8.out`).
     - So the sharper conjecture is that **`S_2(T)` has torsion whenever `T` has at least two Hamiltonian cycles.**
     - A natural route is a sign-reversing involution on the cycle covers in the Jacobi expansion, with the two
       Hamiltonian cycles as the non-cancelling pair.
   - **Not a mod-2 statement.** The torsion is not always 2-primary. Among strong classes with at least two
     Hamiltonian cycles, 8 at order 7 and 225 at order 8 have only odd invariant factors above 1. So a proof cannot
     work over `F_2` alone.
   - **A first test of any proof.** Explain why the coarser one-arc observer, an orientation of the star graph
     (THM-4846 (A)), escapes the law.
2. **Which unique-cycle tournaments are torsion-free?**
   - **Data at `k = 2`.** The unique-Hamiltonian-cycle classes number `1, 1, 3, 8, 21, 55` (orders 3 to 8; Douglas
     1970, the Fibonacci numbers `F(2n − 6)` from order 4, Garey 1972). Their torsion-free members number
     `1, 1, 3, 8, 16, 36`.
   - **OEIS.** A search for `1,3,8,16,36` gave one unrelated hit on 2026-10-29.
   - **Data at `k = 3`.** The same counts hold through order 7.
   - **Where it breaks.** The first torsion among unique-cycle classes appears at order 7, in 5 classes, with
     invariant factors such as `3·12³` and `2·6²`.
   - **Singular or unimodular.** For `k = 2` at orders 7 and 8, a unique-cycle class carries torsion exactly when
     `M_2` is singular, and every nonsingular one is unimodular. Through order 8, `det M_2 ∈ {0, ±1}` for every
     unique-cycle tournament. At order 6 the one singular class is torsion-free.
   - **The nonsingular counts.** These are `1, 1, 3, 7, 16, 36` for orders 3 to 8, which fits `(n + 1) · 2^(n−6)` for
     `5 ≤ n ≤ 8`. Several OEIS sequences share these terms (A045891 among them).
   - **PRE-REGISTERED PREDICTION** (written 2026-10-29 before the order-9 exhaustive run `stream_m2.py 9` reported):
     - at order 9, the 144 unique-cycle classes split into exactly 80 with unimodular `M_2` and 64 with torsion;
     - every class with at least two Hamiltonian cycles has torsion.

     The run records the torsion-free count, which equals the unimodular count unless a singular torsion-free class
     occurs, as at order 6.
   - **Question.** Express the criterion through Douglas's structural decomposition of these tournaments.
3. **The non-star kernel of the second permanental compound.**
   - **Background.** The two-arc unordered observer is `per A[β, α]` on pairs (THM-4846 (B)). Through order 5 its
     kernel is always spanned by vertex stars, vectors supported on the pairs through one vertex.
   - **Failures.** This fails for 12 of 56 classes at order 6, 168 of 456 at order 7 and 2808 of 6880 at order 8.
   - **Question.** Describe the non-star kernel vectors. Candidates are signed 4-cycles of pairs, or relations
     coming from 2-element modules.
4. **The mod-3 ghost formula.**
   - THM-4846 (D) proves `spur_2 = rank(B̄ : ker P̄ → coker P̄)`.
   - **Question.** Give the analogous formula for `spur_3` through the uniserial mod-3 structure of `Z_(3)[S_3]`. Here
     the trivial and sign characters are glued to the invariant flag of the standard sector; see THREE_RUN_GLUE (R3).
   - **Data.** The data for orders up to 7 are in `census_n*_k3.out`.
5. **Is T5 the only boundary killer?**
   - **The property.** T5's plus sector annihilates every lifted `(k+1)`-set boundary for `k = 3`.
   - **Data.** No other class of order 5, 6 or 7 does this for `k = 3`.
   - **Conjecture.** For every `n` and every `2 ≤ k ≤ n − 2`, T5 (with `k = 3`) is the only boundary-killing pair.
   - **A constraint from the star criterion.** For `k = n − 2` and a vertex `w` with `|N⁺(w)| ≥ 3`, a star kernel
     vector needs its coefficient vector in the left kernel of `A[V − w, V]`.
6. **The irreducible blocks of the one-arc observer.**
   - **Known.** `det M_{n−1}(T) = Π_λ det(M_λ)^(f^λ)`, with `M_(n−1) = Aᵀ` and the sign block
     `(−1)^(|x−y|−1) A[x][y]`. The order-6 determinants include `7^10`, `3^30`, `2^42` and `2^40 13^10`.
   - **Question.** Find the other blocks explicitly, for instance as matrices on pairs or flags.
7. **The gcd Gaussian field and A000568.** By THM-4845, `Σ a_n x^n = E[exp(Σ W_p x^p/p)]`.
   - **New results?** Does the field give new congruences, or uniform asymptotics, for A000568 and its family? The
     family includes the self-converse counts A002785 and the k-valued counts of THM-4720. The oriented graphs
     A001174 already fit: the base-3 field reproduces them through order 9 (THM-4845 (B)).
   - **A closed form?** Is there a closed form for the core GF `E[exp(v Σ_{p≥3} W_p x^p/p)]` beyond the average,
     for example through product identities for the divisor field?
8. **Fermat.**
   - **(a) A primality criterion?** At a Fermat prime, `s_3(F_n) ≡ −1 (mod F_n)` (THM-4847 (C3)). Is the converse
     true for Fermat numbers? `F_5` fails at both factors. The next tests need `s_3(r) mod p` at the known factors of
     `F_6, …, F_11`, with `r` up to about `10^7`. That needs EGF inversion mod `p`.
   - **(b) A tournament-native invariant?** Is there one that detects compositeness of `F_n` uniformly?
     - Counting modulo `F_n` is blind (THM-4847 (C1)).
     - The rigidity characterization (THM-4644) applies only at primes.
     - Any such invariant would be a new route to "exactly five Fermat primes".
     - Its absence is consistent with all the barriers of THM-4847 (A).
