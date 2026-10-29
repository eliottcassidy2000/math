---
id: HYP-9445
title: "After THM-4845/4846/4847: (1) a strong tournament with at least two Hamiltonian cycles has torsion in the cokernel of its two- and three-run observers (FINITE-EXACT: k = 2 through order 9, exhaustively, 184 304 classes; k = 3 through order 7; false for the one-arc observer k = n-1); (2) RESOLVED: which unique-Hamiltonian-cycle (Douglas) tournaments have a unimodular two-run observer (det M_2 in {0, +-1} through order 9; unimodular counts 1, 1, 3, 7, 16, 36, 81 for orders 3..9; the pre-registered 80 at order 9 FAILED); (3) the non-star kernel of the second permanental compound from order 6; (4) RESOLVED: the mod-3 ghost formula (THM-4846 (D')); (5) T5 as the only boundary-killing class; (6) the irreducible blocks of the one-arc determinant; (7) consequences of the gcd Gaussian field for A000568; (8) Fermat: is the strong semicomplete test a primality criterion for Fermat numbers, and is there a tournament-native invariant that sees their compositeness"
status: >
  OPEN in the remaining questions (mac-mini-2026-10-29-growth, 2026-10-29).
  - Item (2), unimodularity/singularity classification, is PROVED in session103_words_necklaces_20261029,
    independently audited twice. The full singular Smith-torsion classification remains open.
  - Independently audited twice on 2026-10-29, with corrections applied
    (04-computation/experiments/growth_gf_runs_fermat_20261029/audit/REPORT.md and audit2/REPORT.md, MISTAKE-695).
  - The second audit certified every verdict in (1) and (2) through order 8 without trusting Smith forms
    (determinants, ranks mod p and gcds of maximal minors).
  - The audit recomputed every count through order 8 with independent code, and the order-9 unique-cycle split.
  - It did not re-check the order-9 claim for classes with at least two Hamiltonian cycles, or the random samples.
  - Except for the explicitly marked later PROVED classification in (2), evidence below is FINITE-EXACT, over all isomorphism classes from gentourng with exact Smith
    forms (PARI matsnf), unless marked otherwise.
  (4) RESOLVED by THM-4846 (D'). In (2) a pre-registered prediction FAILED at order 9; see (2).
session: mac-mini-2026-10-29-growth
source: 05-knowledge/results/tournament_growth_gaussian_gf_runs_fermat_20261029.md
related:
  - 01-canon/theorems/THM-4845-a000568-as-a-gaussian-average-gcd-covariance-cycle-index-exponential-and-the-meromorphic-normalized-series-with-poles-at-powers-of-four.md
  - 01-canon/theorems/THM-4846-run-observers-of-tournaments-the-nauru-orientation-the-permanental-compound-t5-stars-and-the-ghost-lemma-behind-its-four-spurious-z6-factors.md
  - 01-canon/theorems/THM-4847-fermat-primes-barriers-to-proving-exactly-five-the-function-field-analogue-and-why-tournament-counts-are-blind-while-three-state-counts-are-pepin.md
---

# HYP-9445 — run-observer torsion and Hamiltonian cycles, the star law, Fermat routes, Gaussian-GF questions

**2026-10-29 continuation, PROVED and independently audited twice:**
[Session 103 all-order classification](../results/session103_words_necklaces_20261029/README.md)
resolves the unimodularity part of item (2). In the actual five-letter classical
unique-cycle grammar, M_2 is singular exactly when the word contains `be` or `ee`.
Avoidance makes every ordered-k path observer permutation-equivalent to a unit
upper-triangular matrix. The forbidden direction has an explicit antisymmetric
star kernel transported from the last forbidden pair. Equivalently, avoiding
classes have their Hamilton cycle as their only directed cycle cover.

The previously proved language series (1+x)/(1-2x-x²+x³) therefore counts
unimodular M_2 classes at every order n>=4 (x marks n-4); 2,065 at order 13 is
now a theorem consequence. The earlier [operation-closed census](../results/operation_closed_synthesis_20261029/README.md)
checked all 4,180 classes through order 12, including counts 182,409,919.
Session 103 includes two mathematical audits and 18 fresh normal/optimized replay
checks. This does not classify all singular Smith forms, remove the singular
torsion-free order-six exception, extend the HC>=2 census, or prove item (1).
The failed 80 forecast and subsequent numerology/predictions below are retained
as historical provenance, not current status.

`M_k(T)` is the ordered `k`-run observer (THM-4846 §0): spanning covers of `T` by `k` directed paths with prescribed
starts and ends, which are the base path of a metagraph tiling cut into `k` runs. "Torsion-free" means every nonzero
invariant factor of `M_k(T)` equals 1.

1. **Hamiltonian cycles force torsion (conjecture).** If `T` is strong and has at least two Hamiltonian cycles, then
   `coker M_k(T)` has nontrivial torsion, for `k = 2` and `k = 3`. Consequently a unimodular `M_2` or `M_3` forces a
   unique Hamiltonian cycle.
   - **Evidence for `k = 2`.** All strong classes with at least two Hamiltonian cycles:
     3 (order 5), 27 (order 6), 332 (order 7), 5953 (order 8) and **177 989 (order 9, exhaustive)**, 184 304 in all.
     The order-9 run (`stream_m2.py 9`, four gentourng slices, `out/stream_m2_n9_part*.out`) found no counterexample.
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
     - For unique-cycle classes `S_2` is torsion-free through order 7 (34 of 34, counting C3 and S4; the second audit
       corrected an earlier "32 of 32", which counted orders 5 to 7 only).
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
   - **RESOLVED subquestion:** the all-order unimodularity/singularity criterion and recurrence are proved in the continuation above. Torsion-freeness is different: singular torsion-free cases exist. The following records the earlier finite evidence and prediction history; it is not a complete all-order Smith classification.
   - **Data at `k = 2`.** The unique-Hamiltonian-cycle classes number `1, 1, 3, 8, 21, 55` (orders 3 to 8).
     Douglas (1970) characterised them, and Garey (1972) counted them as the Fibonacci numbers `F(2n − 6)` from
     order 4. Their torsion-free members number `1, 1, 3, 8, 16, 36`.
   - **OEIS.** A search for `1,3,8,16,36` gave one unrelated hit on 2026-10-29.
   - **Data at `k = 3`.** The same counts hold through order 7.
   - **Where it breaks.** The first torsion among unique-cycle classes appears at order 7, in 5 classes. Their
     non-unit invariant factors are `2²` (twice), `3³·6` and `4²` (twice) at `k = 2`, and include `3·12³` and `2·6²`
     at `k = 3`.
   - **Singular or unimodular.** For `k = 2` at orders 7 and 8, a unique-cycle class carries torsion exactly when
     `M_2` is singular, and every nonsingular one is unimodular. Through order 8, `det M_2 ∈ {0, ±1}` for every
     unique-cycle tournament. At order 6 the one singular class is torsion-free.
   - **The nonsingular counts.** These are `1, 1, 3, 7, 16, 36` for orders 3 to 8, which fits `(n + 1) · 2^(n−6)` for
     `5 ≤ n ≤ 8`. Several OEIS sequences share these terms (A045891 among them).
   - **PRE-REGISTERED PREDICTION** (written 2026-10-29 before the order-9 exhaustive run `stream_m2.py 9` reported,
     commit 5047208384):
     - at order 9, the 144 unique-cycle classes split into exactly 80 with unimodular `M_2` and 64 with torsion;
     - every class with at least two Hamiltonian cycles has torsion.

     The run records the torsion-free count, which equals the unimodular count unless a singular torsion-free class
     occurs, as at order 6.
   - **OUTCOME (order 9, exhaustive):**
     - *The prediction FAILED.* The 144 unique-cycle classes (Douglas's `F(12) = 144`, confirmed) split into 81
       unimodular, 0 singular torsion-free and 63 with torsion. All 63 torsion classes are singular
       (`out/n9_unique_hc_split.out`). So `(n + 1)·2^(n−6)` is REFUTED at order 9.
     - *The second clause held:* every one of the 177 989 classes with at least two Hamiltonian cycles has torsion.
     - *What survives.* For a unique-cycle tournament, `M_2` has torsion iff it is singular (orders 7, 8 and 9), and
       `det M_2 ∈ {0, ±1}` (orders 3 to 9).
     - *The corrected counts* are `1, 1, 3, 7, 16, 36, 81` for orders 3 to 9.
       - From order 4 these coincide with A033303, the expansion of `(1 + x)/(1 − 2x − x² + x³)`, whose recurrence is
         `a(n) = 2a(n−1) + a(n−2) − a(n−3)`. Its growth constant is `1 + 2cos(2π/7) ≈ 2.24698`, the largest root of
         `y³ − 2y² − y + 1`. The audit corrected an earlier statement naming `y³ − y² − 2y + 1`, whose largest root
         is 1.80194. The Garey total grows like `φ²`.
       - Both fits, `(n + 1)·2^(n−6)` and A033303, are NUMEROLOGY, exact fits with no spare data point. The first is
         already refuted.
       - **New prediction for order 10:** 182 unimodular among the 377 unique-cycle classes. The test needs a direct
         generator of unique-cycle tournaments, since the full order-10 census has 9 733 056 classes.
   - **Question.** Express the criterion through Douglas's structural decomposition of these tournaments.
3. **The non-star kernel of the second permanental compound.** (`out/nonstar_kernel_n6.out` lists the 12 order-6
   classes. Two strong classes there have no star vector at all, with kernels of dimension 2 and 3.)
   - **Background.** The two-arc unordered observer is `per A[β, α]` on pairs (THM-4846 (B)). Through order 5 its
     kernel is always spanned by vertex stars, vectors supported on the pairs through one vertex.
   - **Failures.** This fails for 12 of 56 classes at order 6, 168 of 456 at order 7 and 2808 of 6880 at order 8.
   - **Question.** Describe the non-star kernel vectors. Candidates are signed 4-cycles of pairs, or relations
     coming from 2-element modules.
4. **The mod-3 ghost formula — RESOLVED by THM-4846 (D′).**
   - In the basis adapted to the invariant flag, the glued observer reduces mod 3 to two block-lower-triangular
     pieces with diagonals `(P, S, P)` and `(S, P, S)`.
   - Hence `spur_3 = rank I + rank II − rank P̄ − rank S̄ − 2 rank [[P̄, L̄12],[0, S̄]]`. This reduces to one connecting
     rank when `S̄` is invertible.
   - Verified on every class through order 6.
   - **Still open.** A characterization of when the connecting maps are injective on the trivial-sector kernel, as
     they are for T5.
5. **Is T5 the only boundary killer?**
   - **The property.** T5's plus sector annihilates every lifted `(k+1)`-set boundary for `k = 3`.
   - **Data.** No other class of order 5, 6 or 7 does this for `k = 3`. The audit found none for
     `(n, k) = (4–8, 2)`, `(6, 4)`, `(7, 4)` or `(7, 5)` either.
   - **Conjecture.** For every `n` and every `2 ≤ k ≤ n − 2`, T5 (with `k = 3`) is the only boundary-killing pair.
   - **A constraint from the star criterion.** For `k = n − 2` and a vertex `w` with `|N⁺(w)| ≥ 3`, a star kernel
     vector needs its coefficient vector in the left kernel of `A[V − w, V]`.
6. **The irreducible blocks of the one-arc observer.**
   - **Known.** `det M_{n−1}(T) = Π_λ det(M_λ)^(f^λ)`, with `M_(n−1) = Aᵀ` and the sign block
     `(−1)^(|x−y|−1) A[x][y]`. The order-6 determinants include `7^10`, `3^30`, `2^42` and `2^40 13^10`.
   - **Question.** Find the other blocks explicitly, for instance as matrices on pairs or flags.
7. **The gcd Gaussian field and A000568.** By THM-4845, `Σ a_n x^n = E[exp(Σ W_p x^p/p)]`.
   - **New results?** Does the field give new congruences, or uniform asymptotics, for A000568 and its family? The
     family includes the k-valued counts of THM-4720. The oriented graphs A001174 already fit: the base-3 field
     reproduces them through order 9. So do the self-converse counts A002785, through the twisted form through order
     16 (THM-4845 (B)).
   - **A closed form?** Is there a closed form for the core GF `E[exp(v Σ_{p≥3} W_p x^p/p)]` beyond the average,
     for example through product identities for the divisor field?
8. **Fermat.**
   - **(a) A primality criterion?** At a Fermat prime `p ≥ 5`, `s_3(p) ≡ −1 (mod p)` (THM-4847 (C3)). Is the
     converse true for Fermat numbers?
     - A single failing prime factor decides the test.
     - `F_5`, `F_6`, `F_9` and `F_11` fail at known factors (`F_6`, `F_9`, `F_11` by the audit, using series
       inversion mod `p`).
     - For `F_7` and `F_8` the smallest known factors give residues `r ≈ 4.1·10^16` and `2.3·10^14`. That is far
       beyond direct recurrence or series inversion, so those tests need a new idea, such as quasi-periodicity at a
       second modulus.
   - **(b) A tournament-native invariant?** Is there one that detects compositeness of `F_n` uniformly?
     - Counting modulo `F_n` is blind (THM-4847 (C1)).
     - The rigidity characterization (THM-4644) applies only at primes.
     - Any such invariant would be a new route to "exactly five Fermat primes".
     - Its absence is consistent with all the barriers of THM-4847 (A).
