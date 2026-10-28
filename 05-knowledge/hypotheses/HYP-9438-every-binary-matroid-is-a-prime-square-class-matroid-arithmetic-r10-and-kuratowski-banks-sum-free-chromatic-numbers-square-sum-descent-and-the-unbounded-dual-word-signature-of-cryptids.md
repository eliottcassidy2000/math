---
id: HYP-9438
title: "After THM-4840: (1) every binary matroid is the square-class matroid of a set of distinct primes (simple ones unconditionally for rank <= 6; all of them under Schinzel's hypothesis H for the quadratics d k^2 - 1); (2) the least prime banks whose square-class matroid is R10, M*(K5) or M*(K3,3) exist and have largest prime below 10^6; (3) the chromatic number of the line hypergraph of PG(r-1,2), i.e. the least number of sum-free sets partitioning F2^r minus 0, is 2, 3, 3, 4 for r = 2..5 and grows like log_2 r + O(1); (4) square-sum Hamiltonian paths on [n] admit a Berggren-type descent (unique inverse, decreasing height) on their alternating square sequences; (5) the minimal Nerode quotient of guarded square-sigma multiplication is a function of the matroid and the conflict graph alone; (6) a Pi^0_1 cryptid's halting certificate is a dual word of unbounded support, and this separates the cryptids from the machines decided by regular or bounded-Parikh deciders; (7) the duad/syntheme self-duality of the fifteen-prime bank is induced by an arithmetic involution"
status: >
  OPEN (klein-2026-10-28). (1) FINITE-EXACT for simple matroids of rank <= 6 (THM-4840 (C)); the general statement needs
  one prime value of each quadratic d k^2 - 1 (Bunyakovsky/Schinzel-type) and, for parallel elements, infinitely many.
  (2) untested; a search over the squarefree parts of p + 1 for p <= 10^6 decides it. (3) FINITE-EXACT for r <= 5
  (2, 3, 3, 4); the growth guess is from the sizes of maximal sum-free sets (affine hyperplane complements, 2^(r-1)) and
  may be wrong; what refutes it: chi(PG(5,2)) = 4 or 6. (4) a programme in unifying_relations' third theory; no descent
  is known; the onset 15, 23, 25 (THM-4510, aliquot_66 ATLAS section 6) is the hostile any descent must reproduce.
  (5) session85 measured 86, 87, 96 states for the seven-atom bank under three observation sets; the claim is that
  the first number is determined by (M, conflict graph) and that the +1 for F = U + S is a general one-bit law.
  (6) a definition-level claim; THM-4811 and the eth-false reflection give the vocabulary; a proof would need the
  three-level decider ladder stated as theorems about dual supports. (7) untested; the bijection is combinatorial
  (THM-4840 (B)); an arithmetic involution would have to permute the fifteen primes and preserve square sigma.
session: klein-2026-10-28 (codex/session-synthesis-20261028)
source: 05-knowledge/results/unifying_binary_dependence_free_bits_7_15_56_20261028/README.md (sections 3, 4, 7, 8)
related:
  - 01-canon/theorems/THM-4840-tutte-camion-dichotomy-for-arithmetic-square-class-banks-hamming-15-is-the-equal-parity-degree-space-of-k6-rank-6-prime-banks-unit-cayley-tournaments-c3-t5-p7-and-fano-chromatic-numbers.md
  - 05-knowledge/results/session85_arithmetic_codes_free_bits_20261028/ARITHMETIC_FREE_BITS.md
  - 05-knowledge/results/unifying_relations_20261028/THEORY.md
  - 05-knowledge/results/analogy_reduction_atlas_20261025/PROBLEM_CLASSIFICATION.md
---

# HYP-9438 — seven open problems on binary arithmetic dependence and free bits

1. **Prime realization of binary matroids.** For squarefree `d`, let `L(d)` be the least prime of the form `d k² − 1`.
   THM-4840 (C) exhibits `L(d)` for all `d | 30030` (largest `k = 42`). Prediction: `L(d)` exists for every squarefree
   `d` (one prime value of the quadratic `d k² − 1`, a Schinzel-H instance), with `k ≤ C log d`; hence every simple binary
   matroid is a square-class matroid of primes. Parallel elements need infinitely many primes `d k² − 1` for one `d`
   (Bunyakovsky). Decisive test: `L(d)` for all squarefree `d ≤ 10^5`; a `d` with no prime value below `k = 10^4` would
   be the first warning. The loop (prime 3) is the one element that can never be repeated.
2. **The arithmetic `R_{10}` and the arithmetic Kuratowski banks.** The regular-but-not-graphic matroids are those
   containing `R_{10}`, `M*(K_5)` or `M*(K_{3,3})` but no Fano minor (Seymour). Prediction: prime banks realizing each of
   the three exist with all primes below `10^6`, and the smallest prime bank that is regular but not graphic has ten
   elements. Test: enumerate the squarefree parts of `p + 1` for `p ≤ 10^6`, pick ten columns forming `R_{10}` (five
   row primes suffice).
3. **The sum-free chromatic ladder.** `χ(PG(r−1,2)) = 2, 3, 3, 4` for `r = 2, 3, 4, 5` (THM-4840 (D)). Prediction:
   `χ(PG(5,2)) = 4` and `χ(PG(r−1,2)) ≤ ⌈log_2(r+1)⌉ + 1`; refutation: a rank-six search returning 5. The arithmetic
   reading: the least number of colours on a rank-`r` prime bank avoiding monochromatic square-`σ` circuits.
4. **Square-sum paths by descent.** A square-sum Hamiltonian path on `[n]` is a sequence of `n − 1` squares whose
   alternating partial sums, shifted by one integer, enumerate `[n]` (THEORY §D). Prediction: there is a height on these
   sequences and a unique inverse step (deleting an endpoint and re-closing) under which every path descends to one of
   finitely many roots, which would prove the conjecture for all `n ≥ 25`; refutation: a path at `n ≤ 60` with no
   descending step. Related: THM-4510's `3–4–5` window and the Pell entry/exit events.
5. **The Nerode quotient as an invariant.** For a prime/composite atom bank with square-class matroid `M` and conflict
   graph `G`, the minimal deterministic quotient of "multiply by an atom, observe square `σ` and legality" is a function
   `q(M, G)`; adding the observation `F = U + S` adds exactly one state. Test: recompute session85's 86/87/96 for banks
   with the same `(M, G)` and different primes.
6. **The unbounded-dual-word signature.** A `Π⁰₁` cryptid is a machine whose halting predicate is a linear form of
   unbounded support in the Parikh vector of its branch word (Antihydra: `Σ(2v_i − 3) = −6`), while regular and
   bounded-Parikh deciders see only bounded-support forms. Prediction: every bbchallenge cryptid has such a form and no
   decided machine does. Test: write the halting predicates of the catalogued cryptids as dual words.
7. **Arithmetic self-duality.** The outer automorphism of `S_6` matches the fifteen primes with the fifteen matching-type
   square-`σ` triples (THM-4840 (B)). Prediction: no permutation of the fifteen primes realizes it while preserving square
   `σ` (the correlation is not an automorphism), but the bijection is realized by the Frobenius-type action on a chart
   of `F_{16}`. Test: the 20160 collineations of `PG(3,2)` versus the 16 charts.
