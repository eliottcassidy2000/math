---
id: THM-4821
title: "The Ulam finite section P_K of the Syracuse map (px+1)/2^v (p odd; p = 3 is Collatz) on the odd residues modulo 2^K is doubly stochastic and rank one plus nilpotent: the uniform law on a residue class goes to the uniform law on a class of strictly smaller precision (or to the stationary law), so P_K^n has rows the uniform laws on all classes of precision <= K - n, each class of precision m >= 2 repeated binom(K-m-1, n-1) times, and rank 2^(K-1-n); P_K^(K-1) = Pi != P_K^(K-2); the second-largest eigenvalue modulus is exactly 0 for every K; Jordan type of P_K - Pi on the mean-zero space: one block of size K - 1 and 2^(K-2-s) blocks of size s; exact time to equilibrium h with h - 1 ~ Binomial(K-2, 1/2) over the residues (#{h} = 2 C(K-2, h-1), mean K/2), worst-case total variation 1 - 2^-(K-1-n), Dobrushin coefficients 1, ..., 1, 0; singular values of P_K^n: 1 and sqrt(P(Bin(K-j, 1/2) >= n)) with multiplicity 2^(j-2), j = 2..K, so ||(P_K - Pi)^n||^2 = P(Bin(K-2, 1/2) >= n) and the second singular value is sqrt(1 - 2^-(K-2)); the matrices for different odd p are permutation-conjugate (the classical parity-vector conjugacy), so none of this depends on p"
status: >
  PROVED (elementary: the Haar class lemma, the bijection c -> pc+1 from odd to even residues, a span count on the
  binary tree of classes, the subset-counting of intermediate precisions, the operator identity (P^n)* P^n =
  sum_m w_m(n) E_m over nested conditional expectations, and the parity-vector conjugacy). FINITE-EXACT for K = 8..20,
  p = 3:
  - the counting matrix C_K = 2^(K-1) P_K is built in exact integers, with row and column sums 2^(K-1) and K 2^(K-1)
    nonzero entries;
  - the class lemma is checked for every class at every level, and so are the reachable sets with their multiplicities
    binom(K-m-1, n-1), the ranks (exact tree span), the Jordan type, the hitting-time counts, E h = K/2, the worst-case
    TV, the Dobrushin sequence and every tail weight P(Bin(K-j, 1/2) >= n) of the singular spectrum;
  - independently, C^(K-1) = 2^((K-1)(K-2)) J is proved by congruence modulo primes with product above the entry bound
    for K = 8..12 (exact integers at K = 8), and C^(K-2) differs modulo 1048573 (a rigorous certificate P^(K-2) != Pi);
  - the rows equal direct Haar lift counts for K = 3..7 and p = 3, 5, 7, 1, -1 (script section I).
  The conjugacy P_K^(p) = Q P_K^(3) Q^(-1) is checked entry by entry for p = 1, 5, 7, 9, 15, -1, -3, 1023 at K = 8..14 and
  p = 5, -1 at K = 15..20. VERIFIED (float guide only): numpy SVD gives the full singular spectrum of every power for
  K <= 11 and its norm at K = 12. numpy's eigenvalues report spurious |lambda_2| = 0.0029..0.0179 at K = 8..12 (natural
  ordering; ordering-dependent), the pseudospectrum of the Jordan block of size K - 1; the exact value is 0.
  CITED (classical): the parity-vector map conjugates x -> x/2, (px+1)/2 to the one-sided Bernoulli shift under Haar
  measure (Terras 1976; Lagarias 1985, Amer. Math. Monthly 92, 3-23; from memory), whose finite-section shadow this is;
  statement 6 is its corollary. Canon THM-4629 proof 2(c)-(d) and THM-4631 prove the corresponding de Bruijn structure
  for the inverse-tree (T-map) operator: 2(c) the halving-branch averaging E_j, rank E_j^q = 2^(j-q); 2(d) the full
  de Bruijn walk with the single nonzero eigenvalue (u + v)/2. The exact index K of the T-map Ulam matrix is classical
  (de Bruijn A^K = J). The case K = 3 (mod-8 chain, P^2 uniform) is in procgen_repunit_20260924 section 2.1. NOT
  CLAIMED: any statement about integer orbits. Independent audit (2026-10-27): SOUND WITH CORRECTIONS; the corrections
  (citations, the conjugacy form of statement 6, the reproducibility of the lift-count check) are applied, and the
  audit's strengthenings (multiplicities, full singular spectrum) are added.
session: opus-2026-10-27-S65
source: 05-knowledge/results/syracuse_transfer_finite_section_20261027.md
scripts:
  - 04-computation/experiments/syracuse_transfer_finite_section_20261027.py (+ .out; sections A, B, C, C2, D, E, H, I)
related:
  - THM-4822 (the 3-adic / p-adic carrier: Tao's chain, rank one plus nilpotent of index J)
  - THM-4629 (proof 2(c)-(d): the de Bruijn structure of the inverse-tree compression; statement 2: silence of finite truncations)
  - THM-4631 (twisted level factorisation; the 2-adic part contributes only the even density; p enters only through the odd branch)
  - THM-4520, THM-4519 (the Fourier-side level operator and the frequency recursion)
  - 05-knowledge/results/procgen_repunit_20260924_digits_rotation_repunits.md section 2.1 (the case K = 3; the T-map window chain)
  - 05-knowledge/results/graph_theory_as_common_carrier_20261013.md section 5 (the T-map class graph mod 2^j is the de Bruijn graph B(2, j))
---

# THM-4821 -- the Syracuse Ulam finite section mod 2^K is rank one plus nilpotent

**Setting.** K ≥ 2, p odd, S_p(x) = (px + 1)/2^(v₂(px+1)) on the odd 2-adic integers with Haar measure (preserved by
S_p). Ω_K = odd residues mod 2^K, π uniform. P_K(a, b) = Haar{x ≡ a : S_p(x) ≡ b}/Haar{x ≡ a} (mod 2^K). Row a is
uniform on the 2^v residues ≡ (pa + 1)/2^v (mod 2^(K−v)) if v = v₂(pa+1) < K, and uniform on Ω_K for the one residue with
pa + 1 ≡ 0.

**Statement.**
1. *(Class lemma.)* μ_C P_K = μ_(Φ(C)) for every class C = c + 2^m Z (c odd, m ≤ K), where Φ(c, m) =
   ((pc + 1)/2^v mod 2^(m−v), m − v) if v = v₂(pc + 1) < m, and Φ(C) = Ω otherwise.
2. For 1 ≤ n ≤ K − 1, the rows of P_K^n are exactly the uniform laws on all classes of precision 1..K − n.
   - A class of precision m ≥ 2 is the row of exactly C(K − m − 1, n − 1) residues.
   - rank P_K^n = 2^(K−1−n).
3. P_K = Π + N, N^(K−1) = 0 ≠ N^(K−2). The spectrum is {1} with 0 otherwise, so SLEM = 0. On the mean-zero space N has
   one Jordan block of size K − 1 and 2^(K−2−s) blocks of size s (1 ≤ s ≤ K − 2).
4. #{a : h(a) = h} = 2 C(K−2, h−1), where h(a) = min{n : P^n(a, ·) = π}. Also:
   - max_a ‖P^n(a, ·) − π‖_TV = 1 − 2^-(K−1−n);
   - δ(P^n) = 1 for n ≤ K − 2 and 0 at n = K − 1.
5. The singular values of P_K^n are 1 together with √P(Bin(K − j, 1/2) ≥ n), of multiplicity 2^(j−2), for j = 2..K.
   - In particular ‖(P_K − Π)^n‖²_(ℓ²(π)) = P(Bin(K − 2, 1/2) ≥ n), and σ₂(P_K) = √(1 − 2^-(K−2)).
6. For every odd p, P_K^(p) = Q P_K^(3) Q^(−1). Here Q is the permutation of the odd residues induced by
   ψ = Q_3^(−1) ∘ Q_p, and Q_p is the parity-vector map of x ↦ x/2, (px + 1)/2, a tree automorphism of Z/2^K. Hence 1–5,
   and every property invariant under relabelling, do not depend on p.

**Proof.** The results note, section 2, gives the full proof; the outline follows.
- *(1)* x ↦ px + 1 is a Haar-preserving bijection of classes. If the valuation is constant on the class, division by 2^v
  maps it onto a class of precision m − v. Otherwise the odd part of a Haar-uniform element of 2^m Z₂ is Haar-uniform on
  Z₂^×.
- *(2)* c ↦ pc + 1 maps the odd residues mod 2^m bijectively onto the even ones, so Φ hits each class of precision
  2..m − 1 once and the root twice.
  - Induction gives the rows and the multiplicities mult_(n+1)(m′) = Σ_(m>m′) mult_n(m), i.e. C(K − m − 1, n − 1).
  - The span of all classes of precision ≤ K − n is that of the 2^(K−n−1) disjoint classes of precision K − n.
- *(3)* The rank sequence gives the Jordan type.
- *(4)* Each set of intermediate precisions in [2, K − 1] comes from exactly two residues. The total variation of a class
  of precision m from π is 1 − 2^(1−m).
- *(5)* By the level-only multiplicities, ‖P^n f‖² = Σ_m w_m(n)‖E_m f‖² for every f, so (P^n)*P^n = Σ_m w_m(n) E_m.
  - The E_m are nested conditional expectations. On range E_j ⊖ range E_(j−1), of dimension 2^(j−2), this operator is
    the scalar Σ_(m≥j) w_m(n), the fraction of residues still at precision ≥ j after n steps.
  - That fraction is P(Bin(K − j, 1/2) ≥ n) by the subset count of (4). The norm is attained at f = ±1 on 1 and 3 mod 4.
- *(6)* Q_p is a bijection of Z/2^K, compatible across K, for every odd p (Terras). So ψ is a Haar-preserving tree
  automorphism that conjugates T_p to T_3.
  - S_p is the first return of T_p to the odd residues, which in parity coordinates is "shift past the next 1" for every
    p. So ψ conjugates S_p to S_3 and maps classes to classes.

∎
