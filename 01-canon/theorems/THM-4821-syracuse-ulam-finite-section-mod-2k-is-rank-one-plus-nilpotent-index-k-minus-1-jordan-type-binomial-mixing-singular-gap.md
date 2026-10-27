---
id: THM-4821
title: "The Ulam finite section P_K of the Syracuse map (px+1)/2^v (p odd; p = 3 is Collatz) on the odd residues modulo 2^K is doubly stochastic and rank one plus nilpotent: the uniform law on a residue class goes to the uniform law on a class of strictly smaller precision (or to the stationary law), so P_K^n has rows the uniform laws on all classes of precision <= K - n and rank 2^(K-1-n), P_K^(K-1) = Pi != P_K^(K-2); the second-largest eigenvalue modulus is exactly 0 for every K; Jordan type of P_K - Pi on the mean-zero space: one block of size K - 1 and 2^(K-2-s) blocks of size s; exact time to equilibrium h with h - 1 ~ Binomial(K-2, 1/2) over the residues (#{h} = 2 C(K-2, h-1), mean K/2), worst-case total variation 1 - 2^-(K-1-n), Dobrushin coefficients 1, ..., 1, 0; l^2 profile ||(P_K - Pi)^n||^2 = P(Bin(K-2, 1/2) >= n), so the second singular value is sqrt(1 - 2^-(K-2)); none of this depends on p"
status: >
  PROVED (elementary: the Haar class lemma, the bijection c -> pc+1 from odd to even residues, a span count on the
  binary tree of classes, the subset-counting of intermediate precisions, and a conditional-expectation bound for the
  l^2 norm). FINITE-EXACT for K = 8..20, p = 3:
  - the counting matrix C_K = 2^(K-1) P_K is built in exact integers, with row and column sums 2^(K-1) and K 2^(K-1)
    nonzero entries;
  - the class lemma is checked for every class at every level; reachable sets, ranks (exact tree span), Jordan type,
    hitting-time counts, E h = K/2, worst-case TV, Dobrushin sequence and the exact l^2 value at the maximiser;
  - independently, C^(K-1) = 2^((K-1)(K-2)) J is proved by congruence modulo primes with product above the entry bound
    for K = 8..12 (exact integers at K = 8), and C^(K-2) differs modulo 1048573 (a rigorous certificate P^(K-2) != Pi);
  - the rows were checked against direct lift counts for K = 3..7.
  Universality is checked for p in {1, 5, 7, -1} at K = 8..14. VERIFIED (float guide only): numpy SVD of every power for
  K <= 12 agrees with the l^2 profile. numpy's eigenvalues report spurious |lambda_2| = 0.0029..0.0179 at K = 8..12
  (pseudospectrum of the Jordan block of size K - 1); the exact value is 0. CITED (classical, from memory): the 2-adic
  Collatz map is a Bernoulli shift under Haar measure (Lagarias 1985), whose finite-section shadow this is. Canon
  THM-4629 proof 2(c) and THM-4631 prove the corresponding de Bruijn nilpotency for the inverse-tree (T-map) operator. NOT
  CLAIMED: any statement about integer orbits. Independent audit: OWED (pre-audit checkpoint).
session: opus-2026-10-27-S65
source: 05-knowledge/results/syracuse_transfer_finite_section_20261027.md
scripts:
  - 04-computation/experiments/syracuse_transfer_finite_section_20261027.py (+ .out; sections A, B, C, C2, D, E, H)
related:
  - THM-4822 (the 3-adic / p-adic carrier: Tao's chain, rank one plus nilpotent of index J)
  - THM-4629 (proof 2(c): de Bruijn nilpotency of the inverse-tree compression; statement 2: silence of finite truncations)
  - THM-4631 (twisted level factorisation; the 2-adic part contributes only the even density)
  - THM-4520, THM-4519 (the Fourier-side level operator and the frequency recursion)
---

# THM-4821 -- the Syracuse Ulam finite section mod 2^K is rank one plus nilpotent

**Setting.** K ≥ 2, p odd, S_p(x) = (px + 1)/2^(v₂(px+1)) on the odd 2-adic integers with Haar measure (preserved by
S_p). Ω_K = odd residues mod 2^K, π uniform. P_K(a, b) = Haar{x ≡ a : S_p(x) ≡ b}/Haar{x ≡ a} (mod 2^K). Row a is
uniform on the 2^v residues ≡ (pa + 1)/2^v (mod 2^(K−v)) if v = v₂(pa+1) < K, and uniform on Ω_K for the one residue with
pa + 1 ≡ 0.

**Statement.**
1. *(Class lemma.)* μ_C P_K = μ_(Φ(C)) for every class C = c + 2^m Z (c odd, m ≤ K), where Φ(c, m) =
   ((pc + 1)/2^v mod 2^(m−v), m − v) if v = v₂(pc + 1) < m, and Φ(C) = Ω otherwise.
2. For 1 ≤ n ≤ K − 1, the rows of P_K^n are exactly the uniform laws on all classes of precision 1..K − n. The
   multiplicity of a class depends only on its precision, and rank P_K^n = 2^(K−1−n).
3. P_K = Π + N, N^(K−1) = 0 ≠ N^(K−2). The spectrum is {1} with 0 otherwise, so SLEM = 0. On the mean-zero space N has
   one Jordan block of size K − 1 and 2^(K−2−s) blocks of size s (1 ≤ s ≤ K − 2).
4. #{a : h(a) = h} = 2 C(K−2, h−1), where h(a) = min{n : P^n(a, ·) = π}. Also:
   - max_a ‖P^n(a, ·) − π‖_TV = 1 − 2^-(K−1−n);
   - δ(P^n) = 1 for n ≤ K − 2 and 0 at n = K − 1.
5. ‖(P_K − Π)^n‖²_(ℓ²(π)) = P(Bin(K − 2, 1/2) ≥ n). In particular σ₂(P_K) = √(1 − 2^-(K−2)).
6. Statements 1–5 do not depend on p.

**Proof.** The results note, section 2, gives the full proof; the outline follows.
- *(1)* x ↦ px + 1 is a Haar-preserving bijection of classes. If the valuation is constant on the class, division by 2^v
  maps it onto a class of precision m − v. Otherwise the odd part of a Haar-uniform element of 2^m Z₂ is Haar-uniform on
  Z₂^×.
- *(2)* c ↦ pc + 1 maps the odd residues mod 2^m bijectively onto the even ones, so Φ hits each class of precision
  2..m − 1 once and the root twice. Induction gives the rows, and the span of all classes of precision ≤ K − n is that of
  the 2^(K−n−1) disjoint classes of precision K − n.
- *(3)* The rank sequence gives the Jordan type.
- *(4)* Each set of intermediate precisions in [2, K − 1] comes from exactly two residues. The total variation of a class
  of precision m from π is 1 − 2^(1−m).
- *(5)* For f of mean 0, ‖P^n f‖² = Σ_(m≥2) w_m(n)‖E_m f‖², maximised by f = ±1 on 1 and 3 mod 4.

∎
