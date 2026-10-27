---
id: THM-4822
title: "Tao's Syracuse chain x -> (3x+1) 2^-a on Z/3^J (a ~ Geom(1/2)), and its analogue x -> (px+1) 2^-a on Z/p^J for every odd p >= 3, is rank one plus nilpotent of index exactly J: P^n(x_0, .) depends only on x_0 mod p^(J-n), P^J = Pi whose rows are the law of Syrac(Z/3^J) for p = 3, every start needs exactly J steps, and rank P^n = p^(J-n) because the one-step Koopman map from functions mod p^d onto functions mod p^(d-1) is surjective (its rows are distinct rows of the circulant of the folded geometric weights, whose eigenvalues 1/(2 zeta - 1) never vanish); hence the second-largest eigenvalue modulus is exactly 0 for every J; the Jordan type of P - Pi on the mean-zero space of Z/p^J is p - 1 blocks of size J and (p-1)^2 p^(J-1-s) blocks of size s, and on the recurrent class {y mod p in <2>} (r = |<2 mod p>|) it is r - 1 blocks of size J, r(p-2) + 1 of size J - 1 and r(p-1)^2 p^(J-2-s) of size s"
status: >
  PROVED (elementary: the closed form x_n = p^n x_0 2^-A_n + sum_j p^(j-1) 2^-(a_(n-j+1)+...+a_n), the p-adic depth
  filtration of the Koopman operator, invertibility of the circulant with eigenvalues 1/(2 zeta - 1), and its
  restriction to the closed recurrent class). FINITE-EXACT (p = 3):
  - exact rational transition matrices with entries 2^(L-a0)/(2^L - 1), L = 2 3^(J-1), for J = 1..6;
  - equal rows of P^J from every start for J <= 4 and from five starts for J = 5, 6; unequal rows of P^(J-1);
  - Tao's mod-9 stationary law (0, 8, 16, 0, 11, 4, 0, 2, 22)/63 reproduced; THM-4520 (5)'s largest primitive Fourier
    coefficients 0.37792, 0.25224, 0.17700, 0.12927, 0.09611 (J = 2..6) reproduced from the exact stationary laws;
  - modulo q = 1048573: P^J has equal rows from every start for J <= 8, P^(J-1) differs (a rigorous certificate that
    P^(J-1) != Pi), and rank_q P^n = 3^(J-n), on the units 2 3^(J-1-n), for J <= 6. The modular rank is a lower bound and
    meets the structural upper bound, so the ranks are exact.
  FINITE-EXACT (p = 5, 7; J <= 3, modulo q): the same, including the ranks on the recurrent class (r = 4, 3), and
  p = 7, where 2 is not a primitive root (the cosets of <2> matter). The independent audit (python-flint, exact) also
  covered J = 5, 6 from every start, p = 9, 11, 13, 15, 17, 21, 23, 31 at J = 2 (primality is not used), the
  recurrent-class type in 12 cases, and THM-4520's J = 7 value 0.07587. VERIFIED (float guide only): numpy's
  |lambda_2| <= 0.0013 (exact 0); the l^2(pi) norm of P_J - Pi on the units is 0.674, 0.971, 0.996, 0.9995, 0.99995 at
  J = 2..6, so the non-normality seen in THM-4821 is present here too (not proved). CITED: the chain's stationary law
  is Syrac(Z/3^J) (Tao, Forum Math. Pi 10 (2022), e12, Lemma 1.12 and Remark 1.13); statement 1 is Tao's (1.21)/(1.24),
  and its J = 2 case is in the TCPC note section 5.2 and procgen_repunit_20260924 section 2.1. The law of Pi in the form
  sum_j p^(j-1) 2^-(a_1+...+a_j) is Tao's (1.29). CONTRAST (canon): the T-map random chain x -> x/2 or (3x+1)/2 on Z/3^J
  has characteristic polynomial (lambda - 1) prod_d (lambda^(L_d) - 2^-(L_d)), i.e. SLEM exactly 1/2 (THM-4631 with
  j = 0 on the line u = v = z, det(I - zU) = (1 - z) prod_d (1 - (z/2)^(L_d)); checked exactly for J <= 3, by the
  audit for J <= 5). NOT CLAIMED: an exact l^2 profile; anything about integer orbits. Independent audit (2026-10-27):
  SOUND WITH CORRECTIONS. This file's corrections (citations, the scope of the Jordan type) are applied, statement 4 is
  added, and the note's mis-indexed display of the law of Pi is corrected (MISTAKE-686; not in this file).
session: opus-2026-10-27-S65
source: 05-knowledge/results/syracuse_transfer_finite_section_20261027.md
scripts:
  - 04-computation/experiments/syracuse_transfer_finite_section_20261027.py (+ .out; sections F, G, H)
related:
  - THM-4821 (the 2-adic carrier)
  - THM-4520 (the level operator: Gauss phases e(2^k/3^n), circulant eigenvalues 1/(2 e(xi/L) - 1) on |w - 1/3| = 2/3)
  - THM-4519 (the frequency recursion and its hypothesis H1; the cocycle of the level blocks)
  - THM-4631 (twisted level factorisation; doubling cycles; the T-map contrast)
  - 05-knowledge/results/procgen_tcpc_20261001_tournament_clock_prime_collatz.md section 5.2 (the mod-9 chain: rank 2 on the units, P^2 = 1 pi)
  - 05-knowledge/results/procgen_repunit_20260924_digits_rotation_repunits.md section 2.1 (Theorem H: the mod-9 chain and Tao's law)
---

# THM-4822 -- Tao's Syracuse chain mod 3^J is rank one plus nilpotent of index J

**Setting.** p ≥ 3 odd, J ≥ 1, a₁, a₂, … i.i.d. with P(a = k) = 2^-k. P is the chain x ↦ (px + 1)·2^-a on Z/p^J (2 is
invertible). For p = 3 this is Tao's chain, whose stationary law is Syrac(Z/3^J).

**Statement.**
1. x_n = p^n x₀ 2^-(A_n) + Σ_(j=1..n) p^(j−1) 2^-(a_(n−j+1) + … + a_n). So P^n(x₀, ·) depends only on x₀ mod p^(J−n),
   and P^J = Π. The rows of Π are the law of Σ_(j=1..J) p^(j−1) 2^-(a_1 + … + a_j).
2. rank P^n = p^(J−n) for 0 ≤ n ≤ J.
3. P = Π + N with N nilpotent of index exactly J. The spectrum is {1} with 0 otherwise, so SLEM = 0.
   - On the mean-zero space N has p − 1 Jordan blocks of size J and (p − 1)² p^(J−1−s) blocks of size s (1 ≤ s ≤ J − 1).
   - No row of P^(J−1) equals π, so every start needs exactly J steps.
4. The type in 3 counts all of Z/p^J, transient states included. Let R = supp π = {y : y mod p ∈ ⟨2⟩} and
   r = |⟨2 mod p⟩| ≥ 2 (for p = 3, R is the units). R is closed.
   - rank P_R^n = r p^(J−1−n) for n ≤ J − 1, and 1 at n = J.
   - For J ≥ 2, N on the π-mean-zero functions on R has r − 1 Jordan blocks of size J, r(p − 2) + 1 of size J − 1, and
     r(p − 1)² p^(J−2−s) of size s ≤ J − 2.
   - For p = 3 that is 1, 3 and 8·3^(J−2−s). The index is J.

**Proof.**
- *(1)* is induction on n. The law of Π follows by reversing the i.i.d. a's.
- *(2)* Let F_d be the functions of y mod p^d and U the Koopman operator. Then U F_d ⊆ F_(d−1), because
  (px + 1)2^-a mod p^d depends only on x mod p^(d−1).
  - The induced map F_d → F_(d−1) is a matrix with rows x̄ ∈ Z/p^(d−1).
  - Row x̄ lives on the coset (px + 1)⟨2⟩ of the units mod p^d. Along that cycle of length L = ord(2 mod p^d) it is the
    folded geometric weight w(j) = 2^-j/(1 − 2^-L) shifted to px + 1.
  - Rows on different cosets are disjointly supported.
  - Rows on the same coset are distinct shifts of w, since x̄ ↦ px + 1 is injective, i.e. distinct rows of Circ(w).
  - Circ(w) has eigenvalues 1/(2ζ − 1) with ζ^L = 1, never 0, so it is invertible. The rows are independent and the
    map is onto.
  - Therefore U^n(F_J) = F_(J−n) and rank P^n = p^(J−n).
- *(3)* Ranks give the Jordan type. The index is J since rank N^(J−1) = p − 1 > 0 on the mean-zero space.
  - The p rows of P^(J−1) are independent.
  - π = Σ_c π_c row_c, where π_c = π{x ≡ c mod p} > 0 for the r ≥ 2 classes c ∈ ⟨2⟩. So π is no single row.
- *(4)* Every image is ≡ 2^-a mod p, so R is closed and P_R^n is the restriction of P^n.
  - U maps the functions on R of y mod p^d (dimension r p^(d−1), or 1 at d = 0) onto those of y mod p^(d−1): extend
    off R, use the surjection of (2), and restrict.
  - So rank P_R^n = r p^(J−1−n), and the block counts follow from consecutive rank differences.

∎
