---
id: THM-4826
title: "Which grid carries a norm shell: for Z[omega] and the Hurwitz order (residue field F_4 over 2, first grid Z + P, a unit of order 3) exactly one third of every odd-norm shell lies on the first grid and every even-norm element does (Lipschitz 8 sigma(n) of Hurwitz 24 sigma(n)); among the rings Z[(1+sqrt-d)/2] (d = 3 mod 4) the centred grid carries a positive share of every odd shell only for d = 3 (d = 7 mod 8: none; d = 3 mod 8, d > 3: prime shells lie on one grid, first-grid share 1/3 by the cubic ring class field); the primary prime over p = 1 mod 3 lies on the first grid iff 2 is a cubic residue iff p = x^2 + 27y^2, so the primes over Mersenne primes 2^q - 1 (q > 3) lie on the first grid; the units of each grid alone are not 4-designs but their unions (hexagon, 24-cell) are 5-designs, which is what Candido's and Liouville's quartic identities need; Eisenstein fourth powers are 0, u or 9u mod 16 (u in mu_3 = mu_6^4, the merge charges), so a rational n = 5 mod 8 is never a sum of four Eisenstein fourth powers"
status: >
  PROVED (elementary: units acting on (O/P)* = F_4*; the Frobenius on F_4; the binomial expansion of (1 + 2 beta)^4; a
  residue count mod 8). CITED: Jacobi's r_4(n) = 8 sigma(n) (n odd) and Hurwitz's 24 sigma(n); the cubic reciprocity
  supplement for 2 (Ireland and Rosen, Prop. 9.6.1-9.6.2, restated as Thm 4.2 of Dechakulkamjorn-Rungtanapirom,
  arXiv:2208.02454) and Gauss's p = x^2 + 27y^2 (Cox, Primes of the Form x^2 + ny^2, Thm 4.15); the ring class field of
  Z[sqrt-d] and the class number formula for orders (Cox Thm 7.24, Thm 9.4); Candido's identity (Candido 1951; the
  special case is older, Proth 1878) and Liouville's identity with Liouville's bound g(4) <= 53 (now g(4) = 19,
  Balasubramanian-Deshouillers-Dress 1986); spherical designs (Delsarte-Goethals-Seidel 1977). FINITE-EXACT: the split
  law for Z[omega] (n <= 20,000) and the Hurwitz order (n <= 300); d = 7 mod 8 (odd n <= 3000, d <= 47); d = 11, 19,
  43, 67 (every split prime below 3*10^5 on one grid; first-grid shares 0.3325, 0.3334, 0.3343, 0.3336); h(-4d)/h(-d)
  for squarefree d = 3 mod 4 below 200; the cubic-character equivalence for all 39,231 primes p = 1 mod 3 below 10^6
  (0 mismatches; residue share 0.3322); x^2 + 27y^2 representations of 2^q - 1 for q = 5, ..., 127; the design table
  and both identities (sympy); Eisenstein fourth powers mod 16; the n = 5 (mod 8) obstruction and the absence of other
  obstructions modulo 2^8, 3^5, 5^3, 7^2 on 0..200. The representability of the remaining n (6 among them) is HYP-9429.
  NEW TO THE REPOSITORY: the uniform two-grid statement of (1)-(2), the design-completion reading of (4), and (6); the
  rest is classical in grid language. Independent audit: OWED (pre-audit checkpoint).
session: opus-2026-10-27-S66
source: 05-knowledge/results/two_grids_ring_or_ideal_class_13_37_20261027.md
scripts:
  - 04-computation/experiments/two_grids_class_group_20261027.py (+ .out; sections D, E, F)
  - 04-computation/experiments/eisenstein_four_fourth_powers_20261027.c (+ .out)
related:
  - THM-4825 (the second grid is a ring or the order-two class over 2; d = 5, 13, 37)
  - THM-4824 (mac-mini-2026-10-27-slab, concurrent: the zero-sum quartic identity (F), the parity-coset shells (E))
  - THM-4820 (Eisenstein merge debt; the charges 1, omega, omega^2 are the unit fourth powers)
  - HYP-9429 (sums of four Eisenstein fourth powers; 6)
---

# THM-4826 -- which grid carries a norm shell; the cubic character of 2; the quartic designs

**Setting.** O is Z[ω] (ω² + ω + 1 = 0) or the Hurwitz order H = Z⟨i, j, k⟩ ∪ (½(1 + i + j + k) + Z⟨i, j, k⟩).
- 𝔓 is the prime of O over 2: 2Z[ω], respectively (1 + i)H = D₄. In both cases O/𝔓 ≅ F₄.
- The *first grid* is L = Z + 𝔓: Z[√−3], respectively the Lipschitz order Z⟨i, j, k⟩. It has index 2.
- The *centred grid* is O ∖ L.

Norms are the field norm, respectively the reduced norm.

**Statement.**
1. **(The 1 : 2 split.)** For every odd n, exactly one third of the elements of O of norm n lie on the first grid. Every
   element of even norm lies on it.
   - Z[ω]: on odd shells the first grid carries r(n)/3 = 2Σ_(e|n) χ₋₃(e) and the centred grid 2r(n)/3; even shells lie
     on the first grid.
   - H: on odd shells the first grid carries 8σ(n) (Jacobi) and the centred grid 16σ(n), of Hurwitz's 24σ(n); even shells
     lie on the first grid.
2. **(The family Z[(1 + √−d)/2], d ≡ 3 (mod 4), first grid Z[√−d].)**
   - If d ≡ 7 (mod 8), every element of odd norm lies on the first grid.
   - If d ≡ 3 (mod 8) and d > 3, the elements of norm p (p a split prime) lie all on the first grid or all on the centred
     grid; the first iff p = x² + dy². When h(−d) = 1, the first-grid primes have density 1/3 among split primes, since
     the ring class field of Z[√−d] has degree h(−4d)/h(−d) = 3 over K.
   - Hence the centred grid carries a positive share of every non-empty odd shell only for d = 3.
3. **(The cubic character of 2.)** Let p ≡ 1 (mod 3), and let π = a + bω be the primary prime over p (π ≡ 2 mod 3).
   - π lies on the first grid (b even) ⟺ 2 is a cubic residue mod p ⟺ p = x² + 27y².
   - For a Mersenne prime p = 2^q − 1 with q > 3, 2 = (2^m)³ where 3m ≡ 1 (mod q). So π lies on the first grid and
     p = x² + 27y². For example, 31 = 2² + 27·1² and 127 = 10² + 27·1².
4. **(Designs.)** The units of the first grid are {±1}, respectively the 8 Lipschitz units (the 16-cell). The units of the
   centred grid are {±ω, ±ω²}, respectively the 16 units ½(±1 ± i ± j ± k) (the tesseract).
   - Neither set is a spherical 4-design.
   - Their unions, the hexagon and the 24-cell, are 5-designs and not 6-designs. Hence:
     - Σ_(u ∈ O*) ⟨u, x⟩⁴ = (9/4)|x|⁴, respectively 3|x|⁴. In integral form these are Candido's identity
       a⁴ + b⁴ + (a + b)⁴ = 2(a² + ab + b²)² and, through (1 + i)H* = the 24 roots of D₄, Liouville's identity
       6(x₁² + x₂² + x₃² + x₄²)² = Σ_(i<j) [(x_i + x_j)⁴ + (x_i − x_j)⁴].
   - The Gaussian units are a 3-design and not a 4-design.
5. **(Fourth powers mod 16.)** For α ∈ Z[ω]:
   - α⁴ ≡ 0 (mod 16) if α ∈ 2Z[ω];
   - otherwise α⁴ ≡ u or 9u (mod 16) with u = α⁴ (mod 2) ∈ μ₃ = {1, ω, ω²}.

   Also μ₆⁴ = μ₃: the merge charges 1, ω, ω² of THM-4820 are exactly the fourth powers of the six units.
6. **(Four Eisenstein fourth powers.)** A rational integer n ≡ 5 (mod 8) is not a sum α⁴ + β⁴ + γ⁴ + δ⁴ with α, β, γ, δ
   ∈ Z[ω]. In particular 13 and 37 are excluded over Z[ω] as well as over Z.

**Proof.**
- *(1)* Elements of odd norm are units modulo 𝔓 (N𝔓 = 4), so they lie in (O/𝔓)* = F₄*. Elements of even norm lie in
  𝔓 ⊂ L (2 is inert in Q(ω); 𝔓 is the unique prime over 2 in H).
  - L ∖ 𝔓 is the class 1. A unit ζ of order 3 (ω, respectively ½(−1 + i + j + k)) is not ≡ 1 (mod 𝔓), so it acts on F₄* by
    a 3-cycle (by left multiplication in H; 𝔓 is two-sided).
  - The norm-n elements are therefore equidistributed over the three classes. The counts are the classical r(n), 8σ(n)
    and 24σ(n).
- *(2)* If d ≡ 7 (mod 8), then O/2O = F₂ × F₂ and L/2O is the diagonal. An element of odd norm is ≡ (1, 1).
  - If d ≡ 3 (mod 8), then O/2O = F₄ and the units ±1 act trivially. Conjugation acts on F₄ as the Frobenius, fixing 1
    and swapping the other two classes. So ±π and ±π̄ are all in L or all outside it.
  - π ∈ L ⟺ p = x² + dy² ⟺ p splits completely in the ring class field of Z[√−d] (Cox Thm 9.4). Its degree over K is
    h(−4d)/h(−d) = 3 (Cox Thm 7.24, with unit index 1 for d > 3), and Chebotarev gives density 1/3.
- *(3)* The supplementary law gives (2/π)₃ ≡ π^((N2 − 1)/3) = π (mod 2), so 2 is a cubic residue iff π ≡ 1 (mod 2), iff
  b is even (Ireland–Rosen 9.6.1–2). For primary π on the first grid, π = (a − b/2) + (b/2)√−3 with 3 | b/2, so
  p = x² + 27y².
- *(4)* These are polynomial identities in the coordinates, verified exactly. The design property of the unions is the
  4-design identity together with antipodality.
- *(5)* If α ∉ 2Z[ω], then α = u₀(1 + 2β) for a unit u₀, since μ₆ → F₄* is onto. So
  α⁴ ≡ u₀⁴(1 + 8(β + 3β²)) (mod 16), and β + 3β² ≡ Tr_(F₄/F₂)(β mod 2) ∈ {0, 1} (mod 2).
- *(6)* Write α_k⁴ ≡ ε_k·u_k (mod 16) with ε_k ∈ {0, 1, 9} and u_k ∈ μ₃. Let A, B and C be the sets of k with ε_k ≠ 0
  and u_k = 1, ω and ω², respectively. In the basis (1, ω):
  - the ω-coordinate is Σ_B ε − Σ_C ε ≡ |B| − |C| (mod 8), since 1 ≡ 9 (mod 8);
  - rationality forces |B| = |C| =: k, as both are at most 4;
  - the rational coordinate is ≡ |A| − k (mod 8), with |A| + 2k ≤ 4.

  So n mod 8 ∈ {0, 1, 2, 3, 4} ∪ {7, 0, 1} ∪ {6}, which never contains 5. ∎
