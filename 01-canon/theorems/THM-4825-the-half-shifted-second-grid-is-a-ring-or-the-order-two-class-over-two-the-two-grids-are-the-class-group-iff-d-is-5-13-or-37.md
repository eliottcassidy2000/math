---
id: THM-4825
title: "The owner's second grid, Z[sqrt-d] plus its (1/2,1/2)-centred copy (d odd squarefree), completes Z[sqrt-d] to the ring of integers when d = 3 mod 4 (Eisenstein at d = 3) and is half the prime ideal over 2 when d = 1 mod 4, an ideal class of order exactly 2 for d > 1; the two grids are the whole class group iff h(-4d) = 2 iff d in {5, 13, 37}; for those d a split prime is a first-grid norm iff it is 1 mod 4 iff it is a Gaussian norm (and this agreement characterises {1, 5, 13, 37} among squarefree d = 1 mod 4); numbers whose prime factors all ramify sit on one grid according to the parity of their power of 2, so the owner's 5, 10, 25, 50, 100 alternate first, second, first, second, first in Z[sqrt-5]"
status: >
  PROVED (elementary: the ideal (2, 1 + sqrt-d) is the checkerboard of Z[sqrt-d]; its square is (2); norm-2 elements
  exist only for d = 1). CITED: the complete list of imaginary quadratic fields of class number 2 (Baker 1971; Stark,
  Math. Comp. 29 (1975) 289-302; OEIS A014603: -15, -20, -24, -35, -40, -51, -52, -88, -91, -115, -123, -148, -187,
  -232, -235, -267, -403, -427), whose members of the form -4d with d = 1 mod 4 are exactly -20, -52, -148; genus
  theory and the Hilbert class field H(Q(sqrt-d)) = Q(sqrt-d, i) for d in {5, 13, 37} (THM-4566, owner-directive canon;
  the S40 note binary_information_class_groups_largest_anomalies_20261017.md section 2, Cox Thm 3.22); Bauer's theorem
  (Chebotarev) for the converse in (4). FINITE-EXACT: the ring/ideal dichotomy for every odd squarefree d <= 199; h(-4d)
  = 2 exactly at d = 5, 13, 37 among the 40,520 squarefree d = 1 mod 4 below 200,000 (exact reduced-form search); the
  reduced forms (1,0,d), (2,2,(1+d)/2) for d = 5, 13, 37; the three-grid agreement on all 26,839 split primes below
  200,000 for d = 5, 13, 37; first counterexamples for d = 17, 21, 29, 33, 41, 53, 61, 65, 73, 85; the owner's chains.
  NEW TO THE REPOSITORY: the identification of the owner's construction with the ring completion / the order-two class
  over 2, and the resulting characterisation of {5, 13, 37}; the ingredients are classical. Independent audit: OWED
  (pre-audit checkpoint).
session: opus-2026-10-27-S66
source: 05-knowledge/results/two_grids_ring_or_ideal_class_13_37_20261027.md
scripts:
  - 04-computation/experiments/two_grids_class_group_20261027.py (+ .out; sections A, B, C)
related:
  - THM-4826 (which grid carries a norm shell; the cubic character of 2; the quartic designs; four Eisenstein fourth powers)
  - THM-4824 (mac-mini-2026-10-27-slab, concurrent: the Eisenstein grid as the first grid glued with its (1/2,1/2) shift, merge debt by lattice shape)
  - THM-4566 (Hilbert-class-field planes are the idoneal planes; chi = 2 for m = 1, 5, 13, 37, 85)
  - THM-4820 (Eisenstein merge debt; its proof (C) uses the coordinates u = 2a - b of the two grids)
---

# THM-4825 -- the half-shifted second grid is a ring or the order-two class over 2

**Setting.** d ≥ 1 odd and squarefree; K = Q(√−d); L_d = Z[√−d] (the *first grid*: the points u + v√−d, (u, v) ∈ Z²);
G_d = L_d ∪ ((1 + √−d)/2 + L_d), the first grid with its copy centred at (½, ½). So 2G_d = {u + v√−d : u ≡ v (mod 2)}
is the checkerboard of L_d.

**Statement.**
1. **(Ring completion.)** G_d is closed under multiplication iff d ≡ 3 (mod 4). Then G_d = O_K, and L_d is the order of
   conductor 2. Example: d = 3, Z[ω] = Z[√−3] ∪ ((1 + √−3)/2 + Z[√−3]).
2. **(The class over 2.)** If d ≡ 1 (mod 4), then O_K = L_d and 2G_d = 𝔭₂ = (2, 1 + √−d), the prime of O_K over 2.
   - It is ramified: 𝔭₂² = (2).
   - 𝔭₂ is principal iff d = 1, where 𝔭₂ = (1 + i) and G_1 is the square grid turned by 45°.
   - For d > 1, the second grid G_d = ½𝔭₂ represents an ideal class of order exactly 2. Its normalised norm form is
     2N(x + y(1 + √−d)/2) = 2x² + 2xy + ((1 + d)/2)y².
3. **(Two grids = the class group.)** For d ≡ 1 (mod 4): Cl(O_K) = {[L_d], [𝔭₂]} iff h(−4d) = 2 iff d ∈ {5, 13, 37}.
4. **(Three grids agree.)** For d ∈ {5, 13, 37} and an odd prime q ∤ d that splits in K:
   - q = x² + dy² ⟺ q ≡ 1 (mod 4) ⟺ q = a² + b²;
   - otherwise q = 2x² + 2xy + ((1 + d)/2)y².

   Among squarefree d ≡ 1 (mod 4), this agreement between the Gaussian grid and the first grid holds exactly for
   d ∈ {1, 5, 13, 37}.
5. **(One grid per ramified norm.)** If every prime factor of n ramifies in K, then n is the norm of exactly one ideal,
   and n lies on exactly one of the two grids: the first iff the exponent of 2 in n is even (d ∈ {5, 13, 37}). So in
   Z[√−5]:
   - 5, 25, 100 lie on the first grid;
   - 10, 50 lie on the second grid;
   - likewise 13, 52 and 37, 148 (first) and 26 and 74 (second).

**Proof.**
- *(1)* ((1 + √−d)/2)² = (1 − d)/4 + ½√−d. This lies in G_d iff (1 − d)/4 ∈ ½ + Z, iff d ≡ 3 (mod 4). Then
  G_d = Z[(1 + √−d)/2] = O_K. If d ≡ 1 (mod 4), the product leaves G_d, so G_d is not a ring.
- *(2)* The Z-span of 2, 1 + √−d, 2√−d and (1 + √−d)√−d = −d + √−d is the checkerboard, since d is odd. It is
  closed under multiplication by √−d, because (u, v) ↦ (−dv, u) preserves u ≡ v (mod 2). So it is the ideal
  (2, 1 + √−d), of index 2.
  - 𝔭₂² = (4, 2 + 2√−d, (1 − d) + 2√−d) = 2·(2, 1 + √−d, (1 − d)/2 + √−d). The last ideal contains (1 + d)/2, which
    is odd, and 2, so it is O_K. Hence 𝔭₂² = (2).
  - 𝔭₂ is principal iff x² + dy² = 2 is solvable, iff d = 1.
- *(3)* By (2), [𝔭₂] has order 2 for d > 1, so h(−4d) ≥ 2, with equality iff Cl = {1, [𝔭₂]}. For squarefree d ≡ 1
  (mod 4), −4d is a fundamental discriminant. The class-number-2 list contains exactly −20, −52 and −148 of this form.
- *(4)* For prime d ∈ {5, 13, 37}, the genus field K(i) equals H(K) (THM-4566; h = 2). A split prime 𝔮 is principal iff
  it splits in K(i), iff q ≡ 1 (mod 4). The non-principal class is [𝔭₂], whose norm form is the second form.
  - Conversely, suppose the agreement holds for squarefree d ≡ 1 (mod 4). Then the primes splitting completely in H(K)
    and in K(i) coincide, so H(K) = K(i) by Bauer's theorem, hence h ≤ 2 and d ∈ {1, 5, 13, 37}.
- *(5)* Ramified primes have a unique prime ideal above them, so an ideal of norm n = 2^a·d^b is 𝔭₂^a·(√−d)^b. Its class
  is [𝔭₂]^a. ∎
