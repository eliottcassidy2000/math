---
id: THM-5242
title: "The triangle law and the arithmetic of 2 and 3. (A) A cyclic join cyc(A, B, C) of three blocks has H = T(v_A, v_B, v_C), a symmetric trilinear form in the port vectors v_X = (P_1(X), P_2(X), ...) (P_k = ordered k-tuples of disjoint nonempty paths covering X): sum of 3 P_k P_k P_k plus the rotations of P_(k+1) P_k P_k and P_(k+1) P_(k+1) P_k. With C fixed, H = v_A^T J(v_C) v_B, where J(p) is tridiagonal with J_aa = p_(a-1) + 3 p_a + p_(a+1) and J_(a,a+1) = p_a + p_(a+1). With C a point, J = [[3,1],[1,1]], the Z[sqrt -2] inner product. (B) A diagonal join cyc(A, A', *) has H = N(beta_A) = (Q - h)^2 + 2 h^2 = 3 h^2 + 8 HC(A) (mod 16) at odd order and 3 h^2 + 16 c_(n-1)(A) (mod 32) at even order; at odd order HC(A) has the parity of the closure number of A's node in the tiling metagraph. (C) For a prime p = 3 (mod 4) the following are equivalent: p = x^2 + 2y^2; p = 3 (mod 8); halving x -> x/2 reverses the Paley tournament P_p; -1 lies in <2> mod p, i.e. the mod-p tables of 3x+1 and 3x-1 alternate (THM-4792). Every prime diagonal-join count is such a prime. (D) In Q_2^x the closure of the Collatz group <2, 3> is the norm group of Q_2(sqrt -2), and negation (which exchanges 3x+1 and 3x-1) is the non-norm class. Every cycle denominator 2^K - 3^k with K >= 3 is a non-norm, so among the cycles of 3x+1 and 3x-1 only the trivial cycle 1 -> 4 -> 2 -> 1 and its mirror have norm-class denominators."
status: >
  (A) PROVED (path decomposition); FINITE-EXACT check: 400 random triples of blocks of orders 1-4 against subset DP,
  0 failures, and the symmetry H(cyc(A,B,C)) = H(cyc(B,A,C)) on the same 400.
  (B) PROVED from the cycle congruences of chanpak_pair_line_census_20261115 (PAIR_LINE_STRUCTURE.md s.4, via the
  Grinberg-Stanley odd-cycle expansion) and THM-4729 (2); FINITE-EXACT: 300 random tournaments of orders 3-10, 0 failures.
  (C) PROVED (second supplementary law, Fermat-Euler, THM-4792 (c)); FINITE-EXACT: all 332,397 primes p = 3 mod 4 with
  7 <= p < 10^7, 0 exceptions (the representation x^2 + 2y^2 checked directly below 10^6).
  (D) PROVED (local class field theory for a ramified quadratic extension, elementary here); FINITE-EXACT: <3> mod 2^j
  equals the unit norms mod 2^j for 3 <= j <= 16.
  KNOWN ingredients are cited in place. The Z[sqrt -2] inner product itself is THM-5172 (C) (mac-mini-2026-11-14-hosts)
  and THM-5222; the cycle congruences for Q are the pair-line census packet's; the alternation classification is
  THM-4792 (c) with Moree's theorem. New here: the three-block law (A), the mod-16/32 refinement and its metagraph
  reading (B), the four-way equivalence for Paley primes (C) and the 2-adic norm-group reading of the Collatz group (D).
  Not audited.
session: mac-mini-2026-11-15-synthesis (wave 38 of chanpak_waves_20261113)
date: 2026-11-15
source: 05-knowledge/results/chanpak_waves_20261113/README.md (wave 38) and GRAND_SYNTHESIS.md
scripts:
  - 05-knowledge/results/chanpak_waves_20261113/scripts/w38_triangle_form.py (A)
  - 05-knowledge/results/chanpak_waves_20261113/scripts/w38_norm_bridges.py (B), (C), (D), Hasse density, Paley data,
    diagonal-join prime census
  - 05-knowledge/results/chanpak_waves_20261113/scripts/w38_hq.c (exact H and Q = Pi_2 by subset DP)
related:
  - THM-5172 (C), THM-5222 (cyclic join with a point), THM-4792 (two sheets mod odd m), THM-4849 (circulants are
    rational Collatz cycles), THM-4729 / THM-4953 (tiling metagraph), THM-4666 (Paley ranking excess and h(-p)),
    THM-4747 (Paley congruences), HYP-9792
---

# THM-5242 - the triangle law and the arithmetic of 2 and 3

**Notation.**
- For a tournament X and k ≥ 1, P_k(X) is the number of ordered k-tuples of vertex-disjoint nonempty paths covering X.
  So P_1 = H. The port vector is v_X = (P_1(X), P_2(X), …). P_0(X) = 0 for nonempty X.
- Q(X) = Π₂(X) = P_2 + 2H counts ordered pairs of possibly empty disjoint paths covering X.
  In the Grinberg–Stanley odd-cycle expansion f_X(x) = Σ_C 2^|C| x^(n−d(C)), one has H = f_X(1) and Q = f_X(2).
- cyc(A, B, C) is the tournament on A ⊔ B ⊔ C with every arc A → B, B → C and C → A. A *point* is a one-vertex block.

## (A) The triangle law

**Statement.** Call a triple (a, b, c) of positive integers *admissible* if, for some k ≥ 1, it is
- (k, k, k), with multiplicity 3;
- a rotation of (k+1, k, k), with multiplicity 1; or
- a rotation of (k+1, k+1, k), with multiplicity 1.

Then

    H(cyc(A, B, C)) = T(v_A, v_B, v_C) := Σ_admissible mult(a, b, c) · P_a(A) P_b(B) P_c(C).

**Proof.**
1. *Block runs.* Cut a Hamiltonian path of cyc(A, B, C) into maximal runs inside one block. The only arcs between
   blocks are A → B, B → C and C → A. So consecutive runs lie in cyclically consecutive blocks, and the run sequence is
   s, s+1, s+2, … (mod 3), determined by its start block s and its length L.
2. *Free combination.* The runs inside one block X, in path order, form an ordered tuple of disjoint nonempty paths
   covering X. Conversely, any choice of such tuples for the three blocks, with the right counts, joins into a
   Hamiltonian path, because every arc the joins need is present.
3. *Visit counts.* Every block is nonempty, so it is visited at least once.
   - If L = 3k, all counts equal k, for each of the 3 start blocks: multiplicity 3.
   - If L = 3k + 1, the start block gets k + 1: each rotation of (k+1, k, k) comes from exactly one s.
   - If L = 3k + 2, the start block and the next one get k + 1: each rotation of (k+1, k+1, k) comes from exactly one s,
     because in a 3-cycle every pair of blocks is consecutive in exactly one order.

   ∎

**Consequences.**
- *Symmetry.* The admissible set is invariant under all permutations of (a, b, c), so T is a symmetric trilinear form.
  In particular H(cyc(A, B, C)) = H(cyc(B, A, C)): reversing the cyclic order does not change H. The same symmetry
  follows from H(T^op) = H(T) and P_k(X^op) = P_k(X).
- *Jacobi matrix.* Fix C with port vector p. Admissible pairs (a, b) satisfy |a − b| ≤ 1, so
  H = v_A^T J(p) v_B with J(p) tridiagonal:

      J_aa = p_(a−1) + 3 p_a + p_(a+1),      J_(a,a+1) = J_(a+1,a) = p_a + p_(a+1).

  Equivalently, J(p) = Σ_c p_c E_c, where E_c = (e_(c−1) + e_c)(e_(c−1) + e_c)^T + (e_c + e_(c+1))(e_c + e_(c+1))^T
  + e_c e_c^T.
- *Specialisations.*
  - **C a point**: J = [[3, 1], [1, 1]] (det 2), so H = 3h_Ah_B + h_AP_B + P_Ah_B + P_AP_B
    = (h_A + P_A)(h_B + P_B) + 2h_Ah_B. This is Re(β_A β̄_B) with β = (H + P₂) + H√−2, the Z[√−2] inner product of
    THM-5172 (C) and THM-5222.
  - **B and C points**: T(v_A, e_1, e_1) = 3H(A) + P₂(A) = H + Q, the triangle reader of chanpak_pair_line_census.
  - **C an arc** (TT₂, ports (1, 2)): J = [[5, 3, 0], [3, 7, 2], [0, 2, 2]], det 32.
  - **C of order 3 or 4**: the determinants are 11232 (TT₃), 22032 (C₃), and 146792448 to 472799232 at order 4. They
    are not powers of 2 in general.
- *The cubic.* For blocks with P_3 = 0 (orders 1 and 2), T(v, v, v) = 3(P_1 + P_2)(P_1² + P_2²). For example
  H(cyc(TT₂, TT₂, TT₂)) = 3 · 3 · 5 = 45, the maximum at order 6.

**Reading.** Each composition law of tournaments is a multilinear form in port vectors:
- ordinal sums multiply the slot polynomials pointwise (f_(A⇒B) = f_A f_B, a split torus);
- a point-cyclic join is the definite binary form of discriminant −8;
- a three-block cyclic join is a symmetric trilinear Jacobi form.

The definite form comes with positivity Q ≥ 2H (wave 35). This is why cyclic joins carry the anisotropic Z[√−2]
arithmetic and no Pell (split) arithmetic.

## (B) Diagonal joins read the metagraph mod 16

Let A have order n, H(A) = h and Q(A) = Q. Let A′ ≅ A and β_A = (Q − h) + h√−2. Then

    H(cyc(A, A′, ·)) = N(β_A) = (Q − h)² + 2h².

**Statement.**
- odd n: N(β_A) ≡ 3h² + 8·HC(A) (mod 16), where HC counts directed Hamiltonian cycles;
- even n ≥ 4: N(β_A) ≡ 3h² + 16·c_(n−1)(A) (mod 32), where c_(n−1) counts directed (n−1)-cycles.

In both cases N ≡ 3h² ≡ 3 (mod 8), recovering THM-5172 (C), and 3h² ≡ 3 or 11 (mod 16) according to h ≡ ±1 or
±3 (mod 8).

**Proof.**
- *Odd n.* By the pair-line census (PAIR_LINE_STRUCTURE.md §4), Q ≡ 4·HC (mod 8). So Q − h ≡ 4HC − h (mod 8), and
  (Q − h)² ≡ h² − 8HC·h ≡ h² + 8HC (mod 16), since h is odd.
- *Even n.* Q ≡ 8c_(n−1) (mod 16), so (Q − h)² ≡ h² − 16c_(n−1)h ≡ h² + 16c_(n−1) (mod 32). ∎

**The metagraph reading.** In the tiling metagraph (THM-4729 (2)), a class A has n·HC(A)/|Aut A| closed tilings.
|Aut A| is odd for every tournament, so at odd n this closure number has the parity of HC(A). Hence:
- **mod 8**, the norm sees only Rédei's parity. N ≡ 3 because y = H is odd (Rédei) and x = Q − h is odd (P₂ is even,
  since its terms pair off as S ↔ S^c);
- **mod 16**, at odd order, the norm sees whether A's metagraph node has an odd number of closed tilings.

## (C) Paley primes, Collatz tables and x² + 2y²

**Statement.** For a prime p ≡ 3 (mod 4), with Paley tournament P_p on F_p (i → j iff j − i is a nonzero square),
the following are equivalent:
1. p = x² + 2y² for integers x, y;
2. p ≡ 3 (mod 8);
3. the halving map x ↦ x/2 is an anti-automorphism of P_p (it reverses every arc);
4. −1 ∈ ⟨2⟩ mod p, i.e. the edge tables of 3x+1 and 3x−1 mod p alternate under halving (THM-4792 (c)).

When they hold, the half-period t with 2^t ≡ −1 is odd, and multiplication by −2 is an automorphism of P_p of odd
order. When they fail (p ≡ 7 mod 8):
- halving is an automorphism of P_p;
- the tables never alternate;
- p = x² − 2y² instead.

Every prime that is a diagonal-join count is ≡ 3 (mod 8) by (B), so it satisfies 1–4.

**Proof.**
- *Multipliers.* For a ∈ F_p^×, multiplication by a sends the arc x → y to ax → ay. Since a(y − x) is a square iff a
  is, a preserves P_p iff a is a square and reverses it otherwise. Halving and doubling lie in the same square class.
  So 3 ⟺ (2/p) = −1 ⟺ p ≡ 3, 5 (mod 8) (second supplementary law) ⟺ 2 (since p ≡ 3 mod 4).
- *Fermat–Euler.* An odd prime is x² + 2y² iff p ≡ 1, 3 (mod 8). So 1 ⟺ 2 for p ≡ 3 (mod 4).
- *Alternation.* The squares form the subgroup of odd order (p − 1)/2.
  - If 2 is a square, ord_p(2) is odd and −1 ∉ ⟨2⟩.
  - If not, 2^((p−1)/2) = −1 ∈ ⟨2⟩, ord_p(2) = 2t, and t divides the odd number (p − 1)/2.
  - Then (−2)^t = (−1)^t 2^t = 1.

  So 4 ⟺ 3. ∎

**Density (KNOWN).** By THM-4792 (c), the alternating primes are those with ord_p(2) even. By Hasse (Math. Ann. 166,
1966) this set has density 17/24. The count below 10^7 gives 0.70817.

## (D) The Collatz group is the 2-adic norm group of Q₂(√−2)

**Statement.**
1. *Norm group.* In Q₂^×, the closure of the group ⟨2, 3⟩ generated by the two Collatz multipliers equals the norm
   group N(Q₂(√−2)^×) = 2^Z × {u ∈ Z₂^× : u ≡ 1, 3 (mod 8)}. Negation (−1 ≡ 7 mod 8), which conjugates 3x+1 into
   3x−1, represents the nontrivial class of Q₂^×/N. Equivalently, the Hilbert symbol (−1, −2)₂ = −1 while
   (2, −2)₂ = (3, −2)₂ = 1.
2. *Cycle denominators.* For K ≥ 3, every cycle denominator D = 2^K − 3^k satisfies D ≡ −3^k ≡ 5 or 7 (mod 8), so D is
   a 2-adic non-norm. Consider primitive cycles with at least one odd step. One with K ≤ 2 halvings has k ≤ K odd steps,
   and the only such cycles are 1 → 4 → 2 → 1 and −1 → −2 → −1 of 3x+1 and their mirrors under negation in 3x−1. Of
   these, the trivial cycle 1 → 4 → 2 → 1 (D = 1) and its mirror −1 → −4 → −2 → −1 are the only cycles of either map
   whose denominator is a norm. The known 3x−1 cycles {5, 7} (D = −1) and {17, …} (D = 2^11 − 3^7 = −139), and every
   hypothetical nontrivial 3x+1 cycle, have non-norm denominators.
3. *Circulant tournaments.* The circulant-tournament cycles of THM-4849, with denominator 2^q − 3^k and q = 2k + 1,
   are ≡ 5 (mod 8) at Paley primes q ≡ 3 (mod 4) and ≡ 7 (mod 8) at q ≡ 1 (mod 4).

**Proof.**
- *Generators are norms.* 2 = N(√−2) and 3 = N(1 + √−2).
- *Unit norms.* A unit norm x² + 2y² has x odd, so it is ≡ 1 or 3 (mod 8). Conversely, u ≡ 1 (mod 8) is a square,
  and u ≡ 3 (mod 8) is 3 times a square.
- *Index.* So the unit norms are exactly the classes 1, 3 mod 8. The group 2^Z × {1, 3 mod 8} already has index 2 in
  Q₂^×, which is the index of any quadratic norm group, so it is all of N.
- *Closure.* ord_(2^j)(3) = 2^(j−2) for j ≥ 3, and 3^k ≡ 1, 3 (mod 8). So the closure of 3^Z is exactly
  {u ≡ 1, 3 (mod 8)}.
- *Denominators.* For K ≥ 3, D ≡ −3^k (mod 8), which is 5 or 7. With K ≤ 2, each odd step is followed by at least one
  halving, so k ≤ K. Solving the k = 1, 2 cases leaves only the cycles named. ∎

## Dictionary (what (A)–(D) say together)

| | tournament side | Collatz side |
|---|---|---|
| parity forced | H odd (Rédei) | 3x ± 1 even for odd x |
| the number 3 | H(C₃) = 3 = N(1 + √−2) = N(β(point)) | the multiplier 3 = N(1 + √−2) |
| the number 2 | Q(point) = 2 = N(√−2) | the halving 2 = N(√−2) |
| the character χ₋₈ (n ≡ 1, 3 mod 8) | every diagonal join is a norm, ≡ 3 mod 8; prime ones are Paley primes where halving reverses P_p | the closure of ⟨2, 3⟩ is ker χ₋₈ on Q₂^×; 3x+1 ↔ 3x−1 is the non-norm class |
| next 2-adic digit | N(β) mod 16 reads the closure number of A's metagraph node | v₂(3x + 1) is read from x mod 8 |
| non-norms | — | every cycle denominator 2^K − 3^k, K ≥ 3 |
