# Golden lane: eleven squares, base φ, 105/223/233/332/425, and Collatz

Session mac-mini-2026-10-07-golden, golden reader. The repo was read only. Scripts and outputs are in `scratchpad/golden/golden/`.

## 0. Bottom line

**The Collatz conjecture remains OPEN. Nothing here moves it.** What is new is exact structure, plus several numerology verdicts.

1. **The Collatz analogue of Theorem D (PROVED).**
   - For S ≤ F_q^×, the group {x ↦ sx + c} is pair-transitive iff S ∪ −S = F_q^×. It is pair-regular iff S is the squares and q ≡ 3 (mod 4). It is never transitive on a composite Z/N.
   - The least regular primes are 7 for ⟨x+1, 2x⟩, 11 for ⟨x+1, 3x⟩, and 23 for all four groups.
   - Every safe prime q ≥ 11 is regular for ⟨x+1, 3x⟩.
   - The −17 cycle's denominator 139 is regular for ⟨x+1, (3/2)x⟩.
   - THM-4581's relation group ⟨x+1, 3x⟩ is regular on odd-difference pairs of Z/2^K for every K ≥ 3.
2. **The paper's R_n are the golden reading of Collatz's rational T-cycles (PROVED + FINITE-EXACT).** The eleven cells are 0, the −5 cycle (−H) and the 1/13 cycle (H). G_5 is BS(1,3) reduced mod 11.
3. **The merged tree is 27's trunk.**
   - 223 meets 233 by THM-4581's coalescence of the consecutive pair (334, 335).
   - The golden features on the tree are NUMEROLOGY (null p = 0.15).
4. **Base-φ binary readings: killed.** No target can arise, because each contains "11" in binary.
5. **A golden pair chain works verbatim.** Its Lyapunov weight degenerates, and its merge tail is T^(−1/2).

## 1. Sources

**The paper.**
- R_n = Z[φ]/(φ^n − 1) has |R_n| = L_n − 1 − (−1)^n.
- Theorem D: G_n = ⟨x+1, φx⟩ is pair-transitive only for n = 3, 4, 5 (F_4, F_5, F_11), and regular only for n = 5 (order 55). The proof uses the size bound |G_n| ≤ n|R_n| and the fact −1 ∉ H = {1, 3, 4, 5, 9}.
- On O/(11) = F_11², the fifth return acts by (α, β) ↦ (α, −β). The map Ψ_−[α, β] = {α ± β^(−1)} sends its 55 orbits to the 55 pairs of F_11, which split into eleven matchings M_α.
- §7: the shortest modular geodesic, of length 4 log φ, from Q² (trace 3).
- §8: an S_8 octic, already filed as GL(7) with no Collatz content (square11 note).

**Prior repo work.**
- THM-4528: Θ semi-conjugates T to the golden β-map.
- golden_carriers §4: denominators 1, 2, 11 and 76; the 1/13 cycle shares the −5 cycle's ideal (11, φ − 4); φ^18 − 1 = 76φ^9.
- The 2026-10-03 trio (digit_carry, prime_clocks, route_compiler) is the only prior "binary in base φ" work. It typed three readings and computed:
  - 105 = 1001010101.0101001001_φ;
  - the clock O/(105) of order 80;
  - Φ_105 as the first nonflat cyclotomic polynomial;
  - L_23 = 139·461.

  It never read the expansions of primes as integers.
- THM-4489 (223, the Fermat sextic) and the triplet and branch-toll notes (the meetings at 425, the family 223 + 62208t / 233 + 65536t).
- THM-4553: B-invariant tournaments of the PSL(2,7) Borel are only P_7 and its reverse. That is pair-regularity at 7.

## 2. Derivations

### 2.1 Classification (PROVED; brute force FINITE-EXACT)

**Proof.**
- g(x) = sx + c sends {0, 1} to {c, c + s}. It hits {x, y} iff ±(y − x) ∈ S.
- The stabilizer of {0, 1} is {id} together with x ↦ 1 − x when −1 ∈ S. So the action is free iff −1 ∉ S.
- Free plus transitive forces |S| = (q − 1)/2, so S is the squares and q ≡ 3 (mod 4).
- In characteristic 2, x + 1 swaps 0 and 1, so the action is never free. This is the paper's n = 3 case.
- On composite N, a unit multiplier preserves whether a difference is a unit.

**Checks.**
- Orbit brute force on all primes from 5 to 160 and all four groups: 0 mismatches.
- Composites below 130: no transitive case.
- The paper's table for n = 3..12 is reproduced.

| Group | Regular primes | Transitive, stabilizer 2 |
|---|---|---|
| ⟨x+1, 2x⟩ | 7, 23, 47, 71, 79, 103 | 5, 11, 13, 19, 29, … |
| ⟨x+1, 3x⟩ | 11, 23, 47, 59, 71, 83 | 5, 7, 17, 19, 29, … |
| ⟨x+1, 2x, 3x⟩ | 23, 47, 71, 167, 191 | 5, 7, 11, 13, 17, … |
| ⟨x+1, (3/2)x⟩ | 23, 43, 47, 67, **139** | 7, 11, 17, 31, … |

**Consequences.**
- **Safe primes (PROVED).** If q = 2p + 1 ≥ 11, the squares have prime order and 3 is a square, so ⟨x+1, 3x⟩ is regular. ⟨x+1, 2x⟩ is regular iff q ≡ 7 (mod 8).
  - Infinitely many regular primes is CONDITIONAL (Sophie Germain, or Artin-type under GRH).
  - Primes below 10^5 for ⟨2,3⟩: 82.5% transitive, 12.3% regular.
- **Golden versus Collatz.** The golden family has one multiplier, whose order divides n. That bounds |G_n| and leaves finitely many transitive n. Collatz has two independent multipliers and no such bound.
- **Cycle denominators.**
  - 139 = 2^11 − 3^7: ⟨3/2⟩ is regular; the others are transitive.
  - 3299 = 3^9 − 2^14: ⟨3⟩ is regular.
  - 233: transitive for ⟨2,3⟩, ⟨3⟩ and ⟨3/2⟩; not for ⟨2⟩, which has index 8.
  - 1631, 6487 and 7153 are composite. Their factors 7, 23 and 311 are regular, and 499 carries none of the four groups.
- **2-adic version (PROVED).** The closure of ⟨3⟩ in Z_2^× is {u ≡ 1, 3 (mod 8)}. It has index 2 and does not contain −1, which is the paper's structure. So Z_2 ⋊ 3^(Z_2) acts simply transitively on odd-difference pairs.
  - FINITE-EXACT: regular for K = 3..11. When v_2(difference) = j the stabilizer has order 2^j.
- **p-adic lifts (FINITE-EXACT).**
  - ⟨3⟩ fails to be regular modulo 121, because 3^5 ≡ 1 (mod 121) (11 is a base-3 Wieferich prime).
  - Golden φ ↦ 37 (order 55) stays regular modulo 121.
  - 2 and 3 stay regular modulo 529.

### 2.2 Borel subgroups, Galois, THM-4581, and 23

- **G_5 is the reduction of BS(1,3) (PROVED, trivial).** BS(1,3) = ⟨x+1, 3x⟩, the relation group of THM-4581, reduces mod 11 to F_11 ⋊ ⟨3⟩ = G_5, since ⟨3⟩ = ⟨4⟩ = H.
  - Every regular case is the Borel subgroup of PSL(2,p) acting simply transitively on the arcs of the Paley tournament QR_p. THM-4553 is the case p = 7.
- **Galois (DICTIONARY, KNOWN).** The degree-7 and degree-11 actions of PSL(2,p) (point stabilizers S_4 and A_5) restrict on the Borel to the affine action, which is pair-regular.
  - At 5 (point stabilizer A_4) the Borel D_5 has two orbits on pairs.
  - Halving realizes 7 and tripling realizes 11.
  - 23:11 lies in M_23, as 11:5 lies in M_11.
- **QR codes (KNOWN).** ⟨ℓ⟩ = squares mod p iff x^p − 1 = (x − 1)·g·g̃ over F_ℓ. The perfect cases:
  - 7: 2^3 = 1 + 7 (Hamming code);
  - 11: 3^5 = 1 + 2·11 + 4·55 (ternary Golay code);
  - 23: 2^11 = 1 + 23 + 253 + 1771 (binary Golay code).

  The ternary identity gives 3^5 − 1 = 2·11², which is the Wieferich obstruction.
- **THM-4581's index 1 has no mod-p shadow (ANALOGY).**
  - Pair-transitivity of Γ_C mod p fails for about 18% of primes. For example, ⟨−1, 2, 3⟩ is the squares mod 73.
  - The exact shadow is the 2-adic regularity of §2.1. It is STRUCTURAL but says nothing about merging.
- **Verdict on 23 and 253: NUMEROLOGY.** The exact consequences are generic for every ⟨2,3⟩-regular prime:
  - 23 never divides 2^a + 3^b;
  - 23 divides 2^A − 3^r iff A ≡ 8r (mod 11);
  - the image of a word mod 23 is determined by its effect on one pair.

### 2.3 Golden reading of Collatz cycles (PROVED + FINITE-EXACT)

**The map κ_n.** For x with T-period word w of length n, set κ_n(x) = Σ_j w_j φ^(n−1−j) = (φ^n − 1)Θ(x), read in R_n.
- **Equivariance (PROVED).** κ_n(Tx) = φ·κ_n(x), because the digit term dies mod (φ^n − 1).
- **Bijectivity (FINITE-EXACT, n ≤ 22).** κ_n is a bijection from the L_n periodic points onto R_n, except that {0, −1, −2} collapse to 0 when n is even.

| Cycle | Word (n) | 2^A − 3^r | Golden R_n | Pair action |
|---|---|---|---|---|
| 1, 4, 2 | 100 (3) | 1 | F_4 | transitive |
| 1/5 | 1000 (4) | 5 | F_5 | both readings: transitive |
| −5 | 10100 (5) | −1 | F_11, labels −H | **regular** |
| 1/13 | 10000 (5) | 13 | F_11, labels H | golden regular; Z/13 transitive |
| −17 | 18 letters | −139 | O/(76) | Z/139: ⟨3/2⟩ **regular**; golden not transitive |

- **n = 5.** The 1/13 cycle has labels 3, 1, 4, 5, 9, and the −5 cycle has labels 8, 10, 7, 6, 2. M_0 is the in-phase pairing of T^j(1/13) with T^j(−5).
- **n = 10.**
  - The β = 0 cells are the period-5 points (with α = 2·label).
  - T^5 acts as (α, β) ↦ (α, −β) on all 110 primitive period-10 rational points. So the paper's 55 exchanged pairs are these points up to T^5.
  - Ψ_− sends the 10-cycle with denominator 101 = 2^7 − 3^3 onto M_0.
- **n = 18.** R_18 ⊃ F_19². Under φ ↦ 5 the −17 orbit folds to period 9.
- **Reading (DICTIONARY).** Theorem D's transitive cases are the T-periods 3, 4 and 5. Its regular case is the golden image of the −5 cycle together with its twin 1/13.

### 2.4 The merged tree (FINITE-EXACT)

**Words.**
- 83 → 233 has word (1,3,1,1,1,2,2,1,2,1,1), so T^18(332) = 233 with 11 odd steps.
- 233 → 377 → 425 has words (2,1,1,1,2,3,1,1) then (2,1): 15 Terras steps, 10 odd.
- 223 → 425 has word (1,1,1,1,3): 7 Terras steps, 5 odd.

**27's trunk.**
- 332 joins 27's orbit at 94, and 223 joins it at 850.
- All six features lie on 27's trajectory: 47 = L_8, 322 = L_12, 242, 121, 233 and 377.
- 38.5% of n in [200, 500] pass through 9232.

**Pair-chain form.**
- T(223) = 335 and T^9(233) = 334. So the meeting is the start (0, 1) at v = 334, merging at Terras time 6 at 425.
  - 18 of the 64 residues mod 64 merge with v + 1 by time 6 (THM-4581 6(c)).
- In equal-time form, (223, T^8(233) = 668) has (k_0, e_0) = (−1, 1/3), i.e. v = 3u − 1.
- 332 → 233 is a pure clock shift.

**Families.**
- The prior family meets at 425 + 118098t; verified for t < 3000.
- In Terras time (223 + 31104s, 233 + 32768s) suffices, half the period; verified for s < 3000. This is the class v ≡ 334 (mod 64) pulled back along T and T^9.

**2-adic neighbours.**
- 105 ≡ 233 ≡ 425 ≡ 41 and 223 ≡ 31 (mod 64), both residues being points of 27's trajectory.
- 233 = 105 + 2^7 gives a run at k = 0, so T^7(233) − T^7(105) = 3^5. This is the only exact 3^5.

**Base rates (NUMEROLOGY).**
- Fraction of orbits of n in [200, 500] containing each value: 233: 0.30; 377: 0.34; 47: 0.15; 322: 0.18; 121 and 242: 0.18.
- Two random ratio-φ sequences hit the trunk at least 4 times with probability 0.153. 121 and 242 were selected after the fact.
- (233, 377) is the only pair (F_k, F_(k+1)) with F_(k+1) on the orbit of F_k for k = 9..19. The base rate is 0.036 per pair.
- Ellison's 7·233: I found no map to the tree.

**Types.**
- STRUCTURAL: the (334, 335) coalescence and the 2-adic neighbours.
- NUMEROLOGY: F_13, F_14, L_8, L_12, 11² and 3^5 − 1.

### 2.5 105

**Typed roles.**
- Prior work: O/(105), Φ_105, and the base-φ word.
- 105 is one of the 19 HCF planes, Q(√3, √5, √7) with χ ≥ 4 (THM-4566).
- 210 is a Borwein–Choi exception.
- These are class-group facts (DICTIONARY), with no Collatz map.

**Orbit.** The odd values are 105, 79, 119, 179, 269, 101, 19, 29, 11, 17, 13, 5, 1.
- Base rates over n in [80, 130]: 101 appears in 0.12 of orbits, 19 in 0.14, 29 in 0.20, and 11, 17 and 13 in about 0.4.
- 105 meets 233 only on the cycle {1, 2}.
- Verdict: NUMEROLOGY, apart from 105 ≡ 233 (mod 2^7).

### 2.6 Base-φ binary readings

Readings: integer part, full string, and reversed; computed by exact Bergman expansion.
- **The targets cannot arise (PROVED).** 105, 223, 233, 332 and 425 all contain "11" in binary, and normal forms never do. Confirmed for every prime up to 500.
- **Lengths (FINITE-EXACT).** Integer and fractional lengths differ by 0 or 1.
- **The full reading F is ≡ 1 (mod 4).** So its stopping time is always 3.
- **Parity words: killed.** Agreement with p's T-parity word is 0.593, against 0.595 for a random q.
- **L_(2k) maps to 2^(4k) + 1 (PROVED, NUMEROLOGY).** So 47 = L_8 maps to 65537.
- **Primality bias (NUMERICAL).** F is prime more often for prime n: 0.305 against 0.149 for n ≤ 500, shrinking to 0.095 against 0.077 for n ≤ 2·10^4.
  - Mechanism: q | n raises P(q | F(n)), most at q = 11.
  - No Collatz content.

### 2.7 Golden pair chain; the converse

**The map.** On Z_2[φ] (residues in F_4), define G(x) = x/2 on 2O and G(x) = (φx + ρ)/2 otherwise, with ρ ≡ φx (mod 2).
- G is Haar-preserving and 4-to-1.
- **Chain.** The relation u = φ^k v + e follows THM-4581's table, with e staying in O.
  - 0 mismatches against direct orbits over 6·10^4 steps.
  - k is a lazy fair random walk.
- **Weight.** Every normalized branch contracts, by 1/2 or φ/2, in both embeddings. So |f|^θ suffices with no s^|k|, and a.s. merging holds (PROVED at sketch level).
  - BFS absorbs all 624 small states (FINITE-EXACT).
  - √T·q(T) ≈ 3.8, i.e. a T^(−1/2) tail (NUMERICAL).
- **Collatz contrast (PROVED, algebra).** For (mx+1)/2, THM-4581's balance ½((m/2)^θ/s + s/2^θ) = 1 is solvable iff m < 4, the condition for negative Terras drift. So m = 3 is the only nontrivial odd case.

**Converse.** κ_5 identifies the regular F_11 with Collatz's period-5 points. ⟨x+1, 3x⟩ on F_11 is the same group as G_5. Z/139 carries a regular ⟨x+1, (3/2)x⟩.

### 2.8 Numerology ledger

- 7153 = 23·311, both regular (the base rate is 12% per prime).
- 11 = 3^3 − 2^4 has golden reading in R_7, not R_5.
- The 10-cycle's 101 against the 101 in 105's orbit.
- F(47) = 65537 against the period 65536.
- |M_23|/253 = 8!.
- 253 is an idoneal plane.

## 3. Suggested canon filings

1. **THM: pair actions of Collatz affine groups.** The criterion, the table, safe primes, 139 and 3299, the 2-adic regularity, and the Wieferich obstruction at 121. PROVED + FINITE-EXACT; `pair_regularity`.
2. **THM: golden reading of rational Collatz cycles equals the eleven-squares R_n.** Equivariance of κ_n; bijectivity for n ≤ 22; the n = 5 and n = 10 identifications; G_5 as BS(1,3) mod 11. PROVED + FINITE-EXACT + DICTIONARY; inherits golden_carriers §4; `golden_cycle_dictionary`.
3. **Results note: 223/233 is the (334, 335) coalescence on 27's trunk.** The halved family and the numerology verdicts with base rates. FINITE-EXACT + NUMEROLOGY; `merged_tree`, `tree_pairchain`, `golden_null`.
4. **Numerology guard: base-φ binary readings.** PROVED (trivial) + NUMERICAL; `phi_binary_primes`, `phi_binary_bias`.
5. **HYP or result: golden pair chain**, plus the m < 4 weight criterion. Sketch + FINITE-EXACT + NUMERICAL; `golden_pair_chain`.
6. **DICTIONARY (KNOWN):** Paley, Galois degrees 7 and 11, Golay, M_11 and M_23. The Collatz link is NUMEROLOGY.

## 4. Files

All files are in `scratchpad/golden/golden/`, each with a matching `.out`:
- `pair_regularity.py`
- `golden_cycle_dictionary.py`
- `merged_tree.py`
- `golden_null.py`
- `tree_pairchain.py`
- `phi_binary_primes.py`
- `phi_binary_bias.py`
- `golden_pair_chain.py`
