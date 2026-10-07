# PLATONIC lane: {2,3,11}, Platonic solids vs Ellison, integer cycles, mod 18/19 vs idoneal counts

Session mac-mini-2026-10-07-golden. The repo was read only. Collatz remains OPEN.

## 0. Bottom line

1. **Integer-cycle census (FINITE-EXACT).** The Terras map has exactly five integer cycles of period L ≤ 301,993: {0}, {−1}, {1,2}, {−5,−7,−10} and the 11-cycle through −17.
   * The periods are exactly {1,2,3,11} in that range.
   * Two independent programs agree for L ≤ 40.
   * This extends THM-4484's census, which stopped at p ≤ 90.
2. **Platonic solids against log₂3 and Ellison's list: NUMEROLOGY.**
   * The base-rate test gives p = 0.78.
   * 19/12 ↔ 12-TET is one inequality (KNOWN).
   * 12-TET ↔ the icosahedron's 12 vertices is NUMEROLOGY, and no equivariant identification exists: no symmetry of the icosahedron has order 12 (PROVED).
3. **The only honest Ellison link is heuristic.**
   * Ellison's five exceptions are exactly the five shapes with L ≥ 12 of largest naive expected cycle count Lyn(L,k)/|2^L−3^k| (FINITE-EXACT for L ≤ 1500).
   * None of them carries an integer cycle.
   * (27,17) carries the two KNOWN 3x+5 cycles; (19,12) carries eight 3x+23 cycles.
4. **{2,3,11} core: one equation, 3⁵ = 1 + 2·11².** It is at once:
   * Ljunggren's square repunit 11111₃ = 11²;
   * the perfectness of the ternary Golay code;
   * the base-3 Wieferich property of 11;
   * the congruence s* = −1/2 ≡ 11² (mod 3⁵).

   Related facts:
   * ord₁₁(3) = 5 makes G₅ = ⟨x+1, 3x⟩ the Borel subgroup of PSL(2,11). This is the coordinate symmetry group of the Golay code.
   * (3|11) = (5|11) = 1 and 11 ≡ 3 (mod 4) are exactly the hypotheses of the χ ≤ 5 bound for m = 165.

   These are STRUCTURAL inside number theory. Every link to the Collatz periods is NUMEROLOGY.
5. **Idoneal slices (PROVED).**
   * Borwein–Choi exceptions = {1,4} ∪ {idoneal n ≡ 2 (mod 4)}. The proof uses Selling parameters and the genus count.
   * The 19 planes = the squarefree idoneal m ≡ 1 (mod 4) (THM-4566).
   * Both are 2-adic slices of one finite list. The counts 18 and 19 have no bridge to ord₁₉(2) = ord₁₉(3) = 18: NUMEROLOGY.
6. **The 9/4 thread.**
   * The owner's restatement matches THM-4568 as corrected: the identity is PROVED, and "traces back to" is ANALOGY.
   * The point −1/2 (equivalently −1) organizes the odd runs of every cycle (Steiner circuits, KNOWN), but it does not select the periods.
   * One striking post-hoc pattern, typed NUMEROLOGY: at the even predecessors of the odd-run starts of the negative cycles, 2x+1 takes the values −19, −67 and −163, all Heegner numbers.

## 1. What the sources say (two sweeps, about 70 threads)

| thread | claim | repo type |
|---|---|---|
| THM-4553 / sixes_and_sevens | PSL(2,7) > Borel F₂₁ > ⟨x↦2x⟩ = Collatz Frobenius on {1,2,4} | PROVED; DICTIONARY for Collatz |
| level11_dessin_golden_torsion | the order-55 Borel of PSL(2,11) preserves Paley P₁₁; the index-11 A₅ is kept separate | PROVED/CITED |
| legacy `info_music_23.out`, `recurrence_web_23.out`, `moat_boundary_synthesis.out` | "12-TET has h(E₆) = 12 notes"; comma = 3^h(E₆)/2^(h(E₇)+1) | **untyped**; MISTAKE-563 covers the method but not these files |
| THM-4484 | free cycles ⇔ \|2^L−3^k\| = 1; −17 is sporadic (gap 139, 1 of 30 necklaces); census p ≤ 90 | PROVED; FINITE-EXACT |
| THM-4512 | Ellison closure x ∈ {13,14,16,19,27} | PROVED+CITED |
| collatz_lucas_monotile_discrepancy | φ enters 11/7 = L₅/L₄ only through the shared continued-fraction prefix | EXACT |
| procgen_brackets_20260924 | "Ljunggren 11² = (3⁵−1)/2 is REAL but unrelated" | NUMEROLOGY |
| THM-4558 / THM-4566 | Moser plane χ = 4; m = 165 has χ = 5 via the prime over 11 | KNOWN/PROVED |
| mod18_mod19_seven_sixtythree | ord₂₇(2) = ord₁₉(2) = 18; Heegner primes on the tower | PROVED; NUMEROLOGY |

**Gaps in the repo:**
* No file ties polyhedra to Ellison's list or to 1631, 3299, 6487 or 5077565.
* The ternary Golay code [11,6,5]₃ appears nowhere; only [12,6,6] does.

## 2. Derivations and computations

### 2.1 Census (task 2)

**Lemma (PROVED).** Around a cycle, Π_{odd x_i}(3 + 1/x_i) = 2^L. Hence:
* a positive cycle has 3^k < 2^L ≤ 4^k and min x ≤ 1/(2^{L/k} − 3);
* a negative cycle has min|x| ≤ 1/(3 − 2^{L/k}).

**(a) Lyndon-word method (as requested).**
* FKM enumeration of the Lyndon words with k ≥ L/2, keeping d ← 3d + 2^t exact.
* The test for a cycle is (2^L − 3^k) | d.
* Every leaf count matches (1/L)Σμ(e)C(L/e,k/e).
* 3.06·10¹⁰ words for L ≤ 40, in 82 s.
* Hits occur only at L = 1, 2, 3 and 11. The L = 11 hit is the word 00011110111, read from −136.

**(b) Minimal-element scan.**
* Every candidate minimum up to the lemma's bound is iterated until it returns, drops, or reaches L steps.
* For L ≤ 301,993 the bounds are 7.22·10⁹ (positive side, at 125743/79335) and 1.03·10¹⁰ (negative side, at 176251/111202).
* The scan took 40 s. Nothing overflowed and nothing hit the step cap.
* It found exactly the minima 1, −1, −5 and −17.

For L ≤ 40 the worst bounds sit at Ellison shapes: 146.8 at (27,17) and 295.3 at (19,12). The method stops at the convergent 301994/190537, followed by the partial quotient 55, where the bound jumps to 9.8·10¹¹.

**FINITE-EXACT.** The periods for L ≤ 301,993 are exactly {1,2,3,11}.
* 1, 2, 3 come from the free shapes, |2^L − 3^k| = 1 (Catalan).
* 11 comes from the sporadic −17 cycle.

The positive side is much weaker than KNOWN results (Simons–de Weger, Hercher). It is new only as a self-contained check and on the negative side.

### 2.2 Ellison shapes and near-cycles (task 1)

**Ranking (FINITE-EXACT, exact for L ≤ 1500).** Let E(L,k) = Lyn(L,k)/|2^L−3^k| for k ≥ L/2. The top five shapes with L ≥ 12 are:
1. (19,12): 0.371
2. (16,10): 0.077
3. (27,17): 0.062
4. (13,8): 0.061
5. (14,9): 0.043

These are Ellison's list. The next shape is (12,7), at 0.035.

**Reading:**
* This is a DICTIONARY statement: both E and Ellison's criterion measure |2^x − 3^y|/2^x.
* For all L it would follow from an effective irrationality measure of log₂3 (CITED, not checked).
* The expected cycle count summed over the five shapes is 0.61. The observed count is 0.

Reduced denominators of c_w/m over all Lyndon words (FINITE-EXACT):

| shape | gap | denominators realised (count) |
|---|---|---|
| (11,7) | −139 | 1: the −17 cycle |
| (13,8) | 1631 = 7·233 | 233 (15); never 7 |
| (14,9) | −3299 | none |
| (16,10) | 6487 = 13·499 | 499 (41); never 13 |
| (19,12) | −7153 = −23·311 | 23 (8): eight 3x+23 cycles; 311 (117) |
| (27,17) | 5077565 = 5·71·14303 | 5 (2): the 3x+5 cycles 187 and 347 (KNOWN, Lagarias 1990), now located at an Ellison shape; 71 (5); … |

Separately, 7153 = 3¹² − 2¹⁹ is the gap of the Pythagorean comma. That is STRUCTURAL, by definition.

### 2.3 Platonic base rates (task 1)

**Test.**
* Statistic: |P ∩ C(α)|, where C(α) is the set of numerators and denominators ≤ 120 of the convergents and semiconvergents of α.
* Null: α ~ U[1.5, 5/3], which keeps the same continued-fraction prefix; 2·10⁵ draws.

| P | observed (log₂3) | null mean | P(≥ obs) |
|---|---|---|---|
| vertex/edge/face counts and group orders {4,6,8,12,20,24,30,48,60,120} | 2 (8, 12) | 2.02 | 0.78 |
| the same plus the Galois primes {5,7,11} | 5 | 4.47 | 0.62 |
| rotation-group orders {12,24,60} | 1 (12) | 0.25 | 0.26 |

**Verdict: NUMEROLOGY.**
* Galois's exceptional actions (p = 5, 7, 11, with stabilizers A₄, S₄, A₅) are STRUCTURAL and KNOWN. But {5,7,11} ⊂ C(α) holds for 60% of random α.
* **12-TET ↔ icosahedron: PROVED non-equivariant.** I_h = A₅ × C₂ has element orders {1,2,3,5,6,10}, so the fifths cycle (order 12) is not a symmetry. Only A₄, acting regularly on the 12 vertices, matches the count, and A₄ ≇ C₁₂.
* **McKay.** The E₈ exponents (the units mod 30) contain every prime in [7,29], so their overlap with 7, 11, 13, 17, 19 is NUMEROLOGY.

### 2.4 {2,3,11} catalogue (task 2)

| entry | real map? | type |
|---|---|---|
| periods {1,2,3}: \|2^L−3^k\| = 1 | inside Collatz | PROVED (THM-4484) |
| period 11: semiconvergent 11/7, gap −139, 1 necklace of 30 (E = 0.22) | inside Collatz; a Diophantine accident | PROVED/FINITE-EXACT |
| **3⁵ = 1+2·11²**: 11111₃ = 11² (Ljunggren) ⟺ Σ_{i≤2}C(11,i)2^i = 3⁵ (ternary Golay perfect) ⟹ 3⁵ ≡ 1 mod 121 ⟺ −1/2 ≡ 121 mod 3⁵ | yes, the same equation | KNOWN + DICTIONARY |
| ord₁₁(3) = 5 ⟹ ⟨x+1,3x⟩ = Borel(PSL(2,11)) = G₅ = the symmetry of the ternary QR (Golay) code, with the multiplier 3 acting as Frobenius; Gleason–Prange: PSL(2,11) acts on [12,6,6]₃, Aut = 2.M₁₂ | yes | KNOWN |
| p ≡ 3 mod 4 ⟹ the Borel (order C(p,2)) is regular on pairs. p = 7, 11, 23 give the perfect QR codes (binary Hamming, ternary Golay, binary Golay), with Collatz multipliers 2, 3, {2,3} | yes | PROVED (checked)/KNOWN |
| φ ≡ 4, 8 mod 11; 11 = N(4−φ); φ¹⁰ − 1 = 11φ⁵ | yes, in Z[φ] | KNOWN |
| Δ(2,3,11) ↠ PSL(2,11); A₅ = stabilizer in Galois's 11-point action | yes | KNOWN |
| m = 165: F = Q(φ)·Q(√3,√11) (golden field times Moser plane). The prime over 11 used for χ ≤ 5 needs (3\|11) = (5\|11) = 1 (residue degree 1) and 11 ≡ 3 mod 4 (inert in F(i)) | yes | PROVED (checks the hypotheses of THM-4558's Lemma R) |
| Moser spindle: sin²θ = 11/36, i.e. 11 = 4·3 − 1 | in the geometry | PROVED |
| j((1+√−11)/2) = −2¹⁵ | CM (Gross–Zagier) | KNOWN |
| 2¹¹ + 3⁷ = 5·7·11² (a 1-in-11 event: q₁₁(2) ≡ 5) | no | NUMEROLOGY |
| (3⁷−1)/2 = 1093 (Wieferich base 2) against k = 7 | no | NUMEROLOGY |
| 11 as the −17 period against 11 as a Golay length, the PSL(2,11) level, or Heegner | no | NUMEROLOGY |

### 2.5 Idoneal slices (task 3)

**Theorem (PROVED; the explicit list of 18 assumes completeness of the idoneal list, CITED via THM-4566).** n ≠ xy+yz+zx for all x, y, z ≥ 1 ⟺ n ∈ {1,4}, or n ≡ 2 mod 4 and n is idoneal.

*Proof.*
1. Odd n ≥ 3 is represented as 1·y + y·1 + 1. If 4 | n and n ≥ 8, take (2, (n−4)/4, 2).
2. Now let n ≡ 2 mod 4. A triple (x,y,z) ≥ 0 with xy+yz+zx = n is the same thing as the Selling parameters of an obtuse superbase of (x+z)X² − 2zXY + (y+z)Y², a form of determinant n. Selling parameters are a class invariant (Conway).
3. So n is an exception ⟺ every form of determinant n is diagonal.
4. All such forms have odd content g. Each primitive part has discriminant −4n′ with n′ = n/g² ≡ 2 mod 4.
5. The number of reduced diagonal forms is 2^{ω(n′)−1} = 2^r. This equals the number of genera, since μ = r+1.
6. Hence every form is diagonal ⟺ h(−4n′) equals the number of genera ⟺ n′ is idoneal. Exponent ≤ 2 descends from n to n′. ∎

**Why these slices.**
* For n ≡ 1 mod 4 there are twice as many genera as diagonal classes, so rhombic classes always give a representation.
* The planes need −4 to be the 2-adic prime discriminant, i.e. i in the genus field.
* Both slices are therefore 2-adic.

**FINITE-EXACT.**
* There are 101 discriminants of exponent ≤ 2 with |D| ≤ 30000: 65 even and 36 odd, the largest 7392.
* The idoneal counts mod 8 match the lead's.
* A brute force over n ≤ 2·10⁵ finds exactly the 18 exceptions.

**Split-prime lemma (PROVED).** In an order of exponent ≤ 2, a split prime p has 𝔭² = (α) with α ∉ Z, so p² ≥ |D|/4. Consequences:
* for idoneal n > 9, n ≢ 2 mod 3 (3 does not split);
* 2 splits only for |D| ∈ {7, 15};
* checked for all odd p < 200.

So 2 and 3 are non-split in every large exponent-2 field. That is the real "{2,3}" content of the idoneal lists.

### 2.6 Mod 18/19 bridges (task 3): NUMEROLOGY

* **The counts.** 18 = 16 + {1,4} and 19 = 22 − {9,25,45} are slice counts of a 65-element list.
* **Primitive roots.** Among the 10 plane primes > 3, both 2 and 3 are primitive only for 5 and 19. The base rate is 23%.
* **The prime 19.** 19 divides 57 and 133, but every prime ≤ 31 divides some idoneal number.
* **Generic DICTIONARY.** On the planes 57 and 133 the 19-genus character equals (−1)^{ind₂ a}, the parity bit of the mod-19 clock. This holds at any prime where 2 is a primitive root.
* **Residues.** The 101 discriminants show no structure mod 18 or mod 19.
* **253.** 253 = C(23,2) = |Borel(PSL(2,23))| = 11·23 is the product of the two Golay lengths. NUMEROLOGY.

### 2.7 The −1/2 thread (task 4)

**The owner's statement.**
* κ* = 4/3 = 2/(1−s*). Here s* = −1/2 is the fixed point of the a = 2 gadget b ↦ 3b+1, and also the limit of the rescaled fixed points −(a−1)/(2a).
* THM-4555's root collision is 3(−1/2) + 1 = −1/2, reached via f₂(−1) = −1/2 and f_c(−1/2) = −2^{−c−1}. PROVED.
* "Saturated via A004991" is THM-4568 (2). PROVED.
* "Traces back to" is ANALOGY: the same map serves unrelated purposes.
* Caution (MISTAKE-583): the fixed point of the halved map is −1.

**Family (PROVED).**
* f_c(x) = (3x+1)/2^c has fixed point 1/(2^c − 3). For c = 0, 1, 2 this is −1/2, −1, 1.
* The integral members are the free cycles; integrality ⟺ 2^c − 3 = ±1, the root identity 3 = 1 + 2.
* In general a word's fixed point is c_w/(2^L − 3^k).

**Organization of the cycles (KNOWN, Steiner circuits).** In y = x+1 an odd run is y ↦ (3/2)^r y (THM-4556).
* −5 is the 1-circuit −(2²+1) ⇒ −(3²+1) = −2(2²+1). This is Catalan's 3² − 2³ = 1.
* −17 is the 2-circuit with runs −(2⁴+1) ⇒ −(3⁴+1) = −2·41 and −(2³·5+1) ⇒ −(3³·5+1) = −8·17.
* So the fixed point organizes the runs, but it does not pick the periods.

**Expansions.**
* −1 = …111₂ and −1/2 = …111₃.
* The truncations (3^j−1)/2 are 1, 4, 13, 40, 121 = 11², 364, 1093.
* Any link from these to the periods is NUMEROLOGY.

**Values relative to −1/2.**
* In z = 2x+1 the −17 cycle is {−33, −49, −73, −109, −163, −81, −121, −181, −271, −135, −67}.
* The odd-run starts of the negative cycles are 5, 17 and 41. These are Euler lucky numbers, so 4m−1 = 19, 67, 163 are Heegner numbers (Rabinowitsch).
* The chance is about 0.19³ ≈ 0.7%, uncorrected and post hoc, with no mechanism. **NUMEROLOGY; do not file as structure.**

## 3. Suggested canon filings

1. **THM, "Integer cycles of the Terras map to period 301,993"** (FINITE-EXACT, with the PROVED bound lemma).
   * Evidence: `cycles/minscan_mt.c`, `minscan_301993.txt`, `necklace.c`, `necklace_*.txt`, `bounds3.py`.
   * It extends THM-4484's p ≤ 90.
2. **Addendum to THM-4566** (PROVED, list CITED-complete): Borwein–Choi exceptions = {1,4} ∪ idoneal(2 mod 4), plus the split-prime lemma.
   * Evidence: `idoneal.py/.out`.
3. **Remark to THM-4512/THM-4484** (FINITE-EXACT, DICTIONARY): Ellison's five = the top-five E(L,k) shapes, with the rational-cycle table, including 3x+5 at (27,17).
   * Evidence: `heuristic*.py`, `nearcyc*.py`.
4. **Note "3⁵ = 1+2·11²"** (KNOWN + DICTIONARY): Ljunggren, ternary Golay, Wieferich-3, s* mod 3⁵; the p = 7, 11, 23 Borel/QR-code dictionary; the m = 165 prime-over-11 hypotheses.
   * Collatz relevance: NUMEROLOGY.
5. **Numerology register** (MISTAKE-559 rule):
   * 12-TET ↔ icosahedron;
   * the untyped `_23.out` h(E₆)/comma claims;
   * the Platonic base rates;
   * E₈ exponents ↔ Ellison/−17;
   * the counts 18/19 ↔ ord₁₉(2);
   * 2¹¹ + 3⁷ = 5·7·11²;
   * 1093;
   * the Heegner run starts.

## 4. Files (`scratchpad/golden/platonic/`)

* `cycles/`:
  * `necklace.c`, `necklace_1_31.txt`, `necklace_32_40.txt`;
  * `minscan.c`, `minscan_mt.c`, `minscan_{24726,50507,301993}.txt`;
  * `bounds{,2,3}.py`, `bounds3.out`;
  * `heuristic{,2}.py`, `nearcyc{,2}.py`, each with its `.out`.
* `baserate{,2}.py`, `idoneal.py`, `facts.py`, each with its `.out`.
