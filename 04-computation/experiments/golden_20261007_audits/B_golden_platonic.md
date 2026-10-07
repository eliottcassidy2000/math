# Audit B: golden/platonic lane (THM-4592, THM-4591, THM-4566 addendum, HYP-9230 arithmetic)

Independent adversarial audit, 2026-10-07. The repo worktree `math-wt-chessboard-20261006` was read only. All scripts and outputs are in this directory. Every computation below uses my own code, not the session's scripts.

## Bottom line

**Nothing is refuted.** Every computation reproduces, including a full single-threaded re-run of THM-4591's whole census range.

Required corrections:
- one prior-art retyping: the THM-4566 addendum theorem is Borwein–Choi 2000;
- one false side-claim: "exactly when n′ ≡ 2 (mod 4)";
- one wrong proof sentence: THM-4591, the negative-cycle case;
- precision fixes:
  - m < 4 needs s < 1;
  - the n = 10 factor-2 twist;
  - the 3x+23 cycles are negative;
  - THM-4592's title writes T for C;
  - the paper's P is ×5.

---

## 1. THM-4592

### 1.1 Golden reading: CONFIRMED (the title needs qualifying)

- **Equivariance.** Re-derived: κ(w′) = φκ(w) − w₀(φⁿ − 1). Checked for n ≤ 22 against C(x)'s word recomputed from the orbit, not a symbolic rotation.
- **Bijectivity** (`a1_golden.out`: own HNF reduction, own periodic-point recurrence, 2-adic word check). For every n ≤ 22, #classes = |Rₙ| = Lₙ − 1 − (−1)ⁿ. Odd n is bijective; for even n the only non-trivial fibre is {0, −1, −2}. A κ-only extension to n = 23..29 agrees (`a3_bij.out`).
- **|Rₙ| formula.** N(φⁿ − 1) = (−1)ⁿ − Lₙ + 1, valid also at n = 1, 2.
- **Defect (notation).** The title says "T-periodic" and "κₙ(Tx)", but the map is the standard map C, as in the body. THM-4591, from the same session, uses T for the Terras map.
- **Defect (overstatement).** The title also states "a bijection for odd n" without qualification, while the status says FINITE-EXACT n ≤ 22.

### 1.2 The n = 5 table: CONFIRMED

- Labels in C-order:
  - 1/13 → 16/13 → 8/13 → 4/13 → 2/13 has labels 3, 1, 4, 5, 9 = H.
  - −5 → −14 → −7 → −20 → −10 has labels 8, 10, 7, 6, 2 = −H.
- The table's orderings, (3, 9, 5, 4, 1) and (8, 10, 7, 6, 2), are exact.
- φ⁵ − 1 = φ³(4 − φ) = 2 + 5φ, and N(4 − φ) = 11.
- C acts as ×4 on the labels, and 4 has order 5.
- These agree with the paper's Proposition 5.2 (label map = reduction mod 4 − φ) and its Figure 1 colours.
- **Wording defect.** The paper's P is x ↦ 5x (Proposition 6.6; Corollary 6.5 gives labels k ↦ 5ʲk + b). ×4 is the paper's Q. Harmless, since ⟨4⟩ = ⟨5⟩ = H.

### 1.3 G₅: CONFIRMED

Brute force (`a1_golden.out`, `a1b_psl.out`):
- ⟨x+1, 4x⟩ = ⟨x+1, 3x⟩ = {ax + b : a ∈ H}, of order 55.
- The orbit of {0,1} has size 55, with trivial stabilizer.
- The orbit of the arc (0,1) is exactly the Paley QR₁₁ arc set.
- Stab(∞) in PSL(2,11), of order 660, restricted to F₁₁ equals G₅.

### 1.4 n = 10: CONFIRMED WITH CORRECTION

Checks:
- φ¹⁰ − 1 = 11φ⁵, and the HNF is 11·Z².
- There are 110 primitive period-10 points, mapped bijectively onto the 110 cells with β ≠ 0.
- C⁵ acts as (α, β) ↦ (α, −β). I checked this by applying C five times to each rational point, not via equivariance.
- There are 55 C⁵-orbits, and Ψ₋ sends them to 55 distinct pairs.

Two precision problems:
- The β = 0 cells hold 13 points: the ten nonzero period-5 points, plus 0, −1 and −2 in the zero cell.
- **Factor-2 twist.** A period-5 point x sits in cell (2κ₅(x), 0), because doubling the word multiplies by φ⁵ + 1 ≡ 2 mod (4 − φ).
  - Under the paper's projection (α, β) ↦ α (Proposition 5.4), x lies over a base cell of the opposite colour. For example, gold 1/13 lies in (6, 0), over teal cell 6.
  - The session's own script silently divides by 2; the theorem text does not say so.

### 1.5 The m < 4 criterion: CONFIRMED WITH CORRECTION

- The roots are s = 2^θ(1 ± √(1 − (m/4)^θ)).
- They are real iff m ≤ 4.
- **At m = 4 the balance is solvable**, with s = 2^θ > 1. So the claim "solvable iff m < 4" is literally false there.
- What THM-4581 needs is 0 < s < 1: its error series Σ s^(h−1) must converge, and ρ = s·2^(−θ) < 1. That holds for some θ ∈ (0,1) iff m < 4.
  - For m < 4, small θ gives s₋ ≈ 1 − √(θ ln(4/m)) < 1.
  - For m ≥ 4, AM–GM gives a left side > 1 whenever s < 1.
- Numerics (`a2_balance.out`):
  - At m = 3.99, min s₋ = 0.99911.
  - At m = 4, min s₋ = 1.0007.
  - At m = 3, s(θ) < 1 on all of (0,1), consistent with THM-4581's s(1/2) = 0.8966.
- The drift equivalence (m < 4 ⟺ ½ log(m/4) < 0) is right.
- "5x+1 having divergent orbits" is conjectural; no 5x+1 orbit is proved to diverge.

### 1.6 R₁₈ remark: CONFIRMED (its typing can be sharpened)

- φ¹⁸ − 1 = 76φ⁹ = L₉φ⁹.
- The −17 cycle has C-period 18, with 7 odd and 11 even steps.
- 19 splits, with φ ↦ 5 or 15 (orders 9 and 18). |R₁₈| = 76².
- "⊃ O/(19)" should read O/(76) ≅ O/(4) × O/(19).
- The 18/19 meeting is **generic Fermat**: every prime p ≡ ±1 (mod 5) divides φ^(p−1) − 1, so O/(p) is a quotient of R_(p−1). For example, R₁₀ = O/(11).
- So the only specific input is "the −17 cycle has C-period 18". The link is numerology rather than dictionary.

---

## 2. THM-4591

### 2.1 Lemma: CONFIRMED WITH CORRECTION (proof sentence)

- Π_odd(3 + 1/xᵢ) = 2^L is right: multiply xᵢ₊₁/xᵢ around the cycle.
- Both bounds are right.
- Positive cycles also have k ≥ L/2, since each factor is ≤ 4. This justifies the Lyndon pruning.
- **The proof is wrong for negative cycles.** It says "The largest factor is at the minimal |x|".
  - For x < 0, 3 + 1/x = 3 − 1/|x| increases with |x|. So the factor at the minimal |x| is the *smallest*. Example: the −5 cycle has factors 2.8 and 2.857.
  - The negative bound comes from smallest factor ≤ geometric mean, i.e. 3 − 1/m ≤ 2^(L/k).
- The minimal-|x| element is odd, so its factor does occur in the product.

### 2.2 Census: CONFIRMED (independently re-run over the whole range)

**Bounds** (`b1_bounds.out`): I maximised over every L ≤ 301,993, with k = ⌊L/log₂3⌋ and k = ⌈L/log₂3⌉, using 80-digit mpmath.

| side | maximum | attained at | scan limit used |
|---|---|---|---|
| positive | 7,216,102,492.69 | (125743, 79335) | 7,216,102,493 |
| negative | 10,295,871,816.1 | (176251, 111202) | 10,295,871,817 |

- Both scan limits cover the maxima.
- CF(log₂3) = [1;1,1,2,2,3,1,5,2,23,2,2,1,1,55,1,4,…].
- a₁₃ = 1, so no semiconvergent lies between 125743/79335 and 301994/190537.
- 301994/190537 would need a bound of 9.85e11.

**Scan logic (minscan_mt.c): sound.**
- The min-|x| element of a cycle of period L ≤ L_MAX never drops below itself and returns at step L. So "0 capped" certifies every cycle of period ≤ 301,993 is reported (indeed, no other cycle has its minimum in range).
- The conjugate map v ↦ (3v − 1)/2, the u128 overflow guard at 2^120, and necklace.c's uint64 cast (3⁴⁰ < 2⁶⁴) are all sound.
- The Lyndon tally is 30,594,555,229 ≈ 3.06e10 words, all with k ≥ L/2 (justified by the lemma).

**Independent re-runs.**
- `b2_scan_full.c`: signed int128, direct T on negatives, both full ranges. Result: exactly −1, −5, −17 and 1. Capped 0, overflow 0, max 508 steps, which matches their `max_stop_neg`.
- `b2_scan.c` to 10⁷: the same cycles.
- `b3_lyndon.py`: all k, L ≤ 20. It finds only 0, −1, 2, −10 and −136.

### 2.3 Interval table: CONFIRMED

- Gersonides/Størmer: the solutions of |2^L − 3^k| = 1 are (1,0), (1,1), (2,1) and (3,2) only.
- 2¹¹ − 3⁷ = −139, and 2187:2048 is the apotome.
- C(11,7)/11 = 30 necklaces.

### 2.4 Ellison-shape table: CONFIRMED WITH CORRECTION

**Ranking** (`b4_ellison.out`): over **all** k with 12 ≤ L ≤ 1500. The original ranked only k ≥ L/2; the result is the same.

| rank | shape | E |
|---|---|---|
| 1 | (19,12) | 0.37075 |
| 2 | (16,10) | 0.07661 |
| 3 | (27,17) | 0.06154 |
| 4 | (13,8) | 0.06070 |
| 5 | (14,9) | 0.04335 |
| 6 | (12,7) | 0.03457 |
| 7 | (38,24) | 0.03370 |

- The sum of the top five is 0.6130.
- Ellison's inequality fails for 12 ≤ x < 2000 at exactly these five shapes.

**Rational cycles** (`b5_ratcycles.out`): reduced denominators below |gap|, with counts.

| shape | denominators (count) |
|---|---|
| (11,7) | 1 (1) |
| (13,8) | 233 (15) |
| (14,9) | none |
| (16,10) | 499 (41) |
| (19,12) | 23 (8), 311 (117) |
| (27,17) | 5 (2), 71 (5), 355 (16), 14303 (843), 71515 (3415), 1015513 (60499) |

- 2992/5 and 2872/5 give the 3x+5 cycles with minima 187 and 347. Each has period 27 and 17 odd steps.
- **Sign.** 2¹⁹ < 3¹², so the eight 3x+23 cycles at (19,12) lie on the **negative** integers. Equivalently they are positive 3x−23 cycles. Their elements nearest 0 are −2263, −2359, −2743, −2963, −3091, −3415, −3743 and −4819.

### 2.5 Platonic base rate: CONFIRMED

- I re-ran it with my own exact-fraction CF code and a new seed: P(≥ 2) = 0.7777, null mean 2.023 (`d1_baserate.out`).
- A₅ × C₂ has element orders {1, 2, 3, 5, 6, 10}, so no element has order 12.
- The statistic counts numerators and denominators ≤ 120 of convergents **and semiconvergents**. The text says "convergents".

---

## 3. THM-4566 addendum

### 3.1 Characterization: CONFIRMED, but KNOWN

**Brute force** (`c1_xyz.out`, `c1_xyz_1e7.out`):
- The n ≤ 10⁷ that are not xy+yz+zx are exactly the 18.
- An independent idoneal test checked every n ≡ 2 (mod 4) up to 10⁷: n is idoneal iff no non-ambiguous primitive reduced form of discriminant −4n exists. It yields exactly 2, 6, …, 462 (16 values).
- Mismatches with the characterization: 0.

**Prior art.** Borwein–Choi, Exp. Math. 9 (2000) already prove the characterization:
- Theorem 3.1: a squarefree n ≡ 2 (mod 4) is an exception iff −4n is a "disjoint" discriminant, meaning one reduced form per genus.
- Theorem 2.6: 4 and 18 are the only non-squarefree exceptions.
- Lemma 2.2: odd n, and n ≡ 0 (mod 4) with n > 4, are representable.

Since 18 is idoneal, these combine to exactly the addendum's theorem. The Selling-parameter proof is a genuine alternative proof: it treats non-squarefree n by descent. It is not a new theorem.

### 3.2 Proof sketch: CONFIRMED WITH CORRECTION

The logic holds. I checked each step:
- the conorm identity det = xy + yz + zx;
- the conorm multiset is GL₂(Z)-invariant;
- a zero conorm ⟺ the form is diagonal in some basis;
- Cl(−4n) ↠ Cl(−4n/g²).

The written sketch has gaps, all repaired by correction 2 below:
- n′ and r are undefined.
- "is diagonal" should say "is GL₂(Z)-equivalent to a diagonal form".
- The reduction to primitive forms (odd content g) is unstated.
- Gauss's #ambiguous classes = #genera is used silently.
- **False: "exactly when n′ ≡ 2 (mod 4)".** `c3_diag_genera.out` checks all n′ ≤ 3000:

  | n′ | #diagonal = #genera? |
  |---|---|
  | n′ ≡ 2, 3 (mod 4) | always |
  | n′ ≡ 4 (mod 8) | always |
  | n′ ≡ 1 (mod 4) | never, except n′ = 1 |
  | n′ ≡ 0 (mod 8) | never |

  The proof needs only the direction "n′ ≡ 2 (mod 4) ⟹ equal", so the theorem is unaffected.

### 3.3 Split-prime lemma: CONFIRMED (the proof is missing from the text)

- Proof: 𝔭 is invertible, and 𝔭² = (α) because the exponent is ≤ 2. Then α ∉ Z, since 𝔭² ≠ (p) for split p. So α = (a + b√D)/2 with b ≠ 0, and p² = (a² + |D|b²)/4 ≥ |D|/4.
- `c2_exp2.out`:
  - there are 101 discriminants: 65 even, 36 odd, the largest |D| = 7392;
  - no split p has 4p² < |D|;
  - no idoneal n > 9 is ≡ 2 (mod 3);
  - 2 splits only at D = −7 and −15.

### 3.4 Reading: CONFIRMED (dictionary/numerology typing is appropriate)

- 18 = 16 + {1, 4}.
- The 19 planes are the squarefree idoneal m ≡ 1 (mod 4).
- Note that 1 lies in both lists.
- (a|19) = (−1)^(ind₂ a) is the generic identity (a|p) = (−1)^(ind_g a) for any primitive root g.

---

## 4. HYP-9230 arithmetic: CONFIRMED

- The continued fraction and convergents are as stated.
- 17 and 29 are semiconvergent denominators: 27/17 and 46/29 lie between 8/5 and 65/41.
- ‖f log₂3‖ values:

  | f | ‖f log₂3‖ |
  |---|---|
  | 5 | 0.07519 |
  | 12 | 0.01955 |
  | 17 | 0.05564 |
  | 29 | 0.03609 |
  | 41 | 0.01654 |
  | 53 | 0.003013 |
  | 306 | 1.475e-3 |
  | 665 | 6.298e-5 |
  | 15601 | 2.625e-5 |

  All match the HYP.
- Ellison's (16,10) = 2·(8,5), (19,12) and (27,17) correspond to modes 5, 12 and 17. This is valid as a dictionary.
- Time constant for 1054/665: 1/(Cθ²) is 4.58e5 octaves at C = 550 and 9.00e5 at C = 280. This is within "4·10⁵–10⁶".
- "23" is a₉, the partial quotient following 1054/665.
- 2π²·27.9 = 550.7.

---

## Required corrections (exact wording)

1. **THM-4566 addendum, typing.** Replace "Theorem (PROVED; the explicit list relies on …)" with:

   > Theorem (KNOWN: Borwein–Choi, Exp. Math. 9 (2000), Thm 3.1 [squarefree n ≡ 2 mod 4: exception iff −4n has one class per genus], Thm 2.6 [only non-squarefree exceptions 4, 18], Lemma 2.2; re-proved below by an independent Selling-parameter argument. The explicit list relies on …)

2. **THM-4566 addendum, sketch, third bullet.** Replace it with:

   > So n is an exception iff every positive definite integral form of determinant n is GL₂(Z)-equivalent to a diagonal form (a zero Selling parameter). Such a form is g·f′ with odd content g and f′ primitive of discriminant −4n′, n′ = n/g² ≡ 2 (mod 4) (so Gauss-primitive). For n′ ≡ 2 (mod 4) the reduced primitive diagonal forms number 2^r (r = number of odd primes dividing n′), which is the number of genera (μ = r + 1), hence by Gauss the number of ambiguous classes; diagonal classes are ambiguous, so they are exactly the ambiguous classes. Hence all forms are diagonal iff Cl(−4n′) has exponent ≤ 2, i.e. iff n′ is idoneal; exponent ≤ 2 descends from n to n/g² because Cl(−4n) ↠ Cl(−4n/g²). (The count equality also holds for n′ ≡ 3 (mod 4) and n′ ≡ 4 (mod 8).)

3. **THM-4591, proof of 1, second bullet.** Replace it with:

   > All factors are positive and 3 + 1/x is monotone, so the factor at the minimal |x| (an odd element) is the largest factor of a positive cycle and the smallest factor of a negative cycle; comparing with the geometric mean 2^(L/k) gives 3 + 1/m ≥ 2^(L/k), resp. 3 − 1/m ≤ 2^(L/k).

4. **THM-4591, title and the (19,12) row.** Replace "eight 3x+23 cycles" with:

   > eight 3x+23 cycles on the negative integers (equivalently positive 3x−23 cycles, since 2¹⁹ < 3¹²)

5. **THM-4592, title.**
   - Replace "{T-periodic points …}" with "{C-periodic points …}".
   - Replace "κₙ(Tx)" with "κₙ(Cx)".
   - Replace "a bijection for odd n" with "a bijection for odd n (FINITE-EXACT, n ≤ 22)".

6. **THM-4592, statement 3, last bullet.** Replace it with:

   > The paper's cell exchanges S (x ↦ x+1) and P (x ↦ 5x = φ²x; its Prop. 6.6, Cor. 6.5) generate G₅; C acts as the paper's Q (×φ = ×4), and ⟨4⟩ = ⟨5⟩ = H.

7. **THM-4592, statement 5.** Replace "The β = 0 cells are the period-5 points." with:

   > The eleven β = 0 cells are the images of the points of period dividing 5 (the zero cell also receives −1, −2). A period-5 point x sits in cell (2κ₅(x), 0), so under the paper's projection (α, β) ↦ α it lies over base cell 2κ₅(x), of the opposite colour (2 ∉ H; e.g. 1/13, gold, lies over teal cell 6).

8. **THM-4592, remark on the golden pair chain.**
   - Replace "is solvable iff m < 4" with:

     > has a solution with 0 < s < 1 for some 0 < θ < 1 iff m < 4 (at m = 4 the only solution is s = 2^θ > 1)

   - Replace "5x+1 having divergent orbits" with:

     > the conjectured (unproved) divergence of typical 5x+1 orbits

## Optional improvements

- **THM-4592, statement 5.** C⁵ ↦ (α, −β) is PROVED, not just FINITE-EXACT: it follows from equivariance with 4⁵ ≡ 1 and 8⁵ ≡ −1 (mod 11).
- **THM-4592, statement 2 for all odd n.** This is likely provable through the arithmetic (Markov) coding of the golden toral automorphism. As it stands it is FINITE-EXACT for n ≤ 22, and I extended the check to n ≤ 29.
- **THM-4592, colour convention.** Gold versus teal for the two cycles depends on the sign of κ: replacing κ by −κ swaps H and −H. The canonical content is "one 5-cycle ↔ H".
- **THM-4592, R₁₈ remark.** Write O/(76) ≅ O/(4) × O/(19) and note the Fermat genericity (§1.6). Retype the 18/19 link as NUMEROLOGY.
- **THM-4592, golden pair chain.** "√T q ≈ 3.8" rests on 2000 paths, with 48 survivors at T = 25600. That is about ±15% statistical error.
- **THM-4591.**
  - Status (a): add "k ≥ L/2 (sufficient by the lemma)".
  - Statement 3: say the ranking holds over all k.
  - Platonic remark: replace "convergents" with "numerators and denominators ≤ 120 of convergents and semiconvergents".
  - The negative-side novelty was not checked against published 3x−1 verifications.
- **THM-4566 addendum.**
  - Add the one-line split-prime proof (§3.3).
  - Note that 1 lies in both "slices".

## Files

All in this directory, each script with its `.out`: `a1_golden.py`, `a1b_psl.py` (THM-4592), `a2_balance.py` (m < 4), `a3_bij.py` (n = 23..29), `b1_bounds.py` (CF, bounds, HYP arithmetic), `b2_scan.c`, `b2_scan_full.c` (census), `b3_lyndon.py`, `b4_ellison.py`, `b5_ratcycles.py` (Ellison shapes), `c1_xyz.c`, `c2_exp2.py`, `c3_diag_genera.py` (THM-4566), `d1_baserate.py`; text extractions `bc.txt` (Borwein–Choi) and `paper.txt` (eleven squares).
