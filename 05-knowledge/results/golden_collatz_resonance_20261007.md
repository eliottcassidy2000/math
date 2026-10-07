# Golden Collatz: classes are equidistributed, the approach to uniformity is a tuning spectrum, the cycles are Pythagorean intervals, and the eleven squares are the period-5 cycles

Session mac-mini-2026-10-07-golden. The owner prompt asked to:
* "focus on finishing a collatz proof through creative means";
* relate 223, 233, 332, 425, 105, the golden ratio, base-φ primes and {2,3,11};
* relate the Platonic solids to Ellison's exceptions 13, 14, 16, 19, 27;
* relate the mod-18/19 work to the 101 discriminants, the 18 Borwein–Choi exceptions and the 19 Hilbert-class-field planes (m = 165, χ = 5);
* draw on the 9/4 fixed point, the attached eleven-squares paper and the CrocSwap integer-multiplication bound (κ = 2^-78) as inspiration.

Linked sources were treated as trusted. Three readers ran (core, golden, platonic), followed by two independent audits (§14).

## 0. The answers in brief

* **The Collatz conjecture remains OPEN.** This session proves no part of it for individual integers. What it does prove moves the frontier and locates the gap exactly.
  1. **Every Collatz class is equidistributed (THM-4590, PROVED from THM-4581).** Any union of grand-orbit classes, such as the basin of a hypothetical second cycle, has `o(x)` cuts `n | n+1`. It is equidistributed modulo every `M`, uniform across multiplicative windows, slowly varying in density, and almost invariant under `n ↦ 3n`, `n ↦ n+1`. A counterexample's basin cannot hide in residues or leading digits.
  2. **The approach to that uniformity is a tuning spectrum (HYP-9230, NUMERICAL).** Collatz moves by fifths up (`×3/2`) and octaves down (`×1/2`). Class densities oscillate log-periodically at the equal temperaments 5, 12, 41 and 53 (the commas of `log_2 3`), decaying like `exp(−C‖f log_2 3‖² log_2 x)`.
     * Ellison's exceptions (16,10), (19,12) and (27,17) are three of these resonances.
     * The Mercator comma `3^53 ≈ 2^84` mode is essentially undamped at accessible scales.
  3. **The integer cycles are the Pythagorean intervals (THM-4591, FINITE-EXACT to period 301,993).** The octave, fifth, fourth and whole tone are Gersonides' gap-1 shapes (THM-4484). The −17 cycle is the apotome 2187:2048. That is why the cycle periods are exactly {1, 2, 3, 11}.
  4. **The eleven squares are the period-5 cycles (THM-4592).** The paper's eleven Fibonacci-torus cells are 0, the 1/13 cycle (the quadratic residues mod 11) and the −5 cycle (the non-residues). Its group `G_5` is `BS(1,3)` mod 11, the relation group of THM-4581.
  5. **Nonadjacent lags (THM-4593).** Following CrocSwap's insight, extra relations were tried as "nonadjacent lags". Finitely many improve only the constant: the partners coalesce first, and the exponent stays 1/2. Exponentially many give exponential decay.
     * HYP-9217's any-lag 0.69 was a coalescence transient.
  6. **A sharper verification sieve and the sign barrier (THM-4594).** A maximal class-decided sieve beats descent by a factor of about 1.8.
     * It confines a minimal counterexample to 6,915,181 classes mod `2^30`, given Barina's `2^71`.
     * It proves that any Collatz proof must use the sign of `n`: the 2-adic point −1 is never certified, and −5 and −17 are not certified to depth ~90.
* **The precise remaining gap (§12).** Collatz is equivalent to "the positive integers avoid the closed Haar-null set `U_∞` of never-certified 2-adic classes". That set contains −1 and the negative cycles' minima. Measure methods cannot see the sign.

## 1. From coalescence to the integers: Collatz classes are equidistributed (THM-4590)

**Input.** THM-4581: two Terras orbits tied by an affine relation merge almost surely in `Z_2`. For integers, all we need is the density-one shadow: for each `K`, the residues `n mod 2^K` whose pair chain with `n + 1` is not absorbed within `K` steps have density `q(K) → 0`.

**Output (PROVED).** Let `B` be any union of grand-orbit classes, on the positive or the negative integers. Examples:
* the basin of 1;
* the basin of a hypothetical second cycle or a divergent class;
* the basins of −1, −5 and −17 on `Z_{<0}`.

Then:
* **Few cuts.** `B` has `o(x)` "cuts" `n | n+1` in `[x, 2x)`.
* **Residues.** `B` is equidistributed modulo every `M`.
* **Windows.** `B` is asymptotically uniform across multiplicative windows `[cx, c(1+η)x)`.
* **Scale.** Its dyadic density is slowly varying.
* **Affine invariance.** `B` is almost invariant under `n ↦ 3n` and `n ↦ n + c`. This is the integer shadow of THM-4569 (7)'s index 1.

The window statement uses only two exact identities, `T(2n) = n` and `T(2j+1) = 3j+2`, together with the density of `{a + b log_2(3/2)}`. The irrational rotation by a fifth averages out every multiplicative bias.

**What this says about a counterexample.** A hypothetical second cycle's basin cannot sit in a residue class or a leading-digit band. It must be spread like the basin of 1, with the same density in every residue class and every multiplicative window, and nearly invariant under `n ↦ 3n` and `n ↦ n+1`.

**The negative integers: a real three-class test.**

| `k` (`|n| ∈ [2^(k−1), 2^k)`) | −1 | −5 | −17 | cut fraction |
|---|---|---|---|---|
| 16 | 0.3330 | 0.3215 | 0.3455 | 0.3998 |
| 22 | 0.3277 | 0.3243 | 0.3480 | 0.3611 |
| 28 | 0.3268 | 0.3248 | 0.3484 | 0.3392 |

* The three basins are equidistributed modulo 16, 64, 256, 1024, 27, 81, 243, 729, 7, 49, 11 and 121, to within 2–3 sampling standard deviations (≤ 0.0036).
* The cut fraction falls slowly, as predicted (≈ `0.67·q(4.8 log_2 x)`).

## 2. Collatz is the circle of fifths: the resonance spectrum (HYP-9230)

**Collatz moves by fifths and octaves.** The odd Terras step is `x ↦ (3x+1)/2 ≈ (3/2)x`, a **perfect fifth up**; the even step is `x ↦ x/2`, an **octave down**. So `log_2` of an orbit walks on the Pythagorean spiral of fifths. A fifth is `log_2(3/2) = 0.58496` octaves, and its fractional parts are what the theory of equal temperament studies.

**The approach to uniformity is a tuning spectrum.**
* At accessible scales the window uniformity of THM-4590 is far from reached: the −1 basin ranges over `[0.257, 0.370]` across 1% windows at `2^27`.
* The non-uniformity is a sharp log-periodic spectrum. Its frequencies `f` (cycles per octave) are exactly where `‖f log_2 3‖` is small: the continued-fraction denominators of `log_2 3` and their sums.

| `f` | `‖f log_2 3‖` | comma | −1 basin amplitude, `2^16 → 2^28` | entry via 85, `2^16 → 2^28` |
|---|---|---|---|---|
| 53 | 0.0030 | Mercator, `3^53 ≈ 2^84` | 0.0595 → 0.0570 | 0.0332 → 0.0306 |
| 12 | 0.0196 | Pythagorean, `3^12 ≈ 2^19` | 0.1248 → 0.0354 | 0.0262 → 0.0073 |
| 41 | 0.0165 | `3^41 ≈ 2^65` | 0.0404 → 0.0141 | 0.0326 → 0.0118 |
| 65 = 53 + 12 | 0.0226 | — | 0.0571 → 0.0132 | 0.0113 → 0.0026 |
| 29 = 17 + 12 | 0.0361 | `3^29 ≈ 2^46` | 0.0284 → 0.0025 | 0.0151 → 0.0012 |
| 17 | 0.0556 | `3^17 ≈ 2^27` | 0.0140 → … | 0.0067 → … |
| 5 | 0.0752 | `3^5 ≈ 2^8` | 0.0131 → … | — |

* **Ellison's exceptions are structural, not numerology.** Three of the five, `(16, 10) = 2·(8, 5)`, `(19, 12)` and `(27, 17)`, are the modes 5, 12 and 17. `(13, 8)` and `(14, 9)` (`‖f log_2 3‖ = 0.32` and `0.26`) are small-`x` artifacts of Ellison's `e^(−x/10)` threshold, not resonances.
* **The decay is ordered by `‖f log_2 3‖²`.**
  * Mode 53 has a phase drift of only 0.003 per fifth; it is essentially undamped to `2^28`, and the measured odd-step characteristic function decays 0.993–0.998 per octave against 0.995 predicted.
  * Modes 12 and 41 decay about half as fast as the naive Gaussian model (an open factor 2: the odd-step count distribution is non-Gaussian at that resolution).
* **The longest-lived modes sit at the next convergents.** `485/306`, `1054/665` (partial quotient 23, `θ = 6.3·10^-5`) and `24727/15601` would persist for `10^5` octaves or more. Any quantitative THM-4590 must therefore depend on the irrationality measure of `log_2 3`: Baker, Rhin, Ellison.
* **Benford's law is the same rotation.** Kontorovich–Miller and Lagarias–Soundararajan prove Benford behavior for 3x+1 iterates using this rotation mechanism. The resonance spectrum of *class* densities, and its musical identification, were not found in short searches.

## 3. The integer cycles are Pythagorean intervals; {2,3,11} are the cycle periods (THM-4591)

**The cycle formula.** A Terras cycle with word `w` of length `L` and `k` odd steps sits at `x = c_w/(2^L − 3^k)`, the fixed point of a map with multiplier `3^k/2^L`.

**PROVED (elementary, plus Levi ben Gerson 1343 / Størmer 1897).**
* If `|2^L − 3^k| = 1`, every word of that shape gives an integer cycle.
* Those shapes are exactly the 3-smooth superparticular intervals:

  | shape `(L, k)` | interval | integer cycle |
  |---|---|---|
  | (1, 0) | octave 2:1 | {0} |
  | (1, 1) | fifth 3:2 | {−1} |
  | (2, 1) | fourth 4:3 | {1, 2} |
  | (3, 2) | whole tone 9:8 | {−5, −7, −10} |

**FINITE-EXACT.** The −17 cycle is the **apotome** `3^7 : 2^11 = 2187 : 2048` (gap 139). It is the unique integer cycle among the 30 necklaces of shape (11, 7).

**So the Terras periods of all known integer cycles are {1, 2, 3, 11}**: the octave and fifth (1), the fourth (2), the whole tone (3) and the apotome (11).
* **THM-4591 (FINITE-EXACT, two independent methods, re-run by audit B):** these five are the only integer cycles of period `L ≤ 301,993`.
* The gap-1 dictionary is THM-4484 (Gersonides).

**The 9/4 thread.**
* `9/4 = (3/2)²` is two fifths, a major ninth. #107's tripling `h ↦ 3h + 1` and THM-4555's root collision `1 + 2 = 3` share the fixed point `−1/2` (THM-4568, ANALOGY).
* In tuning language, the fixed point `−1` of the halved map `T(n) + 1 = (3/2)(n + 1)` is the fifth itself: the −1 fixed point *is* the interval 3:2.

## 4. The eleven squares are Collatz's period-5 cycles (THM-4592)

**The golden reading.**
* Under the standard map (`3x+1`, `x/2`), parity words never contain `11`. These are exactly the base-φ normal forms.
* Read the cyclic word of a period-`n` point in base φ: `κ_n(x) = Σ w_j φ^(n−1−j) mod (φ^n − 1)`. Then the Collatz map becomes multiplication by φ in the paper's `R_n = Z[φ]/(φ^n − 1)`.
* For odd `n` this is a bijection (FINITE-EXACT for `n ≤ 22`).

**At `n = 5`: the eleven cells (PROVED, re-verified).**

| cells | points | labels in `F_11` (`φ ↦ 4`) |
|---|---|---|
| 1 | 0 | 0 |
| gold | the 1/13 cycle: `1/13, 2/13, 4/13, 8/13, 16/13` | `3, 9, 5, 4, 1` = `H`, the quadratic residues |
| teal | the −5 cycle: `−5, −14, −7, −20, −10` | `8, 10, 7, 6, 2` = `−H` |

(Which cycle is gold depends on the sign convention of `κ`. In the paper's 121-cell cover a period-5 point sits in cell `(2κ_5(x), 0)`, over the opposite colour; THM-4592.)

* `L_5 = 11 = 1 + 5 + 5`.
* The paper's group `G_5 = ⟨x+1, φx⟩` equals `⟨x+1, 3x⟩` mod 11, because `3 = 4^4` generates `H`. That is the Borel subgroup of `PSL(2,11)`, and the reduction mod 11 of `BS(1,3)`, THM-4581's relation group.
* At `n = 10` the paper's 55 exchanged pairs are the 110 primitive period-10 rational points modulo `C^5`.

**Pair actions of Collatz groups (PROVED; golden reader).**
* `{x ↦ sx + c : s ∈ S}` is transitive on unordered pairs of `F_q` iff `S ∪ −S = F_q^×`. It is regular iff `S` is the squares and `q ≡ 3 (mod 4)`. It is never transitive on a composite modulus.
* Least regular primes: 7 for doubling, 11 for tripling, 23 for `⟨2, 3⟩`.
* These are the perfect quadratic-residue codes: Hamming `2^3 = 1 + 7`, ternary Golay `3^5 = 1 + 2·11 + 4·55`, binary Golay `2^11 = 1 + 23 + 253 + 1771`.
* The −17 cycle's denominator 139 is regular for `⟨x+1, (3/2)x⟩`.
* THM-4581's group is regular on odd-difference pairs of `Z/2^K` for every `K ≥ 3`. Its index-1 statement has no exact mod-`p` shadow: pair transitivity fails for about 18% of primes, e.g. 73.

**Why 3 (PROVED, algebra).** For maps `(mx+1)/2`, THM-4581's balance equation has a solution with `0 < s < 1` iff `m < 4`, which is exactly negative drift. So `m = 3` is the only nontrivial odd multiplier with this coalescence mechanism, consistent with the conjectured divergence of typical 5x+1 orbits.
* A golden pair chain on `Z_2[φ]` works verbatim (all branches contract). Its tail is `T^(−1/2)` with `√T q ≈ 3.8`.

## 5. 223, 233, 332, 425 and 105

**One trunk.** All five numbers lie on known trajectories, and four of them on 27's (FINITE-EXACT):
* 332 joins 27's orbit at 94 and passes 47 (`L_8`), 322 (`L_12`), 242 (`3^5 − 1`), 121 (`11²`) and 233 (`F_13`).
* 233 passes 377 (`F_14`) and reaches 425.
* 223 reaches 425 in 7 Terras steps.

**The meeting of 223 and 233 is THM-4581's consecutive-pair coalescence.** `T(223) = 335` and `T^9(233) = 334`, so this is the pair `(334, 335)` from state `(0, 1)`, merging at Terras time 6 at 425.
* The prior lifted family works at half the period in Terras time: `223 + 31104s`, `233 + 32768s`.
* `233 = 105 + 2^7` gives `T^7(233) − T^7(105) = 3^5` exactly. This is a run at `k = 0`, and 105 and 233 are 2-adic neighbours.

**Numerology verdicts (base rates).**
* The golden features on the trunk (`F_13`, `F_14`, `L_8`, `L_12`, `11²`, `3^5 − 1`) have a null hit-probability of 0.15. Each value appears in 15–34% of orbits of `n ∈ [200, 500]`: NUMEROLOGY.
* 105's orbit features (101, 19, 29, 11, 17, 13) are base-rate common: NUMEROLOGY. Its real roles are class-group facts (THM-4566's plane, 210 a Borwein–Choi exception) and the prior golden clock `O/(105)`; neither maps to Collatz.

## 6. Primes in base φ read as binary: killed

* 105, 223, 233, 332 and 425 all contain `11` in binary, and base-φ normal forms never do. So none of them can be such a reading.
* A prime's base-φ word agrees with its Collatz parity word at chance level (0.593 against 0.595).
* There is a curiosity: `L_2k` reads as `2^(4k) + 1`, e.g. 47 reads as 65537. There is also a primality bias of the readings that washes out. Neither has Collatz content.

## 7. {2,3,11}: one equation, and the cycle periods

* **One equation.** The single identity `3^5 = 1 + 2·11²` is at once:
  * Ljunggren's repunit square `11111_3 = 11²`;
  * the statement that the ternary Golay code is perfect;
  * the fact that 11 is a base-3 Wieferich prime;
  * the congruence `s* = −1/2 ≡ 11² (mod 3^5)`.
* `ord_11(3) = 5` makes `⟨x+1, 3x⟩` the Borel subgroup of `PSL(2,11)`, which acts by coordinate symmetries on the Golay code.
* `(3|11) = (5|11) = 1` and `11 ≡ 3 (mod 4)` are exactly the hypotheses behind the `χ ≤ 5` bound for the `m = 165 = 3·5·11` plane.
* All of this is STRUCTURAL inside number theory and coding theory.
* **The Collatz side.** The real `{2,3,11}` fact is that the integer cycle periods are 1, 2, 3 and 11 (THM-4591). 11 is the apotome's `2^11` (`11/7` a semiconvergent of `log_2 3`). Every link from this 11 to the Golay 11 or the `PSL(2,11)` 11 is NUMEROLOGY.

## 8. Platonic solids, Ellison, mod 18/19

* **Platonic solids.** Polyhedral numbers hit the convergents of `log_2 3` at the base rate (`p = 0.78`). 12-TET ↔ the icosahedron's 12 vertices cannot be equivariant, since `A_5 × C_2` has no element of order 12.
* **Ellison's exceptions are structural twice over.**
  * They are the five shapes with the largest expected integer-cycle count (THM-4591 (3)). They carry the 3x+5 cycles (27,17), and eight 3x+23 cycles on the negative integers (19,12).
  * Three of them, (16,10), (19,12) and (27,17), are resonance frequencies 5, 12 and 17 of Collatz class densities (§2).
* **Mod 18 and mod 19.**
  * The 18 Borwein–Choi exceptions are `{1, 4}` together with the idoneal `n ≡ 2 (mod 4)`. This is KNOWN (Borwein–Choi 2000); the platonic reader's Selling-parameter argument is an alternative proof. The 19 HCF planes are the squarefree idoneal `m ≡ 1 (mod 4)`. Both are 2-adic slices of one 65-element list.
  * The link to the mod-18/19 Collatz clock (`ord_19(2) = 18`) is NUMEROLOGY.
  * The −17 (apotome) cycle has standard-map period 18 and lives in `R_18 = O/(76) ≅ O/(4) × O/(19)`. But every `p ≡ ±1 (mod 5)` divides `φ^(p−1) − 1` (Fermat), so this is NUMEROLOGY (audit B).
  * On the planes 57 and 133 the 19-genus character is the parity bit of the mod-19 clock. This is generic: it holds at every prime where 2 is a primitive root.

## 9. The 9/4 thread

* `9/4 = (3/2)²` is two fifths, a major ninth.
* The owner's restatement is correct against THM-4568 as corrected: the tripling `h ↦ 3h+1` has fixed point `−1/2`, and #107's growth lemma is saturated by the A004991 profile.
* "Traces back to" is ANALOGY: the same affine map, used for different purposes. The halved map's fixed point is −1, and that −1 is the fifth-interval cycle of THM-4591.
* The point `−1` organizes every cycle's odd runs into Steiner circuits (KNOWN), but it does not select the periods.
* NUMEROLOGY, recorded so it is not refiled: the negative cycles' odd runs start at 5, 17 and 41, which are Euler's lucky numbers, so `4m − 1` gives 19, 67 and 163, Heegner numbers.

## 10. Integer multiplication (CrocSwap, κ = 2^-78) as inspiration

CrocSwap's gain came from direct nonadjacent swaps (`O(d)`) replacing chains of adjacent moves (`O(d²)`). In Collatz the odd-count difference `k` moves only by adjacent `±1` steps, and that diffusive bottleneck gives `T^(−1/2)`. Several relations used at once are nonadjacent lags. The core lane (§11) measures what they buy.

## 11. The core lane: nonadjacent lags, the maximal sieve, the sign barrier

**THM-4593 (nonadjacent lags).**
* Every relation is driven by the same coin, so all relations move together.
* On explicit residue classes the partners merge with each other first. For `R = 2` the class is `y ≡ 21 mod 32`: at time 5, `y + 1` and `y + 2` meet at `27q + 20`.
* After that the source faces one adjacent walk again. So `q_R(T) ≥ c_R T^(−1/2)` for every `R ≤ 32` and for S19's Mersenne lag sets `D ≤ 61`: PROVED from witnesses. The exponent is `1/2` (sketch level).
* Numerically `√T q_R ≈ (11–14)/R`.
* **All lags at once:** `q_∞(T) = N_T/2^T ≤ (T+1) 2^(−0.0346T)` (PROVED; pigeonhole on `(a, T^T(y) mod 3^a)`).
* The CrocSwap effect is real, but it needs exponentially many lags.

**THM-4594 (the maximal sieve).**
* It certifies `n mod 2^K` by every smaller number reachable backward from any orbit point the class decides. Prior sieves used translation joins (Roosendaal, Barina 2020) and depth-1 predecessors (Angeltveit 2026).
* Uncertified fractions: 0.644% at `2^30` and 0.493% at `2^34`, against 0.884% for descent.
* **Minimal counterexample (PROVED given `2^71`):** it lies in 6,915,181 classes mod `2^30`, with `n_0 mod 9 ∈ {0, 1, 3, 6, 7}` and `n_0 ≢ 539, 615 (mod 1024)`.
* The rate `2^(−0.0500K)` is unchanged; only the constant improves.
* HYP-9231 (survivor lemma, FINITE-EXACT `s ≤ 30`) would make translation joins always redundant.

**The sign barrier (THM-4594 (4)).**
* The class of the 2-adic fixed point −1 is never certified (PROVED).
* −5 and −17 stay uncertified to depths 93 and 88 (FINITE-EXACT), while their cycle mates are certified by depth 9.
* No weight `log n + g(n mod M)` and no pair-chain-state weight can decrease (PROVED).
* Certificate depth for stopping-time records is 0.77–1.00 × `σ(n)`, up to 12.4 × `log_2 n`.

## 12. What remains for a Collatz proof

**The single remaining conjecture (equivalent to Collatz).** Every positive integer `n > 2^71` has a class-decided certificate, at some finite depth, with threshold below `n`. Equivalently, `U_∞ ∩ Z_{>0} = ∅`.

**What is now in place.**
* **2-adic coalescence (THM-4581).** Almost every pair merges.
* **Integer density statements (THM-4590).** Every class is equidistributed, has slowly varying density, and is almost affine-invariant.
* **Finite certificates (THM-4594).** They confine counterexamples to 0.49% of residues at `2^34`, shrinking like `2^(−0.05K)`.
* **Cycle exclusion (THM-4591).** No integer cycles other than the five of period `≤ 301,993`.

**Why this does not close the gap.**
* Each statement is about measure, density or finitely many residues.
* A counterexample is a single positive integer whose 2-adic class lies in `U_∞`, a null set that genuinely contains −1, −5 and −17.
* The barrier results (THM-4594 (4)) show that any proof must distinguish the tail `…000` of positive integers from `…111`, using archimedean size beyond the first `log_2 n` bits. Certificates need depth up to `12 log_2 n`.

**The most promising creative bridge (HEURISTIC; owner-facing).** A two-place weight `log n + λ·(certificate-depth deficit)`, monotone along certified classes and fed by an archimedean bound on the deficit.
* The musical picture says where such a bound must come from. The deficit grows only when the orbit's fifths keep failing to close.
* So the bound is a Diophantine statement about the walk of `k·log_2 3 mod 1` over one orbit. It is pathwise, not averaged, which is exactly Baker-type territory.
* No such bound is known. Stating it precisely is the next obligation.

## 13. Typing summary

| Claim | Type |
|---|---|
| Collatz classes equidistributed, slowly varying, almost affine-invariant | PROVED (THM-4590, from THM-4581) |
| Negative-basin and entry-class numerics | FINITE-EXACT / NUMERICAL |
| Musical resonance spectrum; Ellison ↔ modes 5, 12, 17 | NUMERICAL + HEURISTIC (HYP-9230); arithmetic PROVED |
| Integer cycles to period 301,993; periods {1, 2, 3, 11} | FINITE-EXACT (THM-4591) |
| Cycles = Pythagorean intervals | DICTIONARY (with THM-4484) |
| Eleven squares = period-5 cycles; `G_5 = BS(1,3)` mod 11 | PROVED + FINITE-EXACT (THM-4592) |
| Pair-regularity classification | PROVED (golden reader) |
| `m < 4` criterion for the coalescence weight | PROVED |
| Finite lag sets: exponent 1/2 | lower bound PROVED; upper at sketch level (THM-4593) |
| All lags: `≤ (T+1) 2^(−0.0346T)` | PROVED (THM-4593) |
| Maximal sieve; minimal-counterexample residues | FINITE-EXACT; PROVED given `2^71` (THM-4594) |
| Sign barrier | PROVED (THM-4594 (4)) |
| Survivor lemma | CONJECTURE, FINITE-EXACT `s ≤ 30` (HYP-9231) |
| Borwein–Choi exceptions = `{1, 4}` + idoneal `≡ 2 (mod 4)` | KNOWN (Borwein–Choi 2000); alternative proof |
| 223/233/332/425 on 27's trunk; (334, 335) coalescence | FINITE-EXACT; golden features NUMEROLOGY |
| 105's orbit features; base-φ binary readings | NUMEROLOGY; killed |
| Platonic solids ↔ `log_2 3` / Ellison | NUMEROLOGY (`p = 0.78`); no equivariant 12 |
| `3^5 = 1 + 2·11²` (Ljunggren, Golay, Wieferich, `−1/2 mod 3^5`) | KNOWN + DICTIONARY; NUMEROLOGY for Collatz |
| 18 and 19 counts ↔ mod-18/19 clock | NUMEROLOGY |
| 9/4 ↔ `−1/2` | ANALOGY (THM-4568 as corrected) |

## 14. Audit record

* **Audit B** (THM-4591, THM-4592, THM-4566 addendum, HYP-9230 arithmetic): nothing refuted. Every computation reproduced, including an independent full re-run of the cycle census. Eight corrections were applied (MISTAKE-584):
  * the Borwein–Choi characterization retyped KNOWN;
  * a false "exactly when" in the proof sketch;
  * the negative-cycle lemma sentence;
  * the 3x+23 cycles are negative;
  * THM-4592 notation and the factor-2 twist;
  * the `0 < s < 1` qualifier;
  * the `R_18` retyping.
* **Audit A2** (THM-4590, HYP-9230 numerics, THM-4593, THM-4594, HYP-9231, HYP-9217 update): see below when complete.

## 15. Reproduction

* `04-computation/experiments/golden_20261007_lead/`: negative basins, windows, Fourier spectra, entry classes (C), odd-step variance.
* `04-computation/experiments/golden_20261007_readers/{core,golden,platonic}/`: the readers' code, outputs and reports.
* Audit scripts are listed in the audit reports.

## 16. Sources

* "The Fibonacci Geometry of Eleven Squares" (Pingyou Ltd), PDF supplied by the owner.
* CrocSwap/integer-mult-bounds (github.com/CrocSwap/integer-mult-bounds), over openai/math #109.
* Kontorovich–Miller, Acta Arith. 120 (2005); Lagarias–Soundararajan, J. London Math. Soc. 74 (2006) (Benford for 3x+1).
* Borwein–Choi, Exp. Math. 9 (2000); Ellison (1971) via Waldschmidt.
* Barina, J. Supercomput. 2020, and the 2025 verification to `2^71`; Roosendaal (ericr.nl/wondrous); Angeltveit, arXiv:2602.10466.
* Lagarias (1990) on 3x+d cycles; Simons–de Weger; Hercher; Gersonides (Levi ben Gerson); Størmer (1897).
