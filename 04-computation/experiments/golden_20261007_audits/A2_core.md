# Audit A2 (core): THM-4590, HYP-9230, THM-4593, THM-4594, HYP-9231, HYP-9217 Update 3

Independent adversarial audit, 2026-10-07. All code is my own.
* To match conventions I read the session's basin and Fourier sources (`golden_20261007_lead/*.c`) and the Mersenne chain convention in `qwit.c`.
* I did not read its sieve, chain-simulation, all-lags or survivor code.
* My labelling (stopping-time recursion) and transform (windowless) differ from the session's.
* The repository was not touched. Scripts and outputs are in `audit/A2/`.

## Verdicts

| Claim | Verdict |
|---|---|
| THM-4590 (1)–(5) and the basin table | CONFIRMED (two numerical remarks wrong) |
| HYP-9230 data, f ≤ 70 | CONFIRMED (reproduced exactly) |
| HYP-9230 mode list | CORRECTION: the spectrum continues above 70; f = 106 is the largest −1-basin mode for k ≥ 20 (f ≤ 400; f ≤ 1100 at k = 24) |
| HYP-9230 law `exp(−Cθ²s)`, C ∈ [280, 550] | REFUTED as stated; the exact characteristic root fits all measurable modes, so the "OPEN factor 2" is resolved |
| THM-4593 (1)–(4) and the numbers | CONFIRMED |
| THM-4594 counts, thresholds, maximality, −5/−17 depths | CONFIRMED |
| THM-4594 (3): 539, 615 "beyond the classical sieves" | REFUTED (Angeltveit's odd-even-even rule) |
| THM-4594 gain over Angeltveit (0.825×) | CORRECTION: 0.957× at K = 30 |
| THM-4594 (4) sign barrier | CONFIRMED only for unrefined classes (J = 0) |
| THM-4594 "Equivalently, U_∞ ∩ Z_{>0} = ∅" | REFUTED as an equivalence |
| HYP-9231 | CONFIRMED (FINITE-EXACT to s = 28) |
| HYP-9217 Update 3 | CONFIRMED WITH CORRECTION |

## 1. THM-4590: CONFIRMED

**Proofs.**
* **(1)** Absorption by step K is a function of `n mod 2^K` for every integer, so `|E_K| = q(K)2^K` exactly and the bound `c_B(x) ≤ q(K)x + 3·2^K` is uniform. The hypothesis `|n| > 2^(K+1)` is unnecessary (absorption is an algebraic identity) but harmless.
* **(2)** Correct.
* **(3)** Holds for every fixed η, as a full limit uniform in `c ∈ [1, 2]`, not along subsequences. Each δ needs finitely many pairs (a, b), and the cut bound is uniform in x. The O(δ/η) shift is right, and "slowly varying" follows from η = 1 plus doubling.
* **(4)** Valid by either route:
  * compose exact ×2/÷2 with the density-one ×3 and +1 maps (positive slopes pull density-zero sets back);
  * or use 6(b)'s chains. The routed driver then lies in a fixed coset, and conditional absorption needs "THM-4581 (3) from every admissible start plus the strong Markov property after the forced prefix". The proof omits this sentence.
* **(5)** Valid. For a density-zero endpoint set E, `#{n ∈ [x,2x) : orbit meets E within K steps} ≤ Σ_(s≤K) 2^s |E ∩ [x/2^K, 6^K x]| = o(x)`.
* No circularity.

**Numbers.**
* `negb.c` labels by stopping time, not by path filling. It reproduces every row k ≤ 24 to 4 decimals (densities, cuts, deviations mod 3/4/9) and the k = 28 density 0.3268. There are no new negative cycle minima below 2^28.
* Two remarks are wrong:
  * **Cut heuristic.** `0.67·q(4.8 log2 x)` gives **0.36** at k = 10 and **0.30** at k = 28, not 0.41 and 0.34 (`qT.py`: q(48) = 0.537, q(134) = 0.456, q(20) = 0.601 against the exact 0.598).
  * **Entry classes.** "21 is never an entry" is false: 21 is the entry of the density-zero class {21·2^j}.
* Sanity check of (4) (`affine.c`). At 2^24, label(m) ≠ label(3m) for 44% of m and label(m) ≠ label(m−7) for 52% (random baseline 67%), falling slowly with k. Statement (4) is true, but the convergence is very slow.

## 2. HYP-9230: CONFIRMED WITH CORRECTION

**Reproduction.**
* My windowed W = 2048 DFT (`negb.c`) reproduces the session's −1-basin amplitudes for f ≤ 70 exactly (4 decimals) for k = 16..24.
* A windowless log-measure DFT (`spec.c`, f ≤ 400) agrees with them to ±0.001 for both classes and k = 16..28.

**The spectrum continues above f = 70** (amplitudes in parentheses):

| class, k | top modes |
|---|---|
| −1 basin, 18 | 12 (.101), **106** (.092), 53 (.061), **265** (.060), **94**, 171, 65 (.042), 147, 253 |
| −1 basin, 24 | **106** (.081), 53 (.058), 12 (.054), **265**, **253**, 94, 147, 359, 65 (.022) |
| −1 basin, 28 | **106** (.075), 53 (.057), 12 (.035), 253, 265, 147, 94, 306, 41 (.014) |
| entry 85, 28 | 53 (.031), **106** (.025), **147**, **159**, **200**, 306, 253, 41 (.012), 94, 12 (.007) |

* All the strong modes have small ‖f log2 3‖:
  * 94, 147, 200 and 253 are the semiconvergent denominators between 65/41 and 485/306;
  * 359 is the next one after 306;
  * 106, 159, 265 and 318 are harmonics of 53.
* The qualitative claim stands. The stated dominant list (53, 12, 41, 65, 29, 24, 36, 17, 5) is wrong for k ≥ 20. Modes 5 and 17 are at most 0.0003 at k = 28.
* Extending to f ≤ 1100 at k = 24 (`spec_neg24_hi*.out`), 106 is still the largest mode. Next come 465 (.056), 1077, 771, 412 and 1024, all with ‖f log2 3‖ < 0.016.
* The HYP's open point "A_665 ≠ 0?" is answered numerically: |β̂(24, 665)| = 0.0085 (−1 basin) and 0.0099 (entry 85). That is about 10× my estimated noise floor of ≈ 0.0007 (Bernoulli noise inflated for neighbour correlation).

**The decay law.** The HYP's own mechanism (log2 x walks by +log2(3/2) or −1) gives the transfer `ψ(u) = ½ψ(u + log2(3/2)) + ½ψ(u − 1)`.
* The mode near frequency f is `e^(wu)`, where w(θ) is the root, continued from 0, of `e^(2πiθ + w·log2(3/2)) + e^(−w) = 2`.
* Expanding gives `Re w = −2π²·27.97·θ² + O(θ⁴)`, which is the HYP's Gaussian. But the expansion parameter 2πθ/(1 − log2(3/2)) is already 0.3 at f = 12.
* The exact effective C = −Re w/θ² (`roots2.py`) is 533 (f = 53), 486 (106), 310 (41), 273 (12), 243 (65), 157 (29), 144 (24), 97 (17), 91 (36).
* Least-squares decay per octave over k = 16..28 (`decayfit.py`):

| f | θ | −1 basin | entry 85 | exact root | Gaussian |
|---|---|---|---|---|---|
| 12 | .0196 | .105 | .106 | .105 | .211 |
| 41 | −.0165 | .088 | .084 | .085 | .151 |
| 65 | .0226 | .118 | .125 | .124 | .281 |
| 29 | −.0361 | .206 | — | .204 | .720 |
| 24 | .0391 | .207 | — | .220 | .845 |
| 118 | .0256 | .139 | .147 | .143 | .362 |
| 147 | −.0105 | .046 | .045 | .044 | .061 |
| 265 | .0151 | .073 | .082 | .075 | .125 |

* The phase drift also matches: predicted +0.0327 cycles per octave at f = 12 (measured 0.0324), +0.0354 at f = 65 (measured 0.0357).
* Misfits occur only for modes too slow to measure in 12 octaves (253, 306) and for weak, mixed 159 and 212 in the −1 basin.
* So no constant C exists, and the factor-2 "OPEN" point is an artifact of the quadratic approximation.

**Typing.** OPEN, NUMERICAL and HEURISTIC are honest.

**Ellison's exceptions.** The mapping is correct arithmetic but a common cause (both lists are small ‖y log2 3‖), not evidence. (16, 10) is f = 10 = 2·5, and the dominant modes (53, 106, 253, 306) are not Ellison exceptions. "Structural, not numerology" (results file) overstates.

**Prior art not cited.** Wirsching's predecessor-density program has a known log-periodic correction: Berg–Krüppel 1998 (as an infinite product) and Tavares, arXiv:2608.27617 (Aug 2026, Fourier series, proved non-constant). I found no earlier spectrum of class densities.

## 3. THM-4593: CONFIRMED

* **(1)** Immediate from the chain table.
* **(2) Witnesses** (`t4593.py`, `mers_wit.py`: exact Fraction chains plus direct orbits on random lifts).
  * R = 2: both partners are at (2, 2) at t = 5, with values 3q+2, 27q+20, 27q+20.
  * All 31 witnesses R = 2..32 verify: no absorption at any t ≤ t₀, all partner states equal first at exactly t₀, and the relation holds on direct orbits. R = 23 ends at (0, −54), which needs one departure; c > 0 still holds.
  * Mersenne D ≤ 7: 4 partners equal from t = 1534. D ≤ 61: 31 partners equal from t = 2372. The driver is unmerged, and the forced prefix 1010 is present.
  * The SRW lower bound applies from the witness state: future bits are fresh, k moves only at flips with a fresh sign, and absorption needs k = 0. So "α > 1/2" is refuted rigorously for these P.
  * A witness for P is a witness for every nonempty subset. So R = 32 covers all subsets of {1..32}, and D ≤ 61 covers S19's odd D ≤ 41.
* **(3)** (`allags_check.c`).
  * Brute force over all translates r ≤ 2^(T+2) equals N_T/2^T for T = 4..10.
  * N_T reproduces 0.4961, 0.3748, 0.2941, 0.2355 and 0.1915 (N_24 = 3,213,595).
  * c(w) is injective on equal-weight words.
  * x* = 0.6090898 and the rate is 0.0346156. The bound exceeds 1 for T ≲ 230.
* **(4)** `k_L = O_L(r−1) − L` held for 40/40 partners (L = 30), with no absorption.
* **Numbers** (`qR.py`). √T·q_R at T = 10^5: R = 1: 10.65 ± 1.04; R = 2: 6.01 ± 0.79; R = 4: 3.36 ± 0.15; R = 8: 1.68 ± 0.23. All consistent with the table.

## 4. THM-4594 and HYP-9231

**Counts: CONFIRMED** (`msieve.c`, an own residue-backward search).

| K | descent | maximal |
|---|---|---|
| 16 | 2114 | 1363 |
| 24 | 286,581 | 172,868 |
| 30 | 12,771,274 | 6,915,181 |

* Descent equals OEIS A076227 (fetched).
* K = 32: 41,347,483 / 23,797,887. K = 34: 151,917,636 / 84,720,794 (12 min run). Both equal the THM table.
* `joins.c`:
  * descent + joins = 1720 at K = 16 (Roosendaal);
  * depth-1 = 1856 / 244,392 and depth-1 + joins = 1562 / 206,402 at K = 16 / 24, exactly the THM's columns;
  * join-only classes: 0 for K ≤ 22.
* **Maximality (1)** holds by definition (one affine m ⇒ one backward word).

**(3) Minimal counterexample: CONFIRMED except the novelty example.**
* The maximum threshold of the certificates used to K = 30 (`thr.c`) is 108.01 = 2^6.755, at the first descent s = 27, a = 17 (Ellison's pair). The largest branch threshold is 98.1. So `n₀ mod 2^30 ∈ U_30` holds given 2^71.
* The certificates (27n+7)/32 and (27n+3)/32 are valid (2000 lifts).
* **REFUTED:** "exclusions beyond the classical sieves start at 2^10". Angeltveit (arXiv:2602.10466 §2.4, applied at fixed depth in his Step 1) also has an **odd-even-even rule**: m = T^k(n), followed by ℓ ≥ 1 odd steps and then 2 even steps, joins (m−1)/2.
  * 539 (word 1101111100): m = (T³n − 1)/2, which is 303 for n = 539; its orbit passes through 455.
  * 615: m = (T⁴n − 1)/2.
* Angeltveit's 2-adic rules (descent + path merging + odd-even-even, `msieve2.c`) equal the maximal sieve for every K ≤ 14:

| K | 16 | 20 | 24 | 30 |
|---|---|---|---|---|
| Angeltveit | 1370 | 16,264 | 178,304 | 7,227,826 |
| maximal | 1363 | 15,870 | 172,868 | 6,915,181 |
| ratio | .995 | .976 | .970 | **.957** |

* Adding joins to his rules changes little (177,588 at K = 24).
* The first genuinely new exclusions are mod 2^15: 11247, 12191, 12799, 23743 (`joins3.c`; certificates verified on 300 lifts by `certs.py`). For example, 12799 has m = (6561n + 1953)/8192, the "m ≡ 8 mod 9" depth-2 variant that Angeltveit describes but did not implement.

**(4) Sign barrier: valid only for J = 0.**
* (a) is correct for unrefined classes.
* (b) extends: −5 and −17 are uncertified at every depth ≤ 120 (`cyccert.c`, exhaustive). The cycle mates are certified by depth 9.
* **But the Setting allows refinement mod 3^J, and (3) uses it.** With refinement, all three classes are certified at depth 0 by one backward turn around their cycle (`refined.py`):
  * −1: (2n−1)/3 on n ≡ 2 (mod 3);
  * −5: (8n−5)/9 on n ≡ 4 (mod 9);
  * −17: (2048n−2363)/2187 on n ≡ 2170 (mod 2187).
  * Each fixed point is the cycle point. The ratio 2^L/3^a is below 1 exactly because the cycle is negative.
  * The first two are (3)'s own mod-9 exclusions.
* So "isolates exactly the grand-orbit minima … as it must" and "any proof must use the sign of n" hold only for 2-adic certificates.
* (c) is true, but the period-P sketch controls only multiples of P. Use the fixed point −1 (`n ≡ −1 mod 2^L·M_odd`).
* (d) is correct and trivial.
* (e) is a remark, not a theorem.

**"Equivalently, U_∞ ∩ Z_{>0} = ∅": REFUTED.**
* A positive cycle of length L and weight a has minimum n₀ = c/(2^L − 3^a) with 3^a < 2^L. Its class mod 2^L is therefore descent-certified, with threshold exactly n₀; compare 1 in the class 1 mod 4, threshold 1.
* So cycles avoid U_∞. Collatz ⇒ U_∞ ∩ Z_{>0} = ∅, but not conversely. The threshold form (first sentence) is correct.

**HYP-9231: CONFIRMED.**
* `survcheck.c`: 0 collisions for every s ≤ 28 (3,524,586 survivors; counts equal A076227).
* The reduction from joins to survivor collisions is correct.

**Prior art.**
* Barina 2025 (PDF read): a 2^34 descent-plus-joins sieve and depth-0 3^k sieves up to 3^6.
* The Angeltveit quote about deeper path merging (m ≡ 4, 8 mod 9: "less effective … did not implement it") is accurate.
* The full branch closure is new relative to these sources, but its gain over published fixed-depth rules is about 4%, not 18%.

## 5. HYP-9217 Update 3: CONFIRMED WITH CORRECTION

* S19 measured α = 0.679 on odd D ≤ 61 and 0.782 on odd D ≤ 41. The D ≤ 61 witness covers both (subset property).
* So their true exponent is ≤ 1/2, and "0.69 was a transient" is *proved* for those lag sets.
* But "every finite lag set has α = 1/2" is unproved. α ≥ 1/2 holds for every finite set (sketch level), but α ≤ 1/2 needs a witness. The core REPORT itself calls the general case a CONJECTURE.
* The status line's "OPEN: … alpha (0.66–0.78)" is stale.

## Required corrections

1. **THM-4590 Numbers.** Replace "which gives 0.41 at k = 10 and 0.34 at k = 28" with "which gives 0.36 at k = 10 and 0.30 at k = 28 (q(46) ≈ 0.54, q(132) ≈ 0.46): it reproduces the slow decline, not the level". Fix the results file §1 too.
2. **THM-4590 entry classes.** Replace "(21 is never an entry, since 21 ≡ 0 mod 3.)" with "(21 is the entry only of the density-zero class {21·2^j}: 21 ≡ 0 mod 3 has no odd predecessor.)"
3. **HYP-9230 statement.** Replace the law with: "|β̂(s,f)| = A_f e^(μ(θ_f)s + o(s)), with phase drift 2πν(θ_f) per octave, where μ + 2πiν is the root (continued from 0) of e^(2πiθ + w log2(3/2)) + e^(−w) = 2. For small θ, μ = −552θ² + O(θ⁴); the effective C = −μ/θ² falls from 533 (f = 53) through 273 (f = 12) to 97 (f = 17)."
   * Delete "C between about 280 and 550".
   * Mark the factor-2 discrepancy as "resolved (audit A2)".
4. **HYP-9230 evidence.** Replace the dominant-mode list with: "Spectrum (f ≤ 400) on small ‖f log2 3‖: convergents 53 and 306; semiconvergents 94, 147, 200, 253, 359; harmonics 106, 159, 265, 318; and 12, 41, 65. The −1 basin's largest mode for k ≥ 20 is f = 106 (0.075 at k = 28). For entry via 85 the order at k = 28 is 53, 106, 147, 159, 200, 306, 253, 41, 94, 12."
5. **HYP-9230 Ellison.** Replace the parenthetical with "(Ellison's exceptions (16,10), (19,12), (27,17) and the low modes 10 = 2·5, 12, 17 are the same near-coincidences 2^x ≈ 3^y: a common cause, not evidence)". Cite Berg–Krüppel 1998 and arXiv:2608.27617.
6. **THM-4594 (3).** Replace "Exclusions beyond the classical sieves start at 2^10: n_0 ≢ 539, 615 …" with "n_0 ≢ 539, 615 (mod 1024), already excluded by Angeltveit's odd-even-even rule. Exclusions beyond descent, joins and all of Angeltveit's 2-adic rules start at 2^15: n_0 ≢ 11247, 12191, 12799, 23743 (mod 2^15), e.g. m = (6561n+1953)/8192 for 12799."
7. **THM-4594 title, Setting, (2), Reading, status.** Change "depth-1 predecessors [Angeltveit 2026]" to "depth-1 path merging and the odd-even-even rule [Angeltveit 2026]".
   * Replace "0.825× (depth-1 + joins)" and "0.82× the best classical combination" with "0.957× Angeltveit's 2-adic rules at K = 30 (6,915,181 against 7,227,826)".
   * Recompute the mod-9 comparison (0.409% against 0.490%) with the odd-even-even rule included.
8. **THM-4594 (4) and the title.** Insert "for unrefined classes (J = 0)". Add the three refined certificates of §4 above.
   * Replace "so any proof must use the sign of n" with "so a proof using only 2-adic class-decided certificates must use the sign".
   * Fix the proof of (c); retype (e) as a remark.
9. **THM-4594 "The remaining conjecture".** Delete "Equivalently, U_∞ ∩ Z_{>0} = ∅" and "precisely a statement that…". Replace them with "Collatz implies U_∞ ∩ Z_{>0} = ∅; the converse fails for cycles (a positive cycle's minimum lies in a descent-certified class whose threshold equals it)". Fix the results file §0 and §12 to match.
10. **HYP-9217.** In the status line and Update 3 heading, replace "every finite lag set has alpha = 1/2" with "every nonempty set of odd lags D ≤ 61 has alpha = 1/2 (one witness covers all subsets; lower bound PROVED, upper sketch level); for general finite lag sets, alpha ≥ 1/2 (sketch level) and = 1/2 is a CONJECTURE". Update the stale "alpha (0.66–0.78)" item.

## Optional improvements

* THM-4590: drop `|n| > 2^(K+1)` from proof 1; add the coset sentence to proof 4; quote the slow (4) mismatch shares.
* THM-4593: state the subset property; note that the all-lags bound is vacuous for T ≲ 230.
* THM-4594 (4)(b): extend to "depth 120 (audit A2)".
* HYP-9230: add the frequency-shift prediction ν(θ) as a second testable output.

## Files (`scratchpad/golden/audit/A2/`)

* **Basins and spectra:** `negb.c`, `spec.c`, `spec_top.c`, `roots.py`, `roots2.py`, `decayfit.py`, `affine.c`.
* **THM-4593:** `t4593.py`, `mers_wit.py`, `allags_check.c`, `qT.py`, `qR.py`.
* **THM-4594 and HYP-9231:** `msieve.c`, `msieve2.c`, `joins{,2,3}.c`, `certs.py`, `thr.c`, `cyccert.c`, `refined.py`, `survcheck.c`, `barina2025.txt`.
* Outputs are the matching `.out` files.
