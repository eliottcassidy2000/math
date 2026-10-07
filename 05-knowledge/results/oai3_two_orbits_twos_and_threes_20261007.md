# Two orbits, twos and threes: two Collatz orbits merge almost surely; zeroless powers of two; the idoneal planes; the sharp 9/4

Session mac-mini-2026-10-07-oaimath3, the third openai/math session.

**Owner prompt (abridged).**
* "figure out how the two orbits' exponents correlate over long gaps. consider all the openai math results to be sufficiently correct and verified … just explore heavily for possible connections to leverage, … also pick out a few of your own from that repo …"
* Attached: You, arXiv:2610.07537 (GCH above a strongly compact cardinal); openai/math's free-group-factors isomorphism (#287); "investigate whether every power of two above 2^86 contains a zero"; "consider the importance of Serre's intersection-multiplicity conjecture"; Euler's 65 idoneal numbers; "think of how 9/4 relates to our collatz and dynamical systems work" (#107, `ω ≤ 9/4`).
* A pasted list of openai/math headline claims.

**Method.**
* Six readers covered zeroless powers, number theory, algebra, #107, groups and operator algebras, and the "pick" lane.
* Every claim taken into canon was re-checked by the session itself, with independent code.
* Three adversarial audits followed (§11).
* Owner directive: openai/math results are accepted as correct. Dependent results are typed "PROVED modulo openai/math #N (accepted per owner directive 2026-10-07)".

**Concurrent work on the same prompt.**
* opus-2026-10-07-S20: THM-4565, `two_readers_exponent_correlation_20261007.md`.
* codex-tiling: `collatz_overlap_kernel_integration_20261007.md`, `collatz_orbit_information_openai_synthesis_20261007.md`, `powers_two_decimal_windows_20261007.md`, and others.
* How they fit is in §9. This session's zeroless theorem was renumbered THM-4565 → THM-4580 after a namespace collision (MISTAKE-581).

**Canon from this session.**

| Item | Content | Type |
|---|---|---|
| THM-4581 | Haar coalescence: almost-sure merge; rate `T^(−1/2)` up to logs | PROVED; rate at sketch level |
| THM-4569 | Terras clock: recurrence; box reduction; orbit-equivalence form | PROVED; computer-assisted bounds |
| THM-4564 | two readers on one tape | PROVED; repaired by codex-tiling (§9) |
| THM-4580 | zeroless powers of two: the 5-adic tree; verification to `1.1·10^11` | PROVED + FINITE-EXACT |
| THM-4566 | Hilbert-class-field planes are the idoneal planes | PROVED; complete modulo #003 |
| THM-4567 | THM-1300's Keller map at the prime 2 | PROVED + FINITE-EXACT |
| THM-4568 | the 9/4 growth lemma is exactly sharp | PROVED; finite-size bound modulo #107 |

**Hypotheses.**
* RESOLVED (PROVED): HYP-9213, HYP-9214, HYP-9220.
* PARTLY RESOLVED: HYP-9217 (exponent proved; constant open).
* Formulation refuted: HYP-9218 (by codex-tiling; its target is proved without it).
* OPEN, new: HYP-9219.
* Updates were made to THM-4512, THM-4555, THM-4556, THM-4558, THM-4564 and THM-4569.

---

## 0. The answers in brief

1. **How the two orbits' exponents correlate over long gaps.**
   * In the 2-adic model the two orbits **merge almost surely** (THM-4581). After the merge their exponent streams are identical, shifted by the merge lag.
   * Before the merge they are uncorrelated at long gaps:
     * in Terras time, parities at distinct times are exactly pairwise independent;
     * in odd-step time, exponents couple only at window overlaps, and the measured depth law there is Haar to 0.1% away from small gaps.
   * The probability of no merge by time `T` decays like `T^(−1/2)`. The lower bound is PROVED, and the upper bound holds up to a `(log T)^2` factor at sketch level. The constant (about 11.1 for `y` vs `y + 1`, 16.7 for S19's lag-1 Mersenne pair) obeys a one-big-jump law: `(|k_0| + E[#excursions])·√(4/π)`.
   * Consequences: HYP-9213, HYP-9214 and HYP-9220 are proved, and so is the exponent half of HYP-9217.
2. **Does every `2^n` above `2^86` contain a zero?**
   * OPEN, but verified for `87 ≤ n < 1.1·10^11`. This is an independent confirmation, not a record.
   * New structure: the 5-adic tree, the `9Z_k/2` recursion, growth brackets around 9/2.
   * New obstruction: no finite-digit argument can work.
3. **Serre's intersection-multiplicity conjecture.**
   * No leverage on our fronts. Its open content (positivity over ramified mixed-characteristic bases in dimension `≥ 4`) is far from every intersection we use. It is recorded as a guardrail.
   * The algebra lane produced THM-4567 and a new Collatz collision instead.
4. **The 65 idoneal numbers.**
   * openai/math #003 (no zero with `Re s > 7/8`) excludes Siegel zeros. So the list is complete, modulo #003.
   * Our leverage: the planes whose Gaussian extension is an abelian Hilbert class field are exactly 19 idoneal planes. They include the Moser plane (`χ = 4`) and the 5-chromatic Polymath plane (THM-4566).
5. **How 9/4 relates to our Collatz work.**
   * `9/4 = 3·(3/4)` comes from a shifted tripling `h ↦ 3h + a − 1`. At `a = 2` this is literally the odd Collatz step.
   * Its extremal profile linearizes through the repelling fixed point `−1/2`, the same Schröder chart as the Mersenne chain.
   * 9/4 is exactly the ceiling of #107's argument (THM-4568).
   * The other `(3/2)²` coincidences are NUMEROLOGY.
6. **The pick lane.**
   * You's GCH theorem has an exact DICTIONARY entry: its first-difference map `min(X Δ Y)` is our agreement depth.
   * Of the lane's own picks, #238 (Thorp) shares the martingale-energy identity. #145 and #017 give no leverage.
   * **The lane's pair-chain analysis produced Theorem H, which became THM-4581.**

---

## 1. How the two orbits' exponents correlate over long gaps

**The object.** S19's coupled pair (HYP-9217): `y` 2-adic Haar, `x = 3·2^v y + 1`, compared at lag 1. More generally, any two orbits tied by an affine relation `u = 3^k v + e` with `e ∈ Z[1/3]`.
* Each exponent stream is i.i.d. `Geom(1/2)` on its own.
* The question is their joint law across gaps.

### 1.1 Terras clock: a Markov chain, a random walk, and almost-sure coalescence (THM-4569, THM-4581)

**The chain.**
* Step both orbits with `T(x) = x/2` or `(3x+1)/2`. The relation persists as `u_n = 3^(k_n) v_n + e_n`.
* `(k, e)` is a Markov chain driven by `v`'s fair parity bits (THM-4569 (1)).
* `k` is the difference of the two odd-step counts. It moves `±1` exactly when `e` is odd, the predictable "disagreement" or flip times, and its direction is a fresh fair coin.
* So `k` is a simple random walk on the flip clock, and it is recurrent unconditionally (THM-4569 (2)–(3)). The same chain was found independently by the groups reader and by the pick lane.

**Almost-sure coalescence (THM-4581 (3)).** The missing archimedean step is supplied by one observation and one bookkeeping lemma.
* *Observation.* With `|f| = |e|·3^(−max(k,0))`, a flip toward `k = 0` always multiplies `|f|` by about 3/2, and a flip away by about 1/2. So `|f|^θ s^|k|`, with `s = 2^θ(1 − √(1 − (3/4)^θ))`, is a martingale on flips.
  * A departure from `k = 0` costs a factor 1/2 in `|f|` and a factor `s` in the weight.
  * So along returns to `k = 0`, `E|e|^θ` contracts by `ρ = 1 − √(1 − (3/4)^θ)`, plus a constant: 0.634 at `θ = 1/2`, against 0.59–0.61 measured.
* *Lemma.* Runs of non-flip steps at level `h` cannot be stretched adversarially. Beyond `v_2(3^h − 1) − 1` steps, every continuation costs a fresh fair coin. So the expected time per level is bounded, and the additive errors sum.
* *Conclusion.* `|e|` at returns is tight. Every state at `k = 0` reaches `(0, 0)` with positive probability. Lévy's 0–1 law gives absorption almost surely.
* *Rate (THM-4581 (4); upper bound at sketch level).* `c T^(−1/2) ≤ P(no merge by T) ≤ C T^(−1/2) (log T)^2`.
  * The lower bound: the walk must return to 0.
  * The upper bound: geometric tails for the number of excursions, plus first-passage tails for their lengths.

**Checks (this session).**
* The chain equals direct 2-adic orbits on 450,000 steps, with 0 mismatches.
* The local lemmas hold on all 336,040 (state, bit) pairs with `|k| ≤ 10`.
* Run lengths have mean 2.00 at every level, under the bound `m_h + 1`.
* Survival from `(0, 1)`: `√T q(T) = 11.2, 11.0` at `T = 6400, 25600`.
* Chain-free check: actual big-integer orbits of `y` against `y + 1`, `3y`, `y/3`, `9y + 5`, `y − 7` and `27y + 1/3` merge within 8000 steps in 80–89% of 120 trials each.
* Exact residue census: for every `n mod 2^K` with `K ≤ 20`, chain absorption by step `K` equals an actual integer merge of `n` and `n + 1`, with 0 disagreements. The unmerged fraction is `626933/2^20 = 0.598` at `K = 20`.

**What it gives.**

| Item | Statement | Status |
|---|---|---|
| HYP-9220 | `y` and `y + 1` merge for almost every 2-adic `y` | PROVED |
| THM-4569 (7) | Index `[R_A : R_C] = 1`: up to null sets, the Collatz grand-orbit relation on `Z_2` is the orbit relation of the affine group `Z[1/6] ⋊ ⟨2,3⟩` | PROVED |
| Integers | `n` and `n + 1` (also `n` and `3n`) meet at equal Terras time for a set of `n` of natural density 1; the exceptional residues mod `2^K` have density `≤ C K^(−1/2)(log K)^2` | PROVED (density 1); the bound at sketch level |
| HYP-9213 | Mersenne `σ`-levels are `o(A)`: S19's lag-1 switch is a chain absorption, so `μ_2(S) = 1`, then S18 Prop. 6 | PROVED |
| HYP-9214 | the reset-2 debt resolves with probability tending to 1: `(n, (n−1)/2)` is the chain from `(0, 1)`, then a clopen transfer | PROVED |
| HYP-9217 (1) | almost-sure merging | PROVED |
| HYP-9217 (1) | decay exponent exactly 1/2 (up to logs) | PROVED at sketch level |
| HYP-9217 (1) | the constant `c₁` | OPEN |
| HYP-9217 (2) | any-lag `α ≥ 1/2` | PROVED at sketch level |

**The constant (HEURISTIC + NUMERICAL, THM-4581 (7)).**
* A long survival is one long excursion: `q(T) ~ (|k_0| + E[J])·√(4/(πT))`, where `J` is the number of excursions from 0 before absorption.
* For `(0, 1)`: `E[J] ≈ 9.93` (9.83 in 4000 paths to `T = 2·10^5`, truncation-corrected by audit A), predicting 11.2 against 11.06–11.18 measured to `T = 1.6·10^6`.
* For S19's lag-1 pair (post-prefix state `(2, 1)`): `E[J] ≈ 12.8`, predicting `(2 + 12.8)·1.128 = 16.7`, which is S19's constant (audit A).

**Correlations in this clock (THM-4581 (5)).**
* Parities at distinct times are exactly pairwise independent.
* The equal-time correlation is `1 − 2P(σ_n = 1)`, which tends to 1.

**Prior art.** No prior statement of almost-sure coalescence, or of density-one coalescence of `n` and `n + 1`, was found.
* Garner 1985 studies consecutive heights.
* Burson's withdrawn arXiv:2005.09456 studies first coalescence points numerically.
* Kontorovich–Lagarias arXiv:0910.1944 note that `8n+4` and `8n+5` coalesce.
* This is not a claim of novelty; audit A searched further (§11).

### 1.2 Odd-step clock: where and how the exponents couple (THM-4564; S20's THM-4565; codex-tiling's repairs)

**Coupling sits at window overlaps.**
* Both orbits read the same 2-adic digits of `y`. The `x`-orbit's reading head sits `L` digits behind.
* Each exponent is the agreement length between fresh digits and a past-written 2-adic target.
* Coupling happens only where one orbit re-reads digits the other has read: at **window overlaps** (S20, THM-4565).
* Exact tape alignments `x_t = 3^k y_s + κ`, with `κ` past-measurable, are the `λ = 0` case (THM-4564), about a third of all overlaps.

**The exact law at an overlap.**
* It is ultrametric, with past-written depth `δ` (S20: `M`).
* The covariance is `2 − 6·2^(−δ)`.
* The pair is exactly independent iff the depth is `Geom(1/2)`. Covariance zero alone does not suffice: codex-tiling's mixture with masses 1/3 at 1 and 2/3 at 2 has covariance 0 but is dependent.

**Lockstep continuations.**
* They follow `E' = 3E/2^a − (3^k − 1)`.
* The min-depth rule `δ' = min(δ − a, v_2(3^k − 1))` holds off the equality branch. At equality, cancellation can give a larger depth. codex-tiling's witness: `y_s = 17`, `x_t = 49`, `E' = E = −8`.

**The debt.**
* The normalized debt obeys `ρ' = (3/2^a)ρ − 1 + 2^(−L')`.
* The unforced comparison perpetuity has THM-4554's Moran function as moment function and tail index exactly 1 (`3^j/2^(A_j)` is a mean-one martingale, since `3 = 2² − 1`).
* Applying it to the actual debt needs control of the forcing term.

**Numbers (3000 exact pairs, 20,000-bit sources, 1.66M alignments; S20: 3.7M overlaps).**
* At `L ≥ 8` the depth agrees with Haar to 0.1%.
* Successive fresh depths are consistent with independence (`χ² = 32.6` on 35 dof).
* `|Cov|` is at most 0.0014 at lags up to 160.
* `Var(L)/K = 4.00` up to `K = 256`.
* Detectable deviations sit at `L ≤ 8`, with an even-gap phase pattern. No cutoff is proved.

**HYP-9218 is withdrawn in its original form** (codex-tiling).
* It conditioned the depth on the full joint past, which already determines the depth.
* Its target (recurrence, merging, the `T^(−1/2)` exponent) is now proved without it (THM-4581).
* What remains is fine structure: an invariance principle with variance 4 per odd step, and the exact flip density 1/2 (measured 0.5002).

### 1.3 The answer

**The answer as a table (NUMERICAL).** S19's pair: `y` a random 8000-bit odd integer, `x = 3·2^v y + 1`, 1500 pairs, 1200 odd steps. Script: `oai3_20261007_coalescence/lag_correlation_profile.py`.
* Every one of the 1053 merges has lag 1: `x_s = y_(s+1)` from then on.
* The correlation sits at the merge lag only, and tracks the merged share.

| odd-step window `s` | merged share | `Corr(b_s, a_(s+d))`, `d = −2` | `−1` | `0` | `+1` | `+2` | `+3` |
|---|---|---|---|---|---|---|---|
| [5, 15) | 0.109 | +0.001 | +0.005 | +0.013 | **+0.152** | +0.023 | +0.000 |
| [30, 50) | 0.241 | −0.000 | −0.007 | +0.009 | **+0.262** | +0.002 | +0.001 |
| [100, 150) | 0.391 | −0.001 | −0.003 | +0.005 | **+0.422** | −0.006 | +0.001 |
| [300, 400) | 0.531 | −0.004 | −0.003 | +0.000 | **+0.560** | −0.000 | −0.002 |
| [800, 1000) | 0.664 | −0.000 | +0.002 | −0.002 | **+0.676** | +0.001 | +0.003 |
| [1100, 1195) | 0.697 | −0.001 | +0.004 | −0.000 | **+0.702** | −0.000 | +0.003 |

* Off the merge lag: zero to three decimals.
* At the merge lag: the merged share (which tends to 1, THM-4581), plus a small pre-merge excess. The excess is largest early, from near-merge locking at small `L`.


**Over long gaps the two orbits' exponents are uncorrelated until the orbits merge, and the merge is almost sure.**
* Before the merge, every correlation is carried either by the predictable disagreement bit (Terras clock), which decides *when* the odd-step difference moves but never *which way*, or by ultrametric locking at window overlaps (odd-step clock).
* After the merge, the streams coincide exactly, shifted by the merge lag.
* So the long-gap correlation at the merge lag tends to 1. Its deficit is the non-merge probability, of order `T^(−1/2)` (logarithms aside, upper bound at sketch level), with constant `≈ (|k_0| + E[J])·√(4/π)`.

---

## 2. Does every power of two above `2^86` contain a zero? (THM-4580, HYP-9219)

**Status: OPEN.** None of openai/math's 722 manuscripts concerns decimal digits.

**Verification (FINITE-EXACT).**
* `2^n` contains a 0 for every `87 ≤ n < 1.1·10^11`. For `957 ≤ n < 1.1·10^11` the 0 lies among the last 251 digits. This fails beyond the range: A031142(43) = 181477218727 has its first zero at digit 261.
* This was the zeroless reader's C verifier, validated against Python big integers and checked against exact end states.
* An independent verifier written in this session confirms `n < 2·10^9` and reproduces OEIS A031142's rightmost-zero records 24–38.
* Context:
  * A007377 records a check to `10^10` (Radcliffe 2022).
  * A031142's record table (Griffiths 2012), if complete, already implies the conjecture for `n < 7.88·10^12`, because a zeroless power would set an enormous record.
  * So our run is an independent confirmation, not a new record.
  * codex-tiling checks the witness `n = 103233492954`: its first zero from the right is the 250th digit.

**Structure (PROVED).**
* The trailing `k` digits of `2^n` repeat with period `4·5^(k−1)`: the 5-adic clock of 2, twin of THM-4556's 2-adic clock of 3.
* Each zeroless class has 4 or 5 zeroless lifts, decided by one parity bit. The counts satisfy `Z_(k+1) = (9Z_k + Δ_k)/2`, where `Δ_k` is the trace of an explicit cyclotomic unit.
* `Z_k` equals the number of zeroless `k`-digit multiples of `2^k`. It is computed exactly to `k = 40` (OEIS had 26).
* The growth rate lies in `[4.478, 4.524]` and is numerically exactly 9/2 (HYP-9219: `Z_k ≈ 0.8877·(9/2)^k`).
* `#{n ≤ x : 2^n zeroless} ≪ x^0.938`.
* **Obstruction.** No argument that uses finitely many leading or trailing digits can settle the problem. A proof must couple 5-adic and archimedean information through the middle digits.
  * This is the same adelic shape as Collatz: THM-4564 (6) sets the debt's real size against its 2-adic valuation.
  * It is also the shape of Mahler's 3/2 problem.

**Heuristic.**
* A calibrated model (`P ≈ 1.2034·0.9^(digits)`) predicts 33.9 cases up to 86, against 36 actual.
* It predicted 2.3 more beyond 86, so an empty tail was a 10–15% event a priori.
* Beyond `1.1·10^11` it predicts about `10^(−1.5·10^9)`.

**Link to §1.** Both problems are "two clocks reading one object".
* Collatz: the 2-adic debt against its archimedean size.
* Zeroless powers: trailing digits (the 5-adic clock) against leading digits (rotation by `log_10 2`).
* In both, each clock alone is equidistributed and provably harmless, and the content lies in their coupling.

## 3. How 9/4 relates to our Collatz and dynamics work (THM-4568)

**In the paper (#107).**
* `9/4 = 3·(3/4)`, where 3/4 bounds the mean dot-product exponent `t`.
* It comes from a shifted tripling `P(a, 3h + a − 1) ≥ 3P(a, h)` on the profile of polynomial-multiplication tensors, together with concavity. These force `P(a,a) ≥ a^(4/3)` against the rank bound `(2a−1)^(1/t)`.
* The critical exponent is `4/3 = 2/(1 − s*)`, with `s* = −1/2` the fixed point of `s ↦ 3s + 1`. At `a = 2` the tripling is literally `h ↦ 3h + 1`, the Collatz odd step.

**New and sharp (THM-4568, PROVED).**
* The lemma's hypotheses have a closed-form least profile `σ_m(2M + m − 1)/2`, built from Γ-function ratios. It meets the rank bound exactly at `t = 3/4`.
* So **9/4 is precisely the ceiling of the paper's argument**, and scale-free gadgets of this kind cannot do better.
* Finite-size corollary (modulo #107): `ω ≤ 9/4 + 0.684/ln a`.
  * At `a = 188` this beats 2.371339 (Alman et al., SODA 2025). At `a = 190` it beats 2.371177 (Dupont et al. 2026), the latest pre-#107 bound that #107 cites.
  * The paper's crude bound needs `a ≈ 382,000` for that.

**The Collatz relation.**
* **Exact identity (ANALOGY in use).** The paper's gadgets are the affine maps conjugate to multiplication through the repelling fixed point `−1/2`. In Collatz terms this is the unhalved odd step `3h + 1 + 1/2 = 3(h + 1/2)`, with fixed point `−1/2` (THM-4555). Its halved form `T(n) + 1 = (3/2)(n + 1)`, with fixed point `−1`, drives the Mersenne chain (THM-4556).
* **Instructive difference.**
  * The paper's maps commute and share one fixed point, so a closed-form extremal profile exists.
  * Collatz words have word-dependent fixed points `c_w/(2^L − 3^w)` (the rational cycles), so no simultaneous linearization exists.
  * That is one more way to see why Collatz resists a "least-profile" argument.
  * THM-4581's weight `s^|k|` is the probabilistic substitute: it linearizes only *on average*, along the random walk of `k`.
* **`(3/2)²` recurs.** The two-odd-step multiplier; the Mersenne run's extremal growth; the overshoot bound 3/2 in THM-4564 (6); the tube slope `(9/4)ξ`; and compression by `t = 2/3` in the free group factors, which multiplies `r − 1` by 9/4.
* **Related numbers.** 4/3 (THM-4504's mean offspring), 3/4 (the drift), `4/9 = τ*(2)` (HYP-9210), and now `(3/4)^θ` inside THM-4581's `s(θ)`.
* These coincidences are NUMEROLOGY: they agree only because `3 = 1 + 2` in each setting. They are registered so no later session re-derives them.

## 4. The idoneal numbers and the quasi-Riemann hypothesis (THM-4566)

**The 65 numbers in the prompt are Euler's idoneal numbers.**
* openai/math #003 (zero-free `Re s > 7/8` for every Dirichlet L-function) excludes Siegel zeros. So Tatuzawa's class-number bound has no exception, and Elsenhans–Klüners–Nicolae's search finishes the list. #003's own 11/12 paper says so.
* The session's sieve over negative discriminants up to `2.34·10^14` finds exactly the 101 discriminants of exponent at most 2, among them the 65 idoneal numbers.
* **Corollaries, modulo #003:**
  * the idoneal list is complete;
  * the Borwein–Choi problem is closed: every positive integer except the 18 numbers 1, 2, 4, 6, 10, 18, 22, 30, 42, 58, 70, 78, 102, 130, 190, 210, 330, 462 is `xy + yz + zx` with `x, y, z ≥ 1`;
  * the exponent-4 and exponent-8 field lists (EKN) are complete;
  * class-group exponents tend to infinity effectively.

**Our leverage (THM-4566, PROVED; complete modulo #003).**
* The planes `F²` whose `F(i)` is an abelian Hilbert class field of an imaginary quadratic field are exactly the 19 planes `Q(√p : p | m)²` with `m` squarefree idoneal, `m ≡ 1 (mod 4)`.
* Both record planes of the Hadwiger–Nelson field tower are among them:
  * the Moser plane, `Q(√3,√11)² = H(Q(√−33))`, with `χ = 4`;
  * the plane `Q(√3,√5,√11)² = H(Q(√−165))`, which contains the Polymath field and has `χ = 5`.
* Chromatic numbers:
  * 2 for `m = 1, 5, 13, 37, 85`;
  * 3 for `m = 21, 57, 93, 133, 273` (`m = 133` through an explicit 7-cycle in `Q(√7)²`; cf. Madore);
  * 4 for `m = 33` and `177` (`m = 177` by the generalized spindle `N = 723`, found in audit B);
  * 5 for `m = 165`;
  * 4–5 for 345, `≥ 4` for 105, 357 and 1365, 3–4 for 253, 3–5 for 385;
  * for these planes `χ = 2` iff `m` has no prime factor `≡ 3 (mod 4)`.

**Also.** Ellison's 1971 bound on `|2^x − 3^y|` gives a shorter all-length closure of THM-4512, without the `2^42` bridge or Matveev. Its exceptions were re-checked to `x = 20000`.

## 5. Serre's intersection multiplicity (#193), Hilbert's tenth problem over Q (#004), Kaplansky (#197)

**Serre positivity: no leverage here.**
* Positivity was already known in several cases:
  * equicharacteristic, and power series over any complete DVR (Serre);
  * when both modules are Cohen–Macaulay;
  * for `dim R ≤ 4` (small CM modules plus Roberts/Gillet–Soulé vanishing);
  * the Skalit and KC–Soto Levins cases.
* #193's new content lies in ramified mixed characteristic with `dim R ≥ 5` and a prime quotient of dimension `≥ 3`. (Corrected after audit B: the earlier "dimension at least 4" was not sharp.)
* Every intersection in our fronts is a curve, a complete intersection, or over an unramified base.
* For two "orbit curves" Serre's `χ` is 0.
* Recorded as a guardrail.

**What the lane produced instead.**
* **THM-4567: THM-1300's Keller map at the prime 2.**
  * `F mod 2` is purely inseparable of degree 2.
  * The 2-adic fibre counts depend on `w mod 8`. The image of `Z_2³` has measure exactly 11/32.
  * Merge-partner law `7/16, 3/8, 3/16`.
  * Integer triple collisions for every odd `w` (e.g. `F(1,−1,5) = F(−1,2,8) = F(0,2,−16)`).
* **A gap lemma and a new sporadic collision** (THM-4555 update).
  * The gap lemma makes Collatz collisions decidable length by length.
  * The collision is `(2,2,10,a) ~ (6,3,2,1,a+2)`, with value `8207/2^(13+a)`. It gives a real switch: 53803 and 26901 meet at `U⁵ = 25`.
* **Idoneal labels in #004.** The contact points of #004's height bound are CM points labelled 30, 42, 70, 105, 210, all idoneal. This is forced by the paper's own Riemann–Hurwitz count (PROVED modulo Shimura reciprocity).
* **#197's Fano gadget** is the closed in-neighbourhood design of the Paley tournament `P_7 = T_3`. Every doubly regular tournament of order `n ≡ 7 (mod 8)` gives an even gadget of the same kind: a 2-`(n, (n+1)/2, (n+1)/4)` design with even blocks and intersections, spanning a self-orthogonal `F_2`-code of dimension `(n−1)/2` (proved in audit B). DICTIONARY.

## 6. Free group factors, Thompson's F, amenability (#287, #248, #251, #253, #258, #288)

**No direct transfer.** The groups carrying the Collatz pair relation are solvable, hence amenable: `BS(1,2)`, `BS(1,3)` and `Z[1/6] ⋊ ⟨2,3⟩`. #248 needs breakpoints, and #287 needs freeness.

**What came out of it.**
* Changing the clock from `BS(1,2)` (odd steps) to `BS(1,3)` (Terras steps) produced THM-4569, and through it the frame for THM-4581.
* **The i.i.d. contrast.**
  * A genuine i.i.d. random walk on `BS(1,2)` with the same height law returns to the identity with probability like `exp(−c n^(1/3))`; numerically `2.5·10^(−9)` at `n = 256`, given height 0.
  * The Collatz pair merges at a `T^(−1/2)` rate, because its fibre coordinate contracts and is slaved to the height. THM-4581's weight `s^|k|` makes that slaving quantitative.

**Dictionary.**
* The Collatz orbit relation is hyperfinite of type III_(1/2), with `L(R_C) ≅ R_(1/2)`. Its Deaconu–Renault groupoid, the full-2-shift groupoid, has C*-algebra `O_2` and topological full group Thompson's `V`.
* #287 asks whether adjoining a generator enlarges `L(F_n)`. Our index question asked whether adjoining `y ↦ y + 1` enlarges the Collatz relation. THM-4581 answers it: **no**, the index is 1.
* The Stein groups `F_(2,3)` contain `F`, so they are nonamenable modulo #248.
* Compression by `t = 2/3` sends `r − 1` to `(9/4)(r − 1)` (NUMEROLOGY).

## 7. You's GCH theorem and the pick lane's own picks

**You, arXiv:2610.07537.** If `κ` is strongly compact and GCH holds below `κ`, then GCH holds. This answers Woodin's question.
* The engine is a cover `D` of `j[λ]` of minimal `M`-cardinality.
* Its key step (Claim 3.1.1) is that `j` preserves the first difference `min(X Δ Y)`. A family is then coded by its values on a small set of difference levels, which gives `2^λ ≤ (2^(<λ))^+`.

**DICTIONARY (exact core).**

| You | Collatz / our work |
|---|---|
| `j` preserves `min(X Δ Y)` | the parity map `Q` preserves `v_2(x − y)` (Terras/Lagarias, KNOWN); THM-4560's ruler words; LTE on the Mersenne family `v_2(p_a − p_(a')) = 1 + v_2(3^(a−a') − 1)` |
| coding by difference levels | `u` is coded by `v` plus its flip set |
| a small coding set | a finite flip set, which is equivalent to merging (THM-4581) |

**Erdős 592.**
* Delta-colourings give the classical negative relations.
* With Hajnal 1971, You's theorem gives `(λ⁺)² ↛ ((λ⁺)², 3)²` at every regular `λ` above a strongly compact with GCH below. This is PROVED modulo arXiv:2610.07537 (owner directive).
* Erdős #1169 (`ω_1² ↛ (ω_1², 3)²`) remains open in ZFC.
* Countable 592 is `Π¹_2`, hence absolute: NO LEVERAGE.

**Picks.**
* **#238 (Thorp shuffle, `Θ(d)` mixing).** It uses the same energy identity as `E k_n² = k_0² + E Σσ`; the flip skeleton is an exact simple random walk. Thorp's conditional uniformity fails here, because `v` determines `u`. Merging, not mixing, is the right theorem: ANALOGY.
* **#145 (Rokhlin).** Haar Collatz is Bernoulli: NO LEVERAGE.
* **#017 (`μ(π) = 2`).** It uses exact periods of `exp`, while `log_2 3` needs growing heights: NO LEVERAGE.

## 8. Synthesis: two clocks reading one object

Each front of this session turned on the same structure. One object is read by two clocks. Each clock alone is equidistributed and harmless, and the content lies in the coupling.

| Front | Object | Clock 1 | Clock 2 | Coupling |
|---|---|---|---|---|
| Two orbits | the 2-adic tape of `y` | `y`'s reader | `x`'s reader, `L` digits behind | window overlaps (THM-4564/4565); merging (THM-4581) |
| Terras vs odd-step time | one pair | Terras steps (`BS(1,3)`) | odd steps (`BS(1,2)`) | the flip clock; `k` is an SRW on it |
| Debt | `D = 3Δ + 1 − 2^L` | 2-adic valuation | archimedean size (Kesten perpetuity, tail index 1) | THM-4581's weight `|f|^θ s^|k|` couples them |
| Zeroless powers | `2^n` | trailing digits (5-adic clock, period `4·5^(k−1)`) | leading digits (rotation by `log_10 2`) | the middle digits; no finite-digit argument (THM-4580) |
| GCH (You) | subsets of `λ` | `j`'s first-difference level | `M`-cardinality of covers | minimal covers |
| 9/4 | multiplication tensors | rank (2-adic-like slicing) | profile growth `h ↦ 3h + a − 1` | the fixed point `−1/2` |

**The lesson for Collatz.** THM-4581 succeeded where the odd-step analysis stalled because it found a *clock-weighted* norm, `|f|^θ s^|k|`: the archimedean size is discounted by the position of the other clock. The 2-adic clock (`k`) is a fair random walk, and the archimedean size is multiplied by 3/2 exactly when that walk moves toward the merge. The analogous weight for the decimal problem would discount leading-digit freedom by trailing-digit depth. That is an open suggestion, not a result.

## 9. Concurrent work on the same prompt, and how it fits

* **opus-2026-10-07-S20 (THM-4565).**
  * Two readers at any offset: the later reader is fresh, the other reads `min(fresh, M)`, and `Cov = E[2 − 6·2^(−M); overlap]` at every pair.
  * THM-4564 (2)–(3) is its `λ = 0` case, and THM-4564's title should read "window overlaps".
  * S20's long-gap tables (3.7M overlaps) agree with ours.
* **codex-tiling.**
  * Its overlap-kernel integration audited THM-4564 and HYP-9218: the equality branch of the lockstep rule, the full-past conditioning, and the actual vs comparison perpetuity. It also sharpened the long-gap bound to twice the overlap probability. All corrections are accepted (THM-4564 header, HYP-9218, MISTAKES ledger).
  * Its information-theoretic synthesis shows that a source observation using `M` bits is independent of the whole `X` tail after `v + M` steps. This is the information form of §1.3's "uncorrelated until merge".
  * Its decimal-window note complements THM-4580.
* **This session's additions on top:**
  * the Terras-clock coalescence theorem (THM-4581) with its rate and corollaries;
  * the orbit-equivalence form (THM-4569 (7)).
  * Neither concurrent note claimed merging or recurrence before THM-4581.
* **After THM-4581 was pushed, both concurrent sessions built on it.**
  * S20 (THM-4565, Corollary 6) and codex-tiling (orbit-information synthesis) derive the long-time limit `Cov(a_s, b_(s+d)) → 2·1{d = −1}`, conditional on THM-4581. The table in §1.3 measures the same limit in correlation form: `Corr(b_s, a_(s+1)) →` merged share `→ 1`.
  * codex-tiling also audited THM-4581 independently, with two readers. The core passed. It found the same odd-count and template-clock corrections as audit A and supplied the example `p = 2417`, `q = 805` (§11).

## 10. Typing summary

| Claim | Type |
|---|---|
| Almost-sure coalescence of `u = 3^k v + e` (`e ∈ Z[1/3]`) | PROVED (THM-4581 (3)) |
| `c T^(−1/2) ≤ q(T) ≤ C T^(−1/2)(log T)²` | PROVED at sketch level (THM-4581 (4)) |
| HYP-9213, HYP-9214, HYP-9220; index 1 | PROVED (corollaries) |
| Constant `c₁ = (|k_0| + E[J])·√(4/π)` | HEURISTIC + NUMERICAL |
| Recurrence of the odd-count difference | PROVED (THM-4569) |
| Window-overlap law, covariance kernel | PROVED (THM-4564/4565) |
| Long-gap Haar depth, `Var/K = 4` | NUMERICAL |
| Zeroless: tree, recursion, brackets, obstruction | PROVED; verification FINITE-EXACT |
| Zeroless: growth exactly 9/2 | HYP-9219 (NUMERICAL) |
| Idoneal completeness, HCF planes | PROVED modulo #003 |
| 9/4 ceiling and least profile | PROVED; finite-size bound modulo #107 |
| `(3/2)²` coincidences | NUMEROLOGY |
| GCH dictionary | DICTIONARY; the Erdős relation is PROVED modulo arXiv:2610.07537 + Hajnal |
| Serre | guardrail; no leverage |

## 11. Audit record

Three independent adversarial audits were launched at checkpoint `7439e77cd7`. Corrections are recorded in MISTAKE-583.

**Audit A: THM-4581 and its corollaries.**
* Verdict: **statement 3 (almost-sure coalescence) CORRECT**; statements 1, 2 and 5 correct; statement 4 sound as a sketch; corollaries 6(a)–(f) hold.
* Independent checks:
  * the chain against direct orbits on 660,000 steps over 15 starts;
  * the local lemmas on 700,048 (state, bit) pairs, including the one-step drift at `θ = 1/2`;
  * Green-function visits about 2 per level, return ratio 0.58–0.63;
  * descent from every `(0, e)` with `|e| ≤ 10^4`;
  * GMP simulation of 400,000 paths: `√T q = 10.98–11.18` to `T = 1.6·10^6`, geometric `J` (ratio 0.924);
  * integer checks for 6(c)–(e).
* Corrections applied:
  * odd-step counts differ by `k_0`;
  * the 3-adic excess extension (3');
  * the grand-orbit reading of 6(b), with the explicit reduction;
  * the HYP-9214 sampler gap (odd value 1);
  * template total = merge time + `v_2(merge value)`;
  * `E[J] ≈ 9.93` for `(0,1)`, and `12.8` with constant `(|k_0| + E[J])·√(4/π) = 16.7` for the lag-1 pair;
  * "sketch level" carried downstream;
  * THM-4556's any-lag lower bound withdrawn (any-lag exponent about 0.69).
* Prior art: six searches plus Akin's paper found nothing.

**Codex-tiling's audit of THM-4581 (concurrent, independent).**
* The local Lyapunov bound, run cost, return drift and absorption pass.
* It found the same corrections, plus the clock example `2417/805` and the remark that an all-integer equal-time statement fails (`(1,2)` alternates with `(2,1)`).

**Audit C: THM-4568, THM-4569, THM-4580, HYP-9219.**
* Nothing refuted in substance.
* Confirmed with independent code:
  * the least profile, its admissibility and minimality (with the auditor's off-diagonal proof), the gadget barrier, and the constants;
  * the box certificates: B(3,30) and B(4,100) exactly, and B(9,100) gives 0.586175, 0.385315 and 0.344618;
  * `Z_1..Z_40` by two independent methods;
  * the unit formula, now PROVED for all `m`;
  * every `n ∈ [87, 10^6]` plus 6000 sampled exponents.
* Nine corrections applied:
  * THM-4568: the equality chain, the record labels (Alman et al. SODA 2025; Dupont et al. 2026 beaten at `a = 190`), and the ANALOGY tag with fixed point −1/2 vs −1;
  * THM-4569: `√(2/π) ≈ 0.797`, and `O_2`/`V` attributed to the Deaconu–Renault groupoid;
  * THM-4580: outward brackets, entries 24–42, and `957 ≤ n < 1.1·10^11`;
  * HYP-9219: "Equivalently" and the scope of the sufficient condition.

**Audit B: THM-4566, THM-4567 and the THM-4555/4512/4558 updates.**
* No wrong mathematics in the core claims.
* Confirmed with independent code:
  * the identification (plus the "F Galois ⇒ c central" step);
  * the 19 values, by enumerating all discriminants to `2.1·10^11` (exactly 101 = A003171);
  * completeness modulo #003, against the EKN and #003 texts;
  * every chromatic upper bound, and `κ(3), κ(7), κ(11) = 3, 4, 5` and `κ(19) ≤ 5` by SAT;
  * all of THM-4567: 11/32 exact by level-2 Hensel, so the partner law is now PROVED;
  * the THM-4555 gap lemma, the 8207 collision (symbolically) and the switch 53803/26901;
  * the Ellison citation (verbatim in Waldschmidt) and its exceptions.
* Corrections applied:
  * THM-4566: "no conjugation-fixed prime below 60" was false for 105 and 357; `m = 177` is exactly 4; 345, 357 and 253 are sharpened.
  * THM-4512: the constant `1.11·e^(−0.535j)` did not follow and becomes `e^(0.1)·e^(−0.5346j) ≤ 8.93·10^(−16)`.
  * This note: Borwein–Choi wording, the Serre threshold, the Fano gadget.
  * MISTAKE-583 covers these too.

## 12. Reproduction

* `04-computation/experiments/oai3_20261007_coalescence/`: `coal_check.py` (`exact | runs | drift | surv`), `lemmas_exhaustive.py`, `direct_merge.py`, `excursions_J.py`, with `.out` files.
* `04-computation/experiments/oai3_20261007_readers/{zeroless,nt,alg,mm94,groups,pick}/`: the readers' scripts and outputs.
* `04-computation/experiments/twoorbit_exponents_20261007.py`, `twoorbit_decay_perpetuity_20261007.py`, `twoorbit_conditional_haar_20261007.py` (THM-4564).
* `04-computation/experiments/oai3_20261007_tclock_vi.py` (THM-4569 value iteration).
* `04-computation/experiments/oai3_20261007_zeroless_verify.c` (independent zeroless check, `n < 2·10^9`).
* `04-computation/experiments/oai3_20261007_hcf_planes.py` (THM-4566).

## 13. Sources

* Z. You, *The generalized continuum hypothesis above a strongly compact cardinal*, arXiv:2610.07537.
* openai/math (github.com/openai/math): #003, #004, #017, #107 (Matrix-Multiplication-Nine-Fourths), #145, #193, #197, #238, #248, #287 (An isomorphism of the free group factors), and others, read-only.
* R. Terras (1976); J. Lagarias (1985, 2009 annotated bibliography); D. Bernstein, J. Lagarias, *The 3x+1 conjugacy map* (1996).
* A. Kontorovich, J. Lagarias, arXiv:0910.1944. Burson, arXiv:2005.09456 (withdrawn). L. Garner (1985).
* H. Kesten, *Acta Math.* 131 (1973); C. Goldie, *Ann. Appl. Probab.* 1 (1991); R. Durrett, *Probability: Theory and Examples* (Thm 4.3.1).
* H. Matui (2015) (topological full groups of one-sided shifts).
* J. Madore, arXiv:1509.07023; A.-S. Elsenhans, J. Klüners, F. Nicolae, arXiv:1803.02056; G. Exoo, D. Ismailescu, arXiv:1805.00157; M. Waldschmidt, arXiv:0908.4031; W. J. Ellison (1971).
* OEIS A007377, A181610, A031142, A003171, A000926. Saye (2022).
* Goldberg, *The Ultrapower Axiom and the GCH*; Apter–Dimopoulos–Usuba, arXiv:1901.05313; Baumgartner–Hajnal–Todorčević, arXiv:math/9311207; Garti–Hayut–Shelah, arXiv:2502.16625; erdosproblems.com #592, #1169.
