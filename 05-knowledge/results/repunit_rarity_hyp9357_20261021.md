# Why the repunits 31, 63, 127 are rare Hamiltonian-path counts, and the HYP-9357 test

Session opus-2026-10-21-S51. Independently audited (MISTAKE-657; scripts in `audit_S51/`).

**Owner prompt.** Dig into why repunits like 31, 63, 127 are rare. Find the prior repunit work across the repo. Loop between hypotheses and tests, analogize freely, and run the HYP-9357 test: `H(P_27)` and `H(P_31)` against the rotational circulant.

Canon: [THM-4748](../../01-canon/theorems/THM-4748-two-reversed-arcs-exact-path-counts-the-seven-mod-eight-lines-and-the-mersenne-side-three-arcs-from-transitivity.md), [THM-4749](../../01-canon/theorems/THM-4749-the-repunit-degree-law-phi-2-of-2-is-the-triangle-phi-3-of-2-is-the-gap-and-fermat-products-reverse-the-order.md), [HYP-9359](../hypotheses/HYP-9359-arithmetic-support-law-of-strong-hamiltonian-path-counts-gluings-and-order-three-symmetry-carry-all-divisibility-bias.md).

## 0. The answer in brief

The rarity is not one phenomenon. `2^k − 1` and `2^k + 1` are the two sides of `2^(2k) − 1`, and tournaments treat the sides differently in three independent ways.

1. **Complexity (PROVED, THM-4748).**
   * One reversed arc of span `s` adds exactly `2^(s−1)` paths. So `2^k + 1 = H(NT_(k+2))` is one arc away from transitivity (Moon 1972).
   * `2^k − 1 = 2^k + 1 − 2` needs a subtraction. With at most two reversed arcs the values `≡ 7 (mod 8)` lie only at 15 and on the two lines `3·2^j − 1` (Thābit) and `5·2^j − 1`. So every tournament with feedback arc number `≤ 2` misses `2^k − 1` for `k ≥ 5`.
   * Three arcs suffice for every `k`: the **Fibonacci fan**, an apex losing exactly to `a`, `b`, `a + b` of a chain with `b − a ∈ {1, 2}`, has `H = 2^k − 1`. So both sides live at order `k + 2`, the Fermat side one arc from transitivity and the Mersenne side three.
2. **Arithmetic support (NUMERICAL, HYP-9359).** Strong counts carry a bias toward multiples of 3: `log S(d) ≈ smooth + 0.22[3|d] + 0.06[9|d] + 0.05[5|d]`, and primality is neutral. Eleven predictions committed before their evaluation all passed; their blindness is partial (§3).
   * Almost all of the bias comes from two mechanisms: exact one-vertex gluings (`H = H(A)H(B)`, with `H(C_3) = 3` the cheapest factor) and automorphisms of order 3 (`|Aut|` divides `H`).
   * For odd `k`, `Φ_2(2) = 3` divides `2^k + 1` and not `2^k − 1`. That is the 31/33, 127/129, 511/513, 2047/2049 asymmetry.
3. **Gap shadows (THM-4749; the degree law is KNOWN from Session24).**
   * The two impossible values are themselves repunits: `7 = 111₂ = Φ_3(2)` and `21 = 111₄`.
   * `Φ_3(2) = 7` divides `2^k − 1` iff `3 | k` (also in mac-mini's MSG-2400).
   * 63 is a multiple of 9, but `63 = 3·21 = 9·7`, so neither support is available. Its class count is linear, `N(n, 63) = 2(n − 7)` (KNOWN, Session24), while `N(n, 65)` is quadratic. `511 = 7·73` is blocked the same way.
   * In the large-`n` limit the degree law turns all of this into exact orders of growth. `2^k − 1` is asymptotically rarer for `k = 5, 6, 7, 9, 13, 15`, and commoner for `k = 4, 8, 12, 16`:
     * 15, 255 and 65535 are products of Fermat numbers;
     * 4095 absorbs its 7 into 35.

What remains is a small-order effect, seen at `n = 8, 9` where the rarity was first noticed (S50):
* `S_7(63) = 0`: no strong 7-tournament has 63 paths;
* `S_7(127) = 1`;
* at orders 6–8, `H ≡ 7 (mod 8)` has relative density 0.47–0.78 in the second octave of the spectrum.

63 and 127 stay outliers after every arithmetic control: residuals −1.68 and −0.55 in the full model on `[25, 1098]` (`outliers.py`).

**HYP-9357.** `H(R_27) = 426990256001863046739 > H(P_27)` and `H(R_31) = 20351875726012821643802839 > H(P_31) = 19876729650849384218015865`. Prediction (3) passes at `q = 31`, and the rotational circulant also wins at `q = 27`. This equals mac-mini's independent `Z/2^128` values (`paley_extremal_20261021/R27.out`, `R31.out`).

## 1. Data and scripts

All scripts are in `04-computation/experiments/repunits_hyp9357_20261021/` and print ALL CHECKS PASSED.

| script | content | time |
|---|---|---|
| `hyp9357_rotational_vs_paley.py` | Karp inclusion–exclusion modulo three 58-bit primes, rotation-orbit and weighted contains-0 reductions; validated `n = 13..23` against `circ_counts` and Held–Karp | 9 min |
| `two_arc_formulas.py` | THM-4748 (1)–(5): 9,009 two-arc configurations against an exact DP; the mod-8 lines; the fan; the census `n ≤ 9` | 14 s |
| `three_arc_census.py` | THM-4748 (6): all three-arc configurations, `m ≤ 13` | 16 s |
| `fan_uniqueness.py` | THM-4748 (6): least-complexity realizations of `2^k − 1`, `5 ≤ k ≤ 11` | 2 min |
| `repunit_rarity.py` | the strong inventory to 2100 (723,544 classes); the model; R1–R7 and G1–G4; the gluing/automorphism decomposition; per-order residues; feedback arc number at `n = 8` | 3 min |
| `repunit_degrees.py` | THM-4749: `w(2^k ± 1)` for `k ≤ 40`; exact polynomial checks for `h < 2100` | 3 min |
| `outliers.py` | residuals of the repunits under the full model on `[25, 1098]` and `[301, 1098]` | 17 s |
| `class_residues.py` | `H` mod 3, 8, 16 over the class census `n ≤ 9`, by class and labelled | 9 s |
| `refuted_leads.py` | the leads that failed: ear valuations and residue mixing, `I(Ω, −1)`, feedback-arc shares, the one-apex toy model | 11 s |
| `PREREGISTRATION_out_of_sample.md` | R1–R7 (commit 917fef912a) and G1–G4 (commit be7c76115f) | — |
| `audit_S51/` | the independent audit: `alib51.py`, `c0`–`c9`, `REPORT.md` | — |

* **The strong inventory.**
  * Every strong class with `H ≤ 2100` is found by exact ears from `C_3`. This is complete by Moon's strong-deletion theorem and the insertion lemma, using the code of `audit_S50/b3_ear_census_625.py` with the cap raised.
  * Classes by order: 1, 6, 35, 353, 6008, 146,376, 268,683, 180,284, 88,379, 28,817, 4545, 56 at orders 4–15.
  * Below 625 it reproduces THM-4744's 13,439 / 8,769 / 3,536 / 658.
  * Every odd `3 ≤ d < 2100` except 7, 21 is a strong count.
  * The audit's own ear census and canonical form give the same 723,544 classes, with identical H multisets per order.
* **Class counts.** `N(n, h)` for all `n` and `h < 2100` come by sequences of strong components, and agree with the full census for `n ≤ 9`.

## 2. The Fermat side and the Mersenne side (THM-4748)

* **The run decomposition.** A Hamiltonian path of a near-transitive tournament is a sequence of increasing runs joined by back arcs. A non-endpoint vertex under `c` of the used arcs chooses among `c + 1` runs. So few arcs give sums of `2^x 3^y` terms.
* **Two arcs.** There are exact formulas for the disjoint and touching (product), crossing, nested and shared-endpoint configurations. The crossing case has a ternary term `2^(p+q+1)·3^r` from the doubly covered stretch.
  * Special cases are older: Moon 1972 (one arc, chains), bof on MathOverflow 234029 (the crossing case `p = r = 0`), and El Sahili–Abi Aad 2020 (§9).
* **Mod 8.** In the whole `fa ≤ 2` world, `H ≡ 7 (mod 8)` occurs only at:
  * 15;
  * `3·2^j − 1` (shared endpoint, inner span `L − 1`): 23, 47, 95, 191, …;
  * `5·2^j − 1` (inner span `L − 2`): 39, 79, 159, ….

  So for `k ≥ 5` the Mersenne numbers need three arcs. In the census, the least feedback arc number is 3 for `H = 31, 63, 127, 255` at `n ≤ 9`, against 1 for 33 and 65.
* **The fan.** In the one-apex formula
  `H − 1 = Σ_{t∈A} 2^(t−1) − Σ_{s<t} 2^max(t−s−1,1)`
  (KNOWN, the single-core formula), the in-set `A = {a, b, a+b}` cancels two singletons against two pairs.
  * The result is `H − 1 = 2^(a+b−1) − 2`.
  * The two `n = 8` classes with `H = 63` (THM-344) are the `k = 6` fan and its converse. Their signatures were already in `single-core-cycle-count.md`.
  * For `k = 6, 9, 11` (`H = 63, 511, 2047`) the fan pair is the only least-complexity realization among the strong classes below 2100. For `k = 5, 7, 8, 10` there are others (`fan_uniqueness.py`).
  * mac-mini's survey (MSG-2400) noted the one-apex pattern the same day: minimal apex in-degree 1, 2, 2, 3 for `2^k + 1`, `3·2^k − 1`, `5·2^k − 1`, `2^k − 1`, at `n ≤ 8`.
* **Three arcs.** All configurations up to `m = 13` reach `2^k − 1` for `k = 4..11` and 15. Their `7 (mod 8)` values have unbounded odd parts of `H + 1` (PROVED in THM-4748 (6): a crossing pair with `p = q = 0` and odd `r`, plus a disjoint `C_3`). So the "lines" structure is special to two arcs.
* **Multiplicity.** At `n = 8`, values first reached with one arc (33, 65) sit at log-residual +0.65 above the smooth fit. The correlation of residual with least feedback arc number is −0.19 (`repunit_rarity.py` (O)). Complexity is a real but secondary predictor.

## 3. Arithmetic support (HYP-9359)

* **The fit.** On `[301, 1098]`: `β_3 = +0.239`, `β_9 = +0.088`, `β_5 = +0.065`, `β_prime = −0.002`, `β_blocked = +0.028`, `β_(7 mod 8) = −0.027`.
* **Out of sample on `[1101, 2099]`:**
  * `+0.219, +0.060, +0.054, +0.005, +0.008, −0.015`;
  * `S(2047) = 1599 < S(2049) = 1968`;
  * `S(1365) = 1037`, above the smooth fit;
  * the order-14 `3|d` coefficient is +1.19.
* **The mechanism, out of sample.** These fits use the decomposition model, without the prime and blocked terms.
  * Removing the classes with an exact one-vertex cover lowers `β_3` to 0.135 and kills `β_5` and `β_9`.
  * Removing also the classes with `3 | |Aut|` lowers `β_3` to `0.010 ± 0.005`, about 2 s.e.
* **What is left.** The generic remainder keeps mainly `β(1 mod 8) = +0.017 ± 0.005`: a weak 2-adic attraction to `+1` rather than repulsion from `−1`.
* **How blind the tests were** (the audit's chronology, `audit_S51/c7`):
  * R1–R7 were committed 9.6 s, by local timestamps, before the order-≥10 test data were complete. The remote push came later.
  * The order-9 slice of the test range (16% of its classes) had been explored earlier through the `n = 9` class census.
  * G1–G4 were committed after the test-range `S(d)` and the R verdicts were known, and before the decomposition was run.
* **Census residues** (`class_residues.py`).
  * At `n = 9` the class counts by `H mod 3` are 66,320 / 62,129 / 63,087, while the labelled shares are within 2% of uniform.
  * The labelled relative density of `H ≡ 7 (mod 8)` (share × 4) is 0.25, 0.80, 0.81, 0.94, 0.993 at `n = 5..9` (HYP-9355's decaying pull).

## 4. The large-n degree law (THM-4749)

* **The law.** `N(n, h)` is a polynomial of degree `w(h)` (KNOWN, Session24), where `w(h)` is the longest factorization into strong values.
* **Repunit consequences.**
  * `w(31) = w(63) = w(127) = w(511) = w(8191) = 1`, so these counts are linear in `n`.
  * `w(33) = w(65) = w(129) = 2`, `w(513) = 4`.
  * `w(15) = 2` against `w(17) = 1`; `w(255) = 3` against `w(257) = 1`.
  * `2^(2^j) − 1 = F_0 ⋯ F_(j−1)` gives `w ≥ j` unconditionally, since every Fermat number is `H(NT)`. This is mac-mini's "Fermat monoid" (MSG-2400).
* **Ties.** At `k = 10, 11, 14` the degrees tie.
  * The exact polynomials decide `k = 10` for 1023.
  * They decide `k = 11` for 2049: `2·S(23)·S(89) = 184 < 2·S(683) = 516`.
  * An unbalanced factorization `3 × large` beats a balanced one.
* **The record.** At `n = 14`, for `(2^k − 1, 2^k + 1)` with `k = 5..11`: 34 vs 125, 14 vs 151, 82 vs 424, 1176 vs 455, 946 vs 2427, 4770 vs 3133, 7426 vs 10168.
* **Certification.** The value set `V` (all odd numbers except 7 and 21) is certified below 80,265; the conjecture goes back to bof and Royle, MathOverflow 232751 (2016). So the degrees are certified for every `k` whose prime factors are below 80,265: all `k ≤ 16`, and also `k = 18, 20, 21, 22, 24, 25, 26, 30, 36` below 40.

## 5. Small orders: why n = 8, 9 showed the effect most strongly

* **The decomposition.** For `h > 45`, `N(8, h) = 2·S_7(h) + S_8(h)`. Likewise `N(9, h) = 3·S_7(h) + 2·S_8(h) + S_9(h)`, plus a few products. The S51 wide-net agent noted this, and it reproduces the S50 histograms exactly.
* **The order-7 strong spectrum.** It has holes at 63, 107, 119, 149 in `(45, 150)` (S50). Its thin values include 61, 85, 87, 125, 127, each attained by one self-converse class.
  * Since every strong 7-tournament is an ear of a strong 6-tournament, `S_7(63) = 0` says that no out-set of a new vertex over any strong 6-tournament gives 63.
* **The second octave.** By position `H/f(m)` in the strong spectrum, `H ≡ 7 (mod 8)` has relative density 0.47, 0.53, 0.78, 0.89 in the octave `[2f, 4f)` at orders 6–9. It is about 1 in `[4f, 8f)` (`repunit_rarity.py` (O)). 63 at order 7 (`63/25 = 2.5`) and 127 at order 8 (`2.8`) sit in that octave.
* **Unexplained.** Neither the depletion nor the exact holes (63 at order 7) are explained (open in HYP-9359). The tempting reading through the `fa ≤ 2` lines fails (§8).

## 6. HYP-9357: the rotational circulant against Paley at q = 27, 31

`R_q = Z_q{1, …, (q−1)/2}`. Values of `H·2^(n−1)/n!`:

| `n` | rotational | Paley | `H(R)/H(P)` |
|---|---|---|---|
| 19 | 2.552 | 2.527 | 1.0098 |
| 23 | 2.598 | 2.557 | 1.0159 |
| 27 | 2.632 | 2.579 | 1.0205 |
| 31 | 2.658 | 2.596 | 1.0239 |

The gap widens. HYP-9357 prediction (3) holds at the prime `q = 31`, and the prime power 27 behaves the same way. mac-mini records the values in HYP-9357, its own file; its `Z/2^128` values equal this session's three-prime modular count.

## 7. What the objects represent

The first bullet is PROVED; the rest are ANALOGY unless marked.

* **Two sides of one number.** `2^(2k) − 1 = (2^k − 1)(2^k + 1) = ∏_{d | 2k} Φ_d(2)`.
  * The Fermat side (even `d` not dividing `k`) is natively a single back arc.
  * The Mersenne side (`d | k`) is natively a three-arc cancellation.
  * The cheapest strong tournament is `Φ_2(2) = 3`, and the impossible one is `Φ_3(2) = 7`.
* **The gaps are repunits.** 7 and 21 are the smallest binary and base-4 repunits with no Fermat-side factorization, and `7 = I(K_3, 2)` and `21 = I(P_4, 2)` are exactly the forbidden conflict graphs (THM-200/202). The Jacobsthal–Lucas transfer matrix `[[1,2],[1,0]]` has eigenvalues 2 and −1; the repunits are its "−1 branch" (S51 wide-net, A4).
* **The tiling cube** (mac-mini's Möbius lane). For one-apex tournaments, singleton coefficients are `+2^(s−1)` and apex-pair coefficients are `−2^max(gap,1)`, and the expansion stops at pairs. The fan is the minimal exact cancellation.
* **Collatz echo.**
  * mac-mini's unifying reading: a base-`b` repunit is a truncation of the fixed point `−1/(b−1)` of `x ↦ bx + 1`, and `−1` is the Collatz odd-step fixed point.
  * The Collatz step sends `2^k − 1` to `3·2^(k−1) − 1`, the Thābit line on which the two-arc world does reach `7 (mod 8)`. The second step gives `9·2^(k−2) − 1`, which that world does not reach.
  * This is typed NUMEROLOGY beyond the first step.
* **Context** (CITED, not used). Bar-Noy and Naor (SIAM J. Discrete Math. 3 (1990) 7–20): `H(T)` equals the number of minimal feedback arc sets.

## 8. Refuted or withdrawn this session

* **2-adic inheritance along ears** (lead 1 of the catalog). The ear increments `ΔH` have geometric 2-adic valuations (1: 0.51, 2: 0.25, 3: 0.12, …) at orders 5–7. The child residue mod 8 is close to uniform for every parent residue at orders 6→7 and 7→8: total-variation distance at most 0.09 and 0.023 (`refuted_leads.py` (a)). At 5→6 it is rougher (up to 0.23), but there are only 6 parents. Withdrawn.
* **Contractible independence complexes.** `H ≡ I(Ω, −1) (mod 3)` holds. But `I(Ω, −1)` lies between −162 and −3 for the strong classes of orders 5–8 and is never 0, so it cannot carry the mod-3 bias (`refuted_leads.py` (b)).
* **Low feedback arc number as the source of the 7 (mod 8) depletion.** The share of `H ≡ 7 (mod 8)` by least feedback arc number is not consistently low (`refuted_leads.py` (c)):
  * order 6: 0.462 at `fa = 2`;
  * orders 7–8: 0.11–0.22 at `2 ≤ fa ≤ 4`;
  * order 9: 0.082 at `fa = 2` but 0.332 at `fa = 3`.

  The `fa ≤ 2` lines of THM-4748 (2) are a proved fact that does not explain the depletion.
* **A one-apex toy model of the multiplicities.** The word multiplicities are at most 8 per value; there is no signal (`refuted_leads.py` (d)).
* **Binary weight of `H − 1`.** It is a proxy for trailing ones (mac-mini), and both are secondary to the 3-adic support at `d > 300`.
* **This session's own draft claims, corrected:**
  * before canon: "fa ≤ 2 never gives `7 (mod 8)`" (the Thābit lines are counterexamples) and "odd `k` gives a strict degree gap" (it ties at `k = 11`);
  * after the audit: the remark "subtractions come only from arcs sharing an endpoint" (false outside one-apex tournaments), the THM-4749 title's KNOWN typing of the value set, and the overstated blindness of the pre-registration (MISTAKE-657).
* **Primality as the reason** (S50 audit, a8) — confirmed neutral here, with `β_prime ≈ 0`.

## 9. Prior work

This is a condensed catalog. The full lists come from the S51 catalog agent, mac-mini's MSG-2399 and the S51 audit.

* **Literature.**
  * J. W. Moon, *Publ. Math. Debrecen* 19 (1972) 101–104: the single-arc count `1 + 2^(n−2)` and the chain product.
  * bof, MathOverflow 234029 (2016): two reversed arcs `v_1v_3`, `v_2v_n` give `3·2^(n−4) + 3`, the crossing case `p = r = 0`. Math.SE 1688155 gives the converse case.
  * A. El Sahili and M. Abi Aad, *Discrete Appl. Math.* 277 (2020) 115–126: exact `H` for tournaments obtained from the nearly transitive tournament by reversing further arcs, and tournaments whose `H` is a power of 3.
  * M. Abi Aad, *Discrete Appl. Math.* 394 (2026) 220–228. Its Proposition 1 overlaps the one-apex formula.
  * The full text of neither DAM paper was accessed (ScienceDirect CAPTCHA). Priority for the general two-arc formulas is unsettled.
  * bof and Royle, MathOverflow 232751 (2016): the origin of the "only 7 and 21" conjecture, with every odd value below 80,557 except 7, 21 found on at most 11 vertices.
  * None of these gives the `fa ≤ 2` seven-mod-eight lines, the three-arc Mersenne theorem or the Fibonacci fan for all `k`. The same holds for El Sahili–Ghazo Hanna (arXiv 2101.00713, 2207.11510), the papers citing them, and Moon's *Topics on Tournaments*.
* **The gaps (repo).**
  * THM-1370 (`THM-1370-h-spectrum-omits-7-21-all-n`): 7 and 21 are never `H`.
  * THM-200/343: `H = 7` iff `Ω = K_3`.
  * THM-075/115/202: `H = 21`, `P_4`.
  * THM-517: the cluster reading.
  * HYP-3616/3617: `{7, 21} = {Φ_3(2), 3Φ_3(2)}`.
* **The 63 atoms.**
  * THM-344: two classes at `n = 8`, `Ω = K_31`, one-apex. Their signatures 1001100 and 1100110 are in `05-knowledge/variables/single-core-cycle-count.md`; they are the `k = 6` fan and its converse.
  * HYP-1754 was refuted (MISTAKE-024/050), and MISTAKE-028 refuted "31, 127 forbidden". HYP-1754 lives only in `INDEX-HISTORICAL-THROUGH-2026-07-21.md`.
  * mac-mini's MISTAKE-653 corrected three repunit errata the same day.
* **Digits.** THM-466 (the 2-adic digits of `H` are odd-cycle collection counts), THM-478, THM-1980, THM-791 (at even `n`, self-converse implies `H ≡ 1 mod 4`), HYP-9355.
* **Families.**
  * THM-250: the single arc.
  * `deep_extension_s20al.out` §2 (2026-03-22): all 45 two-arc values at `n = 6`, without a closed form.
  * `arithmetic_seams_20260921_tournaments.md` eq. (12): `H(E_n) = F_r` at `n = 2^r + 2`.
  * The single-core formula and HYP-2198/2199 (one-apex gaps `r = 3, 10`, i.e. `H = 7, 21`).
  * HYP-1762 (one-apex counts of `r = 31, 42, 63`). It is in the historical index only, and its ID collided.
  * Session16's feedback-word spectrum (`H(1^a 0^b) = (2^a − 1)(2^b − 1) + 2`), Session25's modular completeness of one-apex values, Session26's one-apex classification and product law.
  * THM-1355/1360 (Proth/Fermat tables), THM-1425 (Jacobsthal values).
* **Multiplicities.** S50 (`h_spectrum_multiplicities_20261021.md`) and Session24 (the degree law, the 7-adic defect, `N(n, 63) = 2(n − 7)`).
* **Collatz repunits.** THM-4473 (`T^k(2^k − 1) = 3^k − 1`), THM-4555/4556 (the Mersenne line), HYP-9333 (terminal channels 5, 85, 341), and the Mersenne-line exits.

## 10. Open

1. **The order-7 hole at 63.** Find a structural proof or reduction.
2. **The second-octave depletion and the residual 2-adic attraction.** Explain both (HYP-9359 (3)).
3. **A closed form for the support coefficients.** Count exact `C_3` gluings and order-3 symmetric classes per value (HYP-9359 (2)).
4. **Feedback number versus multiplicity.** Is `N(n, h)` controlled by the number of `fa ≤ f` realizations of `h` for a fixed `f`, uniformly in `n`?
5. **Large `k`.** Decide `V`-membership of the large prime factors (the 2016 bof–Royle conjecture) and so remove the remaining conditions of THM-4749.
6. **The DAM papers.** Read El Sahili–Abi Aad 2020 and Abi Aad 2026 in full, to settle priority for the two-arc formulas.

## 11. Audit trail

The independent audit (`audit_S51/`) reproduced every number with code independent of the session's:
* the two-arc formulas for 18,564 configurations (`m ≤ 16`);
* the `fa ≤ 2` lines for every tournament with `m ≤ 16`;
* the fan to `k = 18` by DP and to `k = 3000` by an OCF-free count;
* the census `n ≤ 9`;
* the inventory to 2100, with identical H multisets;
* the regression and decomposition, to 0.0005;
* the Session24 polynomials and leading coefficients;
* all eleven verdicts.

It found two false sentences in a THM-4748 remark, a wrong KNOWN typing in THM-4749's title, overstated pre-registration claims, missing prior art (Moon 1972, bof 2016, El Sahili as co-author, repo THM-250 and the signatures of the `k = 6` fan) and missing credits (mac-mini MSG-2400). All are corrected; see MISTAKE-657.
