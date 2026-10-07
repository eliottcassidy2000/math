# Audit C: THM-4568 (9/4 lemma), THM-4569 (Terras clock), THM-4580 + HYP-9219 (zeroless 2^n)

Auditor: independent adversarial audit, 2026-10-07. Repo was read-only (`math-wt-chessboard-20261006`, HEAD 7b6fc1d53a; the four audited files did not change during the audit).
Scripts and outputs are in this directory (`a1_*` for THM-4568, `a2_*` for THM-4569, `a3_*` for THM-4580). #107 was fetched read-only (`mm107_paper.tex`). OEIS data was fetched with curl (`A0*.json`, `b0*.txt`).

**Bottom line.**
* No statement is REFUTED in substance. Every mathematical core survives.
* There are 9 required corrections. Each is a precision, rounding or labelling error, and several would let a later reader infer something false:
  * an inward-rounded "PROVED" bracket;
  * an "n ≥ 957" claim that is false beyond the verified range;
  * a lower-bound constant rounded up;
  * an equality claim that fails at two steps;
  * a mislabelled prior record.

---

## THM-4568: the 9/4 growth lemma is exactly sharp

I checked #107's Lemma 4.1 (`lem:growth`) hypotheses verbatim. They are: `P: N²→R_{>0}` symmetric, `P(1,b)=b`, concavity for `a≥1, b≥2`, and tripling for `a,h≥1`. These are exactly the hypotheses THM-4568 uses.

**(1) Least element: CONFIRMED.**
* **Admissibility.** Exact rational checks pass (`a1_pmin.py`): (B) for `b<560`; (C) for `a<140, b<560`; (T) for `a,h<140`. There are 0 violations.
  * There is also a short proof. Below the diagonal the increments are `σ_b(2b+a)/(3b)`. They are nonincreasing iff `b+1≤a`, and the increment from `a−1` to `a` equals `σ_a`, which is the slope above the diagonal.
* **Minimality off the diagonal (my proof; the file states only the diagonal LP).**
  * Tripling at a general `h`, plus concavity, gives `Δ_{a,h} ≥ 2P(a,h)/(2h+a−1)`.
  * Starting from the paper's diagonal bound `P(a,a) ≥ σ_a(3a−1)/2`, induction in `h` gives `Δ_{a,h} ≥ σ_a`. Hence `P(a,b) ≥ σ_a(2b+a−1)/2` for `b≥a`.
  * An LP at 8 off-diagonal points matches `P_min` to 1e−13 (`a1_lp_offdiag.out`).
* **Identities.**
  * The generating function `(1+x)(1−x)^{−7/3}` is verified to `a<60`.
  * `9^{m−1}σ_m` = 1, 12, 126, 1260, 12285, …, which equals OEIS A004991 (fetched).
  * The asymptotic constant `3/(2Γ(4/3)) = 1.67977` is confirmed.

**(2) Exact ceiling: CONFIRMED WITH CORRECTION.**
* **What holds.**
  * `P_min ≤ (a+b−1)^{4/3}` holds exactly for `a<140, b<560`, with equality only at (1,1).
  * It is also provable for all `a,b`:
    * Wendel's inequality gives `Γ(m+1/3)/Γ(m) ≤ m^{1/3}`, so `D_a ≤ 1.6798 a^{4/3} < (2a−1)^{4/3}` for `a≥2`.
    * Off the diagonal, the difference `s^{4/3} − σ_a(s−(a−1)/2)` (with `s=a+b−1`) is increasing, because `(4/3)s^{1/3} > 1.12a^{1/3} ≥ σ_a`.
* **What is false.** "Every inequality in the paper's chain is an equality on `P_min`" fails at the last two steps:
  * `(1+1/(3m))³ = 1 + 1/m + 1/(3m²) + 1/(27m³) > 1+1/m`;
  * `(3a−1)/2 > a` for `a≥2`.
  * Consequently `P_min(2,2) = 10/3 > 2^{4/3}`, and `P_min(a,a)/a^{4/3} → 1.68`. The mm94 script's check [2] covers only the steps up to `H_a ≥ Π(1+1/(3m))`.
* **Replacement:** "Every inequality in the chain up to `H_a ≥ Π_{m<a}(1+1/(3m))` is an equality on `P_min`. The two final simplifications, `(1+1/(3m))³ ≥ 1+1/m` and `(3a−1)/2 ≥ a`, are strict. They cost only the constant `3/(2Γ(4/3)) ≈ 1.68`, not the exponent."

**(3) Gadget barrier: CONFIRMED.**
* **The shift bound.** Flattening ranks are characters (#107 says so). Their symmetrized profile is `√(ab(a+b−1))`. Requiring `(mh+c)(mh+c+a−1) ≥ m²h(h+a−1)` as `h→∞` gives exactly `c ≥ (m−1)(a−1)/2`.
* **"Pinned at 9/4".**
  * `P_min` satisfies every minimal-shift gadget for `m≤10` and `a,h<120`, checked exactly with 0 violations.
  * There is a general proof whenever `mh+c ≥ a−1`, which covers all `m≥3`. There `P_min(a,·)` equals its tangent line `σ_a(b+(a−1)/2)`, and the concave `P_min` lies below that line.
  * The only remaining case is `m=2` with `2h+c<a−1`. It is checked numerically for `a<6000`, with a minimum ratio of 1.0583; the asymptotic margin is also 1.058.
* **Diagonal-only gadgets.** The formula `κ=2(m−1)/(u−1)` and the floor 2.0584 (at `m=2, u=2.3723`) are confirmed as an asymptotic derivation backed by the LP.
  * Optional: type this as "derivation + LP", since no admissibility proof of the extremal "fake profile" is written.

**(4) Fixed point: the mathematics is CONFIRMED; the typing needs a CORRECTION.**
* The conjugacy, `s* = −1/2` and `κ* = 2/(1−s*) = 4/3` are all correct. So is `D_a^min = (a − b*)σ_a` with `b* = −(a−1)/2`.
* **Issue 1: two different fixed points.** The bullet "This is the Schröder linearization of the Collatz odd step, `T(n)+1=(3/2)(n+1)`" pairs `s* = −1/2`, the fixed point of `3h+1`, with an identity whose fixed point is `−1`.
* **Issue 2: typing.** The bullet sits under "Statements (PROVED)", while the status line calls it ANALOGY.
* **Replacement:** "(ANALOGY) The unhalved odd step satisfies `3h+1+1/2 = 3(h+1/2)`, with fixed point `−1/2` (THM-4555). In the halved normalization this is `T(n)+1 = (3/2)(n+1)`, with fixed point `−1` (THM-4556)."

**(5) Finite size: CONFIRMED WITH CORRECTION.**
* **Values (mpmath, 60 digits).**

  | quantity | value |
  |---|---|
  | `ω_2` | 2.737468 |
  | `ω_187` | 2.371449 |
  | `ω_188` | 2.371335 |
  | `ω_189` | 2.371223 |
  | `ω_190` | 2.371111 |
  | `ω_595974` | 2.30000000332 |
  | `ω_595975` | 2.29999999719 |
  | `K` | 0.684347890891 |

* **Monotonicity.** `ω_a` is strictly decreasing on `[2,5000]` and at 760 log-spaced points up to `10^41`.
* **The constant `K`.** The second-order term predicts `(ω_a−9/4)ln a ≈ K − 0.266/ln a`. This matches 0.6805 at `a=10^30`.
* **Correction: the prior record is mislabelled.** 2.371339 is Alman–Duan–Vassilevska Williams–Xu–Xu–Zhou, SODA 2025 ([AD25] in #107). #107's own introduction names Dupont et al. 2026, `ω<2.371177`, as the latest preceding bound.
  * **Replacement:** "`ω_188 < 2.371339` (Alman et al., SODA 2025); `ω_190 < 2.371177` (Dupont et al. 2026, the last bound before #107 that #107 cites)."
  * The same fix applies to the source note `oai3_two_orbits_twos_and_threes_20261007.md`, line 237.

**Numerology register: CONFIRMED.** It is correctly typed NUMEROLOGY, and its arithmetic is right (`1/4+(1/3)(9/4)=1`; mean offspring 4/3; drift 3/4).

---

## THM-4569: the Terras clock

**(1) Transitions: CONFIRMED.**
* All four rows were re-derived by hand.
* My Fraction-based chain, with transitions coded directly from the table, was compared with direct 2-adic iteration: 100,000 steps over 250 random relations with `j0∈[−5,5]` and `c0∈Z[1/3]`, 0 failures (`a2_tclock.out`).
* **Minor wording.** "in `Z[1/3]` for rational initial relations" should read "when `c_0 ∈ Z[1/3]`". For example, `c_0 = 1/5` stays in `Z_(2) \ Z[1/3]`.

**(2) Random-walk representation, and (3) Dichotomy: CONFIRMED.**
* **Merges happen at equal time.** Given the past, `y_s` is Haar. On each parity cylinder `T^k` is affine with slope `3^w/2^k ≠ 3^j`. So `T^k(y_s) = 3^{j_s}y_s + c_s` holds at one point only, and similarly `x_s = y_s` with `j_s≠0` is a single point.
* **Terras injectivity** gives `N_s → ∞` on the event that the orbits never merge.
* Optional skipping applies because `ε_s` is predictable.

**(4) Rate lower bound: CONFIRMED WITH CORRECTION.**
* **Lag-1 pair: the constant is valid.**
  * `q_1 = 2·3^{a−2} − 1` is Haar on `5+16Z_2`, with forced parities 1,0,0,0.
  * The chain runs `(1,2) → (1,2) → (1,1) → (2,2) → (2,1)`, so `j_* = 2`.
  * Exact binomials give `√T·P(SRW avoids −2) → 2√(2/π) = 1.59577` (`a2_srw_const.out`). This is at least the claimed 1.59, so the claim holds.
* **The pair `x = y+1`: the constant is rounded up.**
  * `(0,1)` goes to `(1,2)` or `(−1,1/3)`, so `j_* = 1`, and the argument gives only `√(2/π) = 0.79788`.
  * "(0.80 + o(1))" is a lower bound rounded up. **Replacement:** "`(√(2/π) + o(1)) T^{−1/2} ≈ 0.797 T^{−1/2}`".

**(5) Box reduction and certificates: CONFIRMED.**
* **Validity.** `V_n` = P(absorbed at (0,0) within `n` steps before leaving the box) ≤ P(merge), since the coins are exactly fair (Terras) and (0,0) is exactly the merge state.
* **Reproduction with independent code** (state `(j, Fraction c)`, not the reader's `N`-recurrence):

  | box | (0,1) | (2,1) | (1,1) | min B(1,10) |
  |---|---|---|---|---|
  | B(3,30) | 0.457151 | 0.193941 | 0.594232 | 0.130370 |
  | B(4,100) | 0.490248 | 0.243208 | 0.618992 | 0.194891 |

  * Both rows match the reader exactly.
  * A sparse linear solve gives the same values, so the iteration has converged.
* **The headline box.** I also ran a separate vectorized implementation on B(9,100): 11,809,419 states, 3200 sweeps, 121 s, 0.76 GB (`a2_vi_9_100.out`).

  | start | value | claimed |
  |---|---|---|
  | (0,1) | 0.586175 | ≥ 0.5861 |
  | (2,1) | 0.385315 | ≥ 0.3853 |
  | min over B(1,10) | 0.344618 | ≥ 0.3446 |

  * The margins are at least 1.5e−5, and float64 error is about 1e−12.
* **The Lévy 0–1 step** is valid with δ = 0.3446.
* **Density transfer.**
  * `b ↦ 9^b` is a measure isomorphism `Z_2 → 1+8Z_2`, and merge-by-`t` is a finite union of residue classes.
  * With σ = number of odd steps (HYP-9213), an equal-time chain merge gives `σ(2^a−1) = σ(2^{a−1}−1)`. So the lower density is at least 0.3853, as claimed.

**(6) The real place: CONFIRMED with a clarifier.** Rows `(p,ε) = (0,1)` and `(1,0)` carry additive errors `1/(2·3^{j+1})` and `1/(2·3^j)`; rows (0,0) and (1,1) are exact. So the statement holds as `j → +∞`. Add "for `j ≥ 0` (error `≤ 3^{−j}/2`)". Mean multiplier 1 and `E[A]=1` at κ=1 are correct.

**(7) Orbit-equivalence form: CONFIRMED WITH CORRECTION (wording).**
* **What holds.**
  * `R_C ⊆ R_A`.
  * Both are ergodic (Kolmogorov 0–1 for the tail relation), hyperfinite (Dougherty–Jackson–Kechris; amenable `Γ_C`) and type III_{1/2} (RN cocycle in `2^Z` via `|2^a3^b|_2`; the pairs `(x, Tx)` put 2 in the ratio set).
  * The index is a.e. constant.
  * The equivalences hold. For "`y~y+1` ⇒ index 1", write `g(y) = 2^a3^b y + r` and route through `2^{max(0,−a)}`, so that `r ∈ Z[1/3]` and every intermediate point stays in `Z_2`.
* **What is inaccurate.** A measured relation has no "C*-algebra `O_2`" and no "topological full group". Those belong to the Deaconu–Renault groupoid `{(x, m−n, y) : T^m x = T^n y}`, which also carries isotropy at eventually periodic points.
* **Replacement:** "Via Lagarias's conjugacy, the Deaconu–Renault groupoid of `(Z_2,T)` is the full-2-shift groupoid. Its C*-algebra is `O_2` and its topological full group is Thompson's `V` (Matui 2015; also Nekrashevych 2004). Its orbit relation `R_C` is the ergodic hyperfinite III_{1/2} relation, with `L(R_C) ≅ R_{1/2}`."

---

## THM-4580 and HYP-9219: zeroless powers of two

**(1) Lift lemma and (2) bijection: CONFIRMED.**
* **Proof.** `D_j = 2^{n−k}(2^{jT_k}−1)/5^k` is even when `n≥k+1`, and `D_j ≡ 2^{n−k} j u (mod 5)`.
* **Brute force** for `k≤6` checks the LTE valuation, the 5 lifts, and the agreement of the parity bit with `q = x/2^k`. It also checks that the tails are exactly `{2^k | r, 5∤r}` and the count of zeroless classes (`a3_lift_doubling.out`).

**(3) Recursion and unit formula: CONFIRMED, and the formula is in fact PROVED for all `m`.**
* **Derivation.**
  * `Δ_{m−1} = −2^{1−m} Σ_{t odd} Π_{i≤m−2} Σ_a ζ^{ta10^i}`.
  * The factors with `i = m−3, m−2` equal 1, because `w^9 = w` for 4th and 8th roots of unity.
  * `ζ^{10^{m−1}} = −1` absorbs the sign.
* **Numerically** the formula reproduces `Δ_0..Δ_17` exactly (`a3_unit_formula.out`). I suggest replacing "Checked for m ≤ 12" with "PROVED (character sum); checked for m ≤ 18".
* **Cross-check.** The reader's DFS gives `Δ_15 = 2788823425 − 2788823405 = 20`, which equals `2Z_16 − 9Z_15`.

**(4) Data: CONFIRMED.**
* **Independent dense DP for `k ≤ 24`.** Its states are `r_i mod 2^{k−i}`, with no meet-in-the-middle. It equals the OEIS A181610 b-file (26 terms, Yamanouchi) and the repo values. The repo's `Z_25` and `Z_26` also equal OEIS.
* **My own MITM C code with a different split.** With `m = 16` low digits (2.5e10 leaves) and `L ≤ 24`, it reproduces `Z_17..Z_40` exactly, including `Z_40 = 119333906141890097435122400` (`a3_zk_m16.out`, 7 min). A third split, `m = 13`, matches to `Z_37`.
* **The constant.** `(2/9)^40 Z_40 = 0.88769404311514826449`, consistent with c.
* **Growth of Δ.** `|Δ_k|^{1/k}` is about 1.50–1.61 for `k = 31..39`, consistent with "≈ 1.6".

**(5) Brackets: CONFIRMED WITH CORRECTION (rounding direction).**
* `min B_25 = 18978171353195236` and `max B_25 = 24419129796700512` are reproduced by my code (`a3_bl_minmax.out`). L = 25 is also the best L.
* The rigorous intervals are:
  * growth in `[4.47847521, 4.52386053]`;
  * `dim_H` in `[0.93155668, 0.93782166]`.
* The stated `[4.47848, 4.52386]` and `[0.93156, 0.93782]` are rounded **inward** at all four ends, so as written they are not certified.
* **Replacement:** "growth rate (liminf and limsup of `Z_k^{1/k}`) in `[4.47847, 4.52387]`; `dim_H ∈ [0.93155, 0.93783]`". Make the same change in the THM-4580 title and in the HYP-9219 status.
* The bound `x^{0.93783}` is correctly rounded.

**(6) Doubling criterion: CONFIRMED.** A zero arises iff a 5 receives no carry, and the carry is 1 iff the right neighbour is ≥ 5. Brute force over all zeroless `x < 10^6` agrees.

**(7) Verification: CONFIRMED WITH CORRECTIONS.**
* **Direct and sampled checks.**
  * Every `n ∈ [87, 10^6]` was checked directly: 0 failures. Exactly 36 zeroless `n ≤ 999` exist, the largest 86. The records on `[957, 10^6]` equal the combined verifier's 10 records.
  * 3000 random `n` in `[957, 1.1e11)`, plus 3000 targeted `n` near `1.1e11` and near `109171987836`, were checked with `pow(2,n,10^300)`. Every one has a 0 among its last 251 digits.
  * `2^{10^10} mod 10^18 = 374549681787109376`, which matches R1's end state. R2 and R3 already have full end-state checks.
* **Records.**
  * For `103233492954` the first zero from the right is digit 250, so there are 249 zeroless trailing digits. This is consistent with THM-4580 and with the concurrent note `powers_two_decimal_windows_20261007.md`.
  * For `109171987836` it is digit 251, so 250 zeroless digits.
* **Correction (a): the OEIS entry range is off by one.** `109171987836` is A031142(**42**). A031143(41,42) = 249, 250.
  * **Replacement:** "entries 24–42".
* **Correction (b): "For n ≥ 957 the 0 lies among the last 251 digits" is false as stated.** A031142(43) = 181477218727 has 260 zeroless trailing digits; I confirmed this with `pow`. The first zero is digit 261. Later entries have 268, 275 and 308.
  * **Replacement:** "For 957 ≤ n < 1.1·10^11 the 0 lies among the last 251 digits."

**(8) Finite-digit obstruction: CONFIRMED.** Every class has at least 4 lifts, which gives nonemptiness and perfectness. The leading-block density is `log10(1+1/D)`.

**Heuristic: CONFIRMED** as HEURISTIC.
* The reader's script gives 33.93 for `n ≤ 86` (actual count 36), 2.3139 for the tail, Poisson 9.9%, clusters 14.6%, and `10^{−1.515·10^9}`.
* Optional: 33.9 comes from the refined finite-`k` lead×trail model. The displayed formula `τλ·0.9^{D(n)}` gives 34.16.

**KNOWN context: CONFIRMED.**
* A007377: "Checked up to k = 10^10. - _David Radcliffe_, Aug 21 2022".
* The A031142 b-file has 46 terms (Griffiths 2012), the last 7879942137257. A zeroless `2^n` with `n > 86` would itself be a record (leading-zero position `D(n)` exceeds every record position, all ≤ 308), so a complete table implies the conjecture for `n ≤ 7.88e12`.
* openai/math has 722 preprint folders; no title is digit-related.

**HYP-9219: OPEN typing is CONFIRMED, with 2 corrections and the bracket fix.**
* **Numbers check out.**
  * Density sup and inf `(2/9)^L max/min B_L` are 1.1413460468 and 0.8870365582, stable to about 3e−10 over `L = 18..25`.
  * `ρ(m=10) = 2.296`.
* **(i) "Equivalently" is too strong.**
  * A continuous density gives `(2/9)^k Z_k → f(0) = c`. Through a bounded density and the mass-distribution principle it also gives dimension `log_5(9/2)`.
  * It does not give the rate `O(0.4^k)`. Conversely, `Z_k`-asymptotics only concern the ball at 0.
  * **Replacement:** "Consequently `(2/9)^k Z_k → c` and the exponent set has dimension `log_5(9/2)`. Numerically, `Z_k = c(9/2)^k(1+O(0.4^k))`."
* **(ii) The sufficient condition is stated for too much.**
  * `max_t |P_m(ζ^t)| = O(ρ^m)` with `ρ < 9/2` gives `μ̂ ∈ ℓ¹`, hence a continuous density, and hence the dimension.
  * It does not give `O(0.4^k)`: the crude bound yields only `O((2ρ/9)^k)`. That rate would need `|Δ_k| = O(1.8^k)`.
  * State the scope accordingly.

---

## Required corrections

1. **THM-4568 (2):** replace "Every inequality in the paper's chain is an equality on `P_min`" with the wording given under (2) above. The cube step and `(3a−1)/2 ≥ a` are strict.
2. **THM-4568 (5):** "`ω_188 < 2.371339`, the bound that preceded #107" becomes "`ω_188 < 2.371339` (Alman et al., SODA 2025); `ω_190 < 2.371177` (Dupont et al. 2026, the latest pre-#107 bound cited by #107)". Apply the same fix to the source note, line 237.
3. **THM-4568 (4):** tag the Collatz bullet ANALOGY, and state the `−1/2` identity `3h+1+1/2 = 3(h+1/2)` rather than only `T(n)+1 = (3/2)(n+1)`, whose fixed point is `−1`.
4. **THM-4569 (4):** "(0.80 + o(1))" becomes "(√(2/π) + o(1)) ≈ 0.797". Optionally write `2√(2/π) ≈ 1.596` for the lag-1 pair.
5. **THM-4569 (7):** attribute `O_2` and `V` to the Deaconu–Renault groupoid, not to the relation `R_C`. Use the wording given under (7) above.
6. **THM-4580 (5), its title, and the HYP-9219 status:** use the outward-rounded brackets `[4.47847, 4.52387]` and `[0.93155, 0.93783]`.
7. **THM-4580 (7):** "entries 24–41" becomes "entries 24–42".
8. **THM-4580 (7):** "For n ≥ 957" becomes "For 957 ≤ n < 1.1·10^11". The counterexample beyond the range is `n = 181477218727`, whose first zero is digit 261.
9. **HYP-9219:** remove "Equivalently". State that the `ρ < 9/2` condition implies continuity and dimension, but not the `O(0.4^k)` rate.

## Optional improvements

* **THM-4568:** write the off-diagonal minimality proof (1), the Wendel proof of "everywhere" (2) and the tangent-line proof of "pinned" (3); type the diagonal-only κ formula as "derivation + LP".
* **THM-4569:** "(when `c_0 ∈ Z[1/3]`)" in (1); "as `j → +∞`" in (6).
* **THM-4580:** upgrade the unit formula to PROVED; say 33.9 is the refined model; cite this audit's `Z_17..Z_40` (`m = 16` split) and, for THM-4569, the B(9,100) reproduction.
