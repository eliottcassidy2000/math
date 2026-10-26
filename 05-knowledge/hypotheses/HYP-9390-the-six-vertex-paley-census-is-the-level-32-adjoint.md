---
id: HYP-9390
title: "The six-vertex Paley census is the level-32 adjoint. (a) For q = 3 (mod 4), TT6(P_q) = q(q-1)(q^4 - 49q^3 + 828q^2 - 5313q + 10193 - 7(q-15)a_q - 8|b(q)|^2)/2^15. Here a_q is the coefficient (Frobenius trace) of the Apery form eta(2z)^4 eta(4z)^4, and b = LMFDB 32.3.c.a is the non-CM weight-3 chi_-4 newform of HYP-9276, with |b(q)|^2 = q^2 minus the Sym^2 Frobenius trace (at primes, |b(q)|^2 = q^2 - b(q^2)). Given THM-4780 and HYP-9274 (B), at primes p = 3 (mod 4) this is equivalent to HYP-9276 (1) at those primes. (b) The wheel identity S_W4(q) = -4 - a_q holds at prime powers. (c) Seven vertices: the census meets 37 further dense sums (multiplicities listed). TT7 at 72 values of q lies outside the span of the six-vertex motives and of 97 GL(2) newform orbits of 2-power level, so seven vertices carry new arithmetic, isolated exactly as g(q); reduced in HYP-9391 to four explicit Legendre triple products, whose new parts behave numerically like four weight-5 motives of rank >= 4 (conjecturally genus-2 Siegel)."
status: >
  OPEN, CONDITIONAL.
  (a) At primes p = 3 (mod 4) it is PROVED conditional on HYP-9276, by THM-4780 (5), HYP-9274 (B) and THM-4653.
  FINITE-EXACT at 75 values of q: 3, 7, 11, the 67 primes 19..751, and the prime powers 27, 243, 343, 1331, 2187.
  The formula was found before the reduction by an exact span search: 14 columns, unique solution, 53 surplus
  equations. In that search the weight-6 level-8 newform, (2/q) twists, a_q^2 and CM forms of Q(sqrt-2) and
  Q(sqrt-7) all enter with coefficient 0 or fail.
  (b) FINITE-EXACT at q = 27, 243, 343, 1331, 2187 through the census (and earlier, by direct evaluation, at 9, 25, 27, 49, 81,
  121, 125, 169 in wheel_apery_extension_20261014h.md section 4); PROVED for primes (HYP-9274 (B)).
  Audited 2026-10-25 (pythagorean_hourglass_20261025/audit/REPORT.md): all FINITE-EXACT values reproduced with independent code
  and PARI's mf eigenforms; the Sym^2 reading of |b(q)|^2 corrected (sign) and Remark 1's coefficient ring corrected (MISTAKE-670).
  (c) The survivor list and multiplicities are PROVED (exhaustive). TT7 at 72 values of q is NOT in the span of the
  six-vertex motives or of 97 GL(2) newform orbits of 2-power level (weights 2-6), singly or by weight; the genuine
  seven-vertex part g(q) is isolated exactly (audit's reduced form). Its arithmetic is OPEN.
session: mac-mini-2026-10-25-pyth
source: 05-knowledge/results/pythagorean_hourglass_census_tower_20261025.md
depends_on: [THM-4780, THM-4653, HYP-9274, HYP-9276]
scripts:
  - 04-computation/experiments/pythagorean_hourglass_20261025/paley_tt6.py (-> paley_tt6_data.out): TT6(P_p), p = 19..751
  - 04-computation/experiments/pythagorean_hourglass_20261025/tt6_adjoint_search.py (+ .out): exact span search, 32 candidate columns
  - 04-computation/experiments/pythagorean_hourglass_20261025/tt6_search.py (+ .out): CM candidates of weights 3 and 5 (no identity)
  - 04-computation/experiments/pythagorean_hourglass_20261025/paley_prime_powers.py (+ .out): TT5, TT6 at prime powers with Frobenius traces
  - 04-computation/experiments/pythagorean_hourglass_20261025/graph_sums6.py (+ .out): S_(K6-P6) = |b|^2 - q^2 - 7 at 23 primes
  - 04-computation/experiments/pythagorean_hourglass_20261025/census_level_forecast.py (+ .out): part (c)
related:
  - 01-canon/theorems/THM-4780-paley-census-tower-to-six-vertices-wheel-octahedron-pgl2-hamiltonian-path-complement.md
  - 05-knowledge/hypotheses/HYP-9276-k6-minus-p6-terminal-is-the-adjoint-of-a-weight-3-level-32-newform.md
  - 05-knowledge/hypotheses/HYP-9275-paley-runs-one-modular-form-octahedron-k6p6.md
  - 05-knowledge/hypotheses/HYP-9371-four-loop-law-doubly-regular-tournaments-see-the-five-vertex-census.md
---

# HYP-9390 — The six-vertex Paley census is the level-32 adjoint

## (a) The formula

By THM-4780 (5), for every `q ≡ 3 (mod 4)`:
`TT6(P_q) = q(q−1)[(q−3)(q³ − 46q² + 682q − 3239) + 7(q − 15)S_W4(q) − 8S_{K6−P6}(q)]/2^15`.

Substitute two identities:
- the wheel identity `S_W4 = −4 − a_q` (HYP-9274 (B) at primes);
- the terminal identity `S_{K6−P6} = |b(q)|² − q² − 7`. This is THM-4653's `V_q/(q(q−1)) = M_q + 2q − 6` together with HYP-9276's `M_q = |b(q)|² − (q+1)²` at Paley primes.

The result is
```
TT6(P_q) = q(q−1) [ q⁴ − 49q³ + 828q² − 5313q + 10193 − 7(q − 15)·a_q − 8|b(q)|² ] / 2^15 .
```

- The adjoint enters only through `K6 − P6`. Its signed multiplicity `−8` is the Désarménien–Foata value `½·2³A_3(−1)`.
- The Apéry form enters with weight `7(q − 15)`. This is the octahedron's `q − 5` plus `6q − 100` from six-vertex graphs that reduce to the wheel.
- For a prime power `q = p^r` (`r` odd), read `a_q = α^r + β^r` with `α + β = a_p` and `αβ = p³`. Read `|b(q)|² = h_r²` with `h_r = u^r + v^r`, `u + v = h(p)`, `uv = p²` and `b(p) = i·h(p)`. These are the Frobenius traces: `|b(q)|² = q² − S_q`, where `S_q = α_b^{2r} + (α_bβ_b)^r + β_b^{2r}` is HYP-9276's `Sym²` trace (`S_p = b(p²)`, so `|b(p)|² − p² = −b(p²)` at primes). For example `S_27 = −1207` and `|b(27)|² = 1936 = 44²`.
- **Evidence (FINITE-EXACT).**
  - Direct census at the 67 primes `19 ≤ p ≤ 751` with `p ≡ 3 (mod 4)` (`paley_tt6_data.out`).
  - `q = 3, 7, 11`, where `TT6 = 0`.
  - The prime powers `27, 243, 343, 1331, 2187`, where `TT6(P_27) = 0` (ST27).
  - Each value is an independent check of HYP-9276 at that `q`.
- **Search history.** The formula was first found by an exact span search over 67 primes, before the reduction was available (`tt6_adjoint_search.out`).
  - `A(p) + B(p)a_p` alone is not in the span.
  - With `|b|²` added (`deg B ≤ 3`, `deg C ≤ 2`) the solution is unique, with 53 surplus equations.
  - A 32-column fit gives the same solution. It includes `(2/p)`-twisted polynomials, `a_p²`, the weight-6 level-8 newform `8.6.a.a` (`q + 20q³ − 74q⁵ − 24q⁷ + …`) and twisted `|b|²`, and every other coefficient is `0`.
  - CM forms of `Q(√−2)` and `Q(√−7)` of weights 3 and 5 give no identity (`tt6_search.out`). The `Q(√−3)` columns there take the trace from the first representation `4p = a² + 3b²` and so ignore the unit ambiguity: the weight-3 column gives `−13`, `−37` at `p = 7, 19`, against `2`, `26` for the newform `η(2z)³η(6z)³`, and the weight-5 column has the same defect. So that case is untested.

## (b) Prime powers

The wheel identity is PROVED at odd primes (Ahlgren–Ono via HYP-9274 (B)). At `q = 27, 243, 343, 1331, 2187` it holds in Frobenius-trace form: the census values satisfy THM-4780 (6) with `a_q = α^r + β^r` (`paley_prime_powers.out`). The same Frobenius form was already FINITE-EXACT, by direct evaluation of the wheel, at 9, 25, 27, 49, 81, 121, 125, 169 (`wheel_apery_extension_20261014h.md` §4); the census adds 243, 343, 1331, 2187. Here `a_q` is the Frobenius trace, not the Fourier coefficient `a(p^r)` (they differ already at `q = 9`, and `a(27) = 152 ≠ 260`). The natural proof is the same Calabi–Yau threefold argument over `F_q`.

## (c) Seven vertices

THM-4780's lemmas give
`2^21 TT7(P_q) = A7(q) + B7(q)T_W4 + C7(q)T_{K6−P6} + Σ_H m_H T_H`,
with `deg A7 ≤ 7`, `deg B7 ≤ 2` and `deg C7 ≤ 1`.

- The sum runs over the 38 dense seven-vertex graphs that have an even number of edges and no odd-reversal automorphism.
- One of them, the 4-regular graph `K7 − C7` with 14 automorphisms, reduces by THM-4780 (4) to `±(q − 6)·T_{K6−P6}`.
- The other 37 are genuine. Their multiplicities satisfy `|m_H| ∈ {1, 3, 4, 6, 8, 12, 16, 24, 32, 40, 48}`. The signs depend on the labelling of each `H`; in geng's labellings (`census_level_forecast.out`) the signed values are `{+1, +3, ±4, +6, ±8, ±12, ±16, ±24, −32, −40, −48}`.

**Question.** Which motives do the 37 genuine five-variable sums carry, and do they collapse in `TT7`?

**Cheapest decisive test.**
- Compute `TT7(P_q) = q·TT6(T[QR])` for `q ≤ 250`.
- Test exact membership in the span of `q^i`, `q^i a_q` (`i ≤ 2`) and `q^i|b(q)|²` (`i ≤ 1`).
- If the data lie in that span, so does the genuine combination `Σ m_H T_H` (it is not zero: `2^21·TT7(P_7) = 0` while the known part is `−25536`). If not, seven vertices carry new arithmetic.
- A natural first candidate is the weight-3 level-32 Siegel-type lift of HYP-9276's notes.

**Result of the test (same session, after the audit; scripts in `pythagorean_hourglass_20261025/tt7/`).**
- **Data.** `TT7(P_q)` at 72 values: `q = 3, 7, 11` and the 69 primes `19 ≤ p ≤ 811` with `p ≡ 3 (mod 4)`.
  - Counted by a chain counter (`paley_tt7.py`), checked against brute-force enumeration at `q ≤ 43`.
  - `TT7 = 0` for `q ≤ 23`, `TT7(P_31) = 465`, `TT7(P_43) = 7224`. The audit's independent values agree at 31 and 43.
- **Not in the base span.** `2^21 TT7` is **not** in the 13-column span of `q^i` (`i ≤ 7`), `q(q−1)q^i a_q` (`i ≤ 2`) and `q(q−1)q^i |b(q)|²` (`i ≤ 1`). So the genuine combination does not reduce to the six-vertex motives.
- **Adding columns does not help.** Tried:
  - any single Galois orbit among 97 newform orbits of 2-power level `8..256`: weights 2–6, characters `1`, `χ_8`, `χ_{−4}`, `χ_{−8}`, each with its admissible multipliers and, for weights `≤ 3`, also `|a_p|²`;
  - all orbits of weight 3, 5 or 6 jointly;
  - `8.6.a.a` together with `4.6.a.a`;
  - twisted polynomials `(D/q)q^j`, `D ∈ {2, ±3, 5, ±7, ±6}`.
  
  Receipts: `tt7_span_test.out`, `tt7_candidate_fit.out`, `tt7_twist_test.out`, `tt7_level48_test.out`.
- **The genuine part, isolated.** The audit's symbolic engine gives the reduced form exactly:
  - `A7 = q(q−1)(q−3)(q⁴ − 87q³ + 2785q² − 38827q + 198842)`;
  - `B7 = 28q² − 1288q + 13238`;
  - `C7 = 1912 − 64q`.
  
  So `g(q) = (2^21 TT7 − A7 − B7T_W4 − C7T_{K6−P6})/(q(q−1)) = Σ m_H S_H` is known at all 72 values (`tt7_genuine_part.out`):
  - `g(7) = 608`, `g(11) = −23616`, `g(19) = 100224`;
  - always divisible by 32;
  - `|g|/q^{5/2}` up to about 640.

  It is not in the span of `q^i` (`i ≤ 3`), `a_q`, `q·a_q`, `|b|²`, plus any single orbit above, or all orbits of weight 3, 5, 6 or `{5, 6}` jointly.
- **Verdict.** The seven-vertex census carries arithmetic outside every GL(2) candidate tried.
- **Reduced in [HYP-9391](HYP-9391-the-seven-vertex-paley-census-and-four-new-weight-five-motives.md)** (same session, 2026-10-26; OPEN, part (c) NUMERICAL):
  - every genuine sum is known terms plus at most one of four explicit two-parameter Legendre triple products;
  - TT7 gets an exact `O(q²)` formula (all 71 counted values);
  - numerically, the new arithmetic behaves like four irreducible weight-5 motives of rank at least 4, with statistics consistent with USp(4) but not fixing the rank (conjecturally genus-2 Siegel), plus `8.6.a.a` and `4.6.a.a` through the twisted fourth Legendre moment (a known identity: the `k = 6` case of Frechette–Ono–Papanikolas 2004, Prop. 2.1).

## Remarks

1. **Whole censuses.** Every 5-vertex class count of `P_q` lies in `Q[q] + Q[q]·S_W4`, and every 6-vertex class count in `Q[q] + Q[q]·S_W4 + Q·q(q−1)·S_{K6−P6}`, by the same reduction. Denominators divide `2^{C(k,2)}·|Aut τ|`, so `Z[1/2]` is right only for labelled counts: already the 3-cycle count is `q(q²−1)/24`, and the regular 5-tournament count has polynomial part with denominator `5·2^10`.
   - So for `P_q` the five-vertex census entering the four-loop term `X_3` of HYP-9371 is affine in `a_q`.
   - The six-vertex census entering its five-loop part is affine in `(a_q, |b(q)|²)`.
2. **A cross-link to check.** McCarthy–Springfield (arXiv:2311.02135, §2.1) meet "a certain non-CM modular form of weight three and level 32" in the order-4 count of the fourth-power Paley digraph (`q ≡ 5 mod 8`).
   - If it is the same `32.3.c.a`, one 2-dimensional Galois representation appears linearly for `k = 4`, four vertices, and through its adjoint for `k = 2`, six vertices.
   - Not checked here.
3. **Falsifier.** Any `q ≡ 3 (mod 4)` at which the formula of (a) fails. THM-4780 makes this the same event as a failure of HYP-9276 (1) at that prime (or, at prime powers, of (b) or of HYP-9276 (2)).
