---
id: HYP-9274
title: "(A) The Paley Hamiltonian-path ratio R(q) and its open uniform asymptotic tail, with formal c_1 = 3/2. (B) The wheel W_4 is the Apery form at odd primes. (C) Exact run reductions through A10, whose wheel coefficient is -34. (D) Exact weighted degrees eight and ten; degree ten retains the six-vertex K6-minus-P6 core."
status: >
  OPEN, with these supports:
  FINITE-EXACT H(P_27) = 418425945688463432031 and H(P_31) = 19876729650849384218015865 (new; the latter confirms the
  pre-registered HYP-2371 interval [2.58949, 2.60249]: R(31) = 2.5955119).
  DERIVED to first order: c_1 = 3/2 = 5/2 (cherry placement) + 1 (two-cherry repulsion) - 2 (quartic runs).
  FINITE-EXACT ingredients: A_2 = q(q-1); A_3 = A_5 = 0; A_4 = 2q(q-1)(q-3) at 11 primes q <= 83;
  A_6 = q(q-1)(q-3)(5q-29) at 6 primes q <= 47, exact under leave-one-out interpolation; two-cherry ratio
  (q-2)(q-9)/((q-4)(q-5)) at 4 primes.
  PROVED + CITED (B) for every odd prime: Ahlgren--Ono's character moment plus the explicit collision correction
  in wheel_apery_extension_20261014h.md. The 27-prime tests remain finite controls.
  PROVED exact A4/A6/A8 and two-cherry formulas by the terminating reduction in wheel_run_interface_20261014h.md;
  PROVED degree-eight weighted Hamiltonian coefficient with modular contribution -7a_q/((q-2)...(q-7)).
  PROVED terminating reduction through ten edges, with FINITE-EXACT symbolic extraction: A10=P10_single-34U_q;
  the weighted degree-ten coefficient retains 64V_q for V_q=I_(K6-minus-P6). Independent q11/q19 walk controls.
  NUMERICAL (A): with c_1 = 3/2 fixed, q^2(1 - R/e - 3/(2q)) = -3.117, -3.109, -3.090 (q = 27, prime power, never
  fitted), -3.097 for q = 19, 23, 27, 31; the holdout error at q = 31 is 1e-5.
session: opus-2026-10-14-S34
source: 05-knowledge/results/paley_h31_exact_rate_and_apery_wheel_20261014.md
scripts:
  - 04-computation/experiments/paley_H31_exact_20261014.py (+ .out)
  - 04-computation/experiments/paley_R_asymptotic_fit_20261014.py (+ .out)
  - 04-computation/experiments/paley_rate_c1_ingredients_20261014.py (+ .out)
  - 04-computation/experiments/paley_A6_polynomiality_20261014.py (+ .out)
  - 04-computation/experiments/paley_W4_charsum_20261014.py (+ .out), paley_W4_apery_form_20261014.py (+ .out)
  - 04-computation/experiments/paley_K4_charsum_cm_20261014.py, paley_K5_charsum_20261014.py
related: [THM-438, THM-4652, HYP-9273, HYP-2371 (CONFIRMED), OPEN-Q-013]
---

# HYP-9274

## (A) The rate of `R(q) → e`

**Statement.** For prime powers `q ≡ 3 (mod 4)`,

`1 − R(q)/e = (3/2)/q + c_2/q^2 + c_3/q^3 + O(q^(−4))`,

with `c_2 ≈ −3.08 ± 0.03`. The coefficients are the same for primes and prime powers (universality).

**Derivation of `c_1 = 3/2` (first order).** Write `R(q) = E_σ Π_k (1 + χ(d_k))` and expand over step subsets, split into maximal runs. At first order the surviving configurations are separated runs of even length (THM-4652 parts 2–3 kill isolated steps and odd totals).
- **Cherries.** A single cherry (run of length 2) has weight exactly `1/(q − 2)`. Placing `j` separated cherries among the `q − 1` steps gives `C(q − 2j, j)/(q − 2)^j = (1/j!)(1 − (5/2) j(j−1)/q + …)`.
- **Repulsion.** Two cherries on six distinct points correlate by the exact factor `(q−2)(q−9)/((q−4)(q−5)) = 1 − 2/q + O(q^(−2))`.
- **Quartic runs.** A quartic run has weight `A_4/(q)_5 = 2/((q−2)(q−4))` in about `q` positions, total `2/q`.

Summing over `j`: `R = e(1 − 5/(2q) − 1/q + 2/q) + O(q^(−2)) = e(1 − 3/(2q)) + O(q^(−2))`. A rigorous version needs uniform control of the remaining configurations, as in THM-438.

**Pre-registered (falsifiable) predictions:**
- `c_1` fixed at `3/2`: `R(43) ∈ [2.62798, 2.62803]`, `R(47) ∈ [2.63531, 2.63535]`;
- model-agnostic (orders 2–3, padded by the worst holdout error): `R(43) ∈ [2.62729, 2.62868]`, `R(47) ∈ [2.63468, 2.63594]`.

`H(P_43)` is out of reach of subset DPs (`2^41` subsets). The test awaits a different algorithm.

## (B) The wheel `W_4` is the Apéry form

**Statement.** Take the wheel `W_4` (hub `0`; rim `a, b, c, d`; edges hub→rim, `a→c`, `a→d`, `b→c`, `b→d`) and its injective pattern sum over `F_q`, for every prime `q >= 5` (both residue classes mod 4):

`Σ*(W_4) = q·Σ_(distinct nonzero a,b,c,d) χ(abcd(c−a)(d−a)(c−b)(d−b)) = −q(q−1)(a_q + 4)`.

Here `a_q` is the `q`-th coefficient of `f = η(2z)⁴η(4z)⁴`, the unique newform of weight 4 and level 8.

**Status, extended 2026-10-14.** **PROVED + CITED for every odd prime** by
[the exact wheel-to-moment interface](../results/wheel_apery_extension_20261014h.md).
The 27-prime eta-product computations remain finite controls. Over extension
fields the same collision correction is proved, but substituting the Fourier
coefficient `a_(p^r)` is false already at q=9; the proposed Frobenius-trace
extension is separately typed as finite-tested.

**Context (CITED).**
- `f` is the modular form of Beukers' Apéry-number supercongruence for `ζ(3)`, proved by Ahlgren–Ono (2000) through Gaussian hypergeometric `₄F₃(1)`.
- Deligne's theorem (the Riemann hypothesis for varieties over finite fields) gives `|a_q| <= 2q^(3/2)`; the observed `|a_q|/(2q^(3/2))` is at most 0.88.
- Completed proof route: allow collisions, adding `(q-1)(4-q)` to the anchored sum, and use scaling to reduce to `Σ_b χ(b)·(Σ_c χ(c(c−1)(c−b)))²`. Ahlgren--Ono Lemma2.2 and Theorem6 evaluate this moment as `-a_p-p` for every odd prime. The linked packet retains the signs, hub factor and injectivity correction.

## (C) Exact run reductions

- `A_2 = q(q−1)`, `A_4 = 2q(q−1)(q−3)` and `A_6 = q(q−1)(q−3)(5q−29)` (leading coefficients `1, 2, 5` = Catalan, THM-438). The odd `A_L` vanish.
- The two-cherry ratio is `(q−2)(q−9)/((q−4)(q−5))`.
- **PROVED extension:** the terminating injective graph reduction in
  [the run interface](../results/wheel_run_interface_20261014h.md) proves these
  formulas for all odd prime powers q=3mod4 and gives
  `A_8=2q(q-1)(q-3)(q-7)(7q-55)`.
- **PROVED extension, with exact symbolic extraction:**
  [the ten-edge core packet](../results/paley_next_cores_20261014i.md) gives
  `A10=2q(q-1)(q-3)(21q^3-576q^2+5300q-16198)-34U_q`.
  Degree-two elimination exposes the wheel although raw path-quotient parity
  would exclude it. This replaces the former open A10 evaluation target.
  The nonzero modular coefficient alone is not asserted to prove that no
  polynomial could agree on a specified infinite domain. Additional identities
  among higher arithmetic terminals remain **OPEN**.

## (D) Where arithmetic enters `R(q)`

**PROVED finite-degree interface.** Retain an edge weight t in
`R_q(t)=E_sigma product_i(1+t chi(d_i))`. For q=3mod4 prime, q>=11,
the exact coefficient `[t^8]R_q(t)` has modular part
`-7a_q/((q-2)(q-3)...(q-7))`, plus the explicit rational polynomial part in
[the run interface](../results/wheel_run_interface_20261014h.md).
The signed wheel coefficients of runs (6,2), (4,4), (4,2,2), (2,2,2,2)
are -14,-6,128,-960; their placement weights leave -14-3+64-40=7.
For the eight-edge single run the raw wheel is impossible by its four odd
vertices. This degree-eight calculation was independently checked by
Held--Karp at11. Deligne bounds this coefficient's modular part by
`14q^(3/2)/((q-2)...(q-7))=O(q^(-9/2))`.

**PROVED degree-ten extension, with exact symbolic extraction:** for q>=11,
`[t^10]R_q=(P10+(7q^2-173q+1166)U_q+64V_q)/(q)_10`, with the polynomial,
natural orientation of `V_q=I_(K6-minus-P6)`, and independent q11 product-DP
check in [the next-core packet](../results/paley_next_cores_20261014i.md).
A [guarded three-matrix evaluator](../results/pairwise_core_trace_20261014i.md)
independently computes V at7,11,19,31. The terminal classification through ten
edges is complete for these reduction rules; no modular evaluation of V is
claimed.

**OPEN:** the combined contribution of higher edge degrees at t=1, and the
uniform tail needed for the asymptotic coefficient c1. The finite-degree
formula is not an interchange-of-limits theorem for the full Hamiltonian count.
