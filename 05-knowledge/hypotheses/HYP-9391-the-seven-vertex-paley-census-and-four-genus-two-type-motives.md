---
id: HYP-9391
title: "The seven-vertex Paley census. (a) For primes q = 3 (mod 4), TT7(P_q) = q(q-1)[q^5 - 90q^4 + 3116q^3 - 49252q^2 + 322374q - 673743 - 2(14q^2 - 701q + 7254)a_q + 46 c_q - 83 d_q - 64(q - 31)|b(q)|^2 - 64 Y02 - 48 Y05 - 64 Y09 + 96 Y18]/2^21. Here a_q, c_q, d_q are the coefficients of 8.4.a.a, 8.6.a.a, 4.6.a.a, b = 32.3.c.a, and the Y are four explicit two-parameter Legendre triple products, so TT7 costs O(q^2). (b) The chi-twisted fourth moment of the Legendre family is sum_l chi(l) a_l^4 = -2q^2 - 3q - 1 - 3q a_q - c_q - d_q. (c) Modulo the known motives, the triple products carry exactly four new motives. They are irreducible, of weight 5, pairwise not quadratic twists, and not GL(2) (normalised traces exceed 2), with USp(4) Sato-Tate statistics: conjecturally the spinor motives of four genus-2 Siegel eigenforms."
status: >
  OPEN.
  (a) FINITE-EXACT: the formula reproduces all 71 directly counted values TT7(P_q), 7 <= q <= 811. It rests on three
  inputs:
  the exact reduced form of 2^21 TT7 (audit of THM-4780; PROVED structure);
  37 graph-sum identities S_H = (known combination) + (at most one triple product, coefficient +-1), FINITE-EXACT at
  43 primes 7..443 (23 primes 7..199 for the four sums without a degree-3 vertex);
  the wheel (HYP-9274 (B), PROVED for primes), the K6 - P6 terminal (HYP-9276, FINITE-EXACT) and (b).
  (b) FINITE-EXACT at the 76 odd primes 5..397, both classes mod 4. A monodromy count predicts the two weight-6 forms; the
  untwisted companion carries 4.6.a.a only.
  (c) NUMERICAL: 338 primes q = 3 (mod 4) up to 4999 (statistics over the 307 primes in [307, 4999]). Identification OPEN.
  Not yet independently audited.
session: mac-mini-2026-10-25-pyth (continuation, 2026-10-26)
source: 05-knowledge/results/pythagorean_hourglass_census_tower_20261025.md
depends_on: [THM-4780, HYP-9390, HYP-9274, HYP-9276]
scripts:
  - 04-computation/experiments/pythagorean_hourglass_20261025/tt7/motives/dump_genuine.py (-> genuine37.json): the 37 genuine graphs and multiplicities from the audit's engine
  - 04-computation/experiments/pythagorean_hourglass_20261025/tt7/motives/sum7.py, sum7_fast.py: exact O(q^5) and O(q^4) (degree-3 elimination) evaluators, validated against brute force
  - 04-computation/experiments/pythagorean_hourglass_20261025/tt7/motives/run_sums.py (-> sums_7_199.jsonl), run_fast.py (-> fast_*.jsonl): all survivor sums; sum m_H S_H = g(q) checked at every prime
  - 04-computation/experiments/pythagorean_hourglass_20261025/tt7/motives/classify_sums.py, exact_tt7_formula.py (+ .out, g_formula.json): the identities and the formula
  - 04-computation/experiments/pythagorean_hourglass_20261025/tt7/motives/cone_reduce.py, test_cones.py, closed_forms_check.py, run_cheap_Y.py (-> cheapY_*.txt): the triple products
  - 04-computation/experiments/pythagorean_hourglass_20261025/tt7/motives/legendre_moments.py (+ .out): part (b)
  - 04-computation/experiments/pythagorean_hourglass_20261025/tt7/motives/analyze_new.py, moments_new.py, moments_mod8.py, twist_tests.py, induced_test.py (+ .out): part (c)
  - 04-computation/experiments/pythagorean_hourglass_20261025/tt7/motives/identify_new.py, identify_products.py, identify_odd_levels.py, identify_Y.py, candidate_forms_odd.gp (+ .out): negative identifications
related:
  - 01-canon/theorems/THM-4780-paley-census-tower-to-six-vertices-wheel-octahedron-pgl2-hamiltonian-path-complement.md
  - 05-knowledge/hypotheses/HYP-9390-the-six-vertex-paley-census-is-the-level-32-adjoint.md
  - 05-knowledge/hypotheses/HYP-9276-k6-minus-p6-terminal-is-the-adjoint-of-a-weight-3-level-32-newform.md
---

# HYP-9391 — The seven-vertex Paley census and four genus-two-type motives

Notation:
- `q ≡ 3 (mod 4)` is prime and `φ` is the quadratic character.
- `L(λ) = Σ_w φ(w(w−1)(w−λ)) = −a_λ` is the Legendre-family trace.
- `a_q`, `c_q`, `d_q` are the coefficients of `8.4.a.a = η(2z)⁴η(4z)⁴`, `8.6.a.a` and `4.6.a.a = η(2z)^12`.
- `|b(q)|²` is the square modulus of the coefficient of `32.3.c.a` (HYP-9276).
- The sums below run over `s ≠ t` with `s, t ∉ {0, 1}`.

## (a) The census in closed form

The four triple products are:
```
Y02 = − Σ φ(s) φ(t) φ(t−1) L(s)² L(t/s)          Y05 = − Σ φ(s−1) φ(t) φ(t−1) L(s)² L(t/s)
Y09 =   Σ φ(s) φ(s−1) φ(t) L(s) L(t) L(t/s)       Y18 =   Σ φ(s−1) φ(t−1) L(s) L(t) L(t/s)
```
Then
```
2^21 TT7(P_q) / (q(q−1)) = q⁵ − 90q⁴ + 3116q³ − 49252q² + 322374q − 673743 − 2(14q² − 701q + 7254) a_q
                           + 46 c_q − 83 d_q − 64(q − 31)|b(q)|² − 64 Y02 − 48 Y05 − 64 Y09 + 96 Y18 .
```

**How it was found.**
- By THM-4780's calculus (exact reduced form from its audit), `2^21 TT7 = A7 + B7 T_W4 + C7 T_{K6−P6} + Σ_{37} m_H T_H`.
- The 37 genuine sums were computed exactly at 43 primes (`sum7.py`, `sum7_fast.py`, `run_*.py`). Each was then fitted in the span of `1, q, q², q³, a_q, q a_q, c_q, d_q, |b|², Y02, Y05, Y09, Y18`.
- Every one is exact with integer coefficients:
  - 21 lie in the span of the known motives alone;
  - the other 16 are known terms plus exactly one `±Y` (`exact_tt7_formula.out`).
- The cone reductions are what bring in the triple products:
  - summing out three pairwise non-adjacent degree-3 vertices with the `L`-table leaves a four-point base `{0, 1, s, t}`;
  - the "merge" corrections lie in the known span;
  - `Y04 = −Y02` and `Y06 = Y05` hold exactly.

**Check.** The formula gives all 71 counted values of `TT7(P_q)` for `7 ≤ q ≤ 811` (`paley_tt7_data*.out`). It costs `O(q²)`, against `O(q⁵)` or worse for counting.

## (b) The twisted fourth Legendre moment

```
Σ_{λ ≠ 0,1} φ(λ) a_λ⁴ = −2q² − 3q − 1 − 3q·a_q − c_q − d_q        (all 76 odd primes 5..397)
```

- The twisted second moment `Σ φ(λ) a_λ² = −q − 1 − a_q` is the wheel (Ahlgren–Ono).
- Prediction of the weight-5 part:
  - the local system `Sym⁴(H¹) ⊗ L_φ` on `P¹ − {0, 1, ∞}` has Euler characteristic `−5`;
  - it has one invariant (at `λ = 1`), so its parabolic `H¹` has dimension 4 and weight 5;
  - the identity finds that part to be `8.6.a.a ⊕ 4.6.a.a`. (For `Sym²` the same count gives `8.4.a.a`.)
- The identity sits in the Ahlgren / Ahlgren–Ono / Frechette–Ono–Papanikolas framework (Legendre moments and Hecke traces on `Γ₀(4)`, `Γ₀(8)`), and may be a known case; not checked against the papers.
- The untwisted companion is `Σ_{λ≠0,1} a_λ⁴ = 2q³ − 4q² − 9q − 3 − d_q` (76 primes; `untwisted_moment.out`). Only `4.6.a.a` appears, the classical `Γ₀(4)` moment, so `8.6.a.a` enters through the `φ(λ)` twist.
- In the census, the `c_q`, `d_q` terms are these fourth moments. For example `K_{3,4}` (four degree-3 vertices on one untwisted triple) has `S = −171 + 50q + 4q² − q³ − d_q`.

## (c) Four new motives

Here `N` is the new part: the residual of a `q^{−5}`-weighted least-squares fit of each `Y` on the nine known columns, normalised by `q^{5/2}`.
- **Number and irreducibility.**
  - Modulo the known span, `Y02, Y05, Y09, Y18` span exactly four new directions (exact rank).
  - The Gram matrix of their new parts is the identity, to within 0.07, in the basis `m1 = N(Y02)`, `m2 = N(Y09)`, `m3 = N(Y18)`, `m4 = N(Y05) + N(Y18)`.
  - `E[m_i²] = 0.96, 0.93, 0.96, 0.96` over the 307 primes in `[307, 4999]`. So there are four distinct irreducible motives of weight 5, each entering one `Y` (or `Y05 + Y18`) once.
- **Not GL(2).**
  - At the 184 primes `≥ 2000`, `max |m_i| = 2.59, 3.16, 3.10, 3.01`.
  - The share of those primes with `|m_i| > 2` is 3–6.5%; USp(4) predicts 5.0%.
  - A GL(2) motive satisfies `|x| ≤ 2`. Odd weight forces a symplectic motive, hence even rank, so the rank is at least 4.
  - The low-order terms that the fit determines poorly contribute less than about 0.05 there.
- **USp(4) statistics.** Compared with Monte Carlo values:

  | | USp(4) | SU(2) | `m1` | `m2` | `m3` | `m4` |
  |---|---|---|---|---|---|---|
  | `E[x⁴]` | 3.0 | 2.0 | 3.21 | 2.34 | 2.60 | 2.70 |
  | `E[x⁶]` | 14.0 | 5.0 | 17.2 | 9.9 | 10.9 | 11.3 |
  | `P(\|x\| < 0.2)` | 0.169 | 0.127 | 0.178 | 0.153 | 0.143 | 0.153 |

  The standard error of `E[x⁴]` is about 0.5 at 307 primes. A tensor product `SU(2)×SU(2)` would give `E[x⁴] = 4`.
- **No coincidences.**
  - No quadratic twist relates two of them, and none is a self-twist: for characters `(2/p), (±3/p), (5/p), (7/p)`, every `|E[χ m_i m_j]|` (including `i = j`) is at most 0.18, about two standard errors, against 1 for a twist.
  - None is induced: there is no vanishing residue class, and only 3–5% of values have `|x| < 0.05`, against about 50% for an induced motive.
- **Negative identifications** (none removes a new direction; tested at 31–71 primes; `identify_*.py`):
  - 97 newform orbits of 2-power level up to 256 (weights 2–6; characters `1, χ8, χ−4, χ−8`);
  - 260 orbits of levels `2^a·{3, 9, 5, 7, 15, 21, 45, 63}` (weights 2, 4, 6);
  - 329 products and cubes of these;
  - `|b|²` and `q²` times weight-2 orbits;
  - twisted polynomials.
- **Conjecture.** The four motives are the spinor motives of four genus-2 Siegel modular eigenforms with 2-power level (paramodular or `Γ₀`-type). Motivic weight 5 means scalar weight 4, or vector-valued weight `(3, 2)`.
- **Mechanism.**
  - Multiplicative convolution of hypergeometric sheaves concatenates parameters (Katz). So `Σ_s φ(…)L(s)L(t/s)` is, up to twists, the rank-4 family of `₄F₃(½,½,½,½; 1,1,1)`, whose generic fibres are rank-4 symplectic.
  - Each `Y` pairs that family with one more Legendre sheaf. The new motives are pieces of `H¹_!` of the resulting rank-8 local system.
- **Caveat.** The split into "known part + new motives" is statistical. The coefficients of `1, q, q², a_q` in the projection are not determined exactly, so the exact integer traces of the four motives are not yet available. The rank and Sato–Tate evidence uses only large primes.

## Next steps

1. **Exact separation.** Compute the boundary and degenerate contributions of each `Y` exactly, giving exact traces `t_k(q)`.
2. **Euler factors.** Evaluate the `Y` over `F_{q²}` (and `F_{q³}`) for small `q` to get degree-4 Euler factors.
3. **Conductor and weight.** Test the functional equation with PARI `lfun` (degree 4, motivic weight 5) over 2-power conductors.
4. **Comparison.** Check against the LMFDB / Poor–Yuen paramodular data.
5. **Seven vertices at prime powers.** Extend the formula of (a) to prime powers, as THM-4780 does at six vertices.
