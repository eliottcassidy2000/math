---
id: HYP-9274
title: "(A) The Paley Hamiltonian-path ratio R(q) = H(P_q) 2^(q-1)/q! = e(1 - 3/(2q) + c_2/q^2 + ...) with c_1 = 3/2 (derived to first order; fits 1.5016) and c_2 near -3.08, universal over prime powers q = 3 mod 4; pre-registered R(43), R(47). (B) The wheel W_4 pattern sum is the Apery form: Sigma*(W_4) = -q(q-1)(a_q + 4), f = eta(2z)^4 eta(4z)^4. (C) Exact run polynomials A_4 = 2q(q-1)(q-3), A_6 = q(q-1)(q-3)(5q-29). (D) The W_4 channel puts its first arithmetic (Deligne-bounded) term into R(q) only through gluings of two runs (order about a_q/q^6); other non-series-parallel channels are open."
status: >
  OPEN, with these supports:
  FINITE-EXACT H(P_27) = 418425945688463432031 and H(P_31) = 19876729650849384218015865 (new; the latter confirms the
  pre-registered HYP-2371 interval [2.58949, 2.60249]: R(31) = 2.5955119).
  DERIVED to first order: c_1 = 3/2 = 5/2 (cherry placement) + 1 (two-cherry repulsion) - 2 (quartic runs).
  FINITE-EXACT ingredients: A_2 = q(q-1); A_3 = A_5 = 0; A_4 = 2q(q-1)(q-3) at 11 primes q <= 83;
  A_6 = q(q-1)(q-3)(5q-29) at 6 primes q <= 47, exact under leave-one-out interpolation; two-cherry ratio
  (q-2)(q-9)/((q-4)(q-5)) at 4 primes.
  FINITE-EXACT (B) at all 27 primes 5 <= q <= 109.
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

**Status.** FINITE-EXACT at all 27 primes `5 <= q <= 109`. The `a_q` are computed from the eta product in the script, not taken from tables.

**Context (CITED).**
- `f` is the modular form of Beukers' Apéry-number supercongruence for `ζ(3)`, proved by Ahlgren–Ono (2000) through Gaussian hypergeometric `₄F₃(1)`.
- Deligne's theorem (the Riemann hypothesis for varieties over finite fields) gives `|a_q| <= 2q^(3/2)`; the observed `|a_q|/(2q^(3/2))` is at most 0.88.
- Expected proof route: drop distinctness (a polynomial correction) and use scaling to reduce to `Σ_b χ(b)·(Σ_c χ(c(c−1)(c−b)))²`. That is a `χ`-twisted second moment of the Legendre family, which is the Ahlgren–Ono `₄F₃` object.

## (C) Exact run polynomials

- `A_2 = q(q−1)`, `A_4 = 2q(q−1)(q−3)` and `A_6 = q(q−1)(q−3)(5q−29)` (leading coefficients `1, 2, 5` = Catalan, THM-438). The odd `A_L` vanish.
- The two-cherry ratio is `(q−2)(q−9)/((q−4)(q−5))`.
- FINITE-EXACT (11, 6 and 4 primes respectively). Proofs would go by Möbius inversion over set partitions.
- **Conjecture (scoped):** `A_8` and `A_10` are also `q(q−1)(q−3)` times polynomials. From `k = 6` on, a single run can glue into Eulerian non-series-parallel patterns. The octahedron `K_(2,2,2)` (12 edges, all degrees 4) is the first candidate. Those patterns may carry their own Deligne-bounded terms, so polynomiality of every `A_(2k)` is not claimed.

## (D) Where arithmetic enters `R(q)`

**Heuristic.**
- **The `W_4` channel.** No single run glues into the `W_4` pattern. A trail and its quotients have at most two vertices of odd pattern-degree (doubling an edge changes neither the character pattern nor the parity), and `W_4` has four. Two runs can. So the Apéry-form contribution to `R(q)` comes from two-run configurations whose ten points glue to five, with size about `a_q/q^6`, at most `2q^(−4.5)` in relative terms.
- **Other channels.** Single runs of length `>= 12` can glue into other Eulerian non-series-parallel patterns (the octahedron first). Their order has not been estimated, and they could enter earlier or later than the `W_4` channel.

**Prediction (`W_4` channel only).** The `W_4` channel contributes no arithmetic term before relative order `q^(−4.5)`. Its first term is a multiple of `a_q(f)/q^6`. Whether other channels keep `R(q)` rational to that order is open.

**Status.** OPEN. This is the precise place where the Riemann hypothesis (Deligne's theorem) "permeates" the Hamiltonian-path count of Paley tournaments.
