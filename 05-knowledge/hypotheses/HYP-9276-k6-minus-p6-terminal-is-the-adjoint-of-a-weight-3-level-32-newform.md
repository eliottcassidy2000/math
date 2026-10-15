---
id: HYP-9276
title: "The Paley K6 - P6 terminal is the adjoint (symmetric square) of the McCarthy-Papanikolas weight-3 newform of level 32: M_q = -b(q^2) - 2q - 1, V_q = q(q-1)(|b(q)|^2 - q^2 + 7 + 14 chi_(-4)(q))"
status: >
  OPEN (a modular identity). FINITE-EXACT at all 3,243 primes 5 <= q <= 30,000 (both residue classes mod 4), and at
  17 prime powers (9, 25, 27, 49, 81, 121, 125, 169, 243, 289, 343, 361, 625, 729, 1331, 2187, 2401) in the Frobenius form of
  statement 2, including the Paley prime powers 27, 243, 343, 1331 and 2187.
  The reduction it feeds (THM-4653) is PROVED. Corollaries (b)-(c) are CONDITIONAL on this hypothesis and use CITED theorems.
  (d)-(e) are FINITE-EXACT observations.
session: mac-mini-2026-10-14-wilson
source: 05-knowledge/results/paley_k6p6_adjoint_weight3_level32_20261014.md
scripts:
  - 04-computation/experiments/paley_k6p6_20261014/k6p6_fft.py (+ .out, .json): M_q at all primes to 30,000
  - 04-computation/experiments/paley_k6p6_20261014/identify_g.gp, g_coeffs.gp (+ .out), g_trace_a2_half.txt: the newform (gp mf package)
  - 04-computation/experiments/paley_k6p6_20261014/k6p6_identity_check.py (+ .out): the identity, Sato-Tate moments, zeros, 2-adic data
  - 04-computation/experiments/paley_k6p6_20261014/k6p6_prime_powers.py (+ .out): prime powers to 361, exact field arithmetic
  - 04-computation/experiments/paley_k6p6_20261014/k6p6_prime_powers_large.py (+ .out): 625, 729, 1331, 2187, 2401 (O(q^2), log tables)
related: [THM-4653, HYP-9274, HYP-9275, THM-4650 (D19 midline curve of conductor 32), paley_next_cores_20261014i.md]
---

# HYP-9276: the K6 − P6 terminal is a symmetric square

## The newform

Let `g = Σ b(n) qⁿ` be the newform of weight 3 on `Γ₀(32)` with nebentypus `χ₋₄ = (−4/·)`:

`g = q + 4i q³ + 2q⁵ − 8i q⁷ − 7q⁹ − 4i q¹¹ − 14q¹³ + 8i q¹⁵ + 18q¹⁷ − 12i q¹⁹ + …`

- `S₃^new(Γ₀(32), χ₋₄)` has dimension 2. It is spanned by `g` and its complex conjugate `ḡ = g ⊗ χ₋₄`.
- `g` has no CM. Its coefficient field is `Q(i)`.
- `b(p)` is real for `p ≡ 1 (mod 4)` and purely imaginary for `p ≡ 3 (mod 4)`.
- This is the form `f₂` of McCarthy–Papanikolas (CITED, arXiv:1205.1006, IJNT 2015, Lemma 5.1); we recomputed all of the above with gp's `mf` package.

At an odd prime `p`, with Satake roots `α + β = b(p)` and `αβ = χ₋₄(p)p²`,

`b(p²) = α² + αβ + β² = χ₋₄(p)(|b(p)|² − p²)`.

This is the trace of Frobenius on the adjoint (symmetric-square) Galois representation. That representation is a cuspidal `GL(3)` object, because `g` has no CM (Gelbart–Jacquet).

## Statement

Let `M_q` be the Legendre bilinear form of THM-4653, and `ε = χ(−1)` in `F_q`.

1. **(Primes.)** For every prime `q >= 5`,

   `M_q = −b(q²) − 2q − 1`.

   Equivalently, `Σ_x φ(x)φ(1−x) ₂F₁(x) ₃F₂(x) = (b(q²) + 2q + 1)/q³` in Greene's normalisation. By THM-4653,

   `V_q = I_(K6−P6)(q) = q(q−1)(|b(q)|² − q² + 7 + 14χ₋₄(q))`.

   For the Paley primes `q ≡ 3 (mod 4)` this reads `V_q = q(q−1)(|b(q)|² − q² − 7) = −q(q−1)(b(q²) + 7)`. Note that `|b(q)|² = −b(q)²`.

2. **(Prime powers.)** For `q = p^r` with `p` odd,

   `M_q = −(α^(2r) + (αβ)^r + β^(2r)) − 2q − 1`,

   which is the trace of `Frob_p^r` on the adjoint representation. Then `V_q = q(q−1)[ε(S_q + 14) + 7]` with `S_q = α^(2r) + (αβ)^r + β^(2r)`. This covers the Paley prime powers `q = 27, 243, 343, 1331, 2187`.

### Values

| q | b(q) | W_q = V_q/(q(q−1)) |
|---|---|---|
| 7 | −8i | 8 |
| 11 | −4i | −112 |
| 19 | −12i | −224 |
| 23 | 40i | 1064 |
| 27 = 3³ | (Frobenius form, `S₂₇ = −1207`) | 1200 |
| 31 | −32i | 56 |
| 43 | 28i | −1072 |
| 47 | 16i | −1960 |

The recorded values `V₇ = 336`, `V₁₁ = −12320`, `V₁₉ = −76608` and `V₃₁ = 52080` of the source note are reproduced.

## Evidence (FINITE-EXACT)

- **Primes.** All 3,243 primes `5 <= q <= 30,000`, with no failure.
  - `M_q` comes from the O(q log q) FFT evaluator of THM-4653 (3). Its float64 rounding is checked against the exact O(q³) values at all 237 primes below 1500, and every entry is an integer of size at most `q²`.
  - `b(q)` comes from gp's `mfcoefs`.
- **Prime powers.** All 12 prime powers up to 361, with exact `F_q` arithmetic (an irreducible polynomial, `χ` from the squares), and 625, 729, 1331, 2187, 2401 by an O(q²) evaluation with discrete-log tables.

The identity was found by fitting `M_q` against `Q(q)`-combinations of trace data of newform orbits. Products and squares were tried at weights 2–4, then the weight-3 square was isolated. It was then tested at every prime and prime power above as a single exact equation. Nothing was fitted at the tested primes.

## Corollaries

**(a) Size (CONDITIONAL; uses Deligne).** Deligne gives `|b(q)| <= 2q`. Hence `−q² − 7 <= W_q <= 3q² − 7` at Paley primes, and `|V_q| <= 3q³(q−1)`.

**(b) Distribution (CONDITIONAL; uses BLGHT 2011).** Let `x_q = (|b(q)|² − q²)/q² = 4cos²θ_q − 1`. Sato–Tate for non-CM newforms with nebentypus (Barnet-Lamb–Geraghty–Harris–Taylor) gives the `SO(3)` law for `x_q` in each residue class mod 4. Its moments are the Riordan numbers `1, 0, 1, 1, 3, 6, 15, 36, 91`.

| Sample | Moments k = 2..8 |
|---|---|
| 1,632 Paley primes below 30,000 | 1.005, 1.000, 3.002, 5.987, 14.974, 35.937, 90.946 |
| 1,611 primes `q ≡ 1 (mod 4)` | 0.966, 0.902, 2.751, 5.349, 13.304, 31.548, 79.220 |

The second row is within about one standard error. The standard errors are `√((m_(2k) − m_k²)/N)`, for example 0.035 for k = 2 and 0.23 for k = 4. The rows are correlated.

**(c) Independence from the wheel (CONDITIONAL; sketch).** There are no `P, B, C ∈ Z[q]`, not all zero, with `P(q) + B(q)U_q + C(q)V_q = 0` for all large Paley primes. This answers, conditionally, the question left open in `paley_next_cores_20261014i.md` §1.

*Sketch.* Write `U_q = −q(q−1)(a_q + 4)` with `y_q = a_q/q^(3/2)`, the Apéry form being non-CM of weight 4 and level 8.
- A relation with `C ≠ 0` forces `x_q = λ y_q + μ + o(1)` along Paley primes, by boundedness of `x_q` and of `y_q`.
- Sato–Tate in residue classes for `f` and for `Ad g` gives:
  - `E[x] = 0` and `E[x²] = 1`;
  - `E[y²] = 1`.
- Rankin–Selberg for `GL(3) × GL(2)` gives `E[xy] = 0`. This uses `L(s, Ad g × f)` and its twist by `χ₋₄`: both are entire and non-vanishing on `Re s = 1` (Jacquet–Piatetski-Shapiro–Shalika, Shahidi).
- The relation then gives `μ = 0` and `λ² = 1`, but also `λ = E[xy] = 0`. Contradiction.
- `C = 0` forces `y_q` to converge, against `E[y²] = 1` and `E[y] = 0`.

**(d) Zeros (FINITE-EXACT).** At Paley primes, `M_q = (|b(q)| − q − 1)(|b(q)| + q + 1)`. So `M_q = 0` iff `|b(q)| = q + 1`, and then `V_q = 2q(q−1)(q−3)`.
- Below 30,000 this happens exactly at `q = 7, 31, 199, 223, 1151, 2239, 2663, 3391, 3583, 23719, 23743`, all `≡ 7 (mod 8)`.
- `V₇ = 336` and `V₃₁ = 52080` are the first two.
- For `q ≡ 1 (mod 4)`, `M_q = (q−1)² − b(q)² − 2` never vanishes, since a difference of two squares is never 2.

**(e) 2-adic (FINITE-EXACT, sharp minima).** Below 30,000:
- `2⁷ | M_q` for `q ≡ 3 (mod 8)`;
- `2¹⁰ | M_q` for `q ≡ 7 (mod 8)`;
- `v₂(M_q) = 1` occurs for `q ≡ 1 (mod 4)`.

In terms of the form: `|b(q)|² ≡ (q+1)²` modulo `2⁷` and `2¹⁰` respectively, which is a congruence for `g`. Unexplained.

## Literature and the level-32 dictionary (typed)

- **McCarthy–Papanikolas (CITED).**
  - Their Theorem 1.8 gives `b(p) = ₃F₂(χ₄, φ, φ; ε, ε | 1)_p` for `p ≡ 1 (mod 4)` (McCarthy's normalisation, quartic `χ₄`). The proof uses Hijikata's trace formula and Schoof's count of Legendre curves with `Z/2 × Z/8` torsion.
  - They give no hypergeometric expression for `b(p)` when `p ≡ 3 (mod 4)`. That is the Paley case.
  - Statement 1 expresses `|b(q)|²` in both classes through quadratic characters only.
  - We found no prior statement of statement 1 in their paper or in the Dawsey–McCarthy survey (arXiv:2101.06303, Birkhäuser 2021; abstract and table of contents only). Searches for Paley `K₅`/`K₆` counts found only `K₃`/`K₄` (Evans–Pulham–Sheehan; Dawsey–McCarthy).
- **Same level, different weight (STRUCTURAL).**
  - The other factor of the McCarthy–Papanikolas level-8 Siegel eigenform is `f₁ = η(4z)²η(8z)²`. That is the congruent-number curve `y² = x³ − x` of conductor 32.
  - Their Theorem 1.6, conditional on the van Geemen–van Straten conjecture, reads `₄F₃(φ, φ, φ, φ; ε, ε, ε | −1)_p = ±b(p)a(p)` for `p ≡ 1 (mod 4)`, and the value vanishes for `p ≡ 3 (mod 4)`.
  - The Guy D19 midline curve of THM-4650 lies in the isogeny class of `f₁`.
  - So the K6 − P6 terminal (via `f₂`) and the D19 midline (via `f₁`) are the two tensor factors of one Siegel spinor L-function (Okazaki's conjecture, CITED). This is a shared level and a shared Siegel lift. No map between the two problems is claimed.

## Proof routes

1. **Trace formula.**
   - `b(q²) = ½ Tr(T_(q²) | S₃^new(32, χ₋₄))`, because the new space is exactly `{g, ḡ}`.
   - Hijikata's formula expresses this trace through Hurwitz class numbers, as in McCarthy–Papanikolas.
   - `M_q` is a pairing of two Legendre traces with the kernel `χ(a x(x−1)(x−a))`. Grouping pairs `(E_a, E_x)` by isomorphism class would give class-number sums.
2. **Hypergeometric.**
   - Expand `₂F₁` and `₃F₂` in multiplicative characters (Jacobi sums), so that `M_q` becomes a double character sum.
   - Then use Hasse–Davenport and McCarthy's finite-field Whipple transformation, the tool behind their factorisation (1.9).
3. **Geometric.**
   - `M_q` is a quadratic-character sum over a 4-fold, a fibre product of double covers over the `(a, x)`-plane.
   - The inner sum `h(x) = ε Σ_v χ(v(v−x)) t(v)` counts points on an elliptic K3 surface, up to trivial terms: the pull-back of the Legendre family by `s² = w(w − x)`. It has fibres `I₄, 2I₂, 2I₂*` and Picard number at least 19. This is the geometric form of Clausen's identity.
   - Statement 2 then says that the motive of the 4-fold is `Sym²(ρ_g) ⊕ Q(−1)² ⊕ Q(0)`.
