# HYP-2322 — Rising/falling factorials are the factorial shadow of Sₙ/Aₙ: rising@1 = |Sₙ|, falling@1 = 0 (the sign/index-2), exchanged by the antipode σ

**Session:** S644
**Status:** CONFIRMED (formalized + verified against permutation cycle counts)
**Provenance forward:** math-lean `Math/Combinatorics/PochhammerSign.lean` (sorry-free)
**Prompt:** continue the investigation, keeping in mind rising and falling factorials.

---

## 0. The connection

The S642/S643 Galois arc turns on `Sₙ` (all permutations; unsolvable, Abel–Ruffini) vs `Aₙ` (the sign
kernel; the solvable/even half, the `Aut(T) ⊆ Aₙ` of S643). Rising and falling factorials are exactly
these two, made into polynomials, with the arc's antipode `σ : x ↦ −x` exchanging them.

| | polynomial | coefficients | at `x = 1` |
|---|---|---|---|
| **rising** `x^(n) = x(x+1)⋯(x+n−1)` | `ascPochhammer` | unsigned cycle counts `c(n,k) = #{σ : k cycles}` | `Σ c(n,k) = n! = |Sₙ|` |
| **falling** `(x)_n = x(x−1)⋯(x−n+1)` | `descPochhammer` | signed `s(n,k) = (−1)^{n−k} c(n,k)` | `Σ s(n,k) = #even − #odd = 0` (`n≥2`) |

Since `sign σ = (−1)^{n − cyc σ}`, the falling factorial's coefficients are the *signed* counts, so its
value at `1` is `#even − #odd`. And `(1)_n = 1·0·(−1)⋯ = 0` for `n ≥ 2` ⟹ `#even = #odd` ⟹
`|Aₙ| = n!/2`. The **rising** factorial counts *all* of `Sₙ`; the **falling** factorial's vanishing
*is* the index-2 sign structure (the kernel `Aₙ`).

---

## 1. Formalized (math-lean, sorry-free): `Math/Combinatorics/PochhammerSign.lean`

- `ascPochhammer_eval_one : (ascPochhammer ℤ n).eval 1 = n!` — **rising@1 = |Sₙ|** (all permutations,
  the generic Galois group, the Abel–Ruffini side).
- `descPochhammer_eval_one : 2 ≤ n → (descPochhammer ℤ n).eval 1 = 0` — **falling@1 = 0**: the even/odd
  balance, `|Aₙ| = n!/2`, the sign/index-2 (the Feit–Thompson / `Aut(T) ⊆ Aₙ` side, S643).
- `descPochhammer_eq_neg_ascPochhammer : (descPochhammer ℤ n).eval x = (−1)^n (ascPochhammer ℤ n).eval
  (−x)` — **the antipode `σ : x ↦ −x` exchanges rising ↔ falling.**

Verified against actual permutations (`pochhammer_sign_permutations_s644.py`): falling coeffs = signed
Stirling-first `s(n,k)`, rising = unsigned `c(n,k)`, both matching the true cycle-count histograms
through `n = 6`; rising@1 = `n!` and falling@1 = `0` (n≥2) with `#even = #odd = n!/2` through `n = 7`; the
reflection `(x)_n = (−1)^n(−x)^(n)` exact through `n = 5`.

---

## 2. What it ties together

1. **The σ-antipode IS the rising↔falling duality.** The involution that has run through the whole arc
   (S638 the swap/converse, S643 the exiled anti-automorphism, S640 the doubling `v↦−v`) is precisely the
   reflection `x ↦ −x` that turns the rising factorial into the falling one (`descPochhammer_eq_neg_…`).
   So `Sₙ`(rising, all) and `Aₙ`(falling, sign) are *exchanged by σ* — the same σ that tournaments forbid
   as an automorphism (S643). The factorial pair is the polynomial avatar of the `Sₙ`/`Aₙ` split.
2. **The "0 switches on at `n = 2`."** `(x)_n` at `1` is `1` for `n = 0, 1` and `0` for `n ≥ 2` — it
   vanishes *exactly* when the sign map `Sₙ → {±1}` becomes onto (`n ≥ 2`), i.e. when `Aₙ` becomes a
   proper index-2 subgroup. The falling factorial's root at `x = 1` is the *onset of the sign structure*.
3. **Falling factorial = AP-rooted polynomial.** `(x)_n = ∏_{i=0}^{n−1}(x − i)` has roots the arithmetic
   progression `0,1,…,n−1` — the canonical additive chain (S617's AP/collapse family). The discriminant
   of `(x)_n` is `∏_{i<j}(i−j)² = ∏_{k<n} k!²` (a square), so `Gal ⊆ Aₙ` — consistent with §1: the AP
   (most degenerate) polynomial sits inside the alternating/sign-kernel side, like `Aut(T)` (S643).
4. **The FTA off-by-one (HYP-2275/S636) again.** `Δ(x)_n = n·(x)_{n−1}` (the falling factorial is the
   discrete analogue of `xⁿ` for the difference operator): the factorials are the natural basis where the
   `n → n−1` shift is clean — the same off-by-one the arc's roots↔coefficients duality runs on.

---

## 3. New threads
- **The discriminant of `(x)_n` is `∏_{k<n} k!²` (a square)** — connects S643's open thread #3 ("the
  tournament discriminant is always a square, find it"): the AP/falling-factorial polynomial is the model
  case where the discriminant is a perfect square and Galois ⊆ Aₙ. Is the tournament's discriminant a
  product-of-factorials / Vandermonde square?
- **Stirling numbers as the bridge.** The rising→monomial→falling transitions are the (unsigned/signed)
  Stirling numbers; the arc's "counts ↔ spectrum" (Vieta, HYP-2275) is here "cycle-counts ↔ factorial
  basis." Worth formalizing `ascPochhammer` coeff = unsigned Stirling-first if/when Mathlib supports it.
- **q-analogue.** The q-Pochhammer / Gaussian binomials are the cyclotomic deformation (ties the
  q-Krawtchouk, S622); at `q → 1` they degenerate to these. The rising/falling pair at a root of unity
  may meet the Φ₃ / forbidden-`7` resonance.
- **The fiber fraction `(1/2)_{n−2}/(n−2)!`** (CLAUDE.md, the tournament-tiling `√π`/Wallis constant) is a
  *rising factorial of `1/2`* — the half-integer Pochhammer. The `1/2` (the apex/σ-fixed point, n/2) sits
  in the rising factorial governing the tiling density. A direct rising-factorial appearance in the
  tournament metagraph worth chasing.

## 4. Open / handoffs
- Formalize the discriminant of `descPochhammer` = `(∏ k!)²` and hence `Gal((x)_n) ⊆ Aₙ`.
- Connect §3 (fiber fraction = `(1/2)`-rising-factorial) to the covering-depth `p₀` line.
- The q-deformation at `ζ₃` (the forbidden-7 resonance) — does the q-Pochhammer see `Φ₃`?
