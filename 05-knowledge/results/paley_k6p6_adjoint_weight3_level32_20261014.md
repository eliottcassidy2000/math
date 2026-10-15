# The Paley K6 − P6 terminal: an exact Legendre bilinear form, conjecturally the symmetric square of a weight-3 newform of level 32

2026-10-14, mac-mini-2026-10-14-wilson. Audited (audit C, corrections MISTAKE-611).

Status of the parts:
- **PROVED:** the reduction of the terminal to one Legendre bilinear form, and its ₂F₁·₃F₂ form ([THM-4653](../../01-canon/theorems/THM-4653-paley-k6-minus-p6-terminal-reduction-to-a-legendre-bilinear-form.md)).
- **FINITE-EXACT:** the modular evaluation ([HYP-9276](../hypotheses/HYP-9276-k6-minus-p6-terminal-is-the-adjoint-of-a-weight-3-level-32-newform.md)).
  - Session: all 3,243 primes `q <= 30,000` and 17 prime powers.
  - Audit C, with independent evaluators: 10,590 primes up to 299,983 and 22 prime powers.
- **CONDITIONAL:** `Z[q]`-linear independence from the wheel, given HYP-9276 and CITED automorphic theorems.
- **PROVED (audit C):** the parity lemma. `M_q ≡ 2 (mod 8)` for `q ≡ 1 (mod 4)`, and `16 | M_q` for `q ≡ 3 (mod 4)`.

This answers the "specialized arithmetic evaluation of K6 − P6" target. The target came from [paley_next_cores_20261014i.md](paley_next_cores_20261014i.md), which left `V_q` as a formal terminal.

## 1. Result

For every odd prime power `q`, write `V_q = I_(K6−P6)(q)` for the injective six-vertex character sum (source note, eq. (4)). Let `t(a) = Σ_y χ(y(y−1)(y−a))` be the Legendre trace, and `ε = χ(−1)` (a sign here, not Greene's trivial character). Then:

`V_q = q(q−1)(−ε M_q + 7 + 13ε − 2εq)`,
`M_q = −Σ_(a,x) χ(a x (x−1)(x−a)) t(a) t(x) = −q³ Σ_x φ(x)φ(1−x) ₂F₁(x) ₃F₂(x)`.
**(PROVED)**

`M_q = −b(q²) − 2q − 1`, where `g = Σ b(n)qⁿ = q + 4iq³ + 2q⁵ − 8iq⁷ + … ∈ S₃^new(Γ₀(32), (−4/·))`. This is LMFDB 32.3.c.a, `g = η(4z)²η(8z)²(θ(4z)² + i(θ(2z)² − θ(4z)²))`.
**(HYP-9276, FINITE-EXACT; out of sample above `q = 1500`)**

For the Paley primes `q ≡ 3 (mod 4)` this gives

`V_q = q(q−1)(|b(q)|² − q² − 7)`,

so the degree-ten run coefficient becomes the following. This is CONDITIONAL on HYP-9276, for primes `q ≡ 3 (mod 4)`, `q >= 11`:

`[t¹⁰] R_q = (P₁₀(q) − q(q−1)(7q² − 173q + 1166)(a_q(f) + 4) + 64 q(q−1)(|b(q)|² − q² − 7)) / (q)₁₀`.

Here `f = η(2z)⁴η(4z)⁴` is the Apéry form of the wheel (HYP-9274 (B), PROVED + CITED for primes).

The two arithmetic coordinates differ in kind:
- the wheel is a `GL(2)` object: weight 4, level 8, `₄F₃(1)`, Ahlgren–Ono;
- the six-vertex terminal is, conjecturally, a `GL(3)` object: the symmetric square (the adjoint twisted by `χ₋₄`, up to Tate twist) of a non-CM weight-3 form of level 32.

## 2. How the graph sum becomes a form

1. **Affine normalisation.** There are ten edges and `{0, 1}` is a non-edge, so we may set `x₀ = 0`, `x₁ = 1` and drop a factor `q(q−1)`.
2. **Inclusion–exclusion over the remaining non-edges.** The four non-edges `12, 23, 34, 45` leave one main term and three explicit corrections. Each correction is linear in `q`: `3 + 3ε`, then `2 + 5ε − εq` twice.
3. **The main term is a bilinear form in Legendre traces.**
   - The path structure 2–4, 2–5, 3–5 makes vertex 5 a four-point character sum `F₅`.
   - Summing `F₅` against the weight of vertex 3 collapses it to `ε t(x)`.
   - Vertex 4 contributes `t(a)`.
   - Hence `M_q = −Σ_x χ(x(x−1)) t(x) h(x)` with `h = (χt) ∗ χ`, an additive convolution. This gives an `O(q log q)` FFT evaluation.
4. **The hypergeometric reading.** Greene's recursion gives `t = φ(−1)q ₂F₁` and `h = q² ₃F₂` exactly.
   - Greene's recursion `₄F₃(1) = (φ(−1)/q) Σ φ(x)φ(1−x) ₃F₂(x)` defines the wheel's value, and Ahlgren–Ono evaluated it (the Apéry form).
   - The six-vertex terminal is the same beta moment with one more `₂F₁(x)` inside.
   - So the Paley run hierarchy climbs the hypergeometric ladder: `₂F₁` (the Legendre trace inside every terminal), then `₄F₃(1)` (the wheel), then `₂F₁·₃F₂` (K6 − P6).

## 3. Finding and testing the form

- The fits (`fit_*.py`, kept as exploration records) first ruled out the obvious candidates. No fit `W_q = P(q) + c·Tr a_q(F)` exists for any single newform orbit `F` scanned. The scan covered 2-power levels 4..256, weights 2..6 and quadratic 2-power characters (`mf_basis.gp`, `mf_basis2.gp`), with quadratic `P`.
- Traces of squares of weight-3 coefficients (`mf_w3sq.gp`) then isolated `Re(b(q)²)` of the level-32 orbit.
- The identity was then tested as one exact equation at every prime `5 <= q <= 30,000`, and in its Frobenius form at the prime powers 9, 25, 27, 49, 81, 121, 125, 169, 243, 289, 343, 361, 625, 729, 1331, 2187 and 2401. There were no failures.
- On the prime-power side the field arithmetic is exact, and the six-vertex brute force agrees at every `q <= 49`.

The prime-power success matters. It strongly suggests that the identity is motivic, not a numerical accident of the prime field.
- The trace of `Frob_p^r` on `Sym²ρ_g` (`≅ Ad ρ_g ⊗ χ₋₄ ⊗ Q_ℓ(−2)`) is `α^(2r) + (αβ)^r + β^(2r)`. (The pre-audit text wrote `Sym²(ρ_g) ⊗ χ₋₄`, whose trace has the opposite sign for `p ≡ 3 (mod 4)`, `r` odd; MISTAKE-611.)
- This matches `−(M_(p^r) + 2p^r + 1)` at every prime power tested.
- The identity was identified from primes below 1500. Everything above 1500 and every prime power is out of sample.

## 4. Consequences

All are typed in HYP-9276.

- **(a) Deligne window, conditional.** `−q² − 7 <= W_q <= 3q² − 7` at Paley primes.
- **(b) Sato–Tate, conditional.** `(W_q + 7)/q² = 4cos²θ_q − 1` follows the `SO(3)` law.
  - On the 1,632 Paley primes the moments are 1.005, 1.000, 3.002, 5.987, 14.974, 35.937, 90.946.
  - The Riordan numbers are 1, 1, 3, 6, 15, 36, 91.
- **(c) Independence, conditional.** No relation `P + B·U_q + C·V_q = 0` with `P, B, C ∈ Z[q]` holds.
  - This is `Z[q]`-linear independence, a weaker statement than the algebraic independence the source note left open.
  - Audit C corrected the argument: Sato–Tate for `Ad g` and for the Apéry form, plus Deligne's bound.
- **(d) Zeros.** Under HYP-9276, `M_q = 0` iff `|b(q)| = q + 1`. At such primes `V_q = 2q(q−1)(q−3)`, which THM-4653 proves.
  - The recorded `V₇ = 336` and `V₃₁ = 52080` are exactly the first two such primes.
  - All twelve below 100,000 (the twelfth is 97,919) are `≡ 7 (mod 8)`.
  - For `q ≡ 1 (mod 4)`, `M_q ≠ 0` by the parity lemma (PROVED).
- **(e) 2-adic laws.**
  - PROVED: `M_q ≡ 2 (mod 8)` for `q ≡ 1 (mod 4)`, and `16 | M_q` for `q ≡ 3 (mod 4)`. The proof uses full 2-torsion of every Legendre curve.
  - FINITE-EXACT to 100,000, with sharp minima: `2⁷ | M_q` for `q ≡ 3 (mod 8)`, and `2¹⁰ | M_q` for `q ≡ 7 (mod 8)`.

## 5. Literature, and a remark on level 32

- `g` is the form `f₂` of McCarthy–Papanikolas (arXiv:1205.1006).
  - They prove `b(p) = ₃F₂(χ₄, φ, φ; ε, ε | 1)_p` for `p ≡ 1 (mod 4)` (quartic `χ₄`, Hijikata trace formula).
  - For `p ≡ 3 (mod 4)` they prove nothing about `b(p)` itself.
  - Their §6 records a numerical observation over `F_(p²)` (checked for `p < 20`). Combined with their (1.5), it already expresses `|b(p)|²` at Paley primes through `₄F₃(φ,φ,φ,φ; ε,ε,ε | −1)_(p²)`, with quadratic characters (audit C, M2).
  - HYP-9276 is, as far as we know, the first expression over the **prime** field: a `₂F₁·₃F₂` beta moment, arising as a tournament pattern sum.
- In their Siegel-modular picture, the degree-2, level-8 eigenform `F₇` has spinor L-function `L(f₁ ⊗ f₂)`. That is the van Geemen–van Straten conjecture (Math. Comp. 61 (1993)), reported resolved by Okazaki (arXiv:1008.1602); CITED, not re-checked.
  - Here `f₁ = η(4z)²η(8z)²` is the conductor-32 congruent-number curve `y² = x³ − x`.
  - That is the isogeny class of the D19 midline curve in [THM-4650](../../01-canon/theorems/THM-4650-guy-d19-special-loci-carry-no-rational-distance-point-and-the-midline-is-fermats-right-triangle.md).
- **REMARK (coincidence of level).** The congruent-number curve that governs the D19 midline and the form behind the K6 − P6 terminal are both level-32 `Q(i)`-objects, paired in the van Geemen–van Straten Siegel form.
  - Level 32 is ubiquitous for `Q(i)`-CM objects.
  - The pairing runs through an unrelated Siegel form.
  - This is not a structural link between D19 and Paley patterns. (The pre-audit text typed it STRUCTURAL and said the problems "meet"; MISTAKE-611.)
- **Further searches.** McCarthy–Springfield (arXiv:2311.02135) and Li–Zhou (arXiv:2512.19312) treat Paley tournaments and graphs, but no six-vertex pattern.

## 6. Reproduction

All scripts are in `04-computation/experiments/paley_k6p6_20261014/`. They use exact integers, plus float64 FFTs with checked rounding.

| Script | What it does | Time |
|---|---|---|
| `k6p6_values.py` | `V_q` by the O(q³) reduction, with O(q⁴) brute force for `q <= 43` | 6 s |
| `k6p6_proof_steps.py` | Every intermediate identity of THM-4653, at the 24 primes 5..101 | |
| `k6p6_fft.py 30000` | `M_q` at all 3,243 primes | 26 s |
| `g_coeffs.gp`, `identify_g.gp` | The newform; writes `g_trace_a2_half.txt` | |
| `k6p6_identity_check.py` | HYP-9276 at all primes, plus moments, zeros and 2-adic data | |
| `k6p6_prime_powers.py` | Exact `F_(p^r)` arithmetic, THM-4653 and HYP-9276 at 25 fields | |

`fit_k6p6.py`, `fit_M.py`, `fit_products.py`, `fit_M_products.py`, `moments.py` and `mf_*.gp` are the exploratory searches that led to the identification. They are kept as records of the negative fits. The two large basis dumps `mf_basis.out` and `mf_basis2.out` are not committed; rerun the `.gp` files with `gp -q` to regenerate them.

## 7. Next obligations

1. **Prove HYP-9276.** Four routes are listed in the HYP: trace formula, hypergeometric, geometric, and Faltings–Serre–Livné. The most concrete is the **descent form** (HYP-9276 statement 3):

   `q³ Σ_x φ(x)φ(1−x) ₂F₁(x)₃F₂(x) = H_(q²) + χ₋₄(q)q² + 2q + 1`, with `H_Q = ₃F₂(χ₄, φ, φ; ε, ε | 1)^*_Q` over `F_(q²)`.

   - It is checked directly at `q = 5..43`.
   - It contains no modular form. It is equivalent to HYP-9276 given the prime-power extension of McCarthy–Papanikolas Theorem 1.8, which is FINITE-EXACT at 20 prime powers, including the inert squares `p²` with `p ≡ 3 (mod 4)`.
   - Geometrically it suggests that the K6 − P6 4-fold is, up to Tate classes, the Weil restriction from `Q(i)` to `Q` of the `₃F₂(1)` surface: tensor induction of `ρ_g|G_(Q(i))` contains `Sym²ρ_g`.
   - The trace-formula route is closest to a finished proof. The new space is exactly `{g, ḡ}`, so `b(q²)` is half a Hecke trace, as in McCarthy–Papanikolas.
   - What remains is to turn the Legendre pairing into Hurwitz class numbers.
2. **Explain the 2-adic refinements (e).** The factor 16 is proved; the `2⁷` / `2¹⁰` refinements look like a congruence between `g` and a CM or Eisenstein form.
3. **Climb further.** At twelve steps (HYP-9275, opus S35) new terminals may bring `Sym³` or `GL(2) × GL(3)` objects.
   - The ladder `₂F₁ → ₄F₃(1) → ₂F₁·₃F₂` suggests watching for `₃F₂·₃F₂` and `₂F₁²·₃F₂` moments.
