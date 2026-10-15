---
id: HYP-9276
title: "The Paley K6 - P6 terminal is the symmetric square (= adjoint x chi_(-4), up to Tate twist) of the McCarthy-Papanikolas weight-3 newform of level 32 (LMFDB 32.3.c.a): M_q = -b(q^2) - 2q - 1, V_q = q(q-1)(|b(q)|^2 - q^2 + 7 + 14 chi_(-4)(q))"
status: >
  OPEN (a modular identity). FINITE-EXACT: all 3,243 primes 5 <= q <= 30,000 (session); independently (audit C, two new
  evaluators) 10,590 primes (every prime to 100,000 plus 1,000 random primes to 299,983) and 22 prime powers to 2401, with no
  failure; the Frobenius form of statement 2 at 9, 25, 27, 49, 81, 121, 125, 169, 243, 289, 343, 361, 529, 625, 729, 841, 961,
  1331, 1369, 2187, 2197, 2401 (Paley prime powers 27, 243, 343, 1331, 2187). The identity was identified from primes below
  1500 (the exploratory fits); everything above 1500 and every prime power is out of sample. The reduction it feeds (THM-4653)
  is PROVED. PROVED unconditionally here (audit C's parity lemma): M_q == 2 (mod 8) for q == 1 (mod 4) and 16 | M_q for
  q == 3 (mod 4). Corollaries (a)-(c) and the zero criterion in (d) are CONDITIONAL on this hypothesis (with CITED theorems).
  The 2^7 / 2^10 refinements in (e) are FINITE-EXACT to 100,000. Audited (audit C, MISTAKE-611).
session: mac-mini-2026-10-14-wilson
source: 05-knowledge/results/paley_k6p6_adjoint_weight3_level32_20261014.md
scripts:
  - 04-computation/experiments/paley_k6p6_20261014/k6p6_fft.py (+ .out, .json): M_q at all primes to 30,000
  - 04-computation/experiments/paley_k6p6_20261014/identify_g.gp, g_coeffs.gp (+ .out), g_trace_a2_half.txt: the newform (gp mf package)
  - 04-computation/experiments/paley_k6p6_20261014/k6p6_identity_check.py (+ .out): the identity, Sato-Tate moments, zeros, 2-adic data
  - 04-computation/experiments/paley_k6p6_20261014/k6p6_prime_powers.py (+ .out): prime powers to 361, exact field arithmetic
  - 04-computation/experiments/paley_k6p6_20261014/k6p6_prime_powers_large.py (+ .out): 625, 729, 1331, 2187, 2401 (O(q^2), log tables)
  - 04-computation/experiments/audit_C_20261014/ (independent evaluators, the eta-theta construction of g, the parity lemma, MP section 6)
related: [THM-4653, HYP-9274, HYP-9275, THM-4650 (D19 midline curve of conductor 32), paley_next_cores_20261014i.md]
---

# HYP-9276: the K6 − P6 terminal is a symmetric square

## The newform

Let `g = Σ b(n) qⁿ` be the newform of weight 3 on `Γ₀(32)` with nebentypus `χ₋₄ = (−4/·)`, LMFDB label **32.3.c.a**:

`g = q + 4i q³ + 2q⁵ − 8i q⁷ − 7q⁹ − 4i q¹¹ − 14q¹³ + 8i q¹⁵ + 18q¹⁷ − 12i q¹⁹ + …`
`  = η(4z)²η(8z)² · (θ(4z)² + i(θ(2z)² − θ(4z)²))`.

The closed form is from audit C: it is the `4i`-eigenvector of `T₃` on `S₃(Γ₀(32), χ₋₄)`, proved by linear algebra up to the Sturm bound. We do not claim it is new.

- `S₃^new(Γ₀(32), χ₋₄)` has dimension 2. It is spanned by `g` and `ḡ = g ⊗ χ₋₄`.
- `g` has no CM. Its coefficient field is `Q(i)`. Its Sato–Tate group is `SU(2)[C₂]` (inner twist by `χ₋₄`).
- `b(p)` is real for `p ≡ 1 (mod 4)` and purely imaginary for `p ≡ 3 (mod 4)`.
- This is the form `f₂` of McCarthy–Papanikolas (CITED, arXiv:1205.1006, IJNT 2015, Lemma 5.1).

At an odd prime `p`, with Satake roots `α + β = b(p)` and `αβ = χ₋₄(p)p²`,

`b(p²) = α² + αβ + β² = χ₋₄(p)(|b(p)|² − p²)`.

- This is the trace of Frobenius on **`Sym²ρ_g`**, and `Sym²ρ_g ≅ Ad(ρ_g) ⊗ χ₋₄ ⊗ Q_ℓ(−2)`.
- The adjoint `Ad = Sym² ⊗ det⁻¹` has trace `x_p = b(p²)/(χ₋₄(p)p²) = (|b(p)|² − p²)/p²`.
- `Ad g` is a cuspidal `GL(3)` object, because `g` has no CM (Gelbart–Jacquet).

## Statement

Let `M_q` be the Legendre bilinear form of THM-4653. Write `φ(−1)` for the sign of `−1` in `F_q`. (`ε` below is Greene's trivial character only where it appears as a hypergeometric parameter.)

1. **(Primes.)** For every prime `q >= 5`,

   `M_q = −b(q²) − 2q − 1`.

   Equivalently, `Σ_x φ(x)φ(1−x) ₂F₁(x) ₃F₂(x) = (b(q²) + 2q + 1)/q³` in Greene's normalisation. By THM-4653,

   `V_q = I_(K6−P6)(q) = q(q−1)(|b(q)|² − q² + 7 + 14χ₋₄(q))`.

   For the Paley primes `q ≡ 3 (mod 4)` this reads `V_q = q(q−1)(|b(q)|² − q² − 7) = −q(q−1)(b(q²) + 7)`. Note that `|b(q)|² = −b(q)²`.

2. **(Prime powers.)** For `q = p^r` with `p` odd,

   `M_q = −(α^(2r) + (αβ)^r + β^(2r)) − 2q − 1`,

   which is minus the trace of `Frob_p^r` on `Sym²ρ_g`, less `2q + 1`. Then `V_q = q(q−1)[φ(−1)(S_q + 14) + 7]` with `S_q = α^(2r) + (αβ)^r + β^(2r)`.

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

The recorded values `V₇ = 336`, `V₁₁ = −12320`, `V₁₉ = −76608` and `V₃₁ = 52080` of the source note are reproduced. Audit C reproduced them again by fully unnormalised enumeration.

## Evidence (FINITE-EXACT)

**Session.**
- All 3,243 primes `5 <= q <= 30,000`.
- `M_q` comes from the O(q log q) additive FFT evaluator of THM-4653 (3). It agrees with the exact O(q³) values at all 237 primes below 1500. (Audit C notes that integrality of the rounded values proves nothing; its exact evaluators supply the exactness check.)
- `b(q)` comes from gp's `mfcoefs`.
- Prime powers to 2401.

**Audit C (independent code).**
- An exact integer evaluator and a multiplicative-FFT evaluator (worst rounding gap `2.4·10^(−9)`) agree at all 758 common primes.
- The identity holds at 10,590 primes up to 299,983 and at 22 prime powers.
- `b(p)` was built independently from eta and theta products, checked against PARI and LMFDB, and verified Hecke-consistent to `N = 300,000`.

**How it was found.** The identity was found by fitting `W_q` against polynomials and traces of newform orbits (`fit_*.py`, `mf_*.gp`), using primes below 1500. It is out of sample everywhere above 1500 and at every prime power.

## Corollaries

**(a) Size (CONDITIONAL; uses Deligne).** Deligne gives `|b(q)| <= 2q`. Hence `−q² − 7 <= W_q <= 3q² − 7` at Paley primes, and `|V_q| <= 3q³(q−1)`.

**(b) Distribution (CONDITIONAL; uses Sato–Tate for non-CM newforms).** Let `x_q = (|b(q)|² − q²)/q² = 4cos²θ_q − 1`, the trace of `Ad ρ_g` normalised.
- Sato–Tate for non-CM newforms with nebentypus (Barnet-Lamb–Geraghty–Harris–Taylor 2011, PRIMS 47; the version for each value of the nebentypus is in that paper's main theorem as we recall it, not re-checked) gives the `SO(3)` law for `x_q` in each residue class mod 4.
- Its moments are the Riordan numbers `1, 0, 1, 1, 3, 6, 15, 36, 91`.

| Sample | Moments k = 2..8 |
|---|---|
| 1,632 Paley primes below 30,000 | 1.005, 1.000, 3.002, 5.987, 14.974, 35.937, 90.946 |
| 1,611 primes `q ≡ 1 (mod 4)` | 0.966, 0.902, 2.751, 5.349, 13.304, 31.548, 79.220 |

- The second row is within about one standard error (`z ≈ −1.0` to `−1.1`).
- To 100,000 (audit C), the Paley class has `|z| <= 0.48`.

**(c) Linear independence from the wheel (CONDITIONAL).** There are no `P, B, C ∈ Z[q]`, not all zero, with `P(q) + B(q)U_q + C(q)V_q = 0` for all large Paley primes. This is `Z[q]`-linear independence of `1, U_q, V_q`, not algebraic independence.

*Argument (audit C's corrected version).* Write `U_q = −q(q−1)(a_q + 4)` with `y_q = a_q/q^(3/2)`, the Apéry form `f` being non-CM of weight 4. A relation with `C ≠ 0` gives `x_q = R(q) + β(q)y_q`, with `R` rational and `β = B/(C√q)` of half-integral degree.
- If `β → 0` (including `B = 0`), then `x_q` converges along Paley primes.
- If `|β| → ∞`, or if `C = 0 ≠ B`, then `y_q → 0` or `y_q` is unbounded.
- Every case contradicts one of three facts along `p ≡ 3 (mod 4)`:
  1. `x_q` does not converge: `E[x] = 0` and `E[x²] = 1`, from the L-functions of `Ad g ⊗ χ₋₄^j` and `Ad g × Ad g ⊗ χ₋₄^j` (`j = 0, 1`), since `Ad g` is cuspidal and `Ad g ≇ Ad g ⊗ χ₋₄`;
  2. `y_q` does not tend to 0: `E[y²] = 1`, by Rankin–Selberg for `f × f` and `f × (f ⊗ χ₋₄)`;
  3. Deligne's bound.

**(d) Zeros.**
- At Paley primes, `M_q = (|b(q)| − q − 1)(|b(q)| + q + 1)`. So `M_q = 0` iff `|b(q)| = q + 1` (CONDITIONAL on statement 1).
- At such primes `V_q = 2q(q−1)(q−3)`. This is PROVED by THM-4653 alone.
- Zeros below 100,000 (FINITE-EXACT): `q = 7, 31, 199, 223, 1151, 2239, 2663, 3391, 3583, 23719, 23743, 97919`, all `≡ 7 (mod 8)`.
- For `q ≡ 1 (mod 4)`, `M_q` never vanishes. This is PROVED unconditionally by the parity lemma below.

**(e) 2-adic laws.**
- **Parity lemma (PROVED; audit C):** `M_q ≡ 2 (mod 8)` for `q ≡ 1 (mod 4)`, and `M_q ≡ 0 (mod 16)` for `q ≡ 3 (mod 4)`, for every odd prime power `q`.
  - *Proof.* The summand is nonzero only for `a ≠ 0`, `x ∉ {0, 1, a}`.
  - The terms with `a = 1` sum to `−(1 + φ(−1))`.
  - For `a ∉ {0, 1}`, `E_a` has full rational 2-torsion, so `4 | #E_a` and `t(a) ≡ −(q + 1) (mod 4)`.
  - For `q ≡ 1 (mod 4)`, each of the `(q−2)(q−3) ≡ 2 (mod 4)` remaining terms is `≡ 4 (mod 8)`.
  - For `q ≡ 3 (mod 4)`, every term is `≡ 0 (mod 16)`.
  - Checked on 5,304 + 5,308 fields. ∎
- Under statement 1, the `q ≡ 1 (mod 4)` case is equivalent to `b(q) ≡ 2 (mod 4)`.
- **Refinements (FINITE-EXACT to 100,000; sharp minima at q = 11 and q = 23):** `2⁷ | M_q` for `q ≡ 3 (mod 8)`, and `2¹⁰ | M_q` for `q ≡ 7 (mod 8)`. Equivalently `|b(q)|² ≡ (q+1)²` modulo `2⁷`, respectively `2¹⁰`. Partially explained: the factor 16 is the parity lemma.

## Literature (CITED and typed)

- **McCarthy–Papanikolas (arXiv:1205.1006; IJNT 2015).**
  - Their Theorem 1.8 gives `b(p) = ₃F₂(χ₄, φ, φ; ε, ε | 1)_p` for `p ≡ 1 (mod 4)` (McCarthy's normalisation, quartic `χ₄`). The proof uses Hijikata's trace formula and Schoof's count of Legendre curves with `Z/2 × Z/8` torsion.
  - Their Theorem 1.6, conditional on the van Geemen–van Straten conjecture, reads `₄F₃(φ,φ,φ,φ; ε,ε,ε | −1)_p = ξ(p)λ(p)` with `ξ = (2/·)`. Combined with their (1.5), the value is **exactly `a(p)b(p)`** for `p ≡ 1 (mod 4)`, where `a(p)` is the coefficient of `f₁ = η(4z)²η(8z)²`, and it vanishes for `p ≡ 3 (mod 4)`.
  - **Closest prior statement (their §6, numerical observation for p < 20).** They observe `a_p² − 2a_(p²) = ₄F₃(φ,φ,φ,φ; ε,ε,ε | −1)_(p²)` over `F_(p²)`, where `a_p = λ(p)`.
    - Through their (1.5) this equals `(a(p)² − 2p)(b(p)² − 2χ₋₄(p)p²)`.
    - So at Paley primes (`a(p) = 0`) it already gives `|b(p)|² = 2p² + ₄F₃(−1)_(p²)/(2p)`, a quadratic-character expression over `F_(p²)`.
    - Audit C verified this for `p <= 23`.
  - Statement 1 is, as far as we know, the first expression for `|b(q)|²` over the **prime** field `F_q` (a `₂F₁·₃F₂` beta moment), the first with a graph-sum origin, and the first in a prime-power Frobenius form.
- **Siegel side.** The conjecture that the spinor L-function of their level-8 Siegel eigenform `F₇` is `L(f₁ ⊗ f₂)` is the **van Geemen–van Straten conjecture** (Math. Comp. 61 (1993)). It is reported resolved by Okazaki (arXiv:1008.1602); we did not check that `F₇` is covered.
- **Other searches.**
  - Dawsey–McCarthy survey (arXiv:2101.06303; abstract only).
  - Paley clique counts `K₃`/`K₄` (Evans–Pulham–Sheehan; Dawsey–McCarthy).
  - McCarthy–Springfield (arXiv:2311.02135): transitive subtournaments of order 4 in Paley digraphs.
  - Li–Zhou (arXiv:2512.19312): asymptotic induced counts.
  - None of these treats six-vertex patterns or this identity.
- **REMARK (coincidence of level).** `f₁` is the congruent-number curve `y² = x³ − x`, the isogeny class of the D19 midline curve (THM-4650). So `f₁` and `f₂` are both level-32 `Q(i)`-objects tied by the van Geemen–van Straten Siegel form. Level 32 is ubiquitous for `Q(i)`-CM objects, and the link runs through an unrelated Siegel form. This is not a structural connection between D19 and Paley patterns.

## Proof routes

1. **Trace formula.**
   - `b(q²) = ½ Tr(T_(q²) | S₃^new(32, χ₋₄))`, because the new space is exactly `{g, ḡ}`.
   - Hijikata's formula expresses this trace through Hurwitz class numbers, as McCarthy and Papanikolas do for `T_p`.
   - `M_q` pairs two Legendre traces with the kernel `χ(a x (x−1)(x−a))`. Grouping pairs `(E_a, E_x)` by isomorphism class would give class-number sums.
2. **Hypergeometric.**
   - Expand `₂F₁` and `₃F₂` in multiplicative characters, so that `M_q` becomes a double Jacobi-sum series.
   - Compare with McCarthy–Papanikolas's `₄F₃(−1)` over `F_(q²)` (their §6) and with `₃F₂(χ₄, φ, φ | 1)²` for `q ≡ 1 (mod 4)`.
   - Tools: Hasse–Davenport and McCarthy's finite-field Whipple transformation.
3. **Geometric (expected, not established).**
   - `M_q` is a quadratic-character sum over a 4-fold, a fibre product of double covers over the `(a, x)`-plane.
   - The inner sum `h(x) = φ(−1) Σ_v χ(v(v−x)) t(v)` counts points on an elliptic K3 surface, up to trivial terms: the pull-back of the Legendre family by `s² = w(w − x)`.
     - It has fibres `I₄, 2I₂, 2I₂*`.
     - By Shioda–Tate its trivial lattice has rank `2 + (4−1) + 2·(2−1) + 2·(7−1) = 19`, since `I₂*` has 7 components. So its Picard number is at least 19. Audit C read rank 17; we keep 19 and flag it for re-check.
     - This is the geometric form of Clausen's identity.
   - Statement 2 then suggests (not proved) that the relevant motive of the 4-fold is `Sym²(ρ_g) ⊕ Q(−1)² ⊕ Q(0)`.
4. **Faltings–Serre–Livné.**
   - Both sides should be traces of 2-adic Galois representations ramified only at 2.
   - The left side would be a 3-dimensional piece of the cohomology of the local system built from the Legendre system by two middle convolutions with `χ` (Katz's algorithm), after removing Tate classes.
   - Once that piece is shown to be 3-dimensional, unramified outside 2, with the right determinant, a finite comparison of traces would suffice. The check to 300,000 far exceeds any such finite test set.
   - What is missing is the structural input: the rank computation via Katz's algorithm, and a Faltings–Serre criterion for 3-dimensional orthogonal representations.
