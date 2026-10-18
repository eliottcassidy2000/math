# Beukers' mod p^2 congruences for Cooper's s10, s18 and the Hadamard products, and the inert CM primes

2026-10-18, session `cloud-2026-10-18-M1` (temporary identity of a cloud thread working beside opus-2026-10-18-S43; no
`.machine-id`, not a registered agent).

**Typing.**
- PROVED here, given the CITED results: Lemma 1, Theorem 1, Corollary 2, Theorem 3, Theorem 4, Corollary 5.
  - Theorem 1 is the congruence `F(t) ≡ F_p(t) F(t^σ) (mod p^2)` that Beukers (arXiv:2403.03301v3, Remark 1.14(2))
    reports from computer experiments for Cooper's s7, s10, s18 and the two Hadamard products of his Appendix C. His
    Theorem 1.13 does not cover them. For s7 it is Theorem A of
    `cooper_s7_weight_one_root_frobenius_20261018.md`; the proof here is shorter and different.
  - Theorem 4 (inert primes) answers the case Beukers leaves open in Remark 1.17 ("If p does not split in Q(α) we do
    not know if our modular form approach is possible") for every family of his Appendix C except
    `F(1/6, 1/2, 5/6; 1, 1 | 1728t)`, whose `Θ` is not a modular form (his section 6). Single cases were proved
    earlier by other methods (Corollary 5 is one). A brief literature search found no general result of this kind, so
    novelty is not established.
  - Corollary 5 is Z.-W. Sun's conjecture on `Σ_(k<p) A_k (mod p^2)` (arXiv:1101.1946; Conjecture 1.2 in Guo-Zeng,
    J. Number Theory 132 (2012) 1731-1740). C. Wang and Z.-W. Sun proved it with the `p`-adic Gamma function (Nanjing
    Univ. J. Math. Biquarterly 41 (2024) 34-56, as cited in Sun, arXiv:2505.02767). The proof here is different.
- CITED: Beukers, arXiv:2403.03301v3: the polynomial `T_p` of his (6), the degree bound `deg T_(p,m) <= m r`
  (section 5), Proposition 5.1(i)-(iii), Lemma 3.2 with Remark 3.3, Lemma 3.7, the congruence
  `T_(p,p) + T_(p,p+1) F/F(t^σ) ≡ 0 (mod p^3)` from the proof of Theorem 1.13, the arguments of Corollaries 4.3 and 5.3
  and of Theorem 1.16, and the data of Appendix C. Shimura's reciprocity law (Shimura, *Introduction to the arithmetic
  theory of automorphic functions*, 1971, Theorem 6.31). Babbage's congruence `C(2p-1, p-1) ≡ 1 (mod p^2)`, `p >= 3`.
- FINITE-EXACT (`beukers_remark114_power_structure_check_20261018.gp`, output `.out`, 16 s): the identities of section 1
  for `n <= 120`; Lemma 1 and Corollary 2 for `p <= 61`; Theorem 1 for `p < 110` (coefficients up to `t^(3p)`);
  24 Newton polygons; all 118 rational CM values of Beukers' Appendix C for `5 <= p < 400` (4600 inert and 4268 split
  prime-value pairs); Sun's conjecture for `3 <= p < 400`.
- VERIFIED (numerics only, no proof): Theorem 1 also holds at the primes `p | N`; `v_p(F_p(t0)) = 1` at the 38
  ramified prime-value pairs.
- OPEN: the higher-power congruences; the primes `p | N`; ramified primes. Theorem 1 and Theorem 4 fail modulo `p^3`.
- Corrections to Beukers' Appendix C (arXiv v3) are listed in section 6.

## 1. Setting

The five families (Beukers' 9.11-9.15). `F(t) = Σ g_n t^n`, `F_p(t) = Σ_(n<p) g_n t^n`; `t(τ) = q + O(q^2)` is a
Hauptmodul for an Atkin-Lehner extension `Γ` of `Γ0(N)`; `Θ(τ) = F(t(τ))` has weight 2 for `Γ` with a character
trivial on `Γ0(N)`, and its only zero modulo `Γ` is at the pole of `t`, of order `r` in the parameter `1/t`.
`η_m = η(mτ)`.

| family | `g_n` | `N` | `Γ` | `t^a Θ^b` (η-product `E`) | `r = a/b` |
| --- | --- | --- | --- | --- | --- |
| s7 | `Σ_k C(n,k)^2 C(n+k,k) C(2k,n)` | 7 | `Γ0(7)+<7>` | `t^2 Θ^3 = (η1 η7)^6` | 2/3 |
| s10 | `Σ_k C(n,k)^4` | 10 | `Γ0(10)+<2,5>` | `t^3 Θ^4 = (η1 η2 η5 η10)^4` | 3/4 |
| s18 | `Σ_(k<=n/2) (-1)^k C(n,k) C(2k,k) C(2n-2k,n-k) (C(2n-3k-1,n) + C(2n-3k,n))`, `C(-1,0) = 0` | 18 | `Γ0(18)+<2,9,w>` | `t^3 Θ^4 = (η3 η6)^8` | 3/4 |
| Apéry2 ⊙ central | `C(2n,n) Σ_k C(n,k)^2 C(n+k,k)` | 5 | `Γ0(5)+<5>` | `t Θ^2 = (η1 η5)^4` | 1/2 |
| Franel ⊙ central | `C(2n,n) Σ_k C(n,k)^3` | 6 | `Γ0(6)+<2,3>` | `t Θ^2 = (η1 η2 η3 η6)^2` | 1/2 |

In every case `a = b - 1`, and `Σ m e_m = 24a` for `E = Π η_m^(e_m)`. Check I verifies the binomial formulas against
the differential equations, `q(t) ∈ Z[[t]]`, `Θ ∈ Z[[q]]` and `F(t(q)) = Θ(q)` with `Θ := (E/t^a)^(1/b)`.

For `p ∤ N`, Beukers' polynomial is

```
T_p(x) = (1 - (Θ(pτ)/Θ(τ)) x) · Π_(k=0)^(p-1) (p^2 - (Θ((τ+k)/p)/Θ(τ)) x)  =  Σ_m T_(p,m)(t) x^m.
```

Its coefficients are polynomials in `t` with `deg T_(p,m) <= m r`, `T_p ∈ Z[x, t]`, `T_(p,m) ≡ 0 (mod p^3)` for
`m < p`, and `T_(p,p+1) ≡ 1 (mod p)` (Beukers, section 5 and Proposition 5.1(i)-(iii); none of this uses his
Assumption (II)(4), which is what fails for these families). `t^σ(τ) = t(pτ)` is the excellent lift, a power series in
`t` with integer coefficients and `t^σ ≡ t^p (mod t^(p+1))`.

## 2. The leading coefficient is an a-th power

**Lemma 1 (PROVED).** Let `p ∤ N` and `N_p(τ) = t(pτ) Π_k t((τ+k)/p)`. Then
`T_(p,p+1)^b = (-1)^((p+1)b + (p-1)a) (t^(p+1)/N_p)^a`, and

```
T_(p,p+1) = (-1)^(p+1) Q^a,   Q ∈ Z[t],  Q(0) = 1,  Q ≡ 1 (mod p),  deg Q <= ⌊(p+1)/b⌋.
```

*Proof.* From the definition, `T_(p,p+1) = (-1)^(p+1) Θ(pτ) Π_k Θ((τ+k)/p) / Θ(τ)^(p+1)`. Raise to the power `b` and
substitute `Θ^b = E t^(-a)`. Beukers' Lemma 3.2 (with Remark 3.3, all `m | N` prime to `p`) gives
`E(pτ) Π_k E((τ+k)/p) = exp(2πi (p-1) Σ m e_m / 48) E^(p+1) = (-1)^((p-1)a) E^(p+1)`, which is the first claim.
`N_p` is a symmetric function of the `t|γ`, `[γ] ∈ Γ0(N)\M(N,p)`, which `Γ` permutes (the weight-0 case of the proof of
Lemma 3.7). So `N_p` is a modular function for `Γ` with rational coefficients, hence `H := t^(p+1)/N_p ∈ Q(t)`.
Now `H^a = ±T_(p,p+1)^b` is a polynomial. Comparing multiplicities of irreducible factors in `Q[t]`, with
`gcd(a, b) = 1`, shows `H ∈ Q[t]` and `T_(p,p+1) = c Q^a` with `Q(0) = 1`. At `q = 0` every `Θ` factor is 1, so
`c = (-1)^(p+1)`. Gauss's lemma gives `Q ∈ Z[t]`. Since `T_(p,p+1) ≡ 1 (mod p)`, the reduction of `Q^a` is a constant,
so `Q ≡ Q(0) = 1 (mod p)`. Finally `a deg Q = deg T_(p,p+1) <= ⌊(p+1)a/b⌋`. ∎

Check T: for every `p <= 61`, `p ∤ N`, `T_(p,p+1)` is `(-1)^(p+1)` times an `a`-th power of a polynomial
`Q ≡ 1 (mod p)`, and both degree bounds are attained (for s10 at `p = 61`: `deg T_(p,p) = 45`, `deg Q = 15`).

## 3. Beukers' Remark 1.14(2)

**Theorem 1 (PROVED).** For each of the five families and every prime `p ∤ N`,
`F(t) ≡ F_p(t) F(t^σ) (mod p^2)` as power series in `t`.

*Proof.* As in the proof of Beukers' Theorem 1.13, `T_p(Θ(τ)/Θ(pτ)) = 0` and Proposition 5.1(ii) give
`T_(p,p) + T_(p,p+1) Λ ≡ 0 (mod p^3)` with `Λ = F(t)/F(t^σ)`. By Lemma 1, `T_(p,p+1) = (-1)^(p+1)(1 + pB)^a` with
`B ∈ Z[t]`, `deg B <= ⌊(p+1)/b⌋`. Multiply by `(-1)^(p+1)(1 + pB)^(-a) ≡ (-1)^(p+1)(1 - apB) (mod p^2)`:

```
Λ ≡ P(t) := (-1)^p T_(p,p)(t) (1 - a p B(t))   (mod p^2),   deg P <= ⌊pa/b⌋ + ⌊(p+1)/b⌋.
```

With `a = b - 1` this is `p - ⌈p/b⌉ + ⌊(p+1)/b⌋ <= p`. (This degree count is exactly what Beukers' argument for
`r < 1/3` cannot close when `r >= 1/3`. Here the power structure replaces `2⌊(p+1)r⌋` by `⌊(p+1)/b⌋`.) Since
`F(t^σ) ≡ 1 + g_1 t^p (mod t^(p+1))`, comparing the coefficients of `t^0, ..., t^p` in `F ≡ P F(t^σ)` gives
`P ≡ F_p + (g_p - g_1) t^p (mod p^2)`. Lemma 2 gives `g_p ≡ g_1`. ∎

**Lemma 2.** `g_p ≡ g_1 (mod p^2)` for `p ∤ N`.
*Proof.* Use `p | C(p,k)` for `0 < k < p` and `C(2p,p) ≡ 2 (mod p^2)` (Babbage, `p >= 3`). Then s7:
`g_p ≡ C(2p,p)^2 ≡ 4 = g_1`. s10: `g_p ≡ 2 = g_1`. s18: for `0 < k <= (p-1)/2` both `C(p,k)` and `C(2p-2k, p-k)` are
divisible by `p` (a carry in base `p`), so `g_p ≡ C(2p,p)(C(2p-1,p) + C(2p,p)) = (3/2) C(2p,p)^2 ≡ 6 = g_1`.
Apéry2 ⊙ central: `g_p ≡ C(2p,p)(1 + C(2p,p)) ≡ 6 = g_1`. Franel ⊙ central: `g_p ≡ 2 · 2 = 4 = g_1`. The prime `p = 2`
(s7 and Apéry2 ⊙ central only) is direct: `g_2 = 48, 114` and `g_1 = 4, 6`. ∎

**Corollary 2 (PROVED).** `T_(p,p)(t) + T_(p,p+1)(t) F_p(t) ≡ 0 (mod p^2)` in `Z[t]`.

Check F (all `p < 110`, coefficients up to `t^(3p)`): the congruence holds, and it fails modulo `p^3` except at `p = 3`
for s7 and s10, where it holds modulo 27 to the computed precision. It also holds at the primes `p | N` (s7: 7;
s10: 2, 5; s18: 2, 3; Apéry2 ⊙ central: 5; Franel ⊙ central: 2, 3), which Theorem 1 does not cover.

## 4. CM values at split primes

**Theorem 3 (PROVED, given Beukers' Proposition 5.1(ii) and the proof of his Theorem 1.16).** For the five families,
Beukers' Theorems 1.15 and 1.16 hold modulo `p^2`: under their hypotheses (`p ∤ N` split in `Q(α)`, `w_e ∈ Γ` of
determinant `e`, integers `c, d` with `N | c`, `e | d`, `ep = |cα + d|^2`, `cα + d` a `p`-adic unit, `t(α) ∈ Z_p`),

```
F_p(t(α)) ≡ χ(w_e) (cα + d)^2 / e   (mod p^2).
```

*Proof.* Beukers' proof of Theorem 1.16 shows that `ϑ = χ(w_e)(cα + d)^2/e` is a `p`-adic unit root of
`T_p(x, t(α))`. By Proposition 5.1(ii), `0 = T_p(ϑ, t(α)) ≡ ϑ^p (T_(p,p)(t(α)) + T_(p,p+1)(t(α)) ϑ) (mod p^3)`. Since
`T_(p,p+1)(t(α)) ≡ 1 (mod p)`, Corollary 2 gives `F_p(t(α)) ≡ ϑ (mod p^2)`. ∎

In terms of a representation `4ep = X^2 + |D| Y^2`, the right side is `±(X^2 - 2ep)/e (mod p^2)`. Section CM finds a
value of this shape at every split prime `5 <= p < 400` for all 58 rational CM values of the five families (for s7 this
is Corollary F of the s7 note).

## 5. Inert primes

**Theorem 4 (PROVED, given Shimura's reciprocity law).** Let `(Γ, t, Θ, F)` satisfy Beukers' Assumption (II), or be
one of the five families. Let `α ∈ H` be imaginary quadratic with `t(α) ∈ Q`, and let `D` be the discriminant of the
multiplier ring of the lattice `Z + Zα`. If `p ∤ N`, `(D/p) = -1` and `t(α) ∈ Z_p`, then

```
F_p(t(α)) ≡ 0   (mod p^2).
```

*Proof.* Put `K = Q(α)`, `L = Z + Zα`, and let `Θ|g = det(g)(cτ + d)^(-2) Θ(gτ)` be the weight-2 slash action.

*Step 1 (roots).* By the definition of `T_p`, the roots of `T_p(x, t(α))` are `x_γ(α) = pΘ(α)/(Θ|γ)(α)` for
`γ ∈ {δ = (p 0; 0 1), β_k = (1 k; 0 p)}`, the upper triangular representatives of `Γ0(N)\M(N,p)`. That is,
`Θ(α)/Θ(pα)` and `p^2 Θ(α)/Θ((α+k)/p)`. Their product is `(-1)^(p+1) p^(2p)/T_(p,p+1)(t(α))`, and
`T_(p,p+1)(t(α)) ≡ 1 (mod p)` is a unit. So the roots are finite and nonzero, and their valuations add up to `2p`.

*Step 2 (one Galois orbit).* `f = x_δ = Θ(τ)/Θ(pτ)` is a modular function for `Γ0(pN)` with rational Fourier
coefficients. For `σ ∈ Γ0(N)`, `f∘σ = pΘ|σ/(Θ|δσ) = pΘ/(Θ|γ)` where `δσ ∈ Γ0(N)γ`. As `σ` runs over
`Γ0(pN)\Γ0(N) ≅ P^1(F_p)`, this gives all `p+1` functions `x_γ`. Take an idele `s` with `s_ℓ = 1` for `ℓ ≠ p` and
`s_p ∈ O_(K,p)^*`. Since `p ∤ D`, `L ⊗ Z_p` is an `O_(K,p)`-module, so `u = q_α(s)` lies in `GL_2(Ẑ)`. It is trivial
away from `p` and is the matrix of multiplication by `s_p` at `p`. Shimura's reciprocity law gives
`f(α)^[s^(-1), K] = f^u(α)`, and since `u ≡ 1 (mod N)` and `f` has rational coefficients, `f^u` is one of the
`f∘σ`. As `s_p` runs over `O_(K,p)^*`, `u mod p` runs over a non-split Cartan subgroup `C ≅ F_(p^2)^*` of
`GL_2(F_p)` (`p` is inert). `C` meets every Borel subgroup only in the scalars, so it acts transitively on
`GL_2(F_p)/B ≅ P^1(F_p)`. Hence every root `x_γ(α)` is `x_δ(α)^[s^(-1), K]` for some such `s`.

*Step 3 (valuations).* Artin symbols of ideles supported at `p` with unit component lie in the inertia group at `p`
of `K^ab/K`. That group is abelian, so these symbols fix every prime `𝔓 | p` of `K^ab`. Therefore all `p+1` roots have
the same `𝔓`-adic valuation, namely `2p/(p+1)` (normalised by `v(p) = 1`). Then
`T_(p,p)(t(α)) = -T_(p,p+1)(t(α)) Σ_γ x_γ(α)` has valuation `>= 2p/(p+1) > 1`, and it lies in `Z_p`, so its
valuation is `>= 2`. Corollary 2 (or, under Assumption (II), the mod `p^3` congruence in Beukers' proof of Theorem 1.13)
gives `F_p(t(α)) ≡ -T_(p,p)(t(α))/T_(p,p+1)(t(α)) ≡ 0 (mod p^2)`. ∎

**The picture behind it (check N).** The Newton polygon of `T_p(x, t0)` was computed from the full polynomial (power
sums `Σ_k Θ((τ+k)/p)^j = p U_p(Θ^j)`, Newton's identities) for s7, s10, s18, both Hadamard products and Apéry at
24 CM prime pairs with `p <= 13`. It always has the shape that the inertia orbits on `P^1(F_p)` suggest:
- inert: one slope, all `p+1` roots of valuation `2p/(p+1)`, and `T_p(x, t0)` is irreducible over `Q_p`, or the
  square of an irreducible polynomial (when two lines give the same value);
- split: one unit root, one root of valuation 2 (the other eigenline) and an orbit of `p-1` roots of valuation 2;
  over `Q_p` the factors have degrees `1, 1, p-1`, or `1, 1` and the square of a factor of degree `(p-1)/2`;
- ramified: one root of valuation 1 and an orbit of `p` roots of valuation `(2p-1)/p`.

So at a ramified prime `F_p(t0)` should have valuation exactly 1. Section CM confirms `v_p(F_p(t0)) = 1` at all 38
ramified pairs. A proof would need the valuation of the root of the `O`-stable line, which Step 2 does not give.

Theorem 4 fails modulo `p^3`. For Apéry at `t0 = 1`, no inert `p < 400` has `F_p(1) ≡ 0 (mod p^3)`.

**Corollary 5 (Z.-W. Sun's conjecture for Apéry numbers, first proved by C. Wang and Z.-W. Sun; PROVED here again,
given Beukers' Theorem 1.15).** Let
`A_k = Σ_j C(k,j)^2 C(k+j,j)^2`. For every odd prime `p`,

```
Σ_(k<p) A_k ≡ 4x^2 - 2p (mod p^2)  if p ≡ 1, 3 (mod 8), p = x^2 + 2y^2;     ≡ 0 (mod p^2)  if p ≡ 5, 7 (mod 8).
```

*Proof.* Apéry is Beukers' 9.5 (Assumption (II), `N = 6`, `t((2+√-2)/6) = 1`, `D = -8`). Inert `p ≡ 5, 7 (mod 8)`:
Theorem 4. Split `p >= 5`: Theorem 1.15 with `c = 6y`, `d = x - 2y` gives `cα + d = x + y√-2` and
`F_p(1) ≡ (x + y√-2)^2 (mod p^3)`. Modulo `p^2` this is `(x + y√-2)^2 + (x - y√-2)^2 = 4x^2 - 2p`, because
`x - y√-2 = p/(x + y√-2)` is divisible by `p` in that embedding. `p = 3`: `1 + 5 + 73 = 79 ≡ -2 (mod 9)`. ∎

Likewise `Σ_(k<p) (-1)^k A_k ≡ 0 (mod p^2)` for `p ≡ 2 (mod 3)`, `p >= 5` (`t((3+√-3)/6) = -1`). That case was
already proved by Z.-W. Sun (J. Number Theory 132 (2012) 2673-2699, arXiv:1101.1946). Theorem 4 gives the inert half
for every rational CM value listed in Beukers' Appendix C. Some of those inert cases were known before by other
methods (both Apéry sums were); I have not surveyed which others.

## 6. The rational CM values, and corrections to Beukers' Appendix C

Section CM takes every rational CM value `t0 = t(α)` listed in Beukers' 9.1-9.3, 9.5, 9.7-9.15. It checks
`|t(α) - t0| < 10^-40` from the η-quotients and computes `D` from the primitive form of `α`. Then for
`5 <= p < 400` with `p ∤ N den(t0)` it tests:
- inert primes: `F_p(t0) ≡ 0 (mod p^2)`. True for all 4600 pairs.
- split primes: `F_p(t0) ≡ s(X^2 - 2ep)/e` for an Atkin-Lehner `e` and `4ep = X^2 + |D|Y^2`. True for all 4268 pairs.
- ramified primes: `v_p(F_p(t0)) = 1`. True for all 38 pairs.

The arXiv v3 tables contain these misprints. The values below are confirmed both numerically and by the congruences:
- 9.1: the operator is `θ^3 - 8t(2θ+1)^3`, not `4t`. Also `t((1+i)/2) = -1/8`, not `-1/2`.
- 9.2: `t((1+3i)/2) = -1/12288`, not `-1/17088`. `t(√(-11/2)) = 1/2509056`, not `1/9801`. Neither printed value is a
  CM value on the imaginary axis.
- 9.3: `t((1+√-11)/6) = 1/64`, not `1/1024`.
- 9.13: the binomial formula for s18 needs the sign `(-1)^k`. Without it `g_3 = 636`, while the operator gives `564`.
- 9.14: `t(√(-3/5)) = 1/147`, not `1/588`.
- 9.15: `t((1+√(-7/3))/2) = -1/112`, not `-1/116`.
- 9.7 and 9.8 (twisted): the printed values belong to the untwisted sequences. For the printed η-quotient `t(q)` and
  the printed (signed) `g_n` they read `t(α + 1/2) = -t0`.

## 7. Open problems

1. **Higher powers.** The Dwork-type congruences `F_(p^s)(t)/F_(p^(s-1))(t^σ) ≡ F(t)/F(t^σ) (mod p^(2s))` and, for
   s7, `c(np^r) ≡ c(np^(r-1)) (mod p^(2r))` hold numerically (s7 note, section 8) but are not proved. Lemma 1 only
   controls `T_p`, a level-`p` object.
2. **Primes dividing the level.** Theorem 1 holds numerically at every `p | N` (section 3). `T_p` is not defined there.
3. **Ramified primes.** `v_p(F_p(t0)) = 1` (section 5).
4. **Other families.** Lemma 1 needs only that `t^a Θ^b` be an η-product with `gcd(a, b) = 1`, and the degree count in
   Theorem 1 goes through whenever `a <= b - 1`, since then `⌊pa/b⌋ + ⌊(p+1)/b⌋ <= p`. Each new family also needs its
   own Lemma 2. I have not looked for further third-order families of this shape.

## 8. Reproduction

`gp -q 04-computation/experiments/beukers_remark114_power_structure_check_20261018.gp` (PARI/GP, 16 s). It runs
sections I, T, F, N, CM and S; the last line is `ALL (I, T, F, N, CM, S): 1`. Output:
`05-knowledge/results/beukers_remark114_power_structure_check_20261018.out`.
