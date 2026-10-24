# Cooper's s7 through its weight-one root: the Chan-Cooper-Sica conjecture and Beukers' mod p^2 congruence

2026-10-18, session `cloud-2026-10-18-M1` (temporary identity of a cloud thread working beside opus-2026-10-18-S43; no
`.machine-id`, not a registered agent).

**Typing.**
- PROVED here: Theorems A, B, D, E, part (ii) of Theorem C, and Corollary F.
  - Theorem D is Chan-Cooper-Sica's Conjecture 5.4 (Int. J. Number Theory 6 (2010) 89-97) together with O'Brien's
    Conjecture 9.3 (MSc thesis, Massey University, 2016), in both directions. Cooper lists it as open in
    arXiv:2302.00757v2, section 11 (May 2024). I found no later proof.
  - Theorem A is the congruence that Beukers (arXiv:2403.03301v3, Remark 1.14(2)) reports for s7 from computer
    experiments; his Theorem 1.13 does not cover s7.
  - Theorem C(ii) is the "only if" half of O'Brien's Conjecture 9.2 (Lucas congruences).
- CITED: Beukers, arXiv:2403.03301 (the polynomial `T_p`, Propositions 4.2 and 5.1, Lemmas 3.1, 3.6, 3.7, the proof of
  Theorem 1.16); Beukers-Tsai-Ye, Bull. LMS 57 (2025), arXiv:2408.16616, Theorem 1.3 (the "if" half of Theorem C);
  Cooper, Ramanujan J. 29 (2012) for `F(t(q)) = θ^2`; von Staudt-Clausen; Sturm's bound; the valence formula.
- PROVED by a check beyond a Sturm bound: the identities (S4), (S6), (S9), C7, D7, D2, E3 of the script.
- FINITE-EXACT: every proof ingredient for split `11 <= p <= 43` and inert `5 <= p <= 41`; Theorems A and B for
  `p < 110`; Theorems C, D, E for `p < 60`, `n < 3000/p`; the CM table for `5 <= p < 400`.
- VERIFIED (numerics only): A, D and E fail modulo `p^3`; `c(np^2) ≡ c(np) (mod p^4)` in a small range.
- OPEN: `c(np^r) ≡ c(np^(r-1)) (mod p^(2r))` for `r >= 2`; Chan-Cooper-Sica Conjecture 5.5 (level 11).
- PROVED in `beukers_remark114_power_structure_and_inert_cm_20261018.md`: the analogue of Theorem A for s10, s18 and
  the Hadamard products of Beukers' Remark 1.14(2), a second proof of Theorem A, and the inert-prime complement
  `F_p(t0) ≡ 0 (mod p^2)` of Corollary F.

## 1. Setting

```
s7(n) = Σ_k C(n,k)^2 C(n+k,k) C(2k,n),    F(t) = Σ s7(n) t^n = 1 + 4t + 48t^2 + 760t^3 + ...,
G(t) = F(t)^(1/2) = Σ c(n) t^n,    c = 1, 2, 22, 336, 6006, ...    (OEIS A279619),
(n+1)^2 c(n+1) = (26n^2 + 13n + 2) c(n) + 3(3n-1)(3n-2) c(n-1).
```

With `q = e^(2πiτ)` and `χ = (-7/·)`:

```
h = (η(7τ)/η(τ))^4,    t = h/(1 + 13h + 49h^2),    s = (1 - 49h^2)/(1 + 13h + 49h^2),
θ(τ) = Σ_(x,y) q^(x^2+xy+2y^2) = 1 + 2 Σ_n (Σ_(d|n) (d/7)) q^n    ∈ M_1(Γ0(7), χ).
```

- (S2), (S3): `F(t(τ)) = θ^2` (Cooper 2012; Beukers, Appendix C 9.11), hence `G(t(τ)) = θ`.
- (S4) `t θ^3 = (η(τ)η(7τ))^3`; (S5) `s^2 = (1 + t)(1 - 27t) = 1 - 26t - 27t^2`; (S6) `q dt/dq = t θ^2 s`.
- `h` is a Hauptmodul of `X0(7)`. The Fricke involution `W7` acts by `h -> 1/(49h)`, so `t -> t` and `s -> -s`, and
  `θ|W7 = -iθ` (the form `x^2 + xy + 2y^2` is equivalent to its adjoint). So `t` is a Hauptmodul of `Γ0(7)+`.
- The two poles of `t` on `X0(7)` are the elliptic points `e1, e2` of order 3 (the roots of `1 + 13h + 49h^2`), and
  `W7` swaps them. By the valence formula (index 8, weight 1) `θ` vanishes exactly at `e1, e2`, to order `1/3` each.
- (R) The functions on `X0(7)` with poles only at `e1, e2`, of order `<= m` at each, form
  `Q[t]_(<=m) ⊕ s·Q[t]_(<=m-1)`. `W7` acts by `+1` on the first summand and by `-1` on the second. An integral
  `q`-expansion gives integral `t`-coefficients (in the second case divide by `s ∈ 1 + tZ[[t]]`).

`q(t)` is the inverse series of `t(q)`, and `t^σ(t) := t(q(t)^p)` is the excellent Frobenius lift:
`t^σ ∈ t^p Z[[t]]`, `t^σ ≡ t^p (mod p)`, and `ε_p := (t^σ/t^p - 1)/p ∈ tZ[[t]]`. For a series `A`, `A^σ = A(t^σ)` and
`A_p` is its truncation below `t^p`; then `G^σ(t(τ)) = θ(pτ)`. `C` is the Cartier operator
`C(Σ a_n t^n) = Σ a_(pn) t^n`, so `C(A(t)·B(t^p)) = C(A)·B`. Finally `ϑ = t d/dt` and
`d_p := ⌊p/3⌋ + ⌊(p+1)/3⌋ <= p - 1`.

## 2. The weight-one polynomial T_p on Γ0(7)

For a prime `p ≠ 7` take Beukers' polynomial (5) for `Θ = θ`:

```
T_p(x) = (1 - θ(pτ)/θ(τ) · x) · Π_(k=0)^(p-1) (p - χ(p) θ((τ+k)/p)/θ(τ) · x) = Σ_(m=0)^(p+1) T_(p,m) x^m.
```

**Lemma T (PROVED).**
1. If `χ(p)^m = 1`, then `T_(p,m) ∈ Z[t]` has degree `<= ⌊m/3⌋`. If `χ(p)^m = -1`, then `T_(p,m) ∈ s·Z[t]` and
   `T_(p,m)/s` has degree `<= ⌊m/3⌋ - 1`.
2. `T_(p,m) ≡ 0 (mod p^2)` for `m < p`.
3. `T_(p,p+1) ≡ χ(p) (mod p)`.
4. `T_(p,p)(t) + T_(p,p+1)(t) · G(t)/G(t^σ) ≡ 0 (mod p^2)` in `Z_(p)[[t]]`.

*Proof.* The coefficients are modular functions on `Γ0(7)` with integral `q`-expansions, and 2 and 3 hold, exactly as
in Beukers' Proposition 4.2 and Lemma 3.6 (CITED). Part 2 rests on `Σ_k θ(τ)/θ((τ+k)/p) = θ·p·U_p(1/θ) ≡ 0 (mod p)`;
part 3 on `θ(pτ) ≡ θ^p` and `Π_k θ((τ+k)/p) ≡ θ (mod p)`. The coefficient of `x^m` is a sum of products of `m`
quotients `θ|γ/θ`, each holomorphic except for a pole of order `<= 1/3` at `e1` and at `e2`.

What is new is that `t` is a Hauptmodul of `Γ0(7)+`, not of `Γ0(7)`. Conjugation by `W7` sends `γ = [[a,b],[7c,d]]`
of determinant `p` to `[[d,-c],[-7b,a]]`, and `χ(a) = χ(p)χ(d)`. Since `θ|W7` is a multiple of `θ`, this gives
`T_p|W7 (x) = T_p(χ(p)x)`, that is `T_(p,m)|W7 = χ(p)^m T_(p,m)`, and (R) gives part 1. For part 4,
`T_p(θ(τ)/θ(pτ)) = 0`; reduce with part 2 and divide by the unit `(θ(τ)/θ(pτ))^p`. ∎

Section T of the script checks parts 1-4 for every `p <= 43`; all the degree bounds are attained.

## 3. Theorems B and A

**Theorem B (PROVED).** Let `p ≠ 7`.
- If `χ(p) = 1` (this includes `p = 2`): `G ≡ G_p(t) · G(t^σ) (mod p^2)`, and `G_p` has degree `<= d_p` modulo `p^2`.
- If `χ(p) = -1`: `H := G/s` satisfies `H ≡ H_p(t) · G(t^σ) (mod p^2)`, and `H_p` has degree `<= d_p - 1` modulo
  `p^2`. Equivalently `G ≡ s · H_p · G(t^σ)`.

*Proof.* By parts 1 and 3, `T_(p,p+1) = χ(p) + pB` with `B ∈ Z[t]` of degree `<= ⌊(p+1)/3⌋`, and
`(χ(p) + pB)(χ(p) - pB) ≡ 1 (mod p^2)`. Part 4 gives `G/G^σ ≡ -(χ(p) - pB) T_(p,p) (mod p^2)`. If `χ(p) = 1` the
right side is a polynomial of degree `<= d_p <= p - 1`, and since `G^σ ≡ 1 (mod t^p)` it is `G_p`. If `χ(p) = -1`,
then `T_(p,p) = s·α` with `deg α <= ⌊p/3⌋ - 1`, so `H/G^σ ≡ (1 + pB)α`, which has degree `<= d_p - 1` and is `H_p`. ∎

**Theorem A (PROVED).** For every prime `p ≠ 7`: `F(t) ≡ F_p(t) · F(t^σ) (mod p^2)`.

*Proof.* `F/F^σ = (G/G^σ)^2 ≡ T_(p,p)^2 · (1 - 2χ(p)pB) (mod p^2)`. `T_(p,p)^2` is a polynomial in `t` of degree
`<= 2⌊p/3⌋` (in the inert case `s^2 = 1 - 26t - 27t^2`), so `F/F^σ` is congruent to a polynomial of degree
`<= 2⌊p/3⌋ + ⌊(p+1)/3⌋`. This is `<= p - 1` for `p = 2` and every `p >= 5`. At `p = 3`, `T_(3,3) = α0·s` with
`α0 ∈ Z` and `s^2 ≡ 1 - 26t (mod 27)`, so the degree is `<= 2`. In every case the polynomial is `F_p`. ∎

Beukers' Theorem 1.13 needs `Θ = θ^2` to be an η-product or to vanish to order `< 1/3`; here the order is `2/3`.
Working with `θ` on the double cover `X0(7) -> X0(7)+` is what makes the degree count close. Straub (Monatsh. Math.
2024) proves Gessel-Lucas congruences modulo `p^2` for s7, which use the naive lift `t -> t^p` and a formal
derivative; that is a different statement. For s7 itself, `s7(pn) ≡ s7(n) (mod p^3)` (`p >= 5`) is
Osburn-Sahu-Straub (2016).

Numerics (section A/B, `p < 110`, coefficients up to `t^(4p)`): A and B hold modulo `p^2`. A fails modulo `p^3` for
every `p < 110` other than 3 (at `p = 3` it holds modulo 27 for the coefficients computed). The naive lift `t -> t^p`
gives only a congruence modulo `p` for `p >= 5`. The inert degree bound of B is attained.

## 4. Lucas congruences (Theorem C)

**Theorem C.** `c(n)` satisfies the Lucas congruences modulo `p` if and only if `p = 7` or `p ≡ 1, 2, 4 (mod 7)`
(O'Brien, Conjecture 9.2).

(i) *If.* This is Beukers-Tsai-Ye's Theorem 1.3 for `Γ0(7)+`, `Z = θ`, `χ(w7) = -1` (CITED). It also follows here. For
split `p`, reduce Theorem B modulo `p` with `t^σ ≡ t^p`: `G ≡ G_p(t) G(t^p) (mod p)`. For `p = 7`,
`θ^6 (1+t)^2 = θ^6 + 2(η(τ)η(7τ))^3 θ^3 + (η(τ)η(7τ))^6` lies in `M_6(Γ0(7))` and is `≡ 1 ≡ E_6 (mod 7)` (Sturm bound
4). So `θ = θ^7/θ^6 ≡ θ(7τ)(1+t)^2`, that is `G ≡ (1+t)^2 G(t^7) (mod 7)`.

(ii) *Only if* (PROVED). Let `p` be inert and suppose the Lucas congruences hold. Then `G ≡ G_p G(t^p)`, and by
Theorem B `G ≡ s H_p G(t^p) (mod p)`, so `s ≡ G_p/H_p` in `F_p[[t]]`. But `s^2 = (1+t)(1-27t)` is not a square in
`F_p(t)`: it has two distinct roots for `p ≠ 2, 3, 7`, and it is `1 + t` for `p = 3`. ∎ The first failure is at
`n = p` for every inert `p < 60`.

## 5. The Chan-Cooper-Sica conjecture (Theorem D)

**Theorem D (PROVED).** For a prime `p` the following are equivalent:
- (a) `c(pn) ≡ c(n) (mod p^2)` for every `n >= 0`;
- (b) `p = 7` or `p ≡ 1, 2, 4 (mod 7)`.

For `p ≡ 3, 5, 6 (mod 7)` the congruence fails already modulo `p`.

**5.1 Reduction.** Let `p` be split. Since `(1 + pε)^k ≡ 1 + kpε (mod p^2)`,
`G(t^σ) = G(t^p (1 + pε_p)) ≡ G(t^p) + pε_p · (ϑG)(t^p)`. Theorem B and `C(G_p) = c(0) = 1` give

```
C(G) ≡ G + p · C(ε_p G_p) · ϑG    (mod p^2).
```

So (a) holds as soon as `ε_p G_p (mod p)` is a polynomial of degree `< p` with constant term 0.

**5.2 Lemma E.** For `p >= 5`, `Φ_p := Π_n (1 - q^(pn)) / Π_n (1 - q^n)^p` satisfies
`Φ_p^2 · E_(p-1) ≡ 1 (mod p^2)`.

*Proof.* From `-log Π(1 - q^n) = Σ_M σ(M)/M · q^M` one gets `log Φ_p = p Σ_M (Σ_(d|M, p∤d) 1/d) q^M ∈ pZ_p[[q]]`.
By von Staudt-Clausen `pB_(p-1) ≡ -1 (mod p)`, so `E_(p-1) ≡ 1 - 2p Σ_M σ_(p-2)(M) q^M (mod p^2)`, and
`σ_(p-2)(M) ≡ Σ_(d|M, p∤d) 1/d (mod p)`. Hence `2 log Φ_p ≡ -(E_(p-1) - 1) ≡ -log E_(p-1) (mod p^2)`. ∎

**5.3 Proposition F.** For split `p >= 5`, as `q`-series,

```
p ε_p ≡ 3 [ θ^(p-1) G_p(t) - (E_(p-1)(τ) + E_(p-1)(7τ))/2 ]    (mod p^2).
```

*Proof.* By (S4), `t^σ/t^p = t(pτ)/t(τ)^p = Φ_p(q)^3 Φ_p(q^7)^3 (θ(τ)^p/θ(pτ))^3`. Lemma E handles the first two
factors. Theorem B gives `θ(pτ) = G^σ ≡ θ/G_p`, so the third factor is `(θ^(p-1) G_p)^3`. Each of `E_(p-1)(τ)`,
`E_(p-1)(7τ)` and `θ^(p-1) G_p` is `≡ 1 (mod p)`, and `(1 + pY)^e ≡ 1 + epY (mod p^2)` for `p`-integral `e`. ∎

**5.4 Lemma W.** Let `k = p - 1`, `Ψ_p := (E_k(τ) + E_k(7τ))/(2θ^k)` and `λ := (-7)^(k/2) ≡ χ(p) (mod p)`.
- If `χ(p) = 1`: `Ψ_p ≡ A_p(t) (mod p^2)` with `A_p ∈ Z_(p)[t]` of degree `<= ⌊(p-1)/3⌋`.
- If `χ(p) = -1`: `Ψ_p ≡ s · B_p(t) (mod p^2)` with `B_p ∈ Z_(p)[t]` of degree `<= ⌊(p-1)/3⌋ - 1`.

In both cases `A_p(0) ≡ B_p(0) ≡ 1 (mod p^2)`, because `Ψ_p(i∞) = s(0) = 1`.

*Proof.* `Ψ_p` is a modular function on `Γ0(7)` (weight `k`, trivial character), with poles only at `e1, e2`, of order
`<= ⌊k/3⌋`. From `E_k|W7 = 7^(k/2) E_k(7τ)`, `E_k(7τ)|W7 = 7^(-k/2) E_k(τ)` and `θ^k|W7 = (-1)^(k/2) θ^k`:

```
Ψ_p - Ψ_p|W7 = (1 - λ)(E_k(7τ) - λ^(-1) E_k(τ)) / (2θ^k),
Ψ_p + Ψ_p|W7 = (1 + λ)(E_k(7τ) + λ^(-1) E_k(τ)) / (2θ^k).
```

The `q`-expansions are `p`-integral and `E_k ≡ 1 (mod p)`. If `λ ≡ 1`, the difference is `≡ 0 (mod p^2)`, so `Ψ_p` is
congruent to `(Ψ_p + Ψ_p|W7)/2`, which is `W7`-invariant and by (R) a polynomial `A_p(t)`. If `λ ≡ -1`, the sum is
`≡ 0 (mod p^2)`, and `Ψ_p ≡ (Ψ_p - Ψ_p|W7)/2 = s · B_p(t)`. ∎

**5.5 Split p >= 11.** By F and W, `p ε_p ≡ 3 θ^(p-1) (G_p - A_p) (mod p^2)`. Modulo `p` this forces `G_p ≡ A_p`, and
Theorem B gives `θ^(p-1) G_p ≡ θ^p/θ(pτ) ≡ 1`. Hence

```
ε_p G_p ≡ 3 (G_p - A_p)/p    (mod p),
```

a polynomial of degree `<= d_p < p` with constant term `(1 - A_p(0))/p ≡ 0`. Now 5.1 applies. ∎

In passing, `G_p ≡ A_p (mod p)` means `c(n) ≡ 0 (mod p)` for `⌊(p-1)/3⌋ < n < p`, and Theorem B means
`c(n) ≡ 0 (mod p^2)` for `d_p < n < p` (both checked for the split `p < 360`; the first bound is sharp there).

**5.6 p = 2.** `θ^9 θ(2τ)^3 (t(2τ) - t^2 - 2t^3)` lies in `M_12(Γ0(14))` with integral coefficients and is `≡ 0 (mod 4)`
beyond the Sturm bound 24. So `t^σ ≡ t^2 (1 + 2t) (mod 4)` and `ε_2 ≡ t (mod 2)`. With `G_2 = 1 + 2t`,
`ε_2 G_2 ≡ t` and `C(t) = 0`, and 5.1 (which only used Theorem B) gives `c(2n) ≡ c(n) (mod 4)`. ∎

**5.7 p = 7.** Here `7 | N` and there is no `T_7`. Use the weight-three form `f := θ^3 s` and
`c(n) = CT_q(f · t^(-n))`, which follows from `G dt/t = f dq/q` (S6).
- `f ∈ M_3(Γ0(7), χ)` (`θ^3` vanishes simply where `s` has simple poles), and `f = 1 - (7/8) A1 - (49/8) A2` with
  `A1 = Σ_n (Σ_(d|n) χ(d) d^2) q^n` and `A2 = Σ_n (Σ_(d|n) χ(n/d) d^2) q^n` (Sturm bound 2).
- From the divisor sums, `U_7 A1 = A1` and `U_7 A2 = 49 A2`, so `U_7 f - f = -294 A2 ≡ 0 (mod 49)`.
- `θ^27 θ(7τ)^3 (t(7τ) - t^7 - 14t^8 - 28t^9)` lies in `M_30(Γ0(49))` and is `≡ 0 (mod 49)` beyond the Sturm bound
  140. So `t(7τ) ≡ t^7 (1 + 14t + 28t^2)` and `ε_7 ≡ 2t + 4t^2 (mod 7)`.

Then `t^(-7n) = t(7τ)^(-n) (1 + 7ε_7)^n ≡ t(7τ)^(-n) (1 + 7n ε_7) (mod 49)`, and

```
c(7n) ≡ CT(f · t(7τ)^(-n)) + 7n · CT(f ε_7 · t(7τ)^(-n))    (mod 49).
```

The first term is `CT(U_7 f · t^(-n)) ≡ CT(f · t^(-n)) = c(n) (mod 49)`. In the second, modulo 7,
`t(7τ)^(-n) ≡ t^(-7n)` and `CT(f ε_7 t^(-7n)) = [t^(7n)](ε_7 G)`. By Theorem C(i),
`ε_7 G ≡ (2t + 4t^2)(1+t)^2 G(t^7) = (2t + t^2 + 3t^3 + 4t^4) G(t^7) (mod 7)`, which has no `t^(7n)` terms. ∎

**5.8 Inert p.** Suppose `c(pn) ≡ c(n) (mod p)` for all `n`, so `C(G) ≡ G`. By Theorem B `G ≡ s H_p G(t^p)`, so
`C(s H_p) ≡ 1 (mod p)`. In `F_p[[t]]`, with `Q = s^2 = (1+t)(1-27t)`, we have `s = s^(1-p) s^p = Q^(-(p-1)/2) s(t^p)`
and `Q^(-(p-1)/2) = Q^((p+1)/2)/Q(t^p)`. Hence `C(s H_p) = s · P/Q` with `P = C(H_p Q^((p+1)/2)) ∈ F_p[t]`, and
`s = Q/P ∈ F_p(t)`, which is impossible as in C(ii). ∎

**5.9 Numerics.** Section D of the script checks (E), (F), (W) and (P) for split `11 <= p <= 43`, and the identities
used for `p = 2` and `p = 7`. Section O checks Theorem D for `p < 60`, `n < 3000/p`. The congruence holds modulo `p^2`
exactly for `p ∈ {2, 7, 11, 23, 29, 37, 43, 53}`, fails modulo `p` for every inert `p`, and fails modulo `p^3` for
every `p`.

## 6. Inert primes: the twisted congruence (Theorem E)

Let `1/s = Σ u(n) t^n`, where `u(n) = [x^n] (x^2 + 13x + 49)^n` are generalized central trinomial coefficients
(`(n+1) u(n+1) = (26n + 13) u(n) + 27n u(n-1)`), and let `H = G/s = Σ d(n) t^n`, so `d(n) = Σ_k c(k) u(n-k)`. In
modular terms `d(n) = CT_q(θ^3 t^(-n))` and `c(n) = CT_q(θ^3 s t^(-n))`.

**Theorem E (PROVED).** For every prime `p ≡ 3, 5, 6 (mod 7)` and every `n >= 0`, `d(pn) ≡ c(n) (mod p^2)`.

*Proof.* With Theorem B for inert `p`, the reduction 5.1 becomes `C(H) ≡ G + p · C(ε_p H_p) · ϑG (mod p^2)`. For
`p >= 5`, the proof of Proposition F with `θ(pτ) ≡ θ/(s H_p)` gives
`p ε_p ≡ 3 [θ^(p-1) s H_p - (E_(p-1)(τ) + E_(p-1)(7τ))/2]`, and Lemma W gives `Ψ_p ≡ s B_p`. As in 5.5,
`ε_p H_p ≡ 3 (H_p - B_p)/p (mod p)`, a polynomial of degree `<= d_p - 1` with constant term 0. For `p = 3`,
`θ^9 θ(3τ)^3 (t(3τ) - t^3)` lies in `M_12(Γ0(21))` and is `≡ 0 (mod 9)` beyond the Sturm bound 32, so
`ε_3 ≡ 0 (mod 3)`. ∎

So at the split primes and at 7, `c(pn) ≡ c(n)`, and at the inert primes, `d(pn) ≡ c(n)`, in both cases modulo `p^2`
and not modulo `p^3`. Section E checks every ingredient for inert `5 <= p <= 41`, and Theorem E for `p < 60`,
`n < 3000/p`. It also shows that `d(pn) ≡ c(n)` fails modulo `p` at the split primes and modulo 49 at `p = 7`.

## 7. CM values (Corollary F)

**Corollary F (PROVED, given Beukers' Proposition 5.1(i)-(iii) and Lemma 3.7).** Beukers' Theorems 1.15 and 1.16 hold
for s7 modulo `p^2`. Under their hypotheses for `Γ = Γ0(7)+` (`e ∈ {1, 7}`, `χ(w_1) = 1`, `χ(w_7) = -1`),
`F_p(t(α)) ≡ χ(w_e) (cα + d)^2 / e (mod p^2)`.

*Proof.* Proposition 5.1(i)-(iii) for `Θ = θ^2` do not use Assumption (II)(4). By (ii) and Theorem A,
`T_(p,p) + T_(p,p+1) F_p ≡ 0 (mod p^2)` as polynomials. The argument of Corollary 5.3 then gives
`F_p(t(α)) ≡ ϑ_p (mod p^2)` for every `p`-adic unit root `ϑ_p` of `T_p(x, t(α))`, and the proof of Theorem 1.16 shows
that `χ(w_e)(cα + d)^2/e` is one. ∎

At Beukers' nine rational CM values the following rule holds for every split prime `5 <= p < 400` (FINITE-EXACT,
section CM). If `4p = X^2 - DY^2` is solvable, then `F_p(t0) ≡ X^2 - 2p (mod p^2)`; otherwise `28p = X^2 - DY^2` and
`F_p(t0) ≡ 2p - X^2/7`. The two branches are the cases `N(cα + d) = p` and `N(cα + d) = 7p` of Corollary F. I have
not written out `α, c, d` for every entry, so the table itself stays FINITE-EXACT. For example, modulo `p^2`:

```
t0 = -1     (D = -7):    Σ_(n<p) (-1)^n s7(n) ≡ x^2 - 2p,    4p = x^2 + 7y^2   (p ≡ 1, 2, 4 mod 7)
t0 = 1/27   (D = -28):   Σ_(n<p) s7(n)/27^n  ≡ 4u^2 - 2p,    p = u^2 + 7v^2
t0 = 1/2    (D = -12):   Σ_(n<p) s7(n)/2^n   ≡ 4u^2 - 2p,    p = u^2 + 3v^2
t0 = 1/8    (D = -19):   Σ_(n<p) s7(n)/8^n   ≡ x^2 - 2p,     4p = x^2 + 19y^2
t0 = 1/24   (D = -27):   Σ_(n<p) s7(n)/24^n  ≡ x^2 - 2p,     4p = x^2 + 27y^2
```

The values `-1/8, -1/64, 1/125, -1/10648` (`D = -35, -91, -112, -427`, class number 2) use both branches.

## 8. Open problems

1. **Higher powers.** `c(np^2) ≡ c(np) (mod p^4)` holds for `p ∈ {2, 7, 11, 23, 29, 37, 43, 53}` whenever
   `np^2 <= 6000`, and at `p = 2, 7` even modulo `p^5` (VERIFIED). The natural conjecture is
   `c(np^r) ≡ c(np^(r-1)) (mod p^(2r))`, the weight-one analogue of Straub's Theorem 1.1 and of Beukers-Vlasenko's
   Conjecture 7.5 (Dwork crystals III). Iterating section 5 would need Theorem B for `G_(p^r)/G_(p^(r-1))(t^σ)` modulo
   `p^(2r)`.
2. **Resolved in Session53 (2026-10-24).** Theorem A modulo 27 at `p = 3`, for the excellent Frobenius lift: [proof](session53_g27_collatz_cooper_20261024/frobenius-mod27/proof.md) and [independent cubic audit](session53_g27_collatz_cooper_20261024/audit-cubic-lift/REPORT.md).
3. *Resolved.* Cooper's s10 and s18 and the Hadamard products of Beukers' Remark 1.14(2) are proved in
   `beukers_remark114_power_structure_and_inert_cm_20261018.md`, which also proves `F_p(t0) ≡ 0 (mod p^2)` at inert
   primes.
4. **Level 11** (Chan-Cooper-Sica Conjecture 5.5, primes with `(p/11) = 1`). There the weight-one theta series has a
   single simple zero, at a fixed point of `W11` on the genus-one curve `X0(11)`, and the count of Theorem B gives a
   polynomial of degree `<= p`, one more than needed.

## 9. Reproduction

`gp -q 04-computation/experiments/cooper_s7_weight1_frobenius_check_20261018.gp` (PARI/GP, about 3.5 minutes). It runs
sections S, T, A/B, C, CM, O, D and E; every check prints 1, and the output is in
`05-knowledge/results/cooper_s7_weight1_frobenius_check_20261018.out`.
