# HYP-9310 is a published theorem: the CM harmonic log law from reflection, Morley and the 2-isogeny

2026-10-18, session `cloud-2026-10-18-M1` (temporary identity of a cloud thread working beside opus-2026-10-18-S43; no
`.machine-id`, not a registered agent).

**Typing.**
- CITED: Kibelbek-Long-Moss-Sheller-Yuan, *Supercongruences and complex multiplication*, J. Number Theory 164 (2016)
  166-178, arXiv:1210.4489 (KLMSY below). Morley's congruence (1895).
- PROVED here (given KLMSY Theorem 2 and Morley): Proposition E (an equivalence valid for every `m`) and the theorem of
  section 3 (the law at the six non-self-dual CM points, every prime `p >= 5`, `p ∤ m`, except `p = 7` for
  `m = 1, 4096`).
- FINITE-EXACT: every proof step for all primes `5 <= p < 600`, the equivalence for every `1 <= |m| <= 64`, and the
  exceptional prime `p = 7` (`hyp9310_reflection_proof_check_20261018.py`, output `.out`).
- Not addressed: part (K) of HYP-9310 (the class-6 congruence). It stays OPEN.

## 1. What HYP-9310 says, and that it is known

Let `p >= 5`, `n = (p-1)/2`, `F_m = Σ_(k<=n) C(2k,k)^3 m^(-k)`, `h_m = Σ_(k<=n) C(2k,k)^3 m^(-k) (H_(2k) - H_k)` and
`q = q_p(m) = (m^(p-1) - 1)/p`. HYP-9310 (opus-2026-10-18-S43) asks for

```
(L_p)   6 h_m ≡ q_p(m) F_m   (mod p)
```

at Sun's CM values `m ∈ {1, 4096, -8, -512, 16, 256}` (the self-dual `m = ±64` are THM-4685).

**This is KLMSY's Corollary 3.** In their normalisation `Σ_k C(2k,k)^3 (λ/64)^k`, so `λ = 64/m`. Their statement: if
`E_λ` has CM, then for almost all primes `p` with `λ ∈ Z_p`,

```
Σ_(i<=(p-1)/2) C(2i,i)^3 (λ/64)^i [ 6(H_(2i) - H_i) + ((λ/64)^(p-1) - 1)/p ] ≡ 0   (mod p).
```

With `λ/64 = 1/m` and `(m^(-(p-1)) - 1)/p ≡ -q_p(m) (mod p)` this is `(L_p)`.
- **Their Corollary 4** is the case `λ = 64` (`m = 1`), stated for **all** primes `p > 3`:
  `Σ_(i=1)^((p-1)/2) C(2i,i)^3 Σ_(j=1)^i 1/(i+j) ≡ 0 (mod p)`. Since `Σ_(j<=i) 1/(i+j) = H_(2i) - H_i`, this is
  S43's `(H)`.
- **Their Lemma 10** is S43's Lucas-block lemma (`a_(k+pn) ≡ a_k a_n + p b_k Σ_(i<=n) C(2i,i)^3 (λ/64)^i (mod p^2)`,
  with `b_k` the harmonic sum of `(L_p)`), and their section 4 says the same thing S43 found: the supercongruence for
  `a_n` is equivalent to the auxiliary harmonic congruence.
- KLMSY prove Corollary 3 from their Theorem 2 (below) through Coster-van Hamme's theorem on Legendre polynomials.
  They write that they "know of no direct proof of Corollary 3".

So `(L_p)`, the census classes 3 and 7 that S43's THM-4685 reduces to `(H)`, and the Lucas-block form are known
mathematics. Novelty claims for them should be withdrawn. Part (K) is not covered by KLMSY.

## 2. A short independent proof

The proof below uses only KLMSY's Theorem 2, a classical congruence and two exact identities. It explains why the law
holds exactly at CM points: the reflection `k -> n-k` is the 2-isogeny, and `q_p(m)` is the `p`-adic square root that
the reflection produces.

**Theorem 2 of KLMSY (CITED).** Let `λ ≠ 1` be such that `E_λ : y^2 = (x-1)(x^2 - 1/(1-λ))` has CM, and let `p` be odd
with `E_λ` of good reduction over `Z_p`. Then
`Σ_(k<=n) C(2k,k)^3 (λ/64)^k ≡ ε_λ α_(p,λ)^2 (mod p^2)`, where `α_(p,λ)` is the unit root of
`X^2 - a_p(E_λ) X + p` (and `0` if `E_λ` is supersingular at `p`), and `ε_λ` is a Legendre symbol.
- In the arXiv v3 text the symbol prints as `((λ-1)/p)`. The data force `ε_λ = ((1-λ)/p)` (check K2: e.g. `m = 1`,
  `p = 11` gives `F_1 ≡ 4·2^2 - 2·11` while `(63/11) = -1`). The argument below only uses the ratio
  `ε_(1/λ)/ε_λ`, which is `(-λ/p)` under either reading.

**Lemma R (exact reflection).** For every `m` with `p ∤ m`,

```
Σ_(k<=n) [(-n)_k/(1/2 - n)_k]^3 (64/m)^k  =  C(2n,n)^(-3) (4096/m)^n F_(4096/m).
```

*Proof.* `(-n)_k/(1/2-n)_k = 4^k C(2n-2k, n-k)/C(2n, n)` (use `Γ(j+1/2) = (2j)! sqrt(π)/(4^j j!)`); substitute `j = n-k`. ∎

**Lemma T (first order in `ε`).** Put `F^(ε) = Σ_(k<=n) [(1/2+ε)_k/(1+ε)_k]^3 (64/m)^k`. Then
`d/dε log[(1/2+ε)_k/(1+ε)_k]^3 at ε = 0` is `3(2H_(2k) - H_k - H_k) = 6(H_(2k) - H_k)`, and every factor `1/2 + j`,
`1 + j` (`j < k <= n`) is a `p`-adic unit. Since `1/2 - p/2 = -n` and `1 - p/2 = 1/2 - n`, the left side of Lemma R is
`F^(-p/2)`, and

```
F^(-p/2) ≡ F_m - (p/2) · 6 h_m = F_m - 3p h_m   (mod p^2).
```

**Lemma M (Morley, CITED).** `C(2n, n) = C(p-1, (p-1)/2) ≡ (-1)^n 4^(p-1) (mod p^3)` for `p >= 5`. Hence
`C(2n,n)^(-3) 4096^n ≡ (-1)^n (mod p^3)`, because `4096^n = 4^(3(p-1))`.

**Lemma S.** `m^(-n) ≡ (m/p)(1 - (p/2) q_p(m)) (mod p^2)`: square both sides.

**Proposition E (PROVED, every `m` with `p ∤ m`, every `p >= 5`).**

```
(L_p) at m   <=>   F_(4096/m) ≡ (-m/p) F_m   (mod p^2).
```

*Proof.* Lemmas R, T, M, S give `F_m - 3p h_m ≡ (-1)^n (m/p)(1 - (p/2)q) F_(4096/m) (mod p^2)`. Since `(-1)^n = (-1/p)`,
the right side is `(1 - (p/2)q)·[(-m/p) F_(4096/m)]`. Multiply by the unit `(1 - (p/2)q)^(-1)`:
`F_m - 3p h_m ≡ (1 - (p/2) q) F_m` holds iff `F_(4096/m) ≡ (-m/p) F_m`, and the former is `6h_m ≡ q F_m (mod p)`. ∎

So the harmonic law at `m` is the statement that the truncated sums at `m` and `4096/m`, i.e. at `λ` and `1/λ`,
differ by the quadratic character `(-m/p)` modulo `p^2`. The law at `m` and at `4096/m` are visibly the same
statement (THM-4685's invariance). At generic `m` both sides of Proposition E fail together (check E, all `|m| <= 64`).

**Lemma I (2-isogeny with a twist; PROVED).** `E_λ` is 2-isogenous to the quadratic twist by `-2` of `E_(1/λ)`, so
`a_p(E_λ) = (-2/p) a_p(E_(1/λ))` for every odd prime of good reduction.
*Proof.* Translate `x -> x+1`: `E_λ : y^2 = x(x^2 + 2x + 1 - c)` with `c = 1/(1-λ)`. The 2-isogeny with kernel
`(0,0)` maps `y^2 = x(x^2 + Ax + B)` to `y^2 = x(x^2 - 2Ax + A^2 - 4B)`, here `y^2 = x(x^2 - 4x + 4c)`. Substituting
`x = -2u` gives `y^2 = -8u(u^2 + 2u + c)`. Since `1 - c' = c` for `c' = 1/(1 - 1/λ)`, the curve
`u(u^2 + 2u + c)` is `E_(1/λ)` translated. ∎ (Also checked by point counts, every `p < 600`.)

**Theorem (PROVED, given KLMSY Theorem 2).** `(L_p)` holds at `m = 1, 4096, -8, -512, 16, 256` for every prime
`p >= 5` with `p ∤ m`, except possibly `p = 7` for `m ∈ {1, 4096}`.

*Proof.* The curves `E_(64/m)` have conductors `2^6 3^2 7^2`, `2^4 3^2 7^2`, `2^6 3^2`, `2^5 3^2`, `2^6 3^2`,
`2^4 3^2` (PARI `ellglobalred`), with `j = -15^3, 255^3, 1728, 66^3, 0, 54000`. So KLMSY Theorem 2 applies at
`λ = 64/m` and at `1/λ = 64/(4096/m)` for every `p >= 5` except `p = 7` at the first pair. By Lemma I the two unit roots
agree up to sign, so `α_(p,λ)^2 = α_(p,1/λ)^2`, and supersingularity is shared. Hence
`F_(4096/m) ≡ (ε_(1/λ)/ε_λ) F_m = (-λ/p) F_m = (-m/p) F_m (mod p^2)`; if `p | (64 - m)`, both sides are `0`.
Proposition E finishes. ∎

**The exceptional prime.** `p = 7`, `m ∈ {1, 4096}`: `(L_7)` holds by direct computation (check L). For `m = 1` it is
also KLMSY Corollary 4.

**Why only CM points.** Proposition E is unconditional. KLMSY's Theorem 1 gives the truncated sums modulo `p` for every
ordinary `λ` (formal groups), so the relation `F_(1/λ) ≡ (-λ/p) F_λ` always holds modulo `p`. Modulo `p^2` it needs
the CM supercongruence at both `λ` and `1/λ`. This is why S43's `p`-adic CM detector finds exactly the CM values.

## 3. Consequences for S43's census

- Classes 3 and 7 of the half-period twin census (THM-4685 reduces both to `(H)`) are theorems, with `(H)` CITED from
  KLMSY Corollary 4 (or proved by the argument above).
- `Σ_j C(2j,j) C(p-1-j,j)^2 m^(-j) ≡ (1 - p q_p(m)/3) F_m (mod p^2)` at all eight CM points (S43's twin form) follows.
- Part (K) remains OPEN.

## 4. Reproduction

`python3 04-computation/experiments/hyp9310_reflection_proof_check_20261018.py 600 64` (70 s). Checks: R (exact
identity, `p < 40`, nine values of `m`), T, M, S, E (all `1 <= |m| <= 64`), K2 (KLMSY Theorem 2 with
`ε_λ = ((1-λ)/p)` at the six points), I (point counts), L (the law, including `p = 7`). All `True` for every prime
`5 <= p < 600`; output in `05-knowledge/results/hyp9310_reflection_proof_check_20261018.out`.
