---
id: THM-4686
title: "The CM harmonic law for sum C(2k,k)^3 m^-k holds at all eight of Sun's CM values (and hence KLMSY Conjecture 9 at s = 1 for every odd m at these points): 6 sum_{k<=(p-1)/2} C(2k,k)^3 m^-k (H_2k - H_k) == q_p(m) sum_{k<=(p-1)/2} C(2k,k)^3 m^-k (mod p) for m in {1, -8, 16, -64, 64, 256, -512, 4096} and every prime p >= 5 not dividing m. It follows from the mod-p^2 CM supercongruences at both ends of each reflection pair m <-> 4096/m (Kibelbek-Long-Moss-Sheller-Yuan, Theorem 2; van Hamme / Ishikawa / Ahlgren at m = 64) and the elementary mod-p^2 reflection F(z) == (-1)^n z^n [F(1/z) + 3p G(1/z)]. Hence also: the census twins sum_j C(2j,j) C(p-1-j,j)^2 == sum_j C(2j,j)^3 (mod p^2) (class 3) and the class-7 twin; the twin form T_m == (1 - p q_p(m)/3) F_m and the Lucas-block factorization at every CM point"
status: >
  PROVED, given CITED published theorems.
  - Input: the mod-p^2 supercongruences at both points of each pair. These are Kibelbek-Long-Moss-Sheller-Yuan,
    "Supercongruences and complex multiplication", J. Number Theory 164 (2016), arXiv:1210.4489, Theorem 2, for every
    CM lambda != 1: lambda = 64/m in {64, -8, 4, -1, 1/4, -1/8, 1/64}. At lambda = 1 (m = 64) the input is van Hamme
    1986 / Ishikawa 1990 / Ahlgren 2001.
  - Mechanism: THM-4685's Lemma 1 (S6) and Lemma 3, together with H_((p-1)/2) == -2 q_p(2).
  - FINITE-EXACT controls:
    - the reflection (R) mod p^2 and the derived law, for p <= 200 (cm_law_from_known_supercongruences_20261018.py);
    - the law itself at all eight points for p < 2000, and at m = 1, -8, 16 for p < 20000 (HYP-9310's scripts).
  - This upgrades HYP-9310's main statement from OPEN to PROVED. As with THM-4670 (S43 audit M1), the law is a
    corollary of two known supercongruences joined by an elementary identity. The contribution is the identity and the
    observation that it encodes the classical pair.
  - Still open in HYP-9310 / HYP-9311:
    - the "only if" direction of the CM detector;
    - the class-6 congruence (K);
    - the four-family law for C(2k,k)^2 C(3k,k), C(2k,k)^2 C(4k,2k), C(2k,k) C(3k,k) C(6k,3k), where no such reflection
      is available;
    - the mod-p^3 refinement at ordinary primes (HYP-9311(O)).
session: opus-2026-10-18-S43
source: 05-knowledge/results/halfperiod_twin_classes_cm_log_law_20261018.md
scripts:
  - 04-computation/experiments/cm_law_from_known_supercongruences_20261018.py (+ .out, .json)
  - 04-computation/experiments/cm_points_forms_legendre_20261018.py (+ .out, .json): the eight supercongruence values (FINITE-EXACT, p <= 200)
  - 04-computation/experiments/cm_harmonic_log_law_20261018.py, cm_harmonic_log_law_extended_20261018.py: the law directly
  - 04-computation/experiments/klmsy_conjecture9_s1_check_20261018.py (+ .out, .json): the Conjecture-9 (s = 1) corollary, M <= 11, p <= 150
related: [THM-4685, HYP-9310, HYP-9311, THM-4670, Kibelbek-Long-Moss-Sheller-Yuan 2016 Theorem 2, Ahlgren 2001]
---

# THM-4686 — the CM harmonic law from the supercongruences at both ends of a 2-isogeny

## Notation

`p >= 5` is prime, `n = (p-1)/2`, `a_k = C(2k,k)/4^k` and `O_k = Σ_(j<=k) 1/(2j-1)`. For `z ∈ Z_(p)^×` write:
- `F(z) = Σ_(k<=n) a_k^3 z^k`, so that `F(64/m) = Σ_(k<=n) C(2k,k)^3 m^(-k)`;
- `G(w) = Σ_(k<=n) a_k^3 w^k Δ_k`, with `Δ_k = O_(n-k) - O_k`;
- `h(z) = Σ_(k<=n) a_k^3 z^k (H_(2k) - H_k)`.

The law at `m` reads `L(m): 6 h(64/m) ≡ q_p(m) F(64/m) (mod p)`.

## Theorem

`L(m)` holds for every `m ∈ {1, -8, 16, -64, 64, 256, -512, 4096}` and every prime `p >= 5`, `p ∤ m`.

**Caveat at `p = 7`.** For `m ∈ {1, 4096}` and `p = 7`, the CM curve (discriminants `-7`, `-28`) has bad reduction,
which is outside KLMSY Theorem 2. There the statement is verified directly (FINITE-EXACT), as are the input
supercongruences at `p = 7`. Every other pair has good reduction at all `p >= 5`: its bad primes lie in `{2, 3}`
(and `7` for `λ = 64, 1/64`).

## Proof

**Step 1: the reflection mod `p^2`.** THM-4685 (S6) gives `a_(n-k) ≡ (-1)^n a_k (1 + pΔ_k) (mod p^2)`. Substituting
`k -> n-k`:
```
(R)   F(z) ≡ (-1)^n z^n [ F(1/z) + 3p G(1/z) ]   (mod p^2).
```

**Step 2: the input.** For `z ∈ {64, -8, 4}` (`m = 1, -8, 16`) and `1/z` (`m = 4096, -512, 256`), KLMSY Theorem 2 gives
```
F(z) ≡ ε_1 A,     F(1/z) ≡ ε_2 A   (mod p^2),     A = α_(p)^2,
```
- `ε_1 = ((z-1)/p)` and `ε_2 = ((1/z - 1)/p)`.
- `A` is the square of the unit root of the CM curve. The curves at `z` and `1/z` are 2-isogenous, so they have the same
  Frobenius and the same `A`.
- At supersingular `p`, `A = 0`.
- Consistency mod `p`: `ε_1 = (-1)^n (z/p) ε_2`. This holds because `((1/z - 1)/p) = ((1-z)/p)(z/p)`.

The eight values were also confirmed numerically mod `p^2` for `p <= 200`: the forms 7AB, 4C, `η(2τ)^3 η(6τ)^3` and
8A with the stated characters.

**Step 3: solve for `G(1/z)`.** Use `z^n ≡ (z/p)(1 + p q_p(z)/2) (mod p^2)`. Then (R) becomes
`ε_1 A ≡ ε_1 (1 + p q_p(z)/2) A + 3p (-1)^n (z/p) G(1/z) (mod p^2)`, so
```
G(1/z) ≡ -(q_p(z)/6) ε_2 A ≡ -(q_p(z)/6) F(1/z)   (mod p).
```
This also holds at supersingular primes, where both sides vanish.

**Step 4: the law at `1/z`.** By THM-4685 (M1)-(M3), `H_(2k) - H_k ≡ -Δ_k + q_p(2) (mod p)`. So
```
h(1/z) ≡ -G(1/z) + q_p(2) F(1/z) ≡ (q_p(z)/6 + q_p(2)) F(1/z) = (q_p(64z)/6) F(1/z)   (mod p),
```
which is `L(64z)`, i.e. `L(4096/m)`. By THM-4685(ii), `L(4096/m)` is equivalent to `L(m)`. This gives the six
non-self-dual points. The self-dual points `m = ±64` are THM-4685(ii). ∎

## Corollaries

- **Class 3** (census, THM-4685). `Σ_(j<=n) C(2j,j) C(p-1-j,j)^2 ≡ Σ_(j<=n) C(2j,j)^3 (mod p^2)`. Both sides are
  `≡ a(p)` for `η(τ)^3 η(7τ)^3` (KLMSY).
- **Class 7.** `Σ C(2n-2k,n-k)^2 C(n+k,2k) C(2n,n+k) ≡ Σ C(n+k,2k)^3 C(2n,n+k) (mod p^2)`.
- **Twin form at every CM point.** `T_m = Σ_j C(2j,j) C(p-1-j,j)^2 m^(-j) ≡ (1 - p q_p(m)/3) F_m (mod p^2)`.
- **Lucas-block factorization** at every CM point:
  `Σ_(k=p)^(p+n) C(2k,k)^3 m^(-k) ≡ (8/m) Σ_(k<=n) C(2k,k)^3 m^(-k) (mod p^2)`. Hence `F_(p^2) ≡ F_p^2 (mod p^2)`
  (HYP-9310, Lucas-block form).

## Corollary: KLMSY Conjecture 9 at `s = 1`

**Statement.** KLMSY (Conjecture 9) predict, for CM `λ`, ordinary `p`, odd `M` and `s >= 1`,
```
F_((Mp^s - 1)/2)(λ) ≡ ((λ-1)/p) α_(p,λ)^2 F_((Mp^(s-1) - 1)/2)(λ)   (mod p^(2s)),
```
with `F_N(λ) = Σ_(k<=N) ((1/2)_k/k!)^3 λ^k`. They prove `M = s = 1` (their Theorem 2). For each of the eight rational
CM values `λ = 64/m`, the case **`s = 1`, every odd `M`**, follows from the law above.

**Proof.**
- **p-shift lemma** (THM-4685, HYP-9310). For every digit `a >= 0` and carry-free `j <= n`,
  `C(2(ap+j), ap+j) ≡ C(2a,a) C(2j,j)(1 + 2ap(H_(2j) - H_j)) (mod p^2)`.
- **Block sums.** With `t_k = C(2k,k)^3 m^(-k)`, the block with base-`p` digit `a` is
  `Σ_(j<=n) t_(ap+j) ≡ t_a m^(-a(p-1)) (F_m + 6ap h_m)`. By the law this is `≡ t_a F_m (mod p^2)`.
- **Negligible terms.** Terms with `n < j < p` carry, so their cubes are divisible by `p^3`.
- **Summation.** Summing the blocks of `[0, (Mp-1)/2]`, whose last block has `j <= n`, gives
  ```
  (D)   F_((Mp-1)/2) ≡ F_((p-1)/2) · F_((M-1)/2)   (mod p^2).
  ```
- **Conclusion.** Theorem 2 gives `F_((p-1)/2) ≡ ((λ-1)/p) α^2`. At supersingular `p` both sides vanish mod `p^2`.

**Checks** (FINITE-EXACT, `klmsy_conjecture9_s1_check_20261018.py`):
- (D) holds at all eight CM values for `M = 3, 5, 7, 9, 11` and every prime `5 <= p <= 150`.
- (D) fails at 25-30 of about 34 primes for each generic `λ = 64/3, 64/5, 64/7, 32/5`.

**Novelty.** KLMSY state that their paper does not prove the general-`M` truncated-sum case. No later proof was found in
the S43 searches; the auditor is asked to check.

## Remarks

- **The pattern.** For both THM-4670 (audit M1) and this theorem, a half-period twin is the shadow of two classical
  supercongruences joined by an elementary `p`-adic identity:
  - THM-4670: the two guises at one CM point.
  - THM-4686: the two ends of a 2-isogeny.
- **The other families.** For `C(2k,k)^2 C(3k,k)`, `C(2k,k)^2 C(4k,2k)` and `C(2k,k) C(3k,k) C(6k,3k)` no reflection
  `k -> n-k` maps the family to itself, and HYP-9311 remains open.
