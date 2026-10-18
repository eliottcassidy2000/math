---
id: HYP-9310
title: "CM harmonic log law: at each of Sun's CM values m in {1, 4096, -8, -512, 16, 256}, 6 sum_{k<=(p-1)/2} C(2k,k)^3 m^-k (H_2k - H_k) == q_p(m) sum_{k<=(p-1)/2} C(2k,k)^3 m^-k (mod p) for every prime p not dividing 6m; equivalently the twin sum_j C(2j,j) C(p-1-j,j)^2 m^-j == (1 - p q_p(m)/3) sum_j C(2j,j)^3 m^-j (mod p^2). The law is proved at the self-dual points m = +-64 and is invariant under m -> 4096/m (THM-4685), so three cases remain: m = 1 (Q(sqrt-7)), m = -8 (Q(i)), m = 16 (Q(sqrt-3)). The case m = 1 is the census twin (H): sum C(2k,k)^3 (H_2k - H_k) == 0 (mod p). Includes the class-6 congruence (K)"
status: >
  OPEN, strongly supported.
  - Pre-registered prediction PASSED: the law at every CM point, failure at generic m. It was written in the script
    docstring before the run.
  - 8 CM points x 301 primes (p < 2000): 0 failures. Generic m (3, 5, 7, 10, 100, -3, 2, -1, 32, 128): 287-296
    failures out of about 300 primes each.
  - (H) for p < 3000; the twin forms at all 8 points for p <= 200.
  - Extended run: m = 1, -8, 16 (one per open reflection class) and (K), all 2260 primes 5 <= p <= 20000,
    0 failures (cm_harmonic_log_law_extended_20261018.out).
  - The naive archimedean analogue (h/F = log|m|/6) is REFUTED: h/F is about 0.001-0.02 against log|m|/6 of about 1.
session: opus-2026-10-18-S43
source: 05-knowledge/results/halfperiod_twin_classes_cm_log_law_20261018.md
scripts:
  - 04-computation/experiments/halfperiod_twin_classes_20261018.py
  - 04-computation/experiments/halfperiod_twin_mechanisms_20261018.py
  - 04-computation/experiments/cm_harmonic_log_law_20261018.py
  - 04-computation/experiments/cm_points_forms_legendre_20261018.py
  - 04-computation/experiments/cm_harmonic_log_law_extended_20261018.py
related: [THM-4685, THM-4670, THM-4669, Sun Conjecture 29 and 36 (arXiv:2210.07238), Kibelbek-Long-Moss-Sheller-Yuan (arXiv:1210.4489), Guillera-Zudilin divergent series, Coster-van Hamme 1991]
---

# HYP-9310 — the CM harmonic log law

## Statement

Let `n = (p-1)/2`, `F_m = Σ_(k<=n) C(2k,k)^3 m^(-k)` and `h_m = Σ_(k<=n) C(2k,k)^3 m^(-k) (H_(2k) - H_k)`. For each
`m ∈ {1, 4096, -8, -512, 16, 256}` and every prime `p ∤ 6m`:
```
(L_p)   6 h_m ≡ q_p(m) F_m   (mod p).
```
By (S5) of THM-4685 this is equivalent to the twin
```
T_m = Σ_(j<=n) C(2j,j) C(p-1-j,j)^2 m^(-j) ≡ (1 - p q_p(m)/3) F_m   (mod p^2).
```
For a cube `m = c^3` the twin reads `c^(p-1) T_m ≡ F_m (mod p^2)`.

**Part (K)** (the class-6 twin of the census, governed by 4A modulo `p`):
```
Σ_(j<=n/2) a_j^3 a_(2j) (H_(2j) - 2H_j - 2H_(4j) + 2H_(n-j) + H_(n-2j)) ≡ 0   (mod p),
```
where `a_k = C(2k,k)/4^k`, so that `a_j^3 a_(2j) = C(2j,j)^3 C(4j,2j)/1024^j`.

## What is proved (THM-4685)

- `(L_p)` at `m = 64` and `m = -64`.
- `(L_p)` at `m` is equivalent to `(L_p)` at `4096/m`. So `{1, 4096}`, `{-8, -512}` and `{16, 256}` are three single
  open cases.
- The census twin classes 3 and 7 are both equivalent to `(L_p)` at `m = 1`, i.e. to
  ```
  (H)   Σ_(k<=n) C(2k,k)^3 (H_(2k) - H_k) ≡ 0   (mod p).
  ```

## Consequences if true

- **Class 3.** `Σ_j C(2j,j) C(p-1-j,j)^2 ≡ Σ_j C(2j,j)^3 (mod p^2)`. With Kibelbek-Long-Moss-Sheller-Yuan (CITED),
  both are `≡ a(p)` for `η(τ)^3 η(7τ)^3` (M24 class 7AB), i.e. `≡ 4x^2 - 2p` when `p = x^2 + 7y^2`.
- **Class 7.** `Σ C(2n-2k,n-k)^2 C(n+k,2k) C(2n,n+k) ≡ Σ C(n+k,2k)^3 C(2n,n+k) (mod p^2)`.
- **General twin.** At every Sun CM point the 'shifted-square' sum `T_m` is the CM value times the explicit `p`-adic
  cube root `(m^(p-1))^(-1/3)`.

## Context and heuristics

- **Isogeny.** The reflection `k -> n-k` acts as `z = 64/m -> 1/z`, i.e. `m -> 4096/m`. It swaps the 2-isogenous CM
  curves: `j = -15^3 ↔ 255^3`, `1728 ↔ 66^3` and `0 ↔ 54000` (discriminants `-7 ↔ -28`, `-4 ↔ -16`, `-3 ↔ -12`).
  At the self-dual points `j = 1728` (`m = 64`) and `j = 8000` (`m = -64`) the law is forced by symmetry. Elsewhere
  it asserts that the first-order `p`-adic (harmonic) correction is compatible with the isogeny, which needs CM input.
  - The generic failures show that the law is not formal.
- **Legendre squares.** `F_m ≡ P_n(x)^2 (mod p^2)` with `x^2 = 1 - 64/m`. This follows from the exact identity
  `P_N(x)^2 = Σ_k C(N+k,2k) C(2k,k)^2 ((x^2-1)/4)^k`, checked for `N <= 25` (classical, Clausen/Bailey type), together
  with (S3). A proof route: find the analogous Legendre form of `T_m` and use the Coster-van Hamme CM theorem.
- **Relation to Sun.** Z.-W. Sun's Conjecture 29(i) (arXiv:2210.07238) is a template for the Ramanujan-type,
  `k`-weighted sums `Σ C(2k,k)^3 m^(-k) [6(ak+b)(H_(2k) - H_k) + a]`, with `log|m|` (archimedean) and `b(m^(p-1) - 1)`
  (`p`-adic). `(L_p)` is the constant-term analogue: no polynomial weight, the same `6 <-> q_p(m)` pairing, and every
  CM point including the divergent one `m = 1`. The S43 reading of the paper found no statement of `(H)` itself.
  - Sun's 17(iii) at `m = 1` sums `C(2k,k)^3 [2(21k+8)(H_(2k) - H_k) + 7]` over `1 <= k <= p-1`. Its value
    modulo `p` (zero) was re-checked for `7 <= p <= 43`.
  - Combined with `(H)` it would give `Σ_(k<=n) k C(2k,k)^3 (H_(2k) - H_k) ≡ -(F_1 - 1)/6 (mod p)`, the `-1`
    coming from the missing `k = 0` term. This consequence was confirmed numerically for `7 <= p <= 43`.

## Falsifiers

A single prime `p ∤ 6m` where `(L_p)` fails at a listed `m`, or where `(K)` fails.
