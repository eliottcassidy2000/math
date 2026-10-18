---
id: HYP-9292
title: "The zeta(2)-class half-period twin (filed as the zeta(2)-class Clausen twin; no Clausen identity is involved): for every prime p >= 5, with n = (p-1)/2, sum_{k=0}^{n} C(n+k,k)^3 == b(n) (mod p^2), where b(n) = sum_k C(n,k)^2 C(n+k,k) are Apery's zeta(2) numbers (governed at the half-period by the M24 class-4C eta-product eta(4 tau)^6); equivalently sum_{k<=n} a_k^3 O_k == 0 (mod p), a_k = C(2k,k)/4^k, O_k = sum_{j<=k} 1/(2j-1); equivalently the harmonic Dixon congruence H_n D_n == S_n (mod p). RESOLVED: PROVED in S43 as THM-4670"
status: >
  PROVED (resolved the same session): THM-4670 proves it by differentiating Dixon's well-poised summation.
  S43 audit (MAJOR M1): the statement is an immediate corollary of the classical supercongruences.
  - b((p-1)/2) == c(p) and sum_{k<=(p-1)/2} a_k^3 == c(p) (mod p^2) are Beukers' zeta(2) supercongruence in its two
    guises: van Hamme 1986, Ishikawa 1990, Ahlgren 2001 (Mortenson 2008, eqs (1.2), (1.3)).
  - sum C(n+k,k)^3 - b == 4 (sum a_k^3 - b) (mod p^2).
  - The lemma sum a_k^3 O_k == 0 is exactly the agreement of the two guises. The Dixon argument is an elementary proof
    of that agreement, and was not a new congruence.
  - odd n: the closed form S_n = Gamma(1-n/2) Gamma(1+3n/2)/(3n Gamma(1+n/2)^2) is divisible by p;
  - n = 2m: E_n = (D_n/2)(H_m - H_(3m)), and H_(3m) - H_m == 0 (mod 4m+1) by the reflection j <-> p - j.
  FINITE-EXACT: the twin congruence holds for all primes 5 <= p < 1000; both equivalent forms hold for all primes
  p < 400.
  Here D_n = sum_k (-1)^k C(n,k)^3 (Dixon's sum, zero for odd n) and S_n = sum_k (-1)^k C(n,k)^3 H_(n-k).
  The reduction (PROVED, by the p-adic perturbation of THM-4669) uses C(n+k,k) == a_k (1 + p O_k) and
  C(n,k) == (-1)^k a_k (1 - p O_k) (mod p^2).
  Found by an exhaustive search over four-factor binomial summands, clausen_twin_search_20261018.py. In that search,
  every other hit for the zeta(3) class was s_7 of THM-4669 in disguise.
session: opus-2026-10-18-S43
source: 05-knowledge/results/factorization_vertex_supersymmetric_bridge_20261018.md
scripts:
  - 04-computation/experiments/clausen_twin_search_20261018.py (+ .out, .json)
  - 04-computation/experiments/zeta2_twin_dixon_20261018.py (+ .out)
related: [THM-4669, THM-4670, THM-4685, Ahlgren-Ono 2000, van Hamme 1986, Ishikawa 1990, Ahlgren 2001, Mortenson 2008, Stienstra-Beukers 1985, Dixon's identity]
---

# HYP-9292 — the ζ(2)-class half-period twin

## Statement

For every prime `p >= 5`, with `n = (p-1)/2`,
```
Σ_(k=0)^n C(n+k,k)^3  ≡  b(n) = Σ_k C(n,k)^2 C(n+k,k)   (mod p^2).
```

Combined with the known half-period supercongruence `b((p-1)/2) ≡ a(p) (mod p^2)` (van Hamme 1986, Ishikawa 1990,
Ahlgren 2001; conjectured by Stienstra-Beukers 1985; `a` the coefficients of `η(4τ)^6`), this would give
`Σ_k C(n+k,k)^3 ≡ a(p) (mod p^2)`.

## Reduction (PROVED)

Write `ε = p/2`, so that `n + 1 = 1/2 + ε`. Then:
- `C(n+k,k) = (1/2 + ε)_k/k! ≡ a_k (1 + p O_k) (mod p^2)`;
- `C(n,k) = (-1)^k (1/2 - ε)_k/k! ≡ (-1)^k a_k (1 - p O_k) (mod p^2)`.

Hence `Σ C(n+k,k)^3 - b(n) ≡ 4p Σ_k a_k^3 O_k (mod p^2)`.

Modulo `p`, `a_k ≡ (-1)^k C(n,k)` and `O_k ≡ -(H_n - H_(n-k))/2`. So the statement is equivalent to
`H_n D_n ≡ S_n (mod p)`.

## Structure of the defect (FINITE-EXACT, n <= 24)

Let `E_n = S_n - H_n D_n` (exact rationals).
- **Odd `n`.** `D_n = 0`, and the numerator of `S_n` is, up to small primes, the product of the primes in
  `(n + 1, 3n]`:
  - `S_3 = -35/3 = -5·7/3`;
  - `S_5 = 1001/5 = 7·11·13/5`;
  - `S_9 = 5·11·13·17·19·23/63`;
  - `S_23 = -3·5^2·29·31·...·67/7429`.

  So `p = 2n + 1 ∈ (n+1, 3n]` divides `S_n`. This is the signature of a Gamma-ratio closed form, obtained by
  differentiating Dixon's well-poised `3F2(1)` summation at a zero of `1/Γ`.
- **Even `n`.** The numerator of `E_n` has large irregular prime factors, for example `74587` at `n = 10`. It is still
  divisible by `2n + 1` whenever `2n + 1` is prime (`n = 2, 6, 8, 14, 18, 20, 24`).

The odd case therefore looks closed-form provable; the even case needs the full derivative of Dixon's theorem.

## Tests

- **Proof: done (THM-4670).** Differentiate the three-parameter Dixon identity
  `Σ_k (-1)^k C(a+b,a+k) C(b+c,b+k) C(c+a,c+k) = (a+b+c)!/(a! b! c!)` along the diagonal.
- **Literature.** Search for `Σ C(2k,k)^3 H_k/64^k`-type congruences (Z.-W. Sun and coauthors).
