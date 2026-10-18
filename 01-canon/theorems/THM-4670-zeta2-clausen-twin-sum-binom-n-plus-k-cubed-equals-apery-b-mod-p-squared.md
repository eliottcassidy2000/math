---
id: THM-4670
title: "The zeta(2)-class half-period twin: for every prime p >= 5, with n = (p-1)/2, sum_{k=0}^{n} C(n+k,k)^3 == b(n) = sum_k C(n,k)^2 C(n+k,k) (mod p^2) (Apery's zeta(2) numbers; hence == a(p), the coefficients of eta(4 tau)^6, the M24 class-4C eta-product). An immediate corollary of the classical supercongruences of van Hamme, Ishikawa and Ahlgren (audit M1); the content here is an elementary proof, by differentiating Dixon's well-poised summation at a = -n, that their two classical guises agree mod p^2: for odd n, sum_k (-1)^k C(n,k)^3 H_(n-k) = Gamma(1-n/2) Gamma(1+3n/2)/(3n Gamma(1+n/2)^2); for n = 2m, S_n - H_n D_n = (D_n/2)(H_m - H_(3m)) with D_n = (-1)^m (3m)!/m!^3"
status: >
  PROVED (elementary + classical summation). It is a COROLLARY of classical results (S43 audit, MAJOR M1).
  - Step 1 below gives sum C(n+k,k)^3 - b(n) == 4 (sum_{k<=n} a_k^3 - b(n)) (mod p^2).
  - Both b((p-1)/2) == c(p) and sum_{k<=(p-1)/2} a_k^3 == c(p) (mod p^2), with c(p) the coefficients of
    eta(4 tau)^6, are classical. They are the two guises of Beukers' zeta(2) supercongruence, conjectured by
    Stienstra-Beukers 1985 and proved by van Hamme 1986, Ishikawa 1990 and Ahlgren 2001 (Mortenson, Proc. AMS 136
    (2008), eqs (1.2) and (1.3), where the second is "another form" of the first). Together they imply the theorem at
    once.
  - The contribution here is a self-contained elementary proof that the two guises agree mod p^2: the key lemma
    sum a_k^3 O_k == 0 (mod p) is exactly that agreement. This proof may itself be classical.
  - The name "Clausen twin" used at filing was a misnomer: no Clausen identity is used (audit m11). The file name keeps
    the historical slug.
  Ingredients of the elementary proof:
  - the p-adic perturbation C(n+k,k) == a_k (1 + p O_k), C(n,k) == (-1)^k a_k (1 - p O_k) (mod p^2), as in THM-4669;
  - Dixon's theorem for the well-poised 3F2(1) (CITED, classical), differentiated in a at a = -n with b = c = -n;
  - the symmetry j <-> p - j of [m+1, 3m] for p = 4m + 1.
  The modular corollary uses b((p-1)/2) == a(p) (mod p^2) for eta(4 tau)^6. That is CITED to van Hamme 1986, Ishikawa
  1990 and Ahlgren 2001 (conjectured by Stienstra-Beukers 1985, who proved it mod p); the earlier attribution to
  Stienstra-Beukers was corrected (audit m5).
  FINITE-EXACT control:
  - the theorem for all primes 5 <= p < 1200, and in the audit for all 667 primes 5 <= p < 5000;
  - both closed forms for n <= 40 (zeta2_twin_dixon_proof_check_20261018.py).
  Resolves HYP-9292, filed and proved the same session. The ID collision with mac-mini's earlier stub reservation was
  resolved by mac-mini ceding THM-4670 / HYP-9292 (commit 11edd404cd; theirs is now THM-4690 / HYP-9325).
session: opus-2026-10-18-S43
source: 05-knowledge/results/factorization_vertex_supersymmetric_bridge_20261018.md
scripts:
  - 04-computation/experiments/clausen_twin_search_20261018.py (+ .out, .json): discovery
  - 04-computation/experiments/zeta2_twin_dixon_20261018.py (+ .out): the reduction, numerically
  - 04-computation/experiments/zeta2_twin_dixon_proof_check_20261018.py (+ .out): closed forms and theorem
related: [THM-4669, THM-4685, HYP-9292, Dixon's theorem, van Hamme 1986, Ishikawa 1990, Ahlgren 2001, Mortenson 2008, Stienstra-Beukers 1985, Mason 1985]
---

# THM-4670 — the ζ(2)-class half-period twin

## Statement

For every prime `p >= 5`, with `n = (p-1)/2`,
```
Σ_(k=0)^n C(n+k,k)^3  ≡  b(n) = Σ_k C(n,k)^2 C(n+k,k)   (mod p^2).
```
Since `b((p-1)/2) ≡ a(p) (mod p^2)` with `Σ a(m) q^m = η(4τ)^6`, both sides are `≡ a(p) (mod p^2)`.

## Proof

**Step 1: reduction (as in THM-4669).** Let `a_k = C(2k,k)/4^k` and `O_k = Σ_(j=1)^k 1/(2j-1)`.
- Since `n + 1 = 1/2 + p/2` and `-n = 1/2 - p/2`, and every `2j - 1 <= 2n - 1` is a `p`-adic unit,
  ```
  C(n+k,k) ≡ a_k (1 + p O_k),   C(n,k) ≡ (-1)^k a_k (1 - p O_k)   (mod p^2).
  ```
- Hence `Σ C(n+k,k)^3 - b(n) ≡ 4p Σ_(k<=n) a_k^3 O_k (mod p^2)`.
- Modulo `p`: `a_k ≡ (-1)^k C(n,k)`, and `O_k ≡ -(H_n - H_(n-k))/2` (because `2j - 1 ≡ -2(n - j + 1)`).
- So, with `D_n = Σ (-1)^k C(n,k)^3`, `S_n = Σ (-1)^k C(n,k)^3 H_(n-k)` and `E_n = S_n - H_n D_n`,
  ```
  Σ a_k^3 O_k ≡ E_n/2   (mod p),
  ```
  and the theorem is equivalent to `E_n ≡ 0 (mod p)`.

**Step 2: Dixon's theorem with a parameter.** Dixon's theorem states
`3F2(a,b,c;1+a-b,1+a-c;1) = Γ(1+a/2)Γ(1+a-b)Γ(1+a-c)Γ(1+a/2-b-c)/(Γ(1+a)Γ(1+a/2-b)Γ(1+a/2-c)Γ(1+a-b-c))`.
- Take `a = -n + δ` and `b = c = -n`. The series terminates for every `δ`, and the identity holds as an identity of
  rational functions in `δ`.
- Differentiate at `δ = 0`. The left side gives `Σ_k (-1)^k C(n,k)^3 [-(H_n - H_(n-k)) - 2 H_k] = -H_n D_n + S_n - 2T_n`,
  with `T_n = Σ (-1)^k C(n,k)^3 H_k`. Under `k <-> n - k`, `T_n = (-1)^n S_n`.
- **`n` odd.** `D_n = 0`, and the right side vanishes at `δ = 0` through `1/Γ(1 - n + δ) ≈ (n-1)! δ`. So
  ```
  3 S_n = Γ(1 - n/2) Γ(1 + 3n/2) / (n Γ(1 + n/2)^2).
  ```
  - In double factorials, `Γ(1+3n/2)/sqrt(π) = (3n)!!/2^((3n+1)/2)` contains the odd number `p = 2n + 1 <= 3n`
    exactly once.
  - The other factors involve only odd numbers `<= n`, the factor `3n` (`p` divides neither `3` nor `n`), and powers
    of `2`.
  - Hence `v_p(S_n) = 1` and `E_n = S_n ≡ 0 (mod p)`.
- **`n = 2m` even.** The pole of `Γ(1 + a/2)` and the zero of `1/Γ(1 + a)` cancel. By the reflection formula their
  ratio is `2(-1)^m cos(πδ/2) Γ(2m - δ)/Γ(m - δ/2)`, whose `cos` factor has zero derivative. The value at `δ = 0` is
  Dixon's `D_n = (-1)^m (3m)!/m!^3`. The logarithmic derivative (digamma values, with `γ` cancelling) gives
  ```
  -H_n D_n - S_n = D_n [-H_(2m-1) + H_(3m)/2 + H_(m-1)/2 - H_m - H_(2m)],
  ```
  that is, `E_n = (D_n/2)(H_m - H_(3m))`.
  - `D_n` is a `p`-adic integer, since `3m < p = 4m + 1`.
  - `H_(3m) - H_m = Σ_(j=m+1)^(3m) 1/j ≡ 0 (mod p)`, because `j -> p - j` maps `[m+1, 3m]` to itself.
  - Hence `E_n ≡ 0 (mod p)`. ∎

## Remarks

- **Comparison with THM-4669.** THM-4669 (the ζ(3) class, M24 class 4A) needs the exact Clausen identity and an
  antisymmetry. THM-4670 (the ζ(2) class, 4C) needs the derivative of Dixon's summation and a reflection symmetry.
- **The search.** `clausen_twin_search_20261018.py` searched all three- and four-factor binomial summands from a library
  of ten factors. These two are the only half-period twins of the Apéry numbers it found; the other hits are
  rewritings. The full census of twin classes is THM-4685.
- **Failure at `p = 3`.** For `p = 3` (`n = 1`), `S_1 = 1`, and the congruence fails.
- **Audit observations** (not claimed):
  - The congruence holds even mod `p^3` at `p = 7, 149, 241`.
  - At `n = (p^2-1)/2`, `v_p(twin - b) >= 2` for `p = 3, 5, ..., 23`, including `p = 3`, where `r = 1` fails. This
    suggests a Frobenius-point extension like THM-4669's.
