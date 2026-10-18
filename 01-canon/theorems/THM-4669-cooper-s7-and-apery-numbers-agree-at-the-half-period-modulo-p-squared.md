---
id: THM-4669
title: "Cooper's level-7 sporadic sequence and Apery's numbers agree at every Frobenius point modulo p^2: for every odd prime p and every r >= 1, with n = (p^r - 1)/2, s_7(n) == (-1)^n A(n) (mod p^2) (in particular at the half-period (p-1)/2), where A(n) = sum_k C(n,k)^2 C(n+k,k)^2 and s_7(n) = sum_k C(n,k)^2 C(n+k,k) C(2k,n); hence, by Ahlgren-Ono, s_7((p-1)/2) == (-1/p) a(p) (mod p^2) with a(p) the coefficients of eta(2 tau)^4 eta(4 tau)^4, Mason's eta-product of the M24 class 4A (cycle shape 2^4 4^4); and s_7((p-1)/2) == (-1)^((p-1)/2) [x^((p-1)/2)] (2K(sqrt x)/pi)^2 (mod p^2), the half-period coefficient of the squared complete elliptic integral"
status: >
  PROVED (elementary).
  Ingredients:
  - an exact binomial identity;
  - the standard half-period congruence C(n,k)C(n+k,k) == (-1)^k a_k^2 (mod p^2);
  - Clausen's formula with the quadratic transformation, 3F2(1/2,1/2,1/2;1,1;4x(1-x)) = 2F1(1/2,1/2;1;x)^2 (CITED,
    classical);
  - a p-adic perturbation of -1/2 = n - p/2 together with an exact symmetry k <-> n - k.

  The modular corollary uses Ahlgren-Ono (2000), A((p-1)/2) == a(p) (mod p^2) (CITED).

  The same proof works at every Frobenius point n = (p^r - 1)/2 (r >= 2). Steps (C1) and (C2) become p-adic valuation
  statements: v_p(2j-1) <= r-1 and v_p(n-i) <= r-1. This resolves HYP-9291, filed and proved in S43.
  FINITE-EXACT control: the theorem for all 237 primes 5 <= p < 1500 (and p = 3); every proof step for all odd primes
  p < 400 at r = 1, and p-adically at (p, r) = (5, 2), (7, 2), (11, 2), (13, 2), (5, 3), (7, 3)
  (s7_apery_frobenius_proof_check_20261018.py); the statement at 13 Frobenius points with r = 2, 3; the exact identities
  for n <= 40.
  - At r = 1 it never holds modulo p^3: all primes 5 <= p < 250 here, and all odd p < 4000 in the S43 audit.
  - At r >= 2 it fails modulo p^3 except at the non-ordinary prime p = 11 (a(11) = -44). There
    v_11(A) = v_11(s_7) = r and the difference has v_11 = r + 1 (r = 2, 3 here; r = 4 in the audit, n = 7320; audit m9).
  - It fails at n = p - 1 (even mod p) and at n = (3p-1)/2 (mod p^2; it holds mod p only at p = 11), for
    5 <= p <= 300.
  - Every one of these lines is produced by s7_apery_halfperiod_extended_20261018.py, rewritten after audit m10
    because the committed version could not produce the earlier .out.
  Independent audit confirmation:
  - r = 1: all 549 odd primes p < 4000;
  - 55 Frobenius points with n up to 14280 (p = 3 to r = 9, p = 5 to r = 6, p = 7 to r = 5, p = 11 and 13 to r = 4,
    p = 17, 19, 23 at r = 3, every prime p <= 149 at r = 2);
  - the Ahlgren-Ono corollary for all 429 odd primes p < 3000.
  The audit re-derived every step (s43_audit_20261018/REPORT.md).

  Novelty: no prior statement was found in any of the following:
  - Osburn-Sahu-Straub 2016 (Gauss-type s_7(mp) == s_7(m) mod p^3);
  - Gorodetsky 2021/23 (Lucas congruences);
  - Straub 2023;
  - Beukers 2024, including Appendix C, entry 9.11 on s_7;
  - two Z.-H. Sun conjecture papers;
  - Z.-W. Sun's open-conjecture list;
  - OEIS A183204 and A036917.
  These were searched in S43 and by the audit. Novelty is unverified.
session: opus-2026-10-18-S43
source: 05-knowledge/results/factorization_vertex_supersymmetric_bridge_20261018.md
scripts:
  - 04-computation/experiments/apery_like_m24_supercongruences_20261018.py (+ .out at PMAX 300 and 900): discovery in the M24 eta-product battery
  - 04-computation/experiments/apery_like_halfperiod_clusters_20261018.py (+ .out): s_7 == (-1/p) A is the only nontrivial new mod p^2 class
  - 04-computation/experiments/s7_apery_halfperiod_extended_20261018.py (+ .out, .json): p < 1500 and the Frobenius points (p^r - 1)/2
  - 04-computation/experiments/s7_apery_proof_check_20261018.py (+ .out): every step of the proof
  - 04-computation/experiments/swap_lemma_halfperiod_20261018.py (+ .out): the swap C(n+k,k) <-> C(2k,n) fails for nine other weights
  - 04-computation/experiments/s7_apery_frobenius_proof_check_20261018.py (+ .out): the proof steps at Frobenius points (p-adic)
related: [HYP-9291, THM-4652, HYP-9274, THM-4670, THM-4685, Ahlgren-Ono 2000, Kilbourn 2006, Cooper 2012 (s_7), Mason 1985 (M24 eta-products)]
---

# THM-4669 — s_7 and the Apéry numbers at the half-period

## Notation

- `p` is an odd prime, `r >= 1` and `n = (p^r - 1)/2`. For `r = 1` this is the half-period `(p-1)/2`.
- The Apéry numbers are `A(n) = Σ_k C(n,k)^2 C(n+k,k)^2`.
- Cooper's level-7 sporadic sequence is `s_7(n) = Σ_k C(n,k)^2 C(n+k,k) C(2k,n)` (Cooper 2012).
- `a_k = C(2k,k)/4^k = (1/2)_k/k! = (-1)^k C(-1/2, k)`, and `H_m = Σ_(j<=m) 1/j`.

## Statement

1. `s_7(n) ≡ (-1)^n A(n) (mod p^2)`.
2. `s_7(n) ≡ (-1)^n Σ_(k=0)^n a_k^2 a_(n-k)^2 = (-1)^n [x^n] 2F1(1/2,1/2;1;x)^2 (mod p^2)`, and
   `2F1(1/2,1/2;1;x) = (2/π) K(sqrt x)`.
3. (Ahlgren-Ono, CITED; `r = 1`) `A((p-1)/2) ≡ a(p) (mod p^2)`, where `Σ a(m) q^m = η(2τ)^4 η(4τ)^4` is the weight-4 newform of
   level 8. Hence `s_7((p-1)/2) ≡ (-1/p) a(p) (mod p^2)`.
   - `η(2τ)^4 η(4τ)^4` is Mason's eta-product of the M24 class with cycle shape `2^4 4^4`, class 4A.
   - `1/η_g = Tr(g q^(L_0 - 1))` is the trace of that element on the Fock space of 24 chiral bosons. `η_g` itself is a
     supertrace on a fermionic Fock space (corrected after audit m8; see the source note, section 3.3).

## Proof

**(E1) Exact.** `C(n,k) C(2k,n) = C(2k,k) C(k, n-k)`. Both sides equal `(2k)!/(k!(n-k)!(2k-n)!)`. Hence
```
s_7(n) = Σ_k [C(n,k) C(n+k,k)] · C(2k,k) C(k, n-k).
```

**(C1) Frobenius-point congruence.** For `0 <= k <= n`:
```
C(n,k) C(n+k,k) = Π_(j=1)^k (n-j+1)(n+j)/j^2,   (n-j+1)(n+j) = p^(2r)/4 - (j - 1/2)^2,
C(n,k) C(n+k,k) / ((-1)^k a_k^2) = Π_(j=1)^k (1 - p^(2r)/(2j-1)^2).
```
- Since `2j - 1 <= 2n - 1 < p^r`, we have `v_p(2j-1) <= r - 1`. So every factor is `≡ 1 (mod p^2)` in `Z_p`.
- `a_k ∈ Z_p`. Hence `C(n,k) C(n+k,k) ≡ (-1)^k a_k^2 (mod p^2)`.
- For `r = 1` this is the standard half-period congruence.
- Consequently `A(n) ≡ Σ_k a_k^4 (mod p^2)`.
- Combined with (E1), and using `C(2k,k) = 4^k a_k`:
  ```
  s_7(n) ≡ Σ_k (-1)^k 4^k a_k^3 C(k, n-k)   (mod p^2).
  ```

**(E2) Exact, via Clausen.** Since `(-1)^(n-k) C(k, n-k) = [x^n] x^k (1-x)^k`,
```
Σ_k (-1)^(n-k) 4^k a_k^3 C(k, n-k) = [x^n] Σ_k a_k^3 (4x(1-x))^k = [x^n] 3F2(1/2,1/2,1/2;1,1; 4x(1-x)).
```
- Clausen's formula gives `2F1(1/4,1/4;1;z)^2 = 3F2(1/2,1/2,1/2;1,1;z)`.
- The quadratic transformation gives `2F1(1/4,1/4;1;4x(1-x)) = 2F1(1/2,1/2;1;x)`. Both are identities of formal power
  series.
- Together they give `3F2(1/2,1/2,1/2;1,1;4x(1-x)) = 2F1(1/2,1/2;1;x)^2 = (Σ_m a_m^2 x^m)^2`.
- Hence the left side equals `Σ_(k=0)^n a_k^2 a_(n-k)^2`, for every `n`.
- With (C1) this proves statement 2.

**(C2) `p`-adic perturbation.** Put `ε = p^r/2`, so that `-1/2 = n - ε`. For `0 <= m <= n`,
```
a_m = (-1)^m C(n - ε, m) = (-1)^m C(n, m) Π_(i=0)^(m-1) (1 - ε/(n - i))
    ≡ (-1)^m C(n, m) (1 - ε S_m)   (mod p^2),     S_m = H_n - H_(n-m),
```
- Here `1 <= n - i < p^r` gives `v_p(n - i) <= r - 1`, so each `ε/(n - i) ∈ p Z_p`.
- Every product of two or more such terms lies in `p^2 Z_p`.
- `ε S_m` may involve denominators divisible by `p`, but it always lies in `p Z_p`.

**(L) Lemma.** `Σ_k a_k^2 a_(n-k)^2 ≡ Σ_k a_k^4 (mod p^2)`.
- By (C2), the difference is `2ε Σ_k C(n,k)^4 (S_k - S_(n-k)) = 2ε Σ_k C(n,k)^4 (H_k - H_(n-k))`.
- This vanishes exactly: `C(n,k)^4` is symmetric and `H_k - H_(n-k)` antisymmetric under `k <-> n-k`.
- The neglected terms are products of two elements of `p Z_p`, hence lie in `p^2 Z_p`.

**Conclusion.**
```
s_7(n) ≡ (-1)^n Σ_k a_k^2 a_(n-k)^2 ≡ (-1)^n Σ_k a_k^4 ≡ (-1)^n A(n)   (mod p^2).
```
This proves statement 1. Statement 3 follows from Ahlgren-Ono. ∎

## Remarks

- **Why `p^2` and not `p^3`.** (C1) and (C2) discard the `ε^2` terms.
  - At `r = 1` the congruence fails modulo `p^3` for every tested prime.
  - At `r >= 2` the non-ordinary prime `11` is the exception: there both sides have `11`-adic valuation `r`.
- **Not a general swap** (among the ten weights tested). Replacing `C(n+k,k)` by `C(2k,n)` preserves half-period values
  only for the weight `C(n,k)^2 C(n+k,k)`. It fails, even modulo `p`, for the nine other tested weights (swap_lemma_halfperiod_20261018.out). The
  mechanism is the Clausen identity (E2), which is specific to the `1/2`-hypergeometric structure.
- **Frobenius points.** The statement at `n = (p^r - 1)/2`, `r >= 2`, was filed as
  [HYP-9291](../../05-knowledge/hypotheses/HYP-9291-s7-and-apery-agree-at-every-frobenius-point-modulo-p-squared.md).
  It is proved by the same argument with the valuation bounds in (C1) and (C2). Denominators `j` with `p | j` cancel in
  the ratio against `a_k^2`.
- **Context.** The discovery came from testing Apéry-like sequences against the 26 Mason eta-products of M24.
  - The governing form here is the class-4A eta-product, and Apéry's ζ(2) numbers are governed by the class-4C
    eta-product `η(4τ)^6` (Stienstra-Beukers).
  - See factorization_vertex_supersymmetric_bridge_20261018.md for the bridge to vertex algebras and supersymmetric
    QFT.
