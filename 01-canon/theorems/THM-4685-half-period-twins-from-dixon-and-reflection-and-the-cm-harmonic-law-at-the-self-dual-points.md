---
id: THM-4685
title: "Half-period twins from Dixon's lemma and the reflection k -> (p-1)/2 - k. For every prime p >= 5, with n = (p-1)/2 and F = sum_{k<=n} C(2k,k)^3/64^k: (i) sum C(n+k,k)C(2k,k)C(2n-2k,n-k) == (-1)^n 4^n F and sum C(n,k)^3 C(2n,n+k) == sum C(n,k)C(n+k,k)^2 C(2n,n+k) == (-1)^n 16^n F (mod p^2), alongside b(n) == sum C(n+k,k)^3 == F (THM-4670); F == 0 (mod p^2) when p == 3 (mod 4). (ii) The CM harmonic law 6 sum_{k<=n} C(2k,k)^3 m^-k (H_2k - H_k) == q_p(m) sum_{k<=n} C(2k,k)^3 m^-k (mod p) holds at the two self-dual CM points m = 64 and m = -64, and the law at m is equivalent to the law at 4096/m. Equivalently 4^(p-1) sum_j C(2j,j) C(p-1-j,j)^2 (+-64)^-j == sum_j C(2j,j)^3 (+-64)^-j (mod p^2)"
status: >
  PROVED (elementary). Ingredients:
  - the p-adic shadows of the binomials at n = (p-1)/2 (Lemma 1);
  - THM-4670's two Dixon closed forms (n odd: v_p(S_n) = 1; n even: S_n == H_n D_n mod p);
  - the reflection k -> n - k (Lemma 3);
  - the classical H_((p-1)/2) == -2 q_p(2) (mod p).
  Modular values are CITED, not reproved: F == a(p) (mod p^2) for eta(4t)^6 (Ishikawa 1990 / van Hamme / Ahlgren 2001);
  C(n, n/2)^2 == 4^n a(p) (Chowla-Dwork-Evans 1986, combined with the former); sum_{k<=n} C(2k,k)^3 (-64)^-k ==
  (-1/p) a(p) for eta(t)^2 eta(2t) eta(4t) eta(8t)^2. This is Sun's Conjecture 5.2, proved mod p^2 in all cases by
  Kibelbek-Long-Moss-Sheller-Yuan 2016, Theorem 2 (character ((1-lambda)/p); v3 misprints ((lambda-1)/p)).
  Z.-H. Sun's paper covers only some cases (S43b n8).
  FINITE-EXACT controls:
  - the census for 5 <= p <= 200 (halfperiod_twin_classes) and the class forms for 5 <= p <= 300 (corrected, S43b
    m6);
  - every step of Lemmas 1-3 and parts (i)-(ii), re-derived by the S43b referee: mod p^2 for all 428 primes
    5 <= p <= 3000, and the mod-p statements for p <= 6000;
  - the shadow formulas for p < 200;
  - the law at m = +-64 for p < 2000.
  The six remaining CM points are KLMSY 2016 Corollary 3 (published; S43b M1), re-derived in THM-4686; see HYP-9310.
  Calibration (S43 audit M1 on THM-4670): Lemma 2(a), sum a_k^3 O_k == 0 (mod p), is exactly the agreement mod p^2 of
  the two classical guises of Beukers' zeta(2) supercongruence (van Hamme 1986, Ishikawa 1990, Ahlgren 2001).
  - Part (i) is therefore a corollary of those classical results, plus the elementary shadows and reflection. The Dixon
    route makes it self-contained.
  - What is new in (i) is only the list of further binomial sums with these values: 4^n F, 16^n F.
  Novelty not established. The law at m = 64 may already be in the harmonic-congruence literature.
session: opus-2026-10-18-S43
source: 05-knowledge/results/halfperiod_twin_classes_cm_log_law_20261018.md
scripts:
  - 04-computation/experiments/halfperiod_twin_classes_20261018.py (+ .out, .json): the census (seven classes)
  - 04-computation/experiments/halfperiod_twin_class_forms_20261018.py (+ .out, .json): governing eta-products
  - 04-computation/experiments/halfperiod_twin_mechanisms_20261018.py (+ .out, .json): shadows, (H), (K)
  - 04-computation/experiments/cm_harmonic_log_law_20261018.py (+ .out, .json): the law at 8 CM and 10 generic m
  - 04-computation/experiments/cm_points_forms_legendre_20261018.py (+ .out, .json): forms per CM point, Legendre squares, twins
related: [THM-4669, THM-4670, HYP-9310, Dixon's theorem, Chowla-Dwork-Evans 1986, Ahlgren 2001, Sun's Conjecture 29 (arXiv:2210.07238)]
---

# THM-4685 — half-period twins from Dixon's lemma and the reflection

## Notation

- `p >= 5` is prime and `n = (p-1)/2`.
- `a_k = C(2k,k)/4^k`, `O_k = Σ_(j<=k) 1/(2j-1)`, `H_k` is the harmonic number, and `q_p(m) = (m^(p-1) - 1)/p`.
- `Δ_k = O_(n-k) - O_k`.
- For `z = 64/m`: `F(z) = Σ_(k<=n) a_k^3 z^k = Σ_(k<=n) C(2k,k)^3 m^(-k)` and
  `h(z) = Σ_(k<=n) a_k^3 z^k (H_(2k) - H_k)`. `F = F(1)` is the truncated `3F2(1/2,1/2,1/2;1,1;1)`.
- From THM-4670: `D_n = Σ (-1)^k C(n,k)^3`, `S_n = Σ (-1)^k C(n,k)^3 H_(n-k)`, and
  `T_n = Σ (-1)^k C(n,k)^3 H_k = (-1)^n S_n`.

## Lemma 1 (shadows)

For `0 <= k <= n` the following hold modulo `p^2`.

| | congruence (mod `p^2`) | reason |
|---|---|---|
| (S1) | `C(n,k) ≡ (-1)^k a_k (1 - p O_k)` | `n - j + 1 = (p - (2j-1))/2` |
| (S2) | `C(n+k,k) ≡ a_k (1 + p O_k)` | `n + j = (p + (2j-1))/2` |
| (S3) | `C(n+k,2k) ≡ (-1)^k a_k / 4^k` | the factors pair to `Π_(i<=k) (p^2 - (2i-1)^2)/4`, so only `p^2` enters |
| (S4) | `C(2n,j) ≡ (-1)^j (1 - p H_j)` for `0 <= j <= 2n` | `2n = p - 1` |
| (S5) | `C(2n-j, j) ≡ (-1)^j C(2j,j) (1 - p (H_(2j) - H_j))` | the numerator is `Π_(i=j+1)^(2j) (p - i)` |
| (S6) | `a_(n-k) ≡ (-1)^n a_k (1 + p Δ_k)` | (S1) applied to `C(n,k) = C(n,n-k)` |

The following hold modulo `p`.

| | congruence (mod `p`) | reason |
|---|---|---|
| (M1) | `O_k ≡ (H_(n-k) - H_n)/2` | `2j - 1 = p - 2(n-j+1)` |
| (M2) | `H_(p-1-i) ≡ H_i` | |
| (M3) | `H_n ≡ -2 q_p(2)` | classical |
| (M4) | `a_k ≡ (-1)^k C(n,k)` | |

Consequences:
- `Δ_k ≡ (H_k - H_(n-k))/2`;
- `H_(2k) - H_k ≡ (H_(n-k) - H_k - H_n)/2`;
- `H_(n+k) ≡ H_(n-k)`.

## Lemma 2 (Dixon inputs)

THM-4670 gives `D_n = 0` and `S_n ≡ 0` for odd `n`, and `S_n ≡ H_n D_n` for even `n` (mod `p`). Hence, modulo `p`:
- (a) `Σ a_k^3 O_k ≡ (S_n - H_n D_n)/2 ≡ 0`;
- (b) `Σ a_k^3 O_(n-k) ≡ (T_n - H_n D_n)/2 ≡ 0`;
- (c) `Σ a_k^3 H_(n-k) ≡ S_n ≡ H_n F`.

**(d) For odd `n`, `F ≡ 0 (mod p^2)`.**
- By (S6), `F = Σ a_(n-k)^3 ≡ -Σ a_k^3 (1 + 3pΔ_k)`, so `2F ≡ -3p Σ a_k^3 Δ_k`.
- Moreover `Σ a_k^3 Δ_k ≡ (T_n - S_n)/2 = -S_n ≡ 0 (mod p)`.

## Theorem, part (i): three twin classes of `F` (all modulo `p^2`)

**(i-a)** By (S2), (S6) and (b):
```
Σ C(n+k,k) C(2k,k) C(2n-2k,n-k) = 4^n Σ C(n+k,k) a_k a_(n-k) ≡ (-1)^n 4^n Σ a_k^3 (1 + p O_(n-k)) ≡ (-1)^n 4^n F.
```

**(i-b)** By (S1), (S4), (M2), (a), (c) and (M3), using `1 - pH_n ≡ 1 + 2p q_p(2) ≡ 2^(2(p-1)) = 16^n`:
```
Σ C(n,k)^3 C(2n,n+k) ≡ (-1)^n Σ a_k^3 (1 - 3p O_k - p H_(n+k)) ≡ (-1)^n (F - p H_n F) ≡ (-1)^n 16^n F.
```

**(i-c)** By (a), `Σ C(n,k) C(n+k,k)^2 C(2n,n+k) ≡ (-1)^n Σ a_k^3 (1 + p O_k - p H_(n+k))` is congruent to the sum in (i-b).
It also needs no Dixon input. Termwise, `C(n+k,k) C(2n,n+k) ≡ (-1)^n 16^n C(n,k) (mod p^2)` for `0 <= k <= n`, so
`(i-c) ≡ (-1)^n 16^n b(n)` exactly as a congruence of sums (S43b n1).

**(i-d)** This is THM-4670 restated: `b(n) ≡ Σ a_k^3 (1 - pO_k)` and `Σ C(n+k,k)^3 ≡ Σ a_k^3 (1 + 3pO_k)`, and both are
`≡ F` by (a).

**Modular corollary (CITED).** `F ≡ a(p)` for `η(4τ)^6`, M24 class 4C, so the three classes take the values `a(p)`,
`4^n a(p)` and `16^n a(p)`. The sign `(-1)^n` is immaterial because `a(p) = 0` for `p ≡ 3 (mod 4)`, matching (d). The
census partner `C(n,n/2)^2 [n even]` of (i-a) is `≡ 4^n a(p)` by Chowla-Dwork-Evans.

## Lemma 3 (reflection, mod `p`)

```
F(z) ≡ (-1)^n z^n F(1/z),        h(z) ≡ -(-1)^n z^n ( h(1/z) + H_n F(1/z) ).
```
Proof: substitute `k -> n-k`, using (M4) and `H_(2n-2k) - H_(n-k) ≡ H_k - H_(2k) - H_n` (by (M1) and (M2)).

## Theorem, part (ii): the CM harmonic law at the self-dual points

Write `L(m)` for the law `6 h(64/m) ≡ q_p(m) F(64/m) (mod p)`.

**`L(m)` is equivalent to `L(4096/m)`.** Use `q_p(4096/m) ≡ 12 q_p(2) - q_p(m)`. If `L(4096/m)` holds, Lemma 3 and (M3)
give
```
6 h(z) ≡ -(-1)^n z^n (q_p(4096/m) - 12 q_p(2)) F(1/z) = q_p(m) (-1)^n z^n F(1/z) ≡ q_p(m) F(z).
```

**`m = -64` (`z = -1`).** Lemma 3 gives `2 h(-1) ≡ -H_n F(-1) = 2 q_p(2) F(-1)`. That is `L(-64)`, since
`q_p(-64) = q_p(64) ≡ 6 q_p(2)`.

**`m = 64` (`z = 1`).**
- For even `n`, the same argument applies.
- For odd `n`, Lemma 3 is empty for `h`. Instead, `h(1) ≡ (S_n - T_n - H_n D_n)/2 = S_n ≡ 0`, and `F(1) ≡ 0` by (d).

**Twin form.** By (S5),
```
T_m := Σ_(j<=n) C(2j,j) C(2n-j,j)^2 m^(-j) ≡ F_m - 2p h(64/m)   (mod p^2),
```
with `F_m` the mod-`p^2` sum. So `L(m)` is equivalent to `T_m ≡ (1 - p q_p(m)/3) F_m (mod p^2)`. At `m = ±64` this is
`4^(p-1) T_(±64) ≡ F_(±64)`. ∎

## Remarks

**Convention for the census** (S43b m7).
- Twins are taken up to the signs `ε ∈ {1, (-1)^n}`.
  - Allowing `±(2/p)`, `±(-2/p)` also adds `[n even] C(n,n/2) C(3n/2,n/2) ≡ (2/p) a_4C(p)` (a CDE corollary) to class 1.
  - Allowing the scalings `4^(±n)`, `16^(±n)` merges classes 1, 2 and 4 (part (i)), and the census then has 56
    classes.
- The S43b referee reproduced exactly 7 classes by a full pairwise union-find on 76 primes `p <= 400`.

**Relation to KLMSY (S43b M1, n2).**
- The CM harmonic law at the non-self-dual points, and `(H)`, are Kibelbek-Long-Moss-Sheller-Yuan 2016 Corollaries 3
  and 4.
- Part (ii) at `m = -64` (`λ = -1`) is also their Corollary 3, but is proved here by pure reflection without CM input.
- `m = 64` (`λ = 1`) lies outside KLMSY.

**The census.** `halfperiod_twin_classes_20261018.py` screens all 648 distinct three- and four-factor products from a
library of 10 binomials. It finds exactly seven twin classes at `p <= 200`:

| class | members | governing form | status |
|---|---|---|---|
| 1 | `b`, `Σ C(n+k,k)^3` | 4C | THM-4670 |
| 2 | (i-a), `C(n,n/2)^2` | 4C | (i-a) here; partner CITED |
| 3 | `Σ C(2k,k)^3`, `T_1` | 7AB, `η(τ)^3 η(7τ)^3` | `⟺ (H)`, PROVED: KLMSY Cor. 4 (and THM-4686) |
| 4 | (i-b), (i-c) | 4C | here |
| 5 | `A`, `(-1)^n s_7` | 4A | THM-4669 |
| 6 | two `4F3`-type sums | 4A, mod `p` only | `⟺ (K)`, OPEN (HYP-9310) |
| 7 | two sums with `C(2n,n+k)` | 7AB, mod `p` only | `⟺ (H)`, PROVED: KLMSY Cor. 4 (and THM-4686) |

**CM points of `Σ C(2k,k)^3 m^(-k)`.** At Sun's eight values, `j = (m-16)^3/m`. The governing forms are:
- `m = 1, 4096`: 7AB (discriminants `-7`, `-28`);
- `m = -8, 64, -512`: 4C (discriminants `-4`, `-4`, `-16`);
- `m = -64`: 8A (discriminant `-8`);
- `m = 16, 256`: `η(2τ)^3 η(6τ)^3` (discriminants `-3`, `-12`). This is not an M24 cycle shape, so the earlier S43
  remark that only M24 eta-products appear holds only for the integer-sum library, not for `m`-weighted sums.

The reflection `m -> 4096/m` swaps the 2-isogenous CM pairs. The self-dual points `m = ±64`, where the law is forced,
are `j = 1728` and `j = 8000`.
