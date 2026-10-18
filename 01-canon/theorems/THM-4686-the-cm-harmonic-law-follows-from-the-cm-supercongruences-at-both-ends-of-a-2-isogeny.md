---
id: THM-4686
title: "An alternative proof of Kibelbek-Long-Moss-Sheller-Yuan's Corollary 3 (the CM harmonic law 6 sum_{k<=(p-1)/2} C(2k,k)^3 m^-k (H_2k - H_k) == q_p(m) sum_{k<=(p-1)/2} C(2k,k)^3 m^-k mod p at the CM values m = 64/lambda). The proof applies their Theorem 2 at BOTH ends of each 2-isogeny m <-> 4096/m and joins the two values by the elementary mod-p^2 reflection F(z) == (-1)^n z^n [F(1/z) + 3p G(1/z)]. Equivalently: the law at m holds iff F_m == (-m/p) F_(4096/m) mod p^2. Corollaries: census twin classes 3 and 7, the twin form T_m, and the block factorization; KLMSY Conjecture 9 at s = 1 for every odd m at the seven rational CM lambda != 1 (immediate also from KLMSY Lemma 10 + Corollary 3; not stated by KLMSY)"
status: >
  PROVED, given CITED published theorems: KLMSY 2016 Theorem 2 at both ends, and THM-4685's reflection.
  - Calibration (S43b audit, MAJOR M1). The law itself is PUBLISHED: Kibelbek-Long-Moss-Sheller-Yuan, J. Number Theory
    164 (2016) 166-178, arXiv:1210.4489v3, Corollary 3 (lambda != 1). (H) is their Corollary 4, and the block
    equivalence is their Lemma 10 / Section 4. When filed, this theorem claimed to "upgrade an open hypothesis"; that
    was wrong, because the paper had been cited but not read beyond its abstract (MISTAKE-632).
  - What is offered here:
    - an alternative derivation of Corollary 3 from Theorem 2 alone, without KLMSY's Lemma 10 / Coster-van Hamme
      route;
    - explicit prime coverage: all p >= 5, with p = 7 at m = 1, 4096 checked directly;
    - the transparent equivalent L(m) <=> F_m == (-m/p) F_(4096/m) (mod p^2) (S43b n12).
  - The S43b referee also notes (n4) that KLMSY's written proof of Corollary 3, read literally, has a factor -8 slip in
    the Legendre rewriting at (3p-1)/2. Their conclusion stands; THM-4686 is a clean derivation.
  - KLMSY v3 prints the character as ((lambda-1)/p) in Theorems 1-2, Proposition 8 and Conjecture 9. The correct
    character is ((1-lambda)/p) (S43b m1: the printed version fails at 217 of 1148 (lambda, p) pairs; the corrected one
    at none). The proof below uses only the ratio of the two characters, which the correction leaves unchanged.
  - m = 64 (lambda = 1, degenerate for KLMSY's family) is THM-4685(ii), proved by reflection without van Hamme /
    Ishikawa / Ahlgren. m = -64 (lambda = -1) is also THM-4685(ii), with no CM input.
  - FINITE-EXACT controls:
    - the reflection (R) and the law, p <= 200 (session) and p < 1000 (referee);
    - the law at all 8 points for p <= 10^5 (referee);
    - the Conjecture-9 corollary for M <= 11, p <= 150 (session), and for every odd M <= 2p+3, p <= 61 (referee).
session: opus-2026-10-18-S43
source: 05-knowledge/results/halfperiod_twin_classes_cm_log_law_20261018.md
scripts:
  - 04-computation/experiments/cm_law_from_known_supercongruences_20261018.py (+ .out, .json): (R) and the law directly
  - 04-computation/experiments/cm_points_forms_legendre_20261018.py (+ .out, .json): the eight supercongruence values (FINITE-EXACT, p <= 200)
  - 04-computation/experiments/cm_harmonic_log_law_20261018.py, cm_harmonic_log_law_extended_20261018.py: the law directly
  - 04-computation/experiments/klmsy_conjecture9_s1_check_20261018.py (+ .out, .json): the Conjecture-9 (s = 1) corollary, M <= 11, p <= 150
  - 04-computation/experiments/s43b_audit_20261018/ (referee: audit_lucas_klmsy, audit_reflection_equiv, audit_conj9_range)
related: [THM-4685, HYP-9310, HYP-9311, THM-4670, Kibelbek-Long-Moss-Sheller-Yuan 2016 (Thm 2, Cor 3, Cor 4, Lemma 10, Conj 9), Long-Ramakrishna 2016]
---

# THM-4686 — an alternative proof of KLMSY Corollary 3 through the 2-isogeny

## Notation

`p >= 5` is prime, `n = (p-1)/2`, `a_k = C(2k,k)/4^k` and `O_k = Σ_(j<=k) 1/(2j-1)`. For `z ∈ Z_(p)^×` write:
- `F(z) = Σ_(k<=n) a_k^3 z^k`, so that `F(64/m) = Σ_(k<=n) C(2k,k)^3 m^(-k)`;
- `G(w) = Σ_(k<=n) a_k^3 w^k Δ_k`, with `Δ_k = O_(n-k) - O_k`;
- `h(z) = Σ_(k<=n) a_k^3 z^k (H_(2k) - H_k)`.

The law at `m` reads `L(m): 6 h(64/m) ≡ q_p(m) F(64/m) (mod p)`.

## Theorem

**(a)** For every `m` with `p ∤ m`, `L(m) ⟺ F_m ≡ (-m/p) F_(4096/m) (mod p^2)`.

**(b)** `L(m)` holds for every `m ∈ {1, -8, 16, -64, 64, 256, -512, 4096}` and every prime `p >= 5`, `p ∤ m`.
- For `λ = 64/m ≠ 1` this is KLMSY Corollary 3.
- Here it is derived from KLMSY Theorem 2 and (a).

**Caveat at `p = 7`.** For `m ∈ {1, 4096}` and `p = 7`, the CM curve has bad reduction, which is outside KLMSY
Theorem 2. On KLMSY's model `Δ = 64 λ^2/(1-λ)^3`, the bad primes are:
- `{2, 3, 7}` for `λ ∈ {64, 1/64}`;
- `{2, 3}` for `λ ∈ {4, 1/4, -8, -1/8}`;
- `{2}` for `λ = -1`.

There the statement is verified directly (FINITE-EXACT), as are the input supercongruences.

## Proof

**Step 1: the reflection mod `p^2`.** THM-4685 (S6) gives `a_(n-k) ≡ (-1)^n a_k (1 + pΔ_k) (mod p^2)`. Substituting
`k -> n-k`:
```
(R)   F(z) ≡ (-1)^n z^n [ F(1/z) + 3p G(1/z) ]   (mod p^2).
```

**Step 2: (a).**
- Put `z = 64/m` and use `z^n ≡ (z/p)(1 + p q_p(z)/2) (mod p^2)`.
- (R) shows `F(z) ≡ (-1)^n (z/p) F(1/z) (mod p^2)` if and only if `G(1/z) ≡ -(q_p(z)/6) F(1/z) (mod p)`.
- By THM-4685 (M1)-(M3), `H_(2k) - H_k ≡ -Δ_k + q_p(2)`. So the latter is `h(1/z) ≡ (q_p(64z)/6) F(1/z)`, which is
  `L(4096/m)`. By THM-4685(ii) that is equivalent to `L(m)`.
- Finally, `(-1)^n (z/p) = (-1/p)(m/p) = (-m/p)`.

**Step 3: the input (KLMSY Theorem 2, corrected character).** For `z ∈ {64, -8, 4}` (`m = 1, -8, 16`) and `1/z`
(`m = 4096, -512, 256`):
```
F(z) ≡ ((1-z)/p) α^2,     F(1/z) ≡ ((1 - 1/z)/p) α'^2   (mod p^2).
```
- `α^2 = α'^2`, because `E_z` is 2-isogenous over `Q` to the quadratic twist by `-2` of `E_(1/z)`. This was verified as
  `a_p(E_z) = (-2/p) a_p(E_(1/z))` for `7 < p < 1000`, all three pairs (S43b m2).
- `((1 - 1/z)/p) = ((z-1)/p)(z/p) = (-1)^n ((1-z)/p)(z/p)`. So
  `F(z) ≡ (-1)^n (z/p) F(1/z) (mod p^2)`, and (a) gives `L(m)` and `L(4096/m)`.
- At supersingular `p` both values vanish mod `p^2`.
- The self-dual points `m = ±64` are THM-4685(ii). ∎

The eight values were also confirmed numerically mod `p^2`, for `p <= 200` (session) and `p <= 1000` (referee): the
forms 7AB, 4C, `η(2τ)^3 η(6τ)^3` and 8A with the stated characters.

## Corollaries

- **Class 3** (census, THM-4685): `Σ_(j<=n) C(2j,j) C(p-1-j,j)^2 ≡ Σ_(j<=n) C(2j,j)^3 (mod p^2)`. Both sides are
  `≡ a(p)` for `η(τ)^3 η(7τ)^3`. This is equivalent to KLMSY Corollary 4 through (S5).
- **Class 7:** `Σ C(2n-2k,n-k)^2 C(n+k,2k) C(2n,n+k) ≡ Σ C(n+k,2k)^3 C(2n,n+k) (mod p^2)`.
- **Twin form** at every CM point: `T_m = Σ_j C(2j,j) C(p-1-j,j)^2 m^(-j) ≡ (1 - p q_p(m)/3) F_m (mod p^2)`. No prior
  statement was found (S43b n6).
- **Block factorization** (KLMSY Lemma 10 with Corollary 3): `Σ_(k=p)^(p+n) C(2k,k)^3 m^(-k) ≡ (8/m) F_m (mod p^2)`,
  and `F_(p^2) ≡ F_p^2 (mod p^2)`.

## Corollary: KLMSY Conjecture 9 at `s = 1`

**Statement.** KLMSY's Conjecture 9 (character corrected) reads, for CM `λ` of `E_λ: y^2 = (x-1)(x^2 - 1/(1-λ))` with
`λ ∈ Z_p`, ordinary `p`, odd `M` and `s >= 1`:
```
F_((Mp^s - 1)/2)(λ) ≡ ((1-λ)/p) α_(p,λ)^2 F_((Mp^(s-1) - 1)/2)(λ)   (mod p^(2s)).
```
For the seven rational CM `λ ≠ 1` and `s = 1`, every odd `M` holds.

**Proof.**
- The `p`-shift lemma holds for every digit `a >= 0`, since `C(2N,N) = C(2a,a) × unit`.
- By the law, the block with digit `a` is `≡ t_a F_m (mod p^2)`.
- Terms with last digit `n < j < p` carry, and vanish mod `p^3`.
- Hence `(D) F_((Mp-1)/2) ≡ F_((p-1)/2) F_((M-1)/2) (mod p^2)`. Theorem 2 finishes.

**Status.**
- It is also immediate from KLMSY Lemma 10 (`k = (p-1)/2`, `n = (M-1)/2`) together with Corollary 3 and Theorem 2
  (S43b M1). KLMSY do not state it.
- No later proof of Conjecture 9, for any `s`, was found (S43 and S43b searches).
- The statement also holds at supersingular `p`, where both sides vanish mod `p^2`.
- At `λ = 1`, Conjecture 9 is not formulated, but (D) still holds.

**Checks.**
- (D) at all eight CM `m`: `M <= 11`, `p <= 150` (session); every odd `M <= 2p+3`, `p <= 61` (referee: 4224 checks).
- At generic `m`, (D) fails in 1581 of 2028 cases.

## Remarks

- **The pattern.** In THM-4670 (audit M1) and here, a "new" congruence turned out to be the shadow of two classical
  supercongruences joined by an elementary `p`-adic identity:
  - THM-4670: the two guises at one CM point;
  - THM-4686: the two ends of a 2-isogeny.
- **The other families.** For `C(2k,k)^2 C(3k,k)`, `C(2k,k)^2 C(4k,2k)` and `C(2k,k) C(3k,k) C(6k,3k)` no mod-`p`
  reflection `m -> c^2/m` exists (`family_reflection_probe_20261018.out`), so this route does not apply (HYP-9311).
