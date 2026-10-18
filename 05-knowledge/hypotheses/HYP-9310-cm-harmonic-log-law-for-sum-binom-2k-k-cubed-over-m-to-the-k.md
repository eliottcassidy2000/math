---
id: HYP-9310
title: "CM harmonic law for sum C(2k,k)^3 m^-k and the p-adic CM detector. The law 6 sum_{k<=(p-1)/2} C(2k,k)^3 m^-k (H_2k - H_k) == q_p(m) sum_{k<=(p-1)/2} C(2k,k)^3 m^-k (mod p) at Sun's CM values is a PUBLISHED theorem: Kibelbek-Long-Moss-Sheller-Yuan 2016, Corollary 3. The census twin (H) sum C(2k,k)^3 (H_2k - H_k) == 0 (mod p) is their Corollary 4, and the Lucas-block equivalence is their Lemma 10 / Section 4 (S43b audit M1). THM-4686 gives an alternative proof. Still open here: the class-6 congruence (K) and the 'only if' direction of the CM detector (c)"
status: >
  RESOLVED for the law (CITED: KLMSY 2016, J. Number Theory 164, Cor. 3 and Cor. 4; alternative proof THM-4686;
  m = 64 and m = -64 by THM-4685(ii), m = 64 being outside KLMSY). OPEN: part (K). OBSERVATION / OPEN: the
  "only if" direction of part (c).
  - S43b audit, MAJOR M1. When filed, the law and (H) were presented as open and "pre-registered". KLMSY, a paper
    cited here, contains them; it had not been read beyond its abstract, and the first S43 audit's note n9 had flagged
    it (MISTAKE-632).
  - Pre-registration. The predictions for the law and for the detector (P1, P2) were written in script docstrings, but
    script and output were committed together (8056a055ae; cf98ae5bce), so git cannot verify the order. P1 also
    restated behaviour already observed. The predictions are kept as records, not as evidence of foresight.
  - FINITE-EXACT, session runs:
    - the law at all 8 CM points for p < 2000, and at m = 1, -8, 16 for p < 20000;
    - (K) for p < 20000;
    - generic m fail.
  - FINITE-EXACT, S43b referee (independent methods):
    - the law at all 8 points for all 9590 primes p <= 10^5 and the 40 primes in [199500, 200000];
    - (K) likewise;
    - 195 non-CM m with a pass rate of at most 4.9%.
  - The naive archimedean analogue (h/F = log|m|/6) is REFUTED: h/F is about 0.001-0.02.
session: opus-2026-10-18-S43
source: 05-knowledge/results/halfperiod_twin_classes_cm_log_law_20261018.md
scripts:
  - 04-computation/experiments/halfperiod_twin_classes_20261018.py
  - 04-computation/experiments/halfperiod_twin_mechanisms_20261018.py
  - 04-computation/experiments/cm_harmonic_log_law_20261018.py
  - 04-computation/experiments/cm_points_forms_legendre_20261018.py
  - 04-computation/experiments/cm_harmonic_log_law_extended_20261018.py
  - 04-computation/experiments/cm_law_lucas_block_form_20261018.py (Lucas-block form, F_(p^2) == F_p^2)
  - 04-computation/experiments/cm_detector_lucas_block_scan_20261018.py (CM detector, four families, |m| <= 10^5)
  - 04-computation/experiments/s43b_audit_20261018/ (referee: audit_cm_law, audit_lucas_klmsy, audit_detector, ...)
related: [THM-4685, THM-4686, HYP-9311, THM-4670, Kibelbek-Long-Moss-Sheller-Yuan 2016 (Thm 2, Cor 3, Cor 4, Lemma 10), Z.-H. Sun arXiv:1012.3898 (1.7), Beukers-Vlasenko Dwork crystals III (arXiv:2105.14841), Sun Conjecture 29 and 36 (arXiv:2210.07238), Gessel 1982]
---

# HYP-9310 — the CM harmonic law (KLMSY Corollary 3) and the CM detector

## Statement and status

Let `n = (p-1)/2`, `F_m = Σ_(k<=n) C(2k,k)^3 m^(-k)` and `h_m = Σ_(k<=n) C(2k,k)^3 m^(-k) (H_(2k) - H_k)`. For each
of Sun's CM values `m ∈ {1, 4096, -8, -512, 16, 256, -64, 64}` and every prime `p >= 5`, `p ∤ m`:
```
(L_p)   6 h_m ≡ q_p(m) F_m   (mod p).
```

**Status: PROVED, CITED.**
- For `λ = 64/m ≠ 1` this is Kibelbek-Long-Moss-Sheller-Yuan (J. Number Theory 164 (2016) 166-178;
  arXiv:1210.4489v3), Corollary 3: for CM `E_λ` and almost all `p` with `λ ∈ Z_p`,
  `Σ_(i<=(p-1)/2) C(2i,i)^3 (λ/64)^i (6(H_(2i) - H_i) + ((λ/64)^(p-1) - 1)/p) ≡ 0 (mod p)`. Since
  `q_p(1/m) ≡ -q_p(m)`, this is `(L_p)`.
- `(H)`, the case `m = 1`, `Σ_(i<=(p-1)/2) C(2i,i)^3 (H_(2i) - H_i) ≡ 0 (mod p)`, is their Corollary 4 ("for all
  primes `p > 3`").
- `m = 64` (`λ = 1`, where `E_λ` degenerates) lies outside KLMSY. It is THM-4685(ii).
- `m = -64` (`λ = -1`) is KLMSY Corollary 3, and THM-4685(ii) proves it with no CM input.
- THM-4686 gives an independent derivation of Corollary 3: KLMSY Theorem 2 at both ends of each 2-isogeny, plus the
  elementary mod-`p^2` reflection.

**Crisp equivalent** (S43b audit n12; follows from THM-4685 (R)). For every `m` with `p ∤ m`,
```
L(m)   ⟺   F_m ≡ (-m/p) F_(4096/m)   (mod p^2).
```
So the law says exactly that the CM supercongruence values agree across the 2-isogeny.

**Twin form.** By (S5) of THM-4685 the law is equivalent to
```
T_m = Σ_(j<=n) C(2j,j) C(p-1-j,j)^2 m^(-j) ≡ (1 - p q_p(m)/3) F_m   (mod p^2).
```
For a cube `m = c^3` this reads `c^(p-1) T_m ≡ F_m (mod p^2)`. No prior statement of this "shifted-square" twin was
found (S43b n6). At `m = 1` it is equivalent to KLMSY Corollary 4 through (S5).

**Census connection.** The census twin classes 3 and 7 (THM-4685) are equivalent to `(H)`, hence PROVED via KLMSY
Corollary 4.

## Part (K): the class-6 congruence (OPEN)

```
Σ_(j<=n/2) a_j^3 a_(2j) (H_(2j) - 2H_j - 2H_(4j) + 2H_(n-j) + H_(n-2j)) ≡ 0   (mod p),
```
where `a_k = C(2k,k)/4^k`, so that `a_j^3 a_(2j) = C(2j,j)^3 C(4j,2j)/1024^j`.
- It is equivalent to the class-6 twin, with constant 1 (re-derived by the S43b referee).
- It holds for all primes `p <= 10^5` (referee) and `p < 20000` (session).
- No reflection forces it (`k_reflection_probe`, S43), and no literature statement was found.

## Lucas-block form (CITED: KLMSY Lemma 10 and Section 4)

**Lemma.** For every digit `a >= 0` and `0 <= j <= n`,
```
C(2(ap+j), ap+j) ≡ C(2a,a) C(2j,j) (1 + 2ap(H_(2j) - H_j))   (mod p^2).
```
Here `C(2N, N) = C(2a,a) × (p-adic unit)` exactly for every `a >= 0`, so the earlier hypothesis `2a < p` is unnecessary
(S43b m8). Checked for `p <= 43` (session) and for `0 <= a <= 2p` (referee).

**Block form and equivalence.** With `t_k = C(2k,k)^3 m^(-k)`,
`Σ_(j<=n) t_(ap+j) ≡ t_a m^(-a(p-1)) (F_m + 6ap h_m)`. Hence `(L_p)` is equivalent to
```
(B1)   Σ_(k=p)^(p+n) C(2k,k)^3 m^(-k) ≡ (8/m) Σ_(k<=n) C(2k,k)^3 m^(-k)   (mod p^2).
```

**Credit.** This is KLMSY Lemma 10: `a_(k+pn) ≡ a_k a_n + p b_k Σ_(i<=n) i C(2i,i)^3 (λ/64)^i (mod p^2)` for
`(p-1)/2 <= k < p` and any `n`. Their Section 4 states the same equivalence; the mechanism goes back to Gessel 1982 for
the Apéry numbers.

**Consequences.**
- Every block factors.
- `F_(p^2) ≡ F_p^2 (mod p^2)`, with `F_N = Σ_(k<N) t_k`.
- In general `F_(p^2) - F_p^2 ≡ p (6h_m - q_p(m) F_m) Σ_(a<=n) a t_a`.
- **Partial converse (corrected, S43b m9).** At the CM points, `Σ_(a<=n) a t_a ≡ 0 (mod p)` exactly at the
  supersingular primes (`m ≠ 64`). So `F_(p^2) ≡ F_p^2` implies the law only at ordinary primes. At generic `m`, 5 of
  268 tested `(m, p)` rows have `F_(p^2) ≡ F_p^2` with the law false.

**Checks.** `cm_law_lucas_block_form_20261018.py` (session); referee `audit_lucas_klmsy`: 8352 checks of Lemma 10,
0 failures.

## Part (c): the Lucas-block congruence detects CM (OBSERVATION; the "only if" direction is OPEN)

The test is
```
(B1')   Σ_(k=p)^(2p-1) t_k ≡ t_1 Σ_(k<p) t_k   (mod p^2),     t_k = A(k) m^(-k).
```
It holds modulo `p` for every `m`.

**Scan** (`cm_detector_lucas_block_scan_20261018.py`). Over all `0 < |m| <= 10^5`, with filter primes `7..43`, then
primes up to `200`:

| family | survivors |
|---|---|
| A1 = `C(2k,k)^3` | `-512, -64, -8, 1, 16, 64, 256, 4096` (Sun's eight) |
| A2 = `C(2k,k)^2 C(3k,k)` | `-8640, -1728, -192, -27, 8, 64, 108, 216, 1458, 3375` |
| A3 = `C(2k,k)^2 C(4k,2k)` | `-82944, -12288, -3969, -1024, -144, 81, 256, 648, 2304, 20736` |
| A4 = `C(2k,k) C(3k,k) C(6k,3k)` | `-32768, -3375, 1728, 8000, 54000` (the rational CM `j` in range) |

**Independent confirmation** (S43b referee, by a different method).
- The survivor lists were reproduced exactly, with re-confirmation on primes up to `600`.
- A supersingular density `>= 0.30` (primes `101..1200`) occurs exactly at the survivors (`0.49-0.54`). The maximum
  among non-survivors is `0.082-0.112`, so the lists are complete in range relative to the CM signature.
- `τ = j^(-1)(j)` is imaginary quadratic for every survivor. The discriminants are:
  - A1: `-16, -8, -4, -7, -3, -4, -12, -28`;
  - A2: `-75, -51, -27/-3, -15, -8, -11, -12, -24, -48, -60`;
  - A3: `-52, -36, -28/-7, -20, -12/-3, -7, -8, -16/-4, -24, -40`.
- A1 completeness is also algebraic: an integer root of `(m-16)^3 = jm` divides `4096`.

**Context.** For A1 the "if" direction is KLMSY Corollary 3, through the equivalence above. For A2-A4 the law is
equivalent to `(B1')` (S43b m3; see HYP-9311).

**Conjecture (c), OPEN in the "only if" direction.** `(B1')` holds for all large `p` only at CM parameters.

## Context and heuristics

- **Isogeny.** The reflection `k -> n-k` acts as `m -> 4096/m`. It swaps the 2-isogenous CM pairs:
  - discriminants `-7 ↔ -28`, `-4 ↔ -16`, `-3 ↔ -12`;
  - `E_z` is 2-isogenous over `Q` to the `-2` twist of `E_(1/z)`, so `α^2` agrees (S43b m2).
  The self-dual points are `m = ±64`.
- **Legendre squares** (CITED).
  - The identity `P_N(x)^2 = Σ_k C(N+k,2k) C(2k,k)^2 ((x^2-1)/4)^k` is Z.-H. Sun's (1.7) (arXiv:1012.3898).
  - `F ≡ P_((p-1)/2)(sqrt(1-λ))^2 (mod p^2)` is in KLMSY Section 3.
  - Re-verified: the identity for `N <= 40`, and the congruence at the 8 points for `p <= 300`.
- **Relation to Sun.** Z.-W. Sun's Conjecture 29(i) (arXiv:2210.07238) is the `k`-weighted Ramanujan-type template.
  `(L_p)` is its constant-term analogue, KLMSY Corollary 3. Sun does not cite KLMSY (S43b n7).
- **Beukers-Vlasenko** (attribution corrected, S43b m5).
  - Dwork crystals III identifies CM fixed points of the excellent lift only for the `n = 2` hypercubic family.
  - For `n = 3`, which is our 3F2, they conjecture only `t^σ(1/8) = 1/8`, that is `m = 64`.
  - "CM points = fixed points" for `n = 3` is this repository's HEURISTIC.
  - In the coordinate `m = 1/t^2`, `t ∈ Q_p` only when `m` is a square in `Q_p`.

## Falsifiers

- A prime `p ∤ 6m` where `(K)` fails.
- A non-CM parameter satisfying `(B1')` for all large `p`.
