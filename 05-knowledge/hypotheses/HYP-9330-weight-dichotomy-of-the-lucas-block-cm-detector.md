---
id: HYP-9330
title: "Weight dichotomy of the Lucas-block CM detector. Writing S43's test (B1') as V_p(x) == u_1 q_p(m) F_p(x) (mod p), with F_p the truncated period and V_p the first-carry Lucas defect: the split-prime part holds at rational CM points and rational cusps in weight 1 and weight 2 alike; the inert-prime part is automatic in weight 2 (rad F_p divides V_p mod p) and impossible in weight 1 (F_p squarefree, gcd(F_p, V_p) = 1). So the detector returns the rational CM points of weight-2 families but only the rational cusps of weight-1 families; its ordinary-restricted form returns CM points and cusps in both weights; and at inert primes the defect valuation equals the weight."
status: >
  OPEN, with FINITE-EXACT evidence; pre-registered predictions P1-P8 are in
  04-computation/cm_detector_sporadic/PREREGISTRATION.md, each pushed before its run.
  - Part (A), the reformulation, is PROVED below (elementary). It assumes the Lucas congruence
    u_{p+j} == u_1 u_j (mod p), j < p.
    - That congruence is verified here for every family and prime used (p <= 199).
    - It is known for these sequences: Gessel 1982 (Apery), Malik-Straub 2016 (sporadic sequences); products of Lucas
      sequences are Lucas.
  - Parts (B) and (C), the inert-prime laws, are FINITE-EXACT for every 7 <= p <= 199, pre-registered as P7 and P8c.
    - (B), weight 1, holds for Zagier A-F.
    - (C), weight 2, holds for Apery zeta(3), Domb, Cooper s_7, s_10 and C(2n,n) b_n (A274786).
    - Part of P7 failed: "no simple root" in weight 2 is false. The simple roots sit exactly on the Fricke-elliptic
      singular fibres (post hoc, all p <= 199).
  - Part (D), the split-prime slope law, is FINITE-EXACT: 25 primes per point at 12 points, and P6.
    - It is expected from Beukers (arXiv:2403.03301) Thm 1.4 (weight 1, mod p^2) and Thm 1.13 (weight 2), combined
      with t_0^sigma = t_0 at rational CM points. That combination is an inference, not a statement of the paper.
    - S43's HYP-9311 heuristic (truncated canonical coordinate = Teichmüller at CM points) is the same picture.
  - Part (E), the survivor sets, is FINITE-EXACT for 0 < |m| <= 10^5 and p <= 200, pre-registered P1, P2, P6, P8a.
  - Part (F), the valuation law, is FINITE-EXACT for p <= 199 (post hoc, exploratory).
  Literature:
  - Beukers proves the split-prime supercongruences at CM values for weight 1 (Thm 1.7) and weight 2 (Thm 1.15,
    1.16). He states that his method does not reach inert primes (Remarks 1.8(4) and 1.17).
  - S43 (HYP-9311(O)) found the mod-p^3 split / mod-p^2 inert pattern for the four hypergeometric weight-2 families.
  - No statement of (B), (C) or the weight dichotomy was found (Beukers 2024; KLMSY 2016; OEIS A274786).
session: kind-pasteur-2026-10-18-S2
source: 05-knowledge/results/lucas_block_cm_detector_sporadic_ad_hauptmoduln_20261018.md
scripts:
  - 04-computation/cm_detector_sporadic/slope_probe.py: (D), slopes
  - 04-computation/cm_detector_sporadic/ordinary_detector_and_mechanism.py (+ _M100000 .out, .json): (B), (C), (E), P6, P7
  - 04-computation/cm_detector_sporadic/level5_leeyang_family.py: P8 (A274786)
  - 04-computation/cm_detector_sporadic/defect_valuation_probe.py (+ .out, .json): (F)
  - 04-computation/cm_detector_sporadic/sporadic_cm_detector.py, high_density_scan.py: (E), P1, P2
related: [HYP-9310, HYP-9311, THM-4708, HYP-9331, THM-4686, Beukers arXiv:2403.03301, Kibelbek-Long-Moss-Sheller-Yuan 2016, Zagier 2009, Malik-Straub 2016, Olson 1974]
---

# HYP-9330 — Weight dichotomy of the Lucas-block CM detector

## Setting

Let `u_0 = 1, u_1, u_2, ...` be an Apery-like sequence with the Lucas congruence `u_{p+j} == u_1 u_j (mod p)` for
`0 <= j < p`. Put

    F_p(x) = sum_{k<p} u_k x^k,   v_j = (u_{p+j} - u_1 u_j)/p,   V_p(x) = sum_{j<p} v_j x^j   (over F_p),

and, for an integer `m` prime to `p`, `x = 1/m` and `q_p(m) = (m^{p-1} - 1)/p`. S43's detector (HYP-9310(c)) is

    (B1')   Delta_p(m) := sum_{k=p}^{2p-1} u_k m^{-k} - u_1 m^{-1} sum_{k<p} u_k m^{-k} == 0   (mod p^2).

The *weight* is 1 for the 2nd-order (elliptic-curve) families: Zagier's sporadic A-F and the `2F1` families. It is 2
for the 3rd-order symmetric-square (K3) families: Apery zeta(3), Domb, Cooper `s_7`, `s_10` and
`C(2n,n) b_n` = A274786.

## Statement

- **(A) Reformulation, PROVED.** `(B1')  <=>  V_p(1/m) == u_1 q_p(m) F_p(1/m)  (mod p)`. Call
  `V_p/(u_1 q_p F_p)` the slope at an ordinary prime, where `F_p(1/m) != 0`. At a non-ordinary prime, where
  `F_p(1/m) == 0`, (B1') says `V_p(1/m) == 0`.
- **(B) Weight 1, inert primes.** For every prime `p >= 7`, `F_p` is squarefree mod `p` and `gcd(F_p, V_p) = 1`. So
  (B1') fails at every non-ordinary prime of every point.
- **(C) Weight 2, inert primes.**
  - For every prime `p >= 7`, `rad(F_p)` divides `V_p` mod `p`. So the non-ordinary half of (B1') holds automatically,
    at CM and non-CM points alike.
  - `F_p` is a square except at the Fricke-elliptic singular fibres, where it has simple roots exactly when `p` is
    inert in their CM field.
- **(D) Split primes, both weights.** At a rational CM point, and at a rational cusp, the slope is 1 at every
  ordinary prime.
- **(E) Survivors** (`0 < |m| <= 10^5`, `p <= 200`):
  - **weight 2:** the survivors of (B1') are the nonzero rational CM points. For the level-5 and level-7 families
    this is the complete list (THM-4708).
  - **weight 1:** the survivors are exactly the rational cusps:

    | family | rational cusps |
    |---|---|
    | A | {-1, 8} |
    | B | none |
    | C | {1, 9} |
    | D | none |
    | E | {4, 8} |
    | F | {8, 9} |

  - **ordinary-restricted test (B1'_ord):** counting only ordinary primes, the survivors are the CM points together
    with the cusps, in both weights. The weight-1 rational CM points are:

    | family | rational CM points |
    |---|---|
    | A | {-4, 2} |
    | B | {3, 9} |
    | C | {-3, 3} |
    | D | none |
    | E | none |
    | F | {6, 12} |

    These match Beukers' appendix tables.
  - Zagier D has no rational CM point: `Gamma_1(5)`, and CM curves over `Q` have no rational 5-torsion
    (Olson 1974).
- **(F) Valuation law.**
  - At the non-ordinary primes of a CM point, `v_p(Delta_p) = w` (the weight): exactly 1 in weight 1, forced by (B);
    `>= 2` in weight 2, forced by (C), and `= 2` apart from sporadic coincidences of frequency about `1/p`.
  - At ordinary primes of CM points and cusps, `v_p >= 2`.
    - `v_p >= 3` at every CM point of Apery and Domb, matching Beukers' mod-`p^3` Dwork congruence (Thm 1.13) and
      S43's (O).
    - `s_7`, `s_10` and `C(2n,n) b_n` are exactly Beukers' Remark 1.14 exceptions (only mod `p^2`). There `v_p = 2`
      at most CM points. The exceptions are `s_7` at `m = -1, 2` and `s_10` at `m = -2, 1`, where `v_p = 3` at
      every ordinary prime.
    - In weight 1, `v_p = 2` except at `A: m = 2` and `B: m = 9`, where `v_p = 3`.
  - At generic `m`: `v_p = 1` at ordinary primes, apart from coincidences of frequency about `1/p`, and `v_p = w`
    at the (density-0) non-ordinary primes.

## Proof of (A)

1. Since `x^p = m^{-p} = x m^{1-p}` and `m^{p-1} = 1 + p q_p(m)`, we have `x^p == x(1 - p q_p(m)) (mod p^2)`.
2. Hence `sum_{j<p} u_{p+j} x^{p+j} == sum_j (u_1 u_j + p v_j) x^{j+1} (1 - p q_p(m))`, which is
   `== u_1 x F_p(x) + p x V_p(x) - p q_p(m) u_1 x F_p(x)  (mod p^2)`.
3. Subtracting `u_1 x F_p(x)` and dividing by `p x` gives the claim. `□`

## Evidence

| Part | Families, primes | Result |
|---|---|---|
| (B) | Zagier A-F, all 7 <= p <= 199 | squarefree and coprime every time (P7, pre-registered) |
| (C) | Apery, Domb, s_7, s_10, all 7 <= p <= 199 | rad F_p divides V_p every time (P7); simple-root part divides the singular-fibre polynomial every time |
| (C) | A274786, all 7 <= p <= 199 | same (P8c, pre-registered); simple roots exactly at p inert in Q(sqrt-5) (spot-checked), D = -20 at the Fricke points 22 +- 10 sqrt5 |
| (D) | 8 weight-1 and 4 weight-2 CM points, p <= 257 | slope 1 at every ordinary prime (25 per point); non-CM controls give random slopes |
| (E) | all 11 families, \|m\| <= 10^5 | survivors as stated; the density >= 0.3 set coincides with the CM points (high_density_scan.py) |
| (F) | 75 points, p <= 199 | histograms in defect_valuation_probe.out |

## Heuristic mechanism (not a proof)

- **Weight 2.** For a symmetric square, the truncated holomorphic period behaves like `f^2` and the first-carry defect
  like `f g`, where `g` is the log-partner. Both vanish where `f` does. The simple roots at the Fricke-elliptic
  points come from the 2:1 Fricke quotient being ramified there.
- **Weight 1.** `(F_p, V_p)` behaves like the pair `(f, g)`. Their Wronskian is a unit, so they have no common zero.
  Separately, Igusa gives simple zeros of the Hasse invariant.
- **Split primes.** At a rational CM point the excellent Frobenius lift fixes the point (Serre-Tate canonical lift).
  Beukers' congruence `F(t) == F_p(t) F(t^sigma)`, truncated below `2p`, then becomes (B1').

## Falsifiers and next tests

- **Breaking (B) or (C).** A prime `p > 199` at which (B) or (C) fails, for any listed family.
- **Wrong survivor sets.** A weight-1 Lucas family with a non-cusp (B1') survivor, or a weight-2 family with a
  rational CM point that fails (B1').
- **More families.** Extend to the remaining sporadic weight-2 sequences (`s_18`, Almkvist-Zudilin), the
  hypergeometric `2F1` and `3F2` families, and the Hadamard products of Beukers' list.
- **Proof route for (C).** A Clausen-type identity mod `p`, `F_p(x) == A(x) f_p(y(x))^2`, together with the same
  identity for the log-partner.
- **Proof route for (B).** Identify `V_p` with the truncation of the Frobenius log-solution modulo `F_p`, then use the
  Wronskian.
