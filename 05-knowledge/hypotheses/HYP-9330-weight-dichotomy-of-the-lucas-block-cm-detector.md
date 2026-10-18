---
id: HYP-9330
title: "Weight dichotomy of the Lucas-block CM detector. Write S43's test (B1') as V_p(x) == u_1 q_p(m) F_p(x) (mod p), with F_p the truncated period and V_p the first-carry Lucas defect. Then: the split-prime part (slope 1) holds at rational CM points and rational cusps in weight 1 and weight 2 alike; the inert-prime part is automatic in weight 2 (rad F_p divides V_p mod p) and impossible in weight 1 (gcd(F_p, V_p) = 1). So the detector returns the rational CM points of weight-2 families but only the rational cusps of weight-1 families, and its ordinary-restricted form returns CM points and cusps in both weights."
status: >
  OPEN, with FINITE-EXACT evidence; repaired after the independent audit (audit_20261018, MISTAKE-638).
  - Pre-registered predictions are in 04-computation/cm_detector_sporadic/PREREGISTRATION.md, each pushed before its
    run. P1, P3, P3' and P7 partly or wholly failed, as recorded there.
  - P2 is PARTLY FAILED in weight 2 (audit M4): the survivor clause passed, but the control clause failed for Apery
    (m = 100: density 0.130) and Domb (m = 999: 0.110).
  - Part (A), the reformulation, is PROVED below (elementary).
    - It assumes the Lucas congruence u_{p+j} == u_1 u_j (mod p), verified here for every family and prime used.
    - Known for these sequences: Gessel 1982 (Apery); Malik-Straub 2016 and Straub arXiv:2301.12248 (sporadic
      sequences).
    - Meaning of V_p: Straub's Thm 1.3 (all 15 sporadic sequences, which includes Zagier A-F, s_7 and s_10) gives
      A(pn+k) == A(k)A(n) + p n A'(k)A(n) (mod p^2). Hence V_p == u_1 sum_{k<p} A'(k) x^k (mod p): a truncated
      formal-derivative series. A274786 is not on Straub's list.
  - Part (B) is PROVED for Zagier A-F and C(2k,k)^2 at every p >= 7 by independently audited
    THM-4711-legendre-wronskian-and-binomial-cube-lucas-carry-divisibility.md.
    Remaining hypergeometric cases below retain FINITE-EXACT status.
    - Author: Zagier A-F at every 7 <= p <= 199 (pre-registered P7).
    - Referee: all 75 primes 7 <= p <= 397, plus four untested weight-1 families at 59 primes 7 <= p <= 293:
      C(2k,k)^2, (3k)!/k!^3, (4k)!/((2k)!k!^2), (6k)!/((3k)!(2k)!k!).
    - Squarefreeness of F_p is Igusa-type. The coprimality to V_p is the new content.
  - Part (C) is PROVED for Apery zeta(3) and C(2k,k)^3 at every p >= 7 by independently audited
    THM-4711-legendre-wronskian-and-binomial-cube-lucas-carry-divisibility.md.
    For the cube the stronger monic gcd(F_p,V_p) = rad(F_p) holds. Apery uses reflected-companion descent.
    The other weight-2 families remain FINITE-EXACT / OPEN.
    - Author: Apery, Domb, s_7, s_10 (P7) and A274786 (P8c) at every 7 <= p <= 199.
    - Referee: all 75 primes up to 397, plus seven untested weight-2 families at 59 primes up to 293: C(2k,k)^3,
      C(2k,k)^2 C(3k,k), C(2k,k)^2 C(4k,2k), C(2k,k)C(3k,k)C(6k,3k), Almkvist-Zudilin, Cooper s_18, Franel ⊙ central.
  - Part (C), second bullet ("F_p is a square except at the Fricke-elliptic singular fibres, which are simple roots
    exactly at inert p"), is PRIOR ART (audit M2): Caruso, Fürnsinn, Vargas-Montoya and Zudilin, arXiv:2510.23298.
    - PROVED there for Apery (Thm 2), Domb (Thm 3) and Almkvist-Zudilin (Thm 4).
    - OBSERVED there (Table 2) for s_7, s_10, A274786, s_18 and others.
    - Their mechanism f(t(x)) = phi(x) h(x)^2 is the one sketched below.
    - The session's post hoc P7 finding rediscovered this.
  - Part (D), the split-prime slope law, is FINITE-EXACT: slope 1 wherever defined, at all CM points and cusps,
    p <= 257 (author and referee).
    - It is SUGGESTED BY, not implied by, Beukers arXiv:2403.03301 v3, Thm 1.4 / 1.13, together with his
      introductory fixed-point heuristic t_0^sigma = t_0.
    - The truncation has a gap (audit m5). Write t^sigma = t^p(1 + p D(t)); then
      V_p(t_0) = u_1 F_p(t_0) D(t_0) - u_1 t_0^p R(t_0), where R is the degree >= p part of F_p D_{<p}. The slope
      law needs both D(t_0) == q_p(m) and R(t_0) == 0.
  - Part (E), the survivor sets, is FINITE-EXACT.
    - Author: 0 < |m| <= 10^5, p <= 200 (P1, P6, P8a).
    - Referee: unchanged for all 15 families to |m| <= 10^7.
  - Part (F), the valuation law, is FINITE-EXACT (post hoc). It was corrected after audit M3: the earlier claim that
    the inert valuation always equals the weight is false at s_10, m = -2.
  Literature:
  - Beukers (v3) proves the split-prime supercongruences at CM values: Thm 1.7 (weight 1) and Thm 1.15/1.16
    (weight 2). He says his method does not reach inert primes (Rem. 1.8(4), 1.17).
  - S43 (HYP-9311(O)) found the mod-p^3 split / mod-p^2 inert pattern for the hypergeometric weight-2 families.
  - Z.-W. Sun conjectures F_p(1/m) == 0 (mod p^2) at inert p for weight-2 CM parameters (e.g. arXiv:1103.4325
    Conj. 3.12). That is a related but different statement.
  - NEW as far as checked (Beukers, KLMSY, CFVZ, Straub, Sun, OEIS): rad F_p | V_p in weight 2,
    gcd(F_p, V_p) = 1 in weight 1, and the resulting weight dichotomy of the mod-p^2 block detector.
session: kind-pasteur-2026-10-18-S2
source: 05-knowledge/results/lucas_block_cm_detector_sporadic_ad_hauptmoduln_20261018.md
scripts:
  - 04-computation/cm_detector_sporadic/slope_probe.py: (D), slopes
  - 04-computation/cm_detector_sporadic/ordinary_detector_and_mechanism.py (+ _M100000 .out, .json): (B), (C), (E), P6, P7
  - 04-computation/cm_detector_sporadic/level5_leeyang_family.py: P8 (A274786)
  - 04-computation/cm_detector_sporadic/defect_valuation_probe.py (+ .out, .json): (F)
  - 04-computation/cm_detector_sporadic/sporadic_cm_detector.py, high_density_scan.py: (E), P1, P2
  - 04-computation/cm_detector_sporadic/audit_20261018/ (a4, a4b, a4c, a6, a6b, a6c, a7, a7b): the independent referee
related: [HYP-9310, HYP-9311, THM-4708, HYP-9331, THM-4686, Beukers arXiv:2403.03301 v3, Caruso-Fürnsinn-Vargas-Montoya-Zudilin arXiv:2510.23298, Straub arXiv:2301.12248, Z.-W. Sun arXiv:1103.4325, Kibelbek-Long-Moss-Sheller-Yuan 2016, Zagier 2009, Malik-Straub 2016, Olson 1974]
---

# HYP-9330 — Weight dichotomy of the Lucas-block CM detector

## Audited partial closure, 2026-10-18

[THM-4711 — Lucas-carry Wronskians and binomial-cube radical](../../01-canon/theorems/THM-4711-legendre-wronskian-and-binomial-cube-lucas-carry-divisibility.md)
proves (B) for all six Zagier families and central-binomial squares, and
(C) for Apéry zeta(3) and central-binomial cubes, at every prime p>=7. The cube has the
stronger equality monic gcd(F_p,V_p)=rad(F_p). The proof keeps the actual
integer carry, all terminal coefficients, and the singular Euler factor.
It uses a precise Clausen parameter tangent, not differentiation in x.

Apéry is closed by a different route: the reflected Franel companion
is invariant under the quadratic deck map, and a pole bound forces its
descent to a polynomial of degree below p. The actual carry recurrence
then identifies it, including the singular-prefactor roots.

**HYP-9330 remains OPEN:** the remaining weight-2 families, other weight-one
hypergeometric cases, slope law and survivor claims retain their stated
evidence status. Each remaining (C) case needs the correct tangent and full
truncation or an audited companion descent with a sufficient degree bound.

## Setting

Let `u_0 = 1, u_1, u_2, ...` be an Apery-like sequence with the Lucas congruence `u_{p+j} == u_1 u_j (mod p)` for
`0 <= j < p`. Put

    F_p(x) = sum_{k<p} u_k x^k,   v_j = (u_{p+j} - u_1 u_j)/p,   V_p(x) = sum_{j<p} v_j x^j   (over F_p),

and, for an integer `m` prime to `p`, `x = 1/m` and `q_p(m) = (m^{p-1} - 1)/p`. By Straub's Thm 1.3, `V_p` is
`u_1` times the truncated formal derivative series `sum_{k<p} A'(k) x^k`. S43's detector (HYP-9310(c)) is

    (B1')   Delta_p(m) := sum_{k=p}^{2p-1} u_k m^{-k} - u_1 m^{-1} sum_{k<p} u_k m^{-k} == 0   (mod p^2).

The *weight* is 1 for the 2nd-order (elliptic-curve) families: Zagier's sporadic A-F and the `2F1` families. It is 2
for the 3rd-order symmetric-square (K3) families: Apery zeta(3), Domb, Cooper `s_7`, `s_10` and
`C(2n,n) b_n` = A274786.

## Statement

- **(A) Reformulation, PROVED.** `(B1')  <=>  V_p(1/m) == u_1 q_p(m) F_p(1/m)  (mod p)`.
  - At an ordinary prime (`F_p(1/m) != 0`) with `q_p(m) != 0`, call `V_p/(u_1 q_p F_p)` the slope.
  - At a non-ordinary prime, (B1') says `V_p(1/m) == 0`.
  - At `m = +-1` we have `q_p = 0`, so the slope is undefined and (B1') says `V_p(+-1) == 0`.
- **(B) Weight 1, inert primes.** For every prime `p >= 7`, `F_p` is squarefree mod `p` and `gcd(F_p, V_p) = 1`. So
  (B1') fails at every non-ordinary prime of every point.
- **(C) Weight 2, inert primes.**
  - For every prime `p >= 7`, `rad(F_p)` divides `V_p` mod `p`. So the non-ordinary half of (B1') holds
    automatically, at CM and non-CM points alike.
  - (PRIOR ART, CFVZ.) `F_p` is a square, up to a factor vanishing exactly at the Fricke-elliptic singular fibres
    when `p` is inert in their CM field.
- **(D) Split primes, both weights.** At a rational CM point, and at a rational cusp, the slope is 1 at every ordinary
  prime where it is defined.
- **(E) Survivors** (`0 < |m| <= 10^5`, `p <= 200`; referee to `10^7`):
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

  - **ordinary-restricted test (B1'_ord):** the survivors are the CM points together with the cusps, in both
    weights. The weight-1 rational CM points are:

    | family | rational CM points |
    |---|---|
    | A | {-4, 2} |
    | B | {3, 9} |
    | C | {-3, 3} |
    | D | none |
    | E | none |
    | F | {6, 12} |

  - These agree with Beukers' App. B tables for A, C, D, E and F. For B, v3 prints `t = 1/9, 1/3` in its twisted
    variable. That is a sign slip, correct in v1/v2; the true set is {3, 9}.
  - Zagier D has no rational CM point: `Gamma_1(5)`, and CM curves over `Q` have no rational 5-torsion
    (Olson 1974).
- **(F) Valuation law (corrected).**
  - At the non-ordinary primes of a CM point, `v_p(Delta_p) >= w`, the weight. This is forced by (B), which makes it
    exactly 1 in weight 1, and by (C) in weight 2.
  - Equality holds generically. There is one systematic exception: `s_10` at `m = -2`, a CM point by `Z[i]`,
    `tau = (2+i)/5`. There `v_p >= 3` at every non-ordinary prime `p <= 400` (39 primes), and `v_p = 3` at every
    ordinary prime.
  - The other 41 weight-2 CM points show only chance-level excess: 27 have no inert prime with `v >= 3`, 10 have one
    and 4 have two, against an expectation `sum 1/p` of about 0.4-0.75.
  - At ordinary primes of CM points and cusps, `v_p >= 2`.
    - `v_p >= 3` at every CM point of Apery and Domb, matching Beukers' mod-`p^3` Thm 1.13 and S43's (O).
    - `s_7`, `s_10` and `C(2n,n) b_n` are three of our five weight-2 families on Beukers' Rem. 1.14(2) exception list.
      That list also contains `1728t`, `s_18` and Franel ⊙ central. For these three families `v_p = 2` at most CM
      points.
    - The exceptions with `v_p = 3` at every ordinary prime are `s_7` at `m = -1, 2` and `s_10` at `m = -2, 1`. In
      weight 1 they are `A: m = 2` and `B: m = 9`.
  - At generic `m`: `v_p = 1` at ordinary primes, apart from coincidences of frequency about `1/p`, and `v_p = w` at
    the (density-0) non-ordinary primes.

## Proof of (A)

1. Since `x^p = m^{-p} = x m^{1-p}` and `m^{p-1} = 1 + p q_p(m)`, we have `x^p == x(1 - p q_p(m)) (mod p^2)`.
2. Hence `sum_{j<p} u_{p+j} x^{p+j} == sum_j (u_1 u_j + p v_j) x^{j+1} (1 - p q_p(m))`, which is
   `== u_1 x F_p(x) + p x V_p(x) - p q_p(m) u_1 x F_p(x)  (mod p^2)`.
3. Subtracting `u_1 x F_p(x)` and dividing by `p x` gives the claim. `□`

## Evidence

| Part | Families, primes | Result |
|---|---|---|
| (B) | Zagier A-F at every 7 <= p <= 199 (pre-registered P7); referee: 75 primes <= 397, plus 4 hypergeometric weight-1 families at 59 primes <= 293 | squarefree and coprime every time |
| (C) | Apery, Domb, s_7, s_10, A274786 at every 7 <= p <= 199 (P7, P8c); referee: 75 primes <= 397, plus 7 more weight-2 families at 59 primes <= 293 | rad F_p divides V_p every time; simple-root pattern as in CFVZ in every case |
| (D) | 8 weight-1 and 4 weight-2 CM points, p <= 257; referee: all CM points and cusps | slope 1 wherever defined; non-CM controls give random slopes |
| (E) | all 11 families, \|m\| <= 10^5; referee: 15 families to 10^7 | survivors as stated; in weight 2 the density >= 0.3 set coincides with the survivors |
| (F) | 75 points, p <= 199; referee p <= 400 at the inert primes of all 42 weight-2 CM points | as corrected above |

## Mechanism: proved cases and remaining heuristic

- **Weight 2 (the F_p part is CFVZ's).** For a symmetric square, the truncated holomorphic period behaves like `f^2`
  (CFVZ: `f(t(x)) = phi(x) h(x)^2`). By Straub's Thm 1.3 the first-carry defect is the truncated derivative series,
  which behaves like `f g`, with `g` the log-partner. Both vanish where `f` does.
- **Weight 1.** THM-4711 proves the precise Wronskian identity for Legendre and all six Zagier families. Its right side is nonzero at every period root; the remaining hypergeometric cases are still heuristic here.
- **Split primes.** At a rational CM point, Beukers' introductory heuristic has the excellent Frobenius lift fixing
  the point (Serre-Tate canonical lift). Truncating his `F(t) == F_p(t) F(t^sigma)` then gives (B1'), provided
  `R(t_0) == 0` (see the status).

## Falsifiers and next tests

- **Breaking (B) or (C).** A prime `p > 397` at which (B) or (C) fails, for any listed family.
- **Wrong survivor sets.** A weight-1 Lucas family with a non-cusp (B1') survivor, or a weight-2 family with a
  rational CM point that fails (B1').
- **The `s_10`, `m = -2` exception.** Explain why this `Z[i]` point has depth 3 at every prime. A natural guess is
  extra symmetry at `tau = (2+i)/5`.
- **Remaining (B) route.** Extend THM-4711 to the three other tested hypergeometric families; the six Zagier cases and Legendre are closed.
- **Remaining (C) route.** Generalize the audited cube deformation or Apéry reflected-companion descent to Domb, Almkvist–Zudilin and the other families, retaining tangent, pole/degree bounds and simple singular roots.
