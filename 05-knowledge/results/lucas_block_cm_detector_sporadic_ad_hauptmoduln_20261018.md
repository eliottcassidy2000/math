# The Lucas-block CM detector on sporadic families: the weight dichotomy, Fricke Hauptmoduln with an Argyres-Douglas reading, and attractor points

kind-pasteur-2026-10-18-S2. This note builds on opus-2026-10-18-S43 (HYP-9310, HYP-9311 and the bridge note's
§4.3–4.4) and on kind-pasteur S1 (Argyres-Douglas chiral algebras, HYP-9297).

- **Code and data:** `04-computation/cm_detector_sporadic/`.
- **Pre-registration:** every prediction was written into `PREREGISTRATION.md` and pushed before its run: P1–P3 in
  8bf465cb8c, P3'/P4/P5 in 2583387d46, P6/P7 in 4a60d0195a, P8 in 543bc0f3ba and P9 in 1846da553b. Outcomes were
  appended later. The file is append-only; the referee verified the push reflog.
- **Audit:** an independent referee (`audit_20261018/REPORT.md`; MAJOR 4, MINOR 6, NOTE 10) confirmed all the
  mathematics. It corrected attributions and two reporting errors (MISTAKE-638). This version incorporates the
  repairs.

## 1. Summary

1. **Weight dichotomy (HYP-9330).**
   - S43's mod-`p^2` Lucas-block test (B1') detects CM points only in weight 2.
   - On Zagier's six sporadic weight-1 families, its survivors are exactly the rational cusps.
   - The reason is located exactly. (B1') is equivalent to `V_p(x) == u_1 q_p(m) F_p(x) (mod p)`, where `V_p` is the
     first-carry Lucas defect (Straub's Thm 1.3 identifies it as a truncated derivative series).
     - At split primes this slope-1 law holds at CM points in both weights.
     - At inert primes it is automatic in weight 2: `rad F_p | V_p` at all 43 primes `7 <= p <= 199`, in 5 families.
       The referee extended this to 75 primes up to 397 and seven more families.
     - It is impossible in weight 1: `gcd(F_p, V_p) = 1` at all 43 primes, in 6 families, plus 4 more (referee).
   - Restricting the test to ordinary primes therefore makes it a CM detector in weight 1 too (P6, exact).
   - Not new: the fact that `F_p` is a square, except at the singular fibres when they are inert, is due to Caruso,
     Fürnsinn, Vargas-Montoya and Zudilin (2025).
2. **Rational CM values of `T_5A + 16` and `T_7A + 9` (THM-4708).**
   - These Fricke Hauptmoduln parametrize `C(2n,n) b_n` (A274786) and Cooper's `s_7`.
   - Their rational CM values: ten each, all integers.
   - The nonzero ones were tabulated by Beukers (2024, App. C) and, at level 5, by Z.-W. Sun (2011). The completeness
     proof is routine.
   - New: they are exactly the detector's survivors (to `|m| <= 10^7`), so the detector misses no CM point at any
     height.
   - Identification of coordinates: `m_N` is built from the product of all characters of the Virasoro `(2,2K+1)`
     model, the chiral algebra of the `(A_1,A_{2K-2})` Argyres-Douglas theory.
3. **Attractor points (HYP-9331).**
   - In weight 3, a blind pre-registered scan of four CY3 families returns exactly the two known rational rank-2
     attractor points (Hulek-Verrill `-7`; `4F3(1/3,1/3,2/3,2/3)` at `-5832`) plus the singular points.
   - The referee added two out-of-sample confirmations: AESZ 11 at `-432` and AESZ 36 at `-64`.
   - The defect valuation matches the modulus of Beukers' conjectured supercongruences.

## 2. Pre-registered tests

| # | prediction | outcome |
|---|---|---|
| P1 | every family has a nonempty finite survivor set | FAILED for Zagier B and D (no survivors); PASSED for the rest |
| P2 | survivors have density 0.35–0.65, controls < 0.10 | FAILED in weight 1 (survivors have density 0); PARTLY FAILED in weight 2: the survivor clause passed, but the control clause failed for Apery (m = 100: 0.130) and Domb (m = 999: 0.110) (audit M4) |
| P3 | Zagier D survivors are CM | VACUOUS (no survivors) |
| P3' | Zagier D has rational m with density >= 0.3, all CM | FAILED: the set is empty. Post hoc: the family lives on `Gamma_1(5)` and Olson's theorem forbids rational CM fibres; ill-posed (MISTAKE-636) |
| P4 | the 7 nonsingular `s_7` survivors are CM | PASSED (D = -427, -91, -35, -12, -19, -27, -112); post hoc, the singular fibres -1, 27 have D = -7, -28 |
| P5 | `s_7` controls 3, 5, 11, 100 are not CM | PASSED |
| P6 | (B1'_ord) survivors in weight 1 = cusps ∪ CM, as listed | PASSED exactly, all six families |
| P7 | weight 1: `F_p` squarefree and `gcd(F_p, V_p) = 1`; weight 2: `rad F_p \| V_p` and no simple roots | PASSED except "no simple roots" (FAILED: the simple roots are the Fricke-elliptic singular fibres, at the inert primes, out of all 43 primes; that pattern is CFVZ's, audit M2) |
| P8 | Lee-Yang weight-2 family: survivors = density set, all CM; `rad F_p \| V_p`; simple part divides `1-44x-16x^2` | PASSED: 9 survivors, D = -235, -115, -35, -15, -11, -16, -19, -40, -60 |
| P9 | weight 3: survivors include -7 (Hulek-Verrill) and -5832 (`(3k)!^2/k!^6`); every survivor is an attractor or singular point | PASSED |

One label was wrong in the pre-registration. The formula `(t^2+13t+49)(t^2+5t+1)^3/t` gives `j(7tau)`, not
`j(tau)`. This does not affect the test, since `tau` is CM iff `7tau` is.

## 3. The weight dichotomy (HYP-9330)

- **Data, weight 1.**
  - Survivors of (B1') with `|m| <= 10^5`: A {-1, 8}, B {}, C {1, 9}, D {}, E {4, 8}, F {8, 9}. These are the
    rational cusps, all with supersingular density 0.
  - Density >= 0.3 set: A {-4, 2}, B {3, 9}, C {-3, 3}, D {}, E {}, F {6, 12}. These are disjoint from the
    survivors.
  - The density sets agree with Beukers' App. B rational CM tables for A, C, D, E and F. For B, v3 prints the twisted
    variable with the wrong sign (v1/v2 are correct).
- **Data, weight 2.** The survivors equal the density set, for Apery, Domb, `s_7`, `s_10` and A274786. The referee
  confirmed this for 15 families to `|m| <= 10^7`.
- **Slopes.** At every CM point and cusp, `V_p/(u_1 q_p F_p) = 1` at every ordinary prime where it is defined
  (`q_p(+-1) = 0`).
  - At the supersingular primes of weight-1 CM points, `V_p != 0` every time.
  - In weight 2, `V_p == 0` at every supersingular prime, CM or not.
  - Beukers v3 Thm 1.4 / 1.13 with his fixed-point heuristic *suggests* the split-prime law, but there is a
    truncation gap (audit m5).
- **Polynomial form.**
  - Weight 1: `F_p` is squarefree and coprime to `V_p`, every time.
  - Weight 2: `rad F_p` divides `V_p`, every time.
  - The square/singular-fibre structure of `F_p` itself is in CFVZ: Thms 2–4, proved for Apery, Domb and
    Almkvist-Zudilin; Table 2, observed for `s_7`, `s_10`, A274786 and `s_18`.
- **Valuations (corrected).**
  - Inert primes: `v_p(Delta) = 1` in weight 1 (forced) and `>= 2` in weight 2.
    - In weight 2 it is generically 2. The systematic exception is `s_10` at `m = -2` (CM by `Z[i]`), where
      `v_p >= 3` at all 39 inert primes up to 400.
  - Split primes: `v_p >= 2`.
    - It is 3 for Apery and Domb (cf. S43's (O)).
    - It is mostly 2 for `s_7`, `s_10` and A274786. These are three of our five weight-2 families on Beukers'
      Rem. 1.14(2) exception list, which also contains `1728t`, `s_18` and Franel ⊙ central.
  - Generic `m`: 1, apart from coincidences of frequency about `1/p`.

## 4. Argyres-Douglas characters and the Hauptmoduln (THM-4708 parts 1-2)

This is an identification of modular coordinates, not physics. It is exact to `q^59` here (`ad_character_hauptmoduln.py`)
and to `q^300` / `q^1499` by the referee.

| AD theory | chiral algebra | character identity | family |
|---|---|---|---|
| `(A_1,A_2)` | Lee-Yang (2,5) | `(chi_0/chi_{-1/5})^5 = q prod (1-q^n)^{5(n/5)}` | Zagier D mirror map (`Gamma_1(5)`, Beukers 1987) |
| `(A_1,A_2)` | Lee-Yang (2,5) | `(chi_0 chi_{-1/5})^6 = (eta5/eta1)^6`; `x = H/(1+22H+125H^2)` | `C(2n,n) b_n` = A274786 (`Gamma_0(5)+`; Beukers App. C §9.14, Chan-Cooper 2012) |
| `(A_1,A_4)` | (2,7) | `(chi_1 chi_2 chi_3)^2 = (eta7/eta1)^4`; `x = h/(1+13h+49h^2)` | Cooper `s_7` (`Gamma_0(7)+`; Cooper 2012) |
| `(A_1,A_6)` | (2,9) | `chi_1 chi_2 chi_3 chi_4 = (eta9/eta1)^3`; `t_B(q) = -P(-q)` | Zagier B mirror map (`Gamma_0(9)`, twisted) |

- **Which levels qualify.** Every odd `N` is `2K+1`. The operative condition for an eta-quotient Hauptmodul is
  `(N-1) | 24`, which leaves the odd levels 5, 7, 9, 13 and 25. Sporadic Apery-like families are known at 5, 7 and 9.
- **The general product identity.** `prod_s chi_s^{(2,2K+1)} = (eta((2K+1)tau)/eta(tau))^{K-1}` is folklore-level
  (for `K = 2` it is Ramanujan's `G(q)H(q)`).
- **Physics context.** CM `tau` is the necessary half of the Gukov-Vafa rational-CFT condition for `E_tau` (S43 §4.4).
- **Lee-Yang CM points.** They are invisible on `Gamma_1(5)`, where no rational CM point exists (Olson). They appear on
  the Fricke quotient: the weight-2 family A274786 has nine.

## 5. The rational CM values (THM-4708 part 3)

- **Classification.** Rational `m` gives `[Q(j(N tau)):Q] <= 2`, hence `h(D) <= 2`, hence one of 42 orders. Exact
  resultants then give:
  - `m_7`: -10648, -64, -8, -1, 0, 2, 8, 24, 27, 125;
  - `m_5`: -15228, -828, -28, -3, 0, 4, 18, 36, 72, 147.
- **Independent check.** The referee reproduced both lists and all 20 D labels by direct eta-product evaluation at
  300 digits.
- **Prior tabulation.** The nine nonzero values at each level are in Beukers v3 App. C §9.11/§9.14, with `tau`, without
  a completeness claim, and with the misprint 1/588 for 1/147. The level-5 values are also in Z.-W. Sun, arXiv:1103.4325,
  §3 (2011).
- **Completeness.** It is the standard description of rational CM points on `X_0^+(N)` (Galbraith 1999).
- **Integrality.** Chen-Yui (Thompson series at CM points are algebraic integers).
- **Huber-Schultz-Ye 2023.** It claims a complete list of rational and quadratic series at singular values. Its body
  was not accessible, so the comparison is open.

## 6. Weight 3 and attractor points (HYP-9331)

| family | survivors | attractor | source | `v_p` at attractor |
|---|---|---|---|---|
| Hulek-Verrill (AESZ 34) | {-7, 1, 9, 25} | -7 | Candelas-de la Ossa-Elmi-van Straten 2020; Beukers Conj. 1.27, mod `p^2` | 2 |
| `((3k)!/k!^3)^2` | {-5832, 729} | -5832 | Boenisch-Elmi-Kashani-Poor-Klemm 2022, Table 1; Beukers Conj. 1.28, mod `p^3` | 3 |
| AESZ 11 (referee, out of sample) | {-432, 1728} | -432 | Boenisch-Elmi-Kashani-Poor-Klemm Table 1 | 3 |
| AESZ 36 (referee, out of sample) | {-64, 64, 128} | -64 | Boenisch-Elmi-Kashani-Poor-Klemm Table 1 | 2 |

- **Controls.** `C(2k,k)^4` gives {256}; the quintic gives {3125}.
- **Post hoc.** The ordinary-restricted test also passes `C(2k,k)^4` at `z = -1`. Its non-ordinary primes are exactly
  `p == 3 (mod 4)`, and `v_p = 1` there except at `p = 127`. It behaves like a CM-by-`Q(i)` point.

## 7. What is new, as far as checked

- **The inert-prime laws.** `rad F_p | V_p` in weight 2 and `gcd(F_p, V_p) = 1` in weight 1, with the weight dichotomy
  of S43's detector that they explain.
- **The ordinary-restricted detector.** A CM detector in weight 1.
- **Detector completeness.** At levels 5 and 7, the detector survivors equal all nonzero rational CM values.
- **Attractor points.** The detector finds rank-2 attractor points in blind scans: 2 pre-registered, plus 2 by the
  referee out of sample.
- **The Argyres-Douglas packaging.** An identification of coordinates.
- **Not new:**
  - `F_p` = square times singular factor at inert primes: CFVZ 2025.
  - `V_p` as a truncated derivative series: Straub 2023.
  - The level-5/7 CM values: Beukers 2024 App. C; Sun 2011.
  - The parametrizations: Cooper 2012; Chan-Cooper 2012; Beukers.
  - Completeness (routine; Galbraith) and integrality (Chen-Yui).
  - The character product (folklore).
  - The split-prime supercongruences (Beukers; KLMSY).
  - S43's detector and harmonic law (KLMSY Cor. 3).
  - A274786 and Bala's product formula.

## 8. Open problems

1. Prove (B) and (C). For (C), combine CFVZ's factorization with Straub's derivative description of `V_p`. For (B),
   use a Wronskian argument mod `p`.
2. Scan the remaining 21 attractor points of Boenisch-Elmi-Kashani-Poor-Klemm Table 1, and operators without known
   attractors (a blind attractor search).
3. Level 9 `(A_1,A_6)`: the Fricke weight-2 family and its rational CM values; then levels 13 and 25.
4. Explain the extra depth:
   - `v_p = 3` at the ordinary primes of `s_7` `m = -1, 2`, `s_10` `m = -2, 1`, Franel `m = 2` and Zagier B `m = 9`;
   - `v_p >= 3` at all inert primes of `s_10` `m = -2` (CM by `Z[i]`, `tau = (2+i)/5`).

## 9. Reproduction (from `04-computation/cm_detector_sporadic/`)

| script | what it does | run time |
|---|---|---|
| `python sporadic_cm_detector.py 100000` | P1/P2 | minutes |
| `python cm_archimedean.py` | P3', P4, P5 | about 10 min |
| `python high_density_scan.py 100000` | density sets | |
| `python slope_probe.py '<json>'` | slopes | |
| `python ordinary_detector_and_mechanism.py 100000` | P6/P7 | |
| `python level5_leeyang_family.py` | P8 | |
| `python integral_cm_classification.py` | classification | |
| `python ad_character_hauptmoduln.py` | identities | |
| `python defect_valuation_probe.py` | valuations | |
| `python weight3_attractor_scan.py` | P9 | |

The referee's independent scripts are listed in `audit_20261018/REPORT.md` §6.
