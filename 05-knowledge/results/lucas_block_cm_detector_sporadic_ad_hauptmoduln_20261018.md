# The Lucas-block CM detector on sporadic families: the weight dichotomy, Argyres-Douglas Hauptmoduln and attractor points

kind-pasteur-2026-10-18-S2. This session builds on opus-2026-10-18-S43 (HYP-9310, HYP-9311 and the bridge note's
§4.3–4.4) and on kind-pasteur S1 (Argyres-Douglas chiral algebras, HYP-9297). Code and data are in
`04-computation/cm_detector_sporadic/`. Every prediction was written into `PREREGISTRATION.md` and pushed before its
run: P1–P3 in 8bf465cb8c, P3'/P4/P5 in 2583387d46, P6/P7 in 4a60d0195a, P8 in 543bc0f3ba and P9 in 1846da553b.
Their outcomes were appended below the predictions afterwards; the prediction text was never edited.

## 1. Summary

1. **Weight dichotomy (HYP-9330).**
   - S43's mod-`p^2` Lucas-block test (B1') detects CM points only in weight 2.
   - On Zagier's six sporadic weight-1 families, its survivors are exactly the rational cusps.
   - The reason is located exactly. (B1') is equivalent to `V_p(x) == u_1 q_p(m) F_p(x) (mod p)`, where `V_p` is the
     first-carry Lucas defect.
     - At split primes this slope-1 law holds at CM points in both weights.
     - At inert primes it is automatic in weight 2: `rad F_p | V_p` for all 41 primes `7 <= p <= 199`, in 5
       families.
     - It is impossible in weight 1: `F_p` is squarefree and `gcd(F_p, V_p) = 1`, for all 41 primes, in 6
       families.
   - Restricting the test to ordinary primes therefore makes it a CM detector in weight 1 too (P6, exact).
   - At inert primes the defect valuation equals the weight.
2. **Argyres-Douglas Hauptmoduln (THM-4708).**
   - The product of all characters of the Virasoro `(2,2K+1)` model, i.e. the chiral algebra of `(A_1,A_{2K-2})`,
     is `(eta((2K+1)tau)/eta(tau))^{K-1}`.
   - The resulting Fricke Hauptmoduln `m_5 = T_5A + 16` (Lee-Yang) and `m_7 = T_7A + 9` (the `(2,7)` model)
     parametrize `C(2n,n) b_n` (A274786) and Cooper's `s_7`.
   - Their rational CM values are classified exactly: ten each, all integers.
   - The nonzero ones are exactly the detector's survivors, so the detector misses no CM point at any height.
3. **Attractor points (HYP-9331).**
   - In weight 3, a blind pre-registered scan of four CY3 families returns exactly the two known rational rank-2
     attractor points (Hulek-Verrill `-7`; `4F3(1/3,1/3,2/3,2/3)` at `-5832`) plus the singular points.
   - The defect valuation there matches the modulus of Beukers' conjectured supercongruences.

## 2. Pre-registered tests

| # | prediction | outcome |
|---|---|---|
| P1 | every family has a nonempty finite survivor set | FAILED for Zagier B and D (no survivors); PASSED for the rest |
| P2 | survivors have density 0.35–0.65, controls < 0.10 | FAILED in weight 1 (survivors have density 0); PASSED in weight 2 |
| P3 | Zagier D survivors are CM | VACUOUS (no survivors) |
| P3' | Zagier D has rational m with density >= 0.3, all CM | FAILED: the set is empty. Post hoc: the family lives on `Gamma_1(5)`, and Olson's theorem forbids rational CM fibres; ill-posed (MISTAKE-636) |
| P4 | the 7 nonsingular `s_7` survivors are CM | PASSED (D = -427, -91, -35, -12, -19, -27, -112); post hoc, the singular fibres -1, 27 have D = -7, -28 |
| P5 | `s_7` controls 3, 5, 11, 100 are not CM | PASSED |
| P6 | (B1'_ord) survivors in weight 1 = cusps ∪ CM, as listed | PASSED exactly, all six families |
| P7 | weight 1: `F_p` squarefree and `gcd(F_p, V_p) = 1`; weight 2: `rad F_p \| V_p` and no simple roots | PASSED except "no simple roots" (FAILED: the simple roots are the Fricke-elliptic singular fibres) |
| P8 | Lee-Yang weight-2 family: survivors = density set, all CM; `rad F_p \| V_p`; simple part divides `1-44x-16x^2` | PASSED: 9 survivors, D = -235, -115, -35, -15, -11, -16, -19, -40, -60 |
| P9 | weight 3: survivors include -7 (Hulek-Verrill) and -5832 (`(3k)!^2/k!^6`); every survivor is an attractor or singular point | PASSED |

One label was wrong in the pre-registration. The formula `(t^2+13t+49)(t^2+5t+1)^3/t` gives `j(7tau)`, not
`j(tau)`. This does not affect the test, since `tau` is CM iff `7tau` is.

## 3. The weight dichotomy (HYP-9330)

- **Data, weight 1.**
  - Survivors of (B1') with `|m| <= 10^5`: A {-1, 8}, B {}, C {1, 9}, D {}, E {4, 8}, F {8, 9}. These are the
    rational cusps, all with supersingular density 0.
  - Density >= 0.3 set: A {-4, 2}, B {3, 9}, C {-3, 3}, D {}, E {}, F {6, 12}. These are disjoint from the
    survivors and agree with Beukers' tables of rational CM values (arXiv:2403.03301, appendix).
- **Data, weight 2.** The survivors equal the density set, for Apery, Domb, `s_7`, `s_10` and A274786.
- **Slopes** (`slope_probe.py`).
  - At the 8 weight-1 CM points, `V_p/(u_1 q_p F_p) = 1` at every ordinary prime up to 257, about 25 each.
  - At the supersingular primes, `V_p != 0` every time, about 27 each.
  - In weight 2, `V_p == 0` at every supersingular prime, CM or not.
- **Polynomial form** (`ordinary_detector_and_mechanism.py`).
  - Weight 1: `F_p` is squarefree and coprime to `V_p`, every time.
  - Weight 2: `rad F_p` divides `V_p`, every time.
  - The simple roots of `F_p` in weight 2 are exactly the singular fibres that are supersingular at `p`. These are
    the Fricke-elliptic CM points, e.g. `s_7`'s -1 and 27 (D = -7, -28) at `p == 3, 5, 6 (mod 7)`.
- **Valuations** (`defect_valuation_probe.py`).
  - Inert primes: `v_p(Delta) = 1` in weight 1 and `2` in weight 2.
  - Split primes: `v_p >= 2`. It is 3 for Apery and Domb (cf. S43's (O)), and mostly 2 for `s_7`, `s_10` and
    A274786. Those three are Beukers' Remark 1.14 exceptions to the mod-`p^3` Dwork congruence.
  - Generic `m`: 1.
- **Literature calibration.**
  - Beukers proves the split-prime supercongruences at CM values (Thm 1.7 in weight 1; 1.15/1.16 in weight 2) and
    the Dwork congruences `F(t) == F_p(t) F(t^sigma)` (Thm 1.4: mod `p^2` in weight 1; Thm 1.13: mod `p^3` in
    weight 2).
  - Truncated below `2p` at a fixed point `t_0^sigma = t_0`, these give our split-prime slope law. That is an
    inference, not a statement of the paper.
  - Beukers says his method does not reach inert primes (Remarks 1.8(4), 1.17). Parts (B) and (C), and the
    dichotomy they produce, are new as far as checked.

## 4. Argyres-Douglas characters as Hauptmoduln (THM-4708 parts 1-2; `ad_character_hauptmoduln.py`, exact to `q^59`)

| AD theory | chiral algebra | character identity | family |
|---|---|---|---|
| `(A_1,A_2)` | Lee-Yang (2,5) | `(chi_0/chi_{-1/5})^5 = q prod (1-q^n)^{5(n/5)}` | Zagier D mirror map (`Gamma_1(5)`, Beukers 1987) |
| `(A_1,A_2)` | Lee-Yang (2,5) | `(chi_0 chi_{-1/5})^6 = (eta5/eta1)^6`; `x = H/(1+22H+125H^2)` | `C(2n,n) b_n` = A274786 (`Gamma_0(5)+`) |
| `(A_1,A_4)` | (2,7) | `(chi_1 chi_2 chi_3)^2 = (eta7/eta1)^4`; `x = h/(1+13h+49h^2)` | Cooper `s_7` (`Gamma_0(7)+`) |
| `(A_1,A_6)` | (2,9) | `chi_1 chi_2 chi_3 chi_4 = (eta9/eta1)^3`; `t_B(q) = -P(-q)` | Zagier B mirror map (`Gamma_0(9)`, twisted) |

So the three odd genus-0 levels 5, 7 and 9 that carry sporadic Apery-like families are exactly the levels
`2K+1` of the Argyres-Douglas theories `(A_1,A_2)`, `(A_1,A_4)` and `(A_1,A_6)`.
- This is an identification of modular coordinates. No physical mechanism is claimed.
- The CM points of these coordinates are the attractor / rational-CFT points of S43's bridge note §4.4.
- The Lee-Yang CM points are invisible on `Gamma_1(5)`, where no rational CM point exists (Olson). They appear on the
  Fricke quotient: the weight-2 family A274786 has nine.

## 5. The exact classification (THM-4708 part 3; `integral_cm_classification.py`)

- **Method.** Rational `m` gives `[Q(j(N tau)):Q] <= 2`, hence `h(D) <= 2`, hence one of 42 orders. Resultants
  against the Hilbert class polynomials then give all rational CM values:
  - `m_7`: -10648, -64, -8, -1, 0, 2, 8, 24, 27, 125;
  - `m_5`: -15228, -828, -28, -3, 0, 4, 18, 36, 72, 147.
- **Integrality.** All are integers. This is consistent with Chen-Yui: singular values of Thompson series are
  algebraic integers.
- **Literature.** The values themselves are presumably in the Ramanujan-Sato literature (Cooper 2012; Huber,
  Schultz and Ye 2023), e.g. `-10648 = -22^3`, which is `D = -427`. No completeness statement was found.

## 6. Weight 3 and attractor points (HYP-9331; `weight3_attractor_scan.py`)

- **Hulek-Verrill.** Survivors {-7, 1, 9, 25}. The attractor is `t_0 = -1/7`, from Candelas, de la Ossa, Elmi and
  van Straten (2020).
- **`((3k)!/k!^3)^2`.** Survivors {-5832, 729}. The attractor is `t_0 = -1/8`, from Boenisch, Elmi, Kashani-Poor
  and Klemm (2022), Table 1.
- **Controls.** `C(2k,k)^4` gives {256}; the quintic gives {3125}.
- **Valuations.** At ordinary primes: 2 at -7 and 3 at -5832, the moduli of Beukers' Conj. 1.27 / 1.28.
- **Post hoc.** The ordinary-restricted test also passes `C(2k,k)^4` at `z = -1`. Its non-ordinary primes are
  exactly `p == 3 (mod 4)`, so it behaves like a CM-by-`Q(i)` point.

## 7. What is new, as far as checked

- **The inert-prime laws.** (B) and (C), and the weight dichotomy of S43's detector they explain.
- **The ordinary-restricted detector.** A CM detector in weight 1.
- **Completeness.** The complete rational-CM classification for the two Fricke Hauptmoduln, together with the
  detector-completeness statement.
- **The Argyres-Douglas reading.** The level-5/7/9 sporadic Hauptmoduln as Argyres-Douglas character data.
- **Attractor points.** The detector finding rank-2 attractor points in a blind scan.
- **Not new:**
  - the split-prime supercongruences (Beukers 2024; KLMSY 2016);
  - S43's detector and harmonic law (KLMSY Cor. 3);
  - A274786 and its product formula (Bala);
  - the classical j-maps and Hauptmoduln;
  - Cooper's `s_7` parametrization;
  - Beukers' weight-1 CM value tables.

## 8. Open problems

1. Prove (B) and (C), e.g. via a Clausen identity mod `p` for the pair (period, log-partner).
2. Extend HYP-9331 to every AESZ operator with known rational attractors, and to operators without known attractors
   (a blind attractor search).
3. Level 9, `(A_1,A_6)`: the Fricke weight-2 family and its rational CM values; then levels 13 and 25.
4. Explain the extra depth (`v_p = 3`) at `s_7` `m = -1, 2`, `s_10` `m = -2, 1`, Franel `m = 2` and Zagier B `m = 9`.

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
