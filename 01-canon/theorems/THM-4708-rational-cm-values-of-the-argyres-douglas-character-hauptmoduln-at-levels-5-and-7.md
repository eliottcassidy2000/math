---
id: THM-4708
title: "Rational CM values of the Fricke Hauptmoduln T_5A + 16 and T_7A + 9, and S43's Lucas-block detector. m_5 = T_5A + 16 and m_7 = T_7A + 9 parametrize the weight-2 Apery-like sequences C(2n,n) b_n (A274786) and Cooper's s_7. At imaginary quadratic tau, m_N(tau) is rational for exactly ten values each, all integers: m_7 in {-10648, -64, -8, -1, 0, 2, 8, 24, 27, 125} and m_5 in {-15228, -828, -28, -3, 0, 4, 18, 36, 72, 147}. The nine nonzero values at each level were tabulated before (Beukers 2024 App. C; Z.-W. Sun 2011); the completeness proof is routine. They are exactly the survivors of S43's Lucas-block CM detector in 0 < |m| <= 10^7, so the detector misses none. Identification of coordinates (not physics): m_N is built from the product of all characters of the Virasoro (2,2K+1) minimal model, the chiral algebra of the (A_1, A_{2K-2}) Argyres-Douglas theory; that product is (eta((2K+1)tau)/eta(tau))^(K-1)."
status: >
  PROVED (computer-assisted) from cited inputs; credits corrected after the independent audit (audit_20261018,
  MISTAKE-638).
  - Part 1 (the character product) is PROVED. It is elementary and folklore-level, not advertised as new:
    - K = 2 is Ramanujan's G(q)H(q) = (q^5;q^5)_inf/(q;q)_inf, and general K follows in one line from Andrews-Gordon.
    - The identification of the (A_1, A_{2K-2}) chiral algebra with Vir(2,2K+1) is CITED: Beem, Lemos, Liendo,
      Peelaers, Rastelli and van Rees, arXiv:1312.5344; Cordova-Shao, arXiv:1506.00265.
  - Part 2 is CITED.
    - Level 7: Cooper, Ramanujan J. 29 (2012); Beukers arXiv:2403.03301 v3, App. C §9.11; the eta form of m_7 is
      OEIS A282877 (Somos 2017).
    - Level 5: Beukers v3 App. C §9.14 ("Apery2 ⊙ central": operator, mirror map t = h_5/(1+22h_5+125h_5^2),
      -4 Theta = E_{2,5}, group Gamma_0(5)+<5>); Chan-Cooper, Math. Proc. Camb. Phil. Soc. 153 (2012), Thm 3.1/3.5
      (as reported by the referee's sub-agent); Cooper, Ramanujan's Theta Functions (2017), (6.73).
    - A274786: recurrence in OEIS (Kotesovec 2017); product formula U_n = C(2n,n) b_n (Bala 2024).
    - Re-verified exactly: here to O(q^60); by the referee to q^161 (generating functions) and q^1499 (Hauptmoduln
      against the McKay-Thompson b-files).
  - Part 3 is PROVED (computer-assisted), with prior tabulation of the values.
    - Prior tabulation: the nine nonzero values at each level are tabulated, with their tau but WITHOUT a
      completeness claim, in Beukers v3 App. C §9.11 and §9.14. The level-5 values are also in Z.-W. Sun,
      arXiv:1103.4325 v3 (2011), §3, as (m, d).
    - Beukers v3 erratum: t(sqrt(-3/5)) is printed 1/588; the correct value is 1/147, as in v1/v2 and in part 3.
    - Completeness over all CM tau, including m = 0, is routine: it is the standard description of rational CM points
      on X_0^+(N) (Galbraith, Exp. Math. 8 (1999); arXiv:1607.04558 §3 for h = 1). It is proved here by exact
      resultants over the 42 orders of class number <= 2.
    - Those orders are CITED: Heegner, Baker, Stark for h = 1; Baker 1971, Stark 1975, Montgomery-Weinberger 1974
      (and Watkins 2004) for fundamental h = 2. Non-maximal orders follow from the class-number formula
      h(f^2 D_K) = h(D_K) f prod(1 - (D_K/l)/l)/[O_K^x : O^x] (Cox, Thm 7.24), which forces phi(f) <= 6.
      Recomputed here: 13 + 29 orders.
    - Rigour of the computation:
      - The author's Hilbert class polynomials come from 120-digit rounding (margin 6.9e-88).
      - The referee reproduced both lists and all 20 D labels by a different method, with no H_D, no j-maps and no
        resultants: eta-product values at all (N+1) h(D) CM points, at 300 digits, margin 2.1e-234.
      - The four j-maps are exact q-series identities through O(q^142).
    - Integrality of the values is also CITED: Chen-Yui, Thm I / Thm 3.4, as reported by the referee's sub-agent.
      Fundamental Thompson series are algebraic integers at CM points with (a, N) = 1.
    - Huber-Schultz-Ye (Acta Arith. 207, 2023): the abstract claims "a complete list of rational and quadratic
      series corresponding to singular values of the parameters". Its body was not accessible to the session or the
      referee, so how its completeness compares with part 3 (all CM tau, including the split values and m = 0) is not
      settled.
  - Part 4 is FINITE-EXACT.
    - Author: (B1') with filter primes 7..43 and confirm primes 47..200, 0 < |m| <= 10^5; pre-registered P4, P5, P8.
    - Referee: unchanged to |m| <= 10^7, and all 18 nonzero values pass at every prime 7 <= p <= 997.
  New as far as checked: the detector-completeness statement (part 4 together with part 3) and the Argyres-Douglas
  packaging, which is an identification of coordinates.
session: kind-pasteur-2026-10-18-S2
source: 05-knowledge/results/lucas_block_cm_detector_sporadic_ad_hauptmoduln_20261018.md
scripts:
  - 04-computation/cm_detector_sporadic/integral_cm_classification.py (+ .out, .json): part 3 (resultants)
  - 04-computation/cm_detector_sporadic/cm_archimedean.py, level5_leeyang_family.py: j-maps, identities, CM matches
  - 04-computation/cm_detector_sporadic/ad_character_hauptmoduln.py: part 1 and the mirror maps, exact to q^59
  - 04-computation/cm_detector_sporadic/sporadic_cm_detector.py, ordinary_detector_and_mechanism.py: part 4
  - 04-computation/cm_detector_sporadic/audit_20261018/ (a1, a2, a3, a3b, a3c, a4c, a4d, a9): the independent referee
  - 04-computation/cm_detector_sporadic/PREREGISTRATION.md: P4, P5, P8, written and pushed before the runs
related: [HYP-9330, HYP-9331, HYP-9310, HYP-9311, THM-4669, HYP-9297, Beukers arXiv:2403.03301 v3 App. C, Z.-W. Sun arXiv:1103.4325, Cooper Ramanujan J. 29 (2012), Chan-Cooper 2012, OEIS A274786, OEIS A282877, Chen-Yui 1996, Galbraith 1999, Conway-Norton 1979]
---

# THM-4708 — Rational CM values of the Fricke Hauptmoduln `T_5A + 16` and `T_7A + 9`, and S43's Lucas-block detector

## Setting

- For `K >= 2`, let `chi_s = chi^{(2,2K+1)}_{1,s}`, `s = 1, ..., K`, be the characters of the Virasoro `(2,2K+1)`
  minimal model, including the factor `q^{h_{1,s} - c/24}`.
- By Andrews-Gordon, `chi_s = q^{h_{1,s} - c/24} prod_{k >= 1, k != 0, +-s (mod 2K+1)} (1 - q^k)^{-1}`.
- This minimal model is the chiral algebra of the 4d `N = 2` Argyres-Douglas theory `(A_1, A_{2K-2})` (CITED,
  Beem et al. 2013; Cordova-Shao 2015). Its vacuum character is the Schur index; the other characters are module
  characters.
- Write `P_5 = chi_0 chi_{-1/5}` (Lee-Yang, `K = 2`) and `P_7 = chi_1 chi_2 chi_3` (`K = 3`).

## Statement

1. **(Character products; PROVED, folklore-level.)** For every `K >= 2`,

       prod_{s=1}^{K} chi_s(tau) = (eta((2K+1)tau)/eta(tau))^{K-1}.

   In particular `P_5 = eta(5tau)/eta(tau)`, `P_7^2 = (eta(7tau)/eta(tau))^4`, and for `K = 4`
   `chi_1 chi_2 chi_3 chi_4 = (eta(9tau)/eta(tau))^3`.

2. **(Fricke Hauptmoduln and the Apery-like families; CITED.)** Let

       m_5 = P_5^{-6} + 22 + 125 P_5^6 = T_5A + 16 = q^{-1} + 16 + 134q + 760q^2 + ...,
       m_7 = P_7^{-2} + 13 + 49 P_7^2  = T_7A + 9  = q^{-1} + 9 + 51q + 204q^2 + ...,

   the Hauptmoduln of `Gamma_0(5)+` and `Gamma_0(7)+`. With `x = 1/m_N`:

       sum_n C(2n,n) b_n x^n = (5 E_2(5tau) - E_2(tau))/4      (N = 5; b_n the Apery numbers for zeta(2); OEIS A274786),
       sum_n s_7(n) x^n      = (7 E_2(7tau) - E_2(tau))/6      (N = 7; Cooper's s_7).

   The coefficients satisfy:
   - `(n+1)^3 U_{n+1} = 2(2n+1)(11n^2+11n+3) U_n + 4n(4n^2-1) U_{n-1}` for `U_n = C(2n,n) b_n`;
   - `(n+1)^3 s_{n+1} = (2n+1)(13n^2+13n+4) s_n + 3n(9n^2-1) s_{n-1}` for `s_7`.

3. **(Classification; PROVED, computer-assisted.)** For imaginary quadratic `tau`, `m_N(tau)` is rational if and only
   if it lies in the following lists. Every value occurs. `D` is the discriminant of the order of `C/(Z + N tau Z)`.
   All 20 points are Heegner points, so `D` is the same at both Fricke-conjugate points; CM points in general need not
   have this property.

   | `m_7` | -10648 | -64 | -8 | -1 | 0 | 2 | 8 | 24 | 27 | 125 |
   |---|---|---|---|---|---|---|---|---|---|---|
   | `D` | -427 | -91 | -35 | -7 | -3 | -12 | -19 | -27 | -28 | -112 |

   | `m_5` | -15228 | -828 | -28 | -3 | 0 | 4 | 18 | 36 | 72 | 147 |
   |---|---|---|---|---|---|---|---|---|---|---|
   | `D` | -235 | -115 | -35 | -15 | -4 | -11 | -16 | -19 | -40 | -60 |

   - Every rational value is an integer.
   - `m_N = 0` is the elliptic point of order 3 (`N = 7`, `D = -3`) or order 2 (`N = 5`, `D = -4`).
   - `m_7 = -1, 27` and `m_5 = 22 +- 10 sqrt5` are the Fricke-fixed points, i.e. the singular fibres of the families.
     At `N = 5` they are irrational: `D = -20` is excluded because `(sqrt-5)` is principal.
   - **Structure** (referee n3). The 20 discriminants are exactly the `D` with `h(D) = 1` and `N` split or ramified,
     together with the `D` with `h(D) = 2`, `N` ramified, `N` not dividing the conductor, and the prime above `N`
     non-principal. This is the standard picture of rational CM points on `X_0^+(N)`.
   - **Prior tabulation.** The nine nonzero values at each level are in Beukers v3 App. C §9.11 / §9.14, with `tau`
     but without a completeness claim. Note the misprint there: 1/588 should be 1/147. The level-5 values are also in
     Z.-W. Sun (2011) §3.

4. **(The p-adic detector finds all of them; FINITE-EXACT.)** For S43's Lucas-block test (B1'),

       sum_{k=p}^{2p-1} u_k m^{-k} == u_1 m^{-1} sum_{k<p} u_k m^{-k}   (mod p^2),

   applied to `u = U` and to `u = s_7` at every prime `7 <= p <= 200` not dividing `m`, the survivors in
   `0 < |m| <= 10^5` are exactly the nonzero values in part 3. The referee extended this to `|m| <= 10^7`.
   - Since every rational CM value is a listed integer, the detector misses none at any height.
   - The survivors also coincide with the integers of supersingular density `>= 0.3`, the CM signature.
   - All 18 nonzero values pass at every prime `7 <= p <= 997`.

## Proof

**Part 1.** Each nonzero residue class mod `2K+1` is `+-s` for exactly one `s` in `{1, ..., K}`. So in the product of
the `K` Andrews-Gordon products every class is omitted once and kept `K-1` times:

    prod_s prod_{k != 0, +-s} (1-q^k)^{-1} = prod_{(2K+1) does not divide k} (1-q^k)^{-(K-1)}.

The right side is `(prod (1-q^{(2K+1)k}) / prod (1-q^k))^{K-1}`. For the leading power, `c = 1 - 3(2K-1)^2/(2K+1)`
and `h_{1,s} = ((2K+1-2s)^2 - (2K-1)^2)/(8(2K+1))`. Using `sum_{s=1}^{K} (2K+1-2s)^2 = K(2K-1)(2K+1)/3`:

    sum_s (h_{1,s} - c/24) = K(K-1)/12 = (K-1) ((2K+1) - 1)/24,

which is the leading exponent of `(eta((2K+1)tau)/eta(tau))^{K-1}`. Exact checks: here `K = 2, 3, 4` through `q^59`
(`ad_character_hauptmoduln.py`); by the referee `K = 2..8` through `q^300`, against Rocha-Caridi.

**Part 2.** This is CITED (see the status). `q`-expansions identify `m_N` with the McKay-Thompson series 5A and 7A
shifted by 16 and 9. The constants 22 and 13 put the order-2 (resp. order-3) elliptic point at `m = 0`, which is what
makes the recurrences three-term. The `N = 5` identity was found here independently, as the unique one-dimensional
nullspace (the constants 21 and 23 admit no recurrence), before the prior sources were located.

**Part 3.**
1. **Degree bound.** Put `t = P_7^{-2}` (resp. `s = P_5^{-6}`). Then `m_7 = t + 13 + 49/t` and
   `m_5 = s + 22 + 125/s`, so for rational `m` the variable `t` (resp. `s`) is a root of `t^2 + (13 - m)t + 49`
   (resp. `s^2 + (22 - m)s + 125`) and has degree `<= 2` over `Q`.
2. **The j-maps** (classical; exact q-series identities through `O(q^142)`):

       j(7tau) = (t^2+13t+49)(t^2+5t+1)^3/t = m((m-8)t - 48)^3,     j(tau) = (t^2+13t+49)(t^2+245t+2401)^3/t^7,
       j(5tau) = (s^2+10s+5)^3/s = ((m-12)s - 120)^3/s,            j(tau) = (s^2+250s+3125)^3/s^5.

   Hence both `j(tau)` and `j(N tau)` lie in a field of degree `<= 2`.
3. **Class number.** If `tau` is imaginary quadratic, `J = j(N tau)` is a singular modulus of some order of
   discriminant `D`, and `[Q(J):Q] = h(D)`. So `h(D) <= 2`.
4. **The finite list.** There are exactly 13 orders with `h = 1` and 29 with `h = 2` (CITED; recomputed for
   `|D| <= 2000`).
5. **Resultants.** For each such `D`, `R_D(m) = Res_J(Res_t(tJ - F_N(t), t^2 + (a_N - m)t + b_N), H_D(J))` is an
   integer polynomial in `m`. Its rational roots are exactly the rational `m` admitting a point with `J` of
   discriminant `D`.
   - A resultant over `J = j(N tau)` alone suffices, because the second root of the quadratic is the Fricke image,
     where `j(N tau') = j(tau)`.
   - **Converse.** For each such root, the Hauptmodul takes every value on the upper half plane, so some `tau` has
     `t(tau) = t`. Then `j(N tau) = J` is CM, and `tau` is imaginary quadratic.
   - **Result.** Computing all 42 resultants and keeping the linear factors gives exactly the two lists. No
     non-integral rational root occurs.
   - **Non-Heegner orbits.** Orders with `N` dividing the conductor give irrational values, e.g.
     `m_5 = -108 +- 48 sqrt5` (`D(5tau) = -75`) and `m_7 = -112 +- 24 sqrt21` (`D(7tau) = -147`).
6. **Independent checks.**
   - Session: Newton inversion of `j` in the fundamental domain returns `tau_Q` of a reduced form of the stated `D`,
     to `1e-40`. At level 7 the eta quotient at a coset point reproduces the root to `< 2e-78`. Controls
     `m = 3, 5, 11, 100` match no form with `|D| <= 4000`.
   - Referee: the eta-product method of the status line, at 300 digits.

**Part 4.** This is a finite computation: `sporadic_cm_detector.py` and `ordinary_detector_and_mechanism.py` for
`s_7`, and `level5_leeyang_family.py` for `U`. Both were pre-registered (P4, P8) and pushed before the runs. The
referee's independent scan (`a4c`, `a4d`) extends it.

## Reading (identification of coordinates; ANALOGY, no mechanism claimed)

- **Two coordinates on one set.** An integer `m` survives the mod-`p^2` Lucas-block congruence of these Apery-like
  coefficients exactly when the character invariant `m_N(tau)` takes the value `m` at a CM point `tau`.
  - In S43's bridge-note language (§4.4), CM `tau` is the necessary half of the Gukov-Vafa rational-CFT condition
    for the torus `E_tau`. They also need the Kähler modulus `rho` to be CM in the same field.
  - The fibres of the weight-2 families are K3 surfaces whose transcendental motive is `Sym^2 H^1(E_tau)`.
- **Which levels carry such a Hauptmodul.** Every odd `N` is `2K+1`. The operative condition for an eta-quotient
  Hauptmodul `(eta(tau)/eta(N tau))^{24/(N-1)}` is `(N-1) | 24`, which leaves the odd levels 5, 7, 9, 13 and 25.
  Among these, sporadic Apery-like families are known at 5, 7 and 9.
- **Level 9.** Zagier's case B has mirror map `t_B(q) = -P_9(-q)` with
  `P_9 = chi_1 chi_2 chi_3 chi_4 = (eta(9tau)/eta(tau))^3`, exact through `q^59` (also verified by the referee).
- **Level 5, weight 1.** Zagier's case D has mirror map `(chi_0/chi_{-1/5})^5` (Beukers 1987).
  - D has no rational CM points: the parameter is a `Gamma_1(5)` Hauptmodul, and CM curves over `Q` have no rational
    point of order 5 (Olson 1974).
  - The Lee-Yang CM points become visible only on the Fricke quotient `Gamma_0(5)+`, i.e. in part 3.

## Falsifiers and extensions

- **Inputs.** Part 3 would fail only if the cited list of class-number-`<= 2` orders were incomplete.
- **Level 9 (`(A_1, A_6)`).** The weight-2 Fricke family and its rational CM values. Prediction: the detector
  survivors equal the rational CM values.
- **Levels 13 and 25.**
