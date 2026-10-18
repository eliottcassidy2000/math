---
id: THM-4708
title: "Rational CM values of the Argyres-Douglas character Hauptmoduln at levels 5 and 7. The product of all K characters of the Virasoro (2,2K+1) minimal model (the chiral algebra of the (A_1,A_{2K-2}) Argyres-Douglas theory) is (eta((2K+1)tau)/eta(tau))^(K-1). Hence m_5 = P_5^-6 + 22 + 125 P_5^6 (P_5 = Lee-Yang product, = T_5A + 16) and m_7 = P_7^-2 + 13 + 49 P_7^2 (P_7 = (2,7) product, = T_7A + 9) are the Fricke Hauptmoduln whose reciprocals parametrize the weight-2 Apery-like sequences C(2n,n) b_n (A274786) and Cooper's s_7. At imaginary quadratic tau, m_N(tau) is rational for exactly ten values each, all integers: m_7 in {-10648, -64, -8, -1, 0, 2, 8, 24, 27, 125} and m_5 in {-15228, -828, -28, -3, 0, 4, 18, 36, 72, 147}, with discriminants as tabulated. The nonzero ones are exactly the survivors of S43's Lucas-block CM detector in 0 < |m| <= 10^5."
status: >
  PROVED, computer-assisted, from cited classical inputs.
  - Part 1 is elementary: Andrews-Gordon products plus a conformal-weight sum, proved below for all K.
  - Part 2 is CITED for s_7 (Cooper 2012; re-verified to q^39). For A274786 it is a q-series identity checked
    exactly to O(q^60) and numerically to 1e-83, together with its recurrence. The product formula
    U_n = C(2n,n) b_n is Bala's (OEIS A274786), re-verified for n <= 520.
  - Part 3, the classification, is an exact resultant computation over the 42 orders of class number <= 2.
    - The list of those orders is CITED: Heegner-Baker-Stark (h = 1); Baker, Stark, Montgomery-Weinberger and Arno
      for h = 2 including non-maximal orders. It is recomputed here by counting forms for |D| <= 2000.
    - The Hilbert class polynomials are obtained by 120-digit evaluation and rounding, with worst margin 6.9e-88.
      This is the only floating-point step.
    - The level-N j-maps are CITED (classical; e.g. Maier, "On rationally parametrized modular equations", 2009)
      and re-verified to 1e-80.
  - Part 4 is FINITE-EXACT: (B1') scans with filter primes 7..43 and confirm primes 47..200, plus the archimedean
    CM matches P4 and P8b of the pre-registration file.
  Novelty, checked against Beukers arXiv:2403.03301 (appendix tables), OEIS A274786, Cooper 2012 and Chen-Yui 1996
  (by reference only):
  - individual CM values of these Hauptmoduln are expected to appear in the Ramanujan-Sato literature
    (Cooper 2012 at level 7; Huber-Schultz-Ye 2023 at level 5) and in Chen-Yui's singular values of Thompson series;
  - no statement of completeness was located;
  - the Argyres-Douglas reading is new as far as checked.
session: kind-pasteur-2026-10-18-S2
source: 05-knowledge/results/lucas_block_cm_detector_sporadic_ad_hauptmoduln_20261018.md
scripts:
  - 04-computation/cm_detector_sporadic/integral_cm_classification.py (+ .out, .json): part 3
  - 04-computation/cm_detector_sporadic/cm_archimedean.py, level5_leeyang_family.py: j-maps, identities, CM matches
  - 04-computation/cm_detector_sporadic/ad_character_hauptmoduln.py: part 1 and the mirror maps, exact to q^59
  - 04-computation/cm_detector_sporadic/sporadic_cm_detector.py, ordinary_detector_and_mechanism.py: part 4
  - 04-computation/cm_detector_sporadic/PREREGISTRATION.md: P4, P5, P8, written and pushed before the runs
related: [HYP-9330, HYP-9331, HYP-9310, HYP-9311, THM-4669, HYP-9297, Cooper Ramanujan J. 29 (2012), Beukers arXiv:2403.03301, OEIS A274786, Chen-Yui 1996, Conway-Norton 1979]
---

# THM-4708 — Rational CM values of the Argyres-Douglas character Hauptmoduln at levels 5 and 7

## Setting

- For `K >= 2`, let `chi_s = chi^{(2,2K+1)}_{1,s}`, `s = 1, ..., K`, be the characters of the Virasoro `(2,2K+1)`
  minimal model, including the factor `q^{h_{1,s} - c/24}`.
- This minimal model is the chiral algebra of the 4d `N = 2` Argyres-Douglas theory `(A_1, A_{2K-2})`. Its vacuum
  character is the Schur index; the other characters are module characters.
- By Andrews-Gordon, `chi_s = q^{h_{1,s} - c/24} prod_{k >= 1, k != 0, +-s (mod 2K+1)} (1 - q^k)^{-1}`.
- Write `P_5 = chi_0 chi_{-1/5}` (Lee-Yang, `K = 2`) and `P_7 = chi_1 chi_2 chi_3` (`K = 3`).

## Statement

1. **(Character products.)** For every `K >= 2`:

       prod_{s=1}^{K} chi_s(tau) = (eta((2K+1)tau)/eta(tau))^{K-1}.

   In particular `P_5 = eta(5tau)/eta(tau)`, `P_7^2 = (eta(7tau)/eta(tau))^4`, and for `K = 4`
   `chi_1 chi_2 chi_3 chi_4 = (eta(9tau)/eta(tau))^3`.

2. **(Fricke Hauptmoduln and the Apery-like families.)** Let

       m_5 = P_5^{-6} + 22 + 125 P_5^6 = T_5A + 16 = q^{-1} + 16 + 134q + 760q^2 + ...,
       m_7 = P_7^{-2} + 13 + 49 P_7^2  = T_7A + 9  = q^{-1} + 9 + 51q + 204q^2 + ...,

   the Hauptmoduln of `Gamma_0(5)+` and `Gamma_0(7)+`. With `x = 1/m_N`:

       sum_n C(2n,n) b_n x^n = (5 E_2(5tau) - E_2(tau))/4      (N = 5; b_n the Apery numbers for zeta(2); OEIS A274786),
       sum_n s_7(n) x^n      = (7 E_2(7tau) - E_2(tau))/6      (N = 7; Cooper's s_7).

   The coefficients satisfy:
   - `(n+1)^3 U_{n+1} = 2(2n+1)(11n^2+11n+3) U_n + 4n(4n^2-1) U_{n-1}` for `U_n = C(2n,n) b_n`;
   - `(n+1)^3 s_{n+1} = (2n+1)(13n^2+13n+4) s_n + 3n(9n^2-1) s_{n-1}` for `s_7`.

3. **(Classification.)** For imaginary quadratic `tau`, `m_N(tau)` is rational if and only if it lies in the
   following lists. Every value occurs. `D` is the discriminant of the order of `C/(Z + N tau Z)`; it is the same at
   both Fricke-conjugate points.

   | `m_7` | -10648 | -64 | -8 | -1 | 0 | 2 | 8 | 24 | 27 | 125 |
   |---|---|---|---|---|---|---|---|---|---|---|
   | `D` | -427 | -91 | -35 | -7 | -3 | -12 | -19 | -27 | -28 | -112 |

   | `m_5` | -15228 | -828 | -28 | -3 | 0 | 4 | 18 | 36 | 72 | 147 |
   |---|---|---|---|---|---|---|---|---|---|---|
   | `D` | -235 | -115 | -35 | -15 | -4 | -11 | -16 | -19 | -40 | -60 |

   - Every rational value is an integer.
   - `m_N = 0` is the elliptic point of order 3 (`N = 7`, `D = -3`) or order 2 (`N = 5`, `D = -4`).
   - `m_7 = -1, 27` and `m_5 = 22 +- 10 sqrt5` are the Fricke-fixed points, i.e. the singular fibres of the
     families. At `N = 5` the Fricke-fixed points are irrational and so do not appear.

4. **(The p-adic detector finds all of them; FINITE-EXACT.)** For S43's Lucas-block test (B1'),

       sum_{k=p}^{2p-1} u_k m^{-k} == u_1 m^{-1} sum_{k<p} u_k m^{-k}   (mod p^2),

   applied to `u = U` and to `u = s_7` at every prime `7 <= p <= 200` not dividing `m`, the survivors in
   `0 < |m| <= 10^5` are exactly the nonzero values in part 3. Since every rational CM value is a listed integer, the
   detector misses none at any height. The survivors also coincide with the integers of "supersingular density"
   `>= 0.3`, the CM signature.

## Proof

**Part 1.** Each nonzero residue class mod `2K+1` is `+-s` for exactly one `s` in `{1, ..., K}`. So in the product of
the `K` Andrews-Gordon products every class is omitted once and kept `K-1` times:

    prod_s prod_{k != 0, +-s} (1-q^k)^{-1} = prod_{(2K+1) does not divide k} (1-q^k)^{-(K-1)}.

The right side is `(prod (1-q^{(2K+1)k}) / prod (1-q^k))^{K-1}`. For the leading power, `c = 1 - 3(2K-1)^2/(2K+1)`
and `h_{1,s} = ((2K+1-2s)^2 - (2K-1)^2)/(8(2K+1))`. Using `sum_{s=1}^{K} (2K+1-2s)^2 = K(2K-1)(2K+1)/3`:

    sum_s (h_{1,s} - c/24) = K(K-1)/12 = (K-1) ((2K+1) - 1)/24,

which is the leading exponent of `(eta((2K+1)tau)/eta(tau))^{K-1}`. Exact q-series check for `K = 2, 3, 4` through
`q^59`: `ad_character_hauptmoduln.py`.

**Part 2.**
- The Hauptmodul property is classical. `t = (eta(tau)/eta(N tau))^{24/(N-1)}` is a Hauptmodul of `Gamma_0(N)` for
  `N = 5, 7`.
- The Fricke involution acts by `t -> N^{12/(N-1)}/t`, so `t + N^{12/(N-1)}/t + const` is a Hauptmodul of
  `Gamma_0(N)+`. The constants 22 and 13 put the order-2 (resp. order-3) elliptic point at `m = 0`, which is what
  makes the recurrences three-term.
- `q`-expansions identify `m_N` with the McKay-Thompson series 5A and 7A, shifted by 16 and 9.
- The `s_7` identity is Cooper (2012). The `N = 5` identity and its recurrence were checked exactly to `O(q^60)`
  (the recurrence was found as the unique one-dimensional nullspace; the constants 21 and 23 give none), and the
  generating function was checked numerically to `1e-83` at `tau = 0.1 + 1.2i`. `U_n = C(2n,n) b_n` holds for
  `n <= 520`; it is Bala's formula in A274786.

**Part 3.**
1. **Degree bound.** Put `t = P_7^{-2}` (resp. `s = P_5^{-6}`). Then `m_7 = t + 13 + 49/t` and
   `m_5 = s + 22 + 125/s`, so for rational `m` the variable `t` (resp. `s`) is a root of `t^2 + (13 - m)t + 49`
   (resp. `s^2 + (22 - m)s + 125`) and has degree `<= 2` over `Q`.
2. **The j-maps** (classical; numerically `< 1e-79` at three random `tau`):

       j(7tau) = (t^2+13t+49)(t^2+5t+1)^3/t = m((m-8)t - 48)^3,     j(tau) = (t^2+13t+49)(t^2+245t+2401)^3/t^7,
       j(5tau) = (s^2+10s+5)^3/s = ((m-12)s - 120)^3/s,            j(tau) = (s^2+250s+3125)^3/s^5.

   Hence `J = j(N tau)` lies in a field of degree `<= 2`.
3. **Class number.** If `tau` is imaginary quadratic, `J` is a singular modulus of some order of discriminant `D`,
   and `[Q(J):Q] = h(D)`. So `h(D) <= 2`.
4. **The finite list.** There are exactly 13 orders with `h = 1` and 29 with `h = 2` (CITED; recomputed for
   `|D| <= 2000`).
5. **Resultants.** For each such `D`, `R_D(m) = Res_J(Res_t(tJ - F_N(t), t^2 + (a_N - m)t + b_N), H_D(J))` is an
   integer polynomial in `m`. Its rational roots are exactly the rational `m` admitting a point with `J` of
   discriminant `D`.
   - **Converse.** For each such root, the Hauptmodul takes every value on the upper half plane, so some `tau` has
     `t(tau) = t`. Then `j(N tau) = J` is CM, and `tau` is imaginary quadratic.
   - **Result.** Computing all 42 resultants and keeping the linear factors gives exactly the two lists. No
     non-integral rational root occurs.
   - **Numerics.** `H_D` was rounded from 120-digit products of `j(tau_Q)` over the reduced forms, with margin
     `<= 6.9e-88`. The resultants themselves are exact (sympy).
6. **Independent check.** For every listed value except `m = 0`, both roots, solving `j(tau) = J` by Newton in the
   fundamental domain returns `tau` equal to `tau_Q` of a reduced form of the stated `D`, to `1e-40`. At level 7 the
   eta quotient `(eta(7tau')/eta(tau'))^4` at the matching one of the 8 coset points reproduces the root `h_0` to
   `< 2e-78` (`cm_archimedean.py`; not done at level 5). Controls `m = 3, 5, 11, 100` match no form with
   `|D| <= 4000`, at either level.

**Part 4.** This is a finite computation: `sporadic_cm_detector.py` and `ordinary_detector_and_mechanism.py` for
`s_7`, and `level5_leeyang_family.py` for `U`. Both were pre-registered (P4, P8) and pushed before the runs.

## Reading

- **A dictionary of two coordinates.** An integer `m` survives the mod-`p^2` Lucas-block congruence of the
  Apery-like coefficients exactly when the Argyres-Douglas character invariant `m_N(tau)` takes the value `m` at a
  CM point `tau`.
  - In S43's bridge-note language (section 4.4, Moore; Gukov-Vafa), these are the attractor / rational-CFT points of
    the torus `E_tau`. The fibres of the weight-2 families are K3 surfaces whose transcendental motive is
    `Sym^2 H^1(E_tau)`.
  - Here the coordinate on that moduli space is built from Schur-index data of the `(A_1, A_2)` and `(A_1, A_4)`
    theories.
  - This is an identification of coordinates, not a physical mechanism.
- **Level 9.** The same characters give Zagier's sporadic weight-1 case B: its mirror map is
  `t_B(q) = -P_9(-q)` with `P_9 = chi_1 chi_2 chi_3 chi_4 = (eta(9tau)/eta(tau))^3`, exact through `q^59`.
- **Level 5, weight 1.** Zagier's case D, the Apery numbers for zeta(2), has mirror map `(chi_0/chi_{-1/5})^5`
  (Beukers 1987).
  - D has no rational CM points at all.
  - The parameter is a `Gamma_1(5)` Hauptmodul, and `j(5tau)` is the j-invariant of the Tate normal form with a
    rational point of order 5. CM curves over `Q` have no such point (Olson 1974).
  - The Lee-Yang CM points become visible only on the Fricke quotient `Gamma_0(5)+`, i.e. in part 3.

## Falsifiers and extensions

- **A gap in part 3.** Part 3 would fail if the list of class-number-`<= 2` orders were incomplete. That would
  contradict published theorems.
- **Level 9 (`(A_1, A_6)`).** The weight-2 family `P_9`-Fricke and its rational CM values. Prediction: the detector
  survivors equal the rational CM values.
- **Levels 13 and 25.** Other genus-0 levels `2K+1` (13, 25) are natural next tests; `(A_1, A_10)` gives
  `(eta(13tau)/eta(tau))^5`.
