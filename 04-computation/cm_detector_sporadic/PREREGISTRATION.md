# Pre-registration (kind-pasteur-2026-10-18-S2), written BEFORE any scan was run

Extension of S43's p-adic CM detector (cm_detector_lucas_block_scan_20261018.py, HYP-9310/9311) from hypergeometric
to sporadic Apery-like families. Test (B1'), with t_k = u_k m^{-k}, at primes 7 <= p <= 43 (filter), then
47 <= p <= 200 (confirm):
    sum_{k=p}^{2p-1} t_k == t_1 sum_{k<p} t_k (mod p^2).

Families:
- zeta(2)-type (Zagier's sporadic cases), from the recurrence
  (n+1)^2 u_{n+1} = (a n^2 + a n + b) u_n - c n^2 u_{n-1}, with u_0 = 1, u_1 = b:
  A = (7,2,-8), B = (9,3,27), C = (10,3,9), D = (11,3,-1) (Apery zeta(2)), E = (12,4,32), F = (17,6,72).
  The integrality of u_n for n <= 400 is used as a check on the parameters.
- zeta(3)-type:
  - Apery A(n) = sum C(n,k)^2 C(n+k,k)^2;
  - Domb D(n) = sum C(n,k)^2 C(2k,k) C(2n-2k,n-k);
  - Cooper s_7(n) = sum C(n,k)^2 C(n+k,k) C(2k,n);
  - s_10(n) = sum C(n,k)^4.

Predictions:
- (P1) Every family has a nonempty, finite survivor set in 0 < |m| <= 10^5.
- (P2) Every survivor has "supersingular density" in [0.35, 0.65]. This density is the fraction of primes
  7 <= p <= 500, p not dividing m, with sum_{k<p} u_k m^{-k} == 0 (mod p), i.e. a CM signature. Non-survivor controls
  (m in {2, 3, 5, 7, 11, 100, -100, 999}) have density < 0.10.
- (P3) For Zagier D (level 5), every survivor m corresponds to a CM point. That is, the solution tau of
  t(tau) = 1/m, with t = q prod (1-q^n)^{5(n/5)} = (chi_0/chi_{-1/5})^5 the Lee-Yang character ratio, is a quadratic
  irrationality. Test: find integers (a, b, c) with |a|, |b|, |c| <= 10^4 such that a tau^2 + b tau + c = 0 to 1e-25.

## Outcome of P1-P3 (recorded after the runs, M = 2000 and M = 100000; predictions above unchanged)

- P1 FAILED for Zagier B and D (no survivors), and P2 FAILED for every zeta(2)-type family. Their survivors are
  exactly the rational singular fibres of the recurrence, with density 0:
  - A = Franel: {-1, 8};
  - C: {1, 9};
  - E: {4, 8};
  - F: {8, 9}.

  B and D have irrational singular fibres and no survivors. Non-survivors such as m = 2 (A) and m = 3 (B, C) have
  density 0.51. The same holds for the four weight-1 hypergeometric families: the only survivor is the singular
  point (16, 27, 64, 432), and the CM parameters (e.g. 8, 32, -16 for C(2k,k)^2) have density 0.53 but fail the
  test.
- P1 and P2 PASSED for all zeta(3)-type families (Apery, Domb, s7, s10): survivors are nonempty and finite, with
  density 0.45-0.54. They include the rational singular fibres (Domb 4, 16; s7 -1, 27; s10 -4, 16).
- P3 is VACUOUS: Zagier D has no survivors.

## Second pre-registration (written BEFORE running these tests)

- (P4) Cooper s7 (level 7). Every block-test survivor m with |m| <= 10^5 that is not a singular fibre (m not in
  {-1, 27}) is a CM point. Test: solve 49h^2 + (13-m)h + 1 = 0 (x = 1/m = h/(1+13h+49h^2), with
  h = (eta(7tau)/eta(tau))^4 = (chi_1 chi_2 chi_3)^2 the product of the (A_1,A_4) Schur characters squared). For each
  root h, j(tau) = (t^2+13t+49)(t^2+5t+1)^3/t with t = 1/h must equal j(tau_Q) for a reduced primitive form Q of
  discriminant -4000 <= D < 0 (|difference| < 1e-20 |j|).
- (P3') Zagier D (level 5). The set S_D = {m : 0 < |m| <= 3000, D-density >= 0.30 over primes 7..500} is
  nonempty. For each m in S_D, j at the point with t(tau) = (chi_0/chi_{-1/5})^5 = 1/m is a CM j-invariant (same
  test, discriminant >= -4000). The j-map for Gamma_1(5) is obtained numerically from the t-expansion; tau is found
  by Newton iteration.
- (P5) Control. Non-CM rational parameters (e.g. m = 3, 5, 11, 100 for s7) fail the j-match.

## Outcome of P3', P4, P5 (recorded after running cm_archimedean.py; predictions above unchanged)

Step 0 (identities, relative error < 1e-78 at three random tau):
- (chi_1 chi_2 chi_3)^2 = (eta(7tau)/eta(tau))^4 for the Virasoro (2,7) = (A_1,A_4) characters.
- (chi_0/chi_{-1/5})^5 = q prod (1-q^n)^{5(n/5)} and chi_0 chi_{-1/5} = eta(5tau)/eta(tau) for Lee-Yang.
- Level 7: CORRECTION OF A LABEL. The pre-registered (t^2+13t+49)(t^2+5t+1)^3/t is j(7tau), not j(tau);
  j(tau) = (t^2+13t+49)(t^2+245t+2401)^3/t^7. Since tau is CM iff 7tau is, the test itself is unaffected.
- Level 5: j(tau) = -(t^4-228t^3+494t^2+228t+1)^3/(t(t^2+11t-1)^5) and
  j(5tau) = -(t^4+12t^3+14t^2-12t+1)^3/(t^5(t^2+11t-1)), with t = r^5.

- P4 PASSED. All 7 non-singular s7 survivors are CM points: for both roots h, the pre-registered j equals
  j(tau_Q) to < 1e-75 relative error. The discriminants (reduced forms, root 1 / root 2) are:

  | m | D | forms |
  |---|---|---|
  | -10648 | -427 | (1,-1,107) / (7,7,17) |
  | -64 | -91 | (1,-1,23) / (5,-3,5) |
  | -8 | -35 | (1,-1,9) / (3,1,3) |
  | 2 | -12 | (1,0,3) |
  | 8 | -19 | (1,+-1,5) |
  | 24 | -27 | (1,+-1,7) |
  | 125 | -112 | (4,0,7) / (1,0,28) |

  - Hauptmodul cross-check: among the 8 Gamma_0(7)-classes of tau' with 7tau' ~ tau_Q, one has
    (eta(7tau')/eta(tau'))^4 = h_0 to < 2e-78, for every root.
  - Post hoc, not pre-registered: the singular fibres -1 and 27 are CM too, with D = -7 and D = -28. They are the
    two fixed points h = -+1/7 of the Fricke involution h -> 1/(49h), i.e. the elliptic points of Gamma_0(7)+.
- P5 PASSED. For m = 3, 5, 11, 100, no reduced form with |D| <= 4000 matches either root.
- P3' FAILED. S_D is EMPTY: every 0 < |m| <= 3000 has D-density < 0.2. The density histogram over the bins
  [0,.05), [.05,.1), [.1,.2), [.2,1] is 3318, 2642, 40, 0.
  - Post-hoc explanation: t = r^5 is a Hauptmodul for Gamma_1(5), and j(5tau) above is the j-invariant of the Tate
    normal form with s = -t.
  - So a rational non-cuspidal t gives an elliptic curve over Q with a rational point of order 5.
  - No CM elliptic curve over Q has one (Olson 1974: torsion is 0, Z/2, Z/3, Z/4, Z/6 or Z/2 x Z/2).
  - Hence family D has no rational CM fibres at all, and P3 and P3' were ill-posed. This is not evidence about the
    detector.

## Third pre-registration (written BEFORE running the scans below)

Context, already observed in slope_probe.py, exploratory, 8 weight-1 and 4 weight-2 points:
- (B1') <=> V(x) == q_p(m) u_1 F_p(x) (mod p), where V(x) = sum_{j<p} v_j x^j and v_j = (u_{p+j} - u_1 u_j)/p.
- At the weight-1 CM points (A: -4, 2; B: 3, 9; C: -3, 3; F: 6, 12):
  - the slope V/(u_1 q_p F_p) is 1 at every ordinary prime 7 <= p <= 257 (25 primes each);
  - V(x) != 0 at every supersingular prime (27 each).
- In weight 2, V(x) == 0 at every supersingular prime seen, at CM points and at non-CM points alike.

Predictions:
- (P6) Ordinary-restricted detector (B1'_ord). Same scan as P1, but a prime counts only if p does not divide m and
  F_p(1/m) != 0 (mod p). For the six Zagier families, the survivors in 0 < |m| <= 10^5 (p <= 200) are exactly the
  rational singular fibres together with the density >= 0.3 set:
  - A {-4, -1, 2, 8};
  - B {3, 9};
  - C {-3, 1, 3, 9};
  - D {};
  - E {4, 8};
  - F {6, 8, 9, 12}.
  For the four weight-2 families, (B1'_ord) and (B1') have the same survivors.
- (P7) Polynomial form of the mechanism, for every prime 7 <= p <= 199:
  - In each weight-2 family (Apery, Domb, s7, s10), V vanishes at every root of F_p in the algebraic closure
    of F_p (the radical of F_p divides V mod p), and F_p has no simple root other than possibly x = 0.
  - In each weight-1 family (A..F), F_p is squarefree mod p and gcd(F_p, V) = 1.

## Outcome of P6 and P7 (recorded after running ordinary_detector_and_mechanism.py, M = 10^5; predictions above unchanged)

- P6 PASSED, exactly as predicted. The (B1'_ord) survivors are:
  - A {-4, -1, 2, 8};
  - B {3, 9};
  - C {-3, 1, 3, 9};
  - D {};
  - E {4, 8};
  - F {6, 8, 9, 12}.

  For the four weight-2 families, (B1'_ord) and (B1') have identical survivors.
- P7, weight-1 half PASSED: for all six Zagier families and every 7 <= p <= 199, F_p is squarefree mod p and
  gcd(F_p, V) = 1.
- P7, weight-2 half SPLIT.
  - "rad(F_p) | V" PASSED: V vanishes at every root of F_p in the algebraic closure of F_p, for all four families
    and every 7 <= p <= 199.
  - "no simple root" FAILED: F_p has simple roots at 21 (Apery), 22 (Domb), 21 (s7) and 35 (s10) of the 41 primes.
  - Post hoc, checked for all p <= 199: the simple-root part of F_p always divides the singular-fibre polynomial.
    The singular fibres are x^2-34x+1 (Apery), (1-4x)(1-16x) (Domb), (1+x)(1-27x) (s7) and (1+4x)(1-16x) (s10).
  - So the simple roots are exactly the singular fibres that are supersingular at p. For s7 these are the primes
    p == 3, 5, 6 (mod 7), inert in Q(sqrt-7), consistent with D = -7 and D = -28.
  - Reading: F_p is a square pulled back along the 2:1 Fricke quotient. It is simple only at the branch points, the
    elliptic CM points of Gamma_0(N)+.

## Fourth pre-registration (written BEFORE running): the Lee-Yang weight-2 family on Gamma_0(5)+

Context (exploratory, scratch): let H = (eta(5tau)/eta(tau))^6 = (chi_0 chi_{-1/5})^6, the Lee-Yang product to the
sixth power, and x = H/(1 + 22H + 125H^2), which is Fricke-invariant (H -> 1/(125H)). Expanding
Z = (5E_2(5tau) - E_2(tau))/4 in x gives U = 1, 6, 114, 2940, 87570, ..., with exact recurrence
(found by nullspace from 50 terms; a = 21 and a = 23 admit none):

    (n+1)^3 U_{n+1} = 2(2n+1)(11n^2+11n+3) U_n + 4n(4n^2-1) U_{n-1}.

Predictions, for the recurrence-defined U (integrality checked to n = 520):
- (P8a) The (B1') survivors in 0 < |m| <= 10^5 are nonempty and equal the density >= 0.3 set. (B1'_ord) gives
  the same set.
- (P8b) Every survivor is a CM point. Test: both roots s of s^2 + (22-m)s + 125 = 0, with s = 1/H. The value
  j(5tau) = (s^2+10s+5)^3/s (to be verified numerically first, together with j(tau) = (s^2+250s+3125)^3/s^5) must
  match j(tau_Q) for a reduced form with -4000 <= D < 0. Controls m = 3, 5, 11, 100 must fail.
- (P8c) For 7 <= p <= 199, rad(F_p) | V, and the simple-root part of F_p divides 1 - 44x - 16x^2, the singular
  polynomial with m = 22 +- 10 sqrt5 (D = -20 expected there).

## Outcome of P8 (recorded after running level5_leeyang_family.py; predictions above unchanged)

Step 0 (relative errors < 3e-80):
- (chi_0 chi_{-1/5})^6 = H;
- j(tau) = (s^2+250s+3125)^3/s^5 and j(5tau) = (s^2+10s+5)^3/s;
- sum_{n<400} U_n x^n = (5E_2(5tau)-E_2(tau))/4 at tau = 0.1+1.2i.

- P8a PASSED. The (B1') survivors, the (B1'_ord) survivors and the density >= 0.3 set all coincide (0 < |m| <= 10^5):
  {-15228, -828, -28, -3, 4, 18, 36, 72, 147}, with densities 0.47-0.54.
- P8b PASSED. Every survivor is CM, with both roots matching a reduced form; controls m = 3, 5, 11, 100 match none.

  | m | D | forms |
  |---|---|---|
  | -15228 | -235 | (5,5,13) / (1,-1,59) |
  | -828 | -115 | (5,5,7) / (1,-1,29) |
  | -28 | -35 | |
  | -3 | -15 | |
  | 4 | -11 | |
  | 18 | -16 | |
  | 36 | -19 | |
  | 72 | -40 | |
  | 147 | -60 | |

- P8c PASSED. For every 7 <= p <= 199, rad(F_p) | V, and the simple-root part of F_p divides 1 - 44x - 16x^2.
  Simple roots occur exactly at p = 11, 13, 17, 19, 31, 37, 53, 59, ...; these are the primes inert in Q(sqrt-5)
  (spot-checked), consistent with D = -20 at m = 22 +- 10 sqrt5.
