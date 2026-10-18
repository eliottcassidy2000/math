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
