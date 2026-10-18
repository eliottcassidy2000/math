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
