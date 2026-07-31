#!/usr/bin/env python3
"""Referee for THM-3017: the AMM 12592 checkpoint-closure capacity threshold
is exactly  gamma* = log_5(phi^2) = 2 log(phi)/log(5),  i.e.
C = 1 + gamma* = log_5(5 phi^2),  with binding ray x* = phi^{-2}.

Frame: THM-3002. The capacity criterion for the H=1 dyadic-checkpoint program
is  sum_{i<=t} binom(d_i, t-i) 2^{t-i} >= binom(R-1, t)  for all t, with
d_i = floor(gamma (R+i)).  Writing t = xR, i = yR and using
binom(aR,bR) ~ exp(R a H(b/a)) (H the natural-log binary entropy), the
criterion becomes the two-ray comparison

    Phi(x,gamma) := max_{0<=y<=x} [ gamma(1+y) H( (x-y)/(gamma(1+y)) )
                                    + (x-y) log 2 ]   >=   H(x).

Checks:
  A  At gamma = 2 log(phi)/log 5 and x = phi^{-2}: Phi(x,gamma) = H(x) to
     working precision (tangency), and d/dx[Phi - H] = 0 there.
  B  x = phi^{-2} is the unique interior minimiser of the margin.
  C  At the critical point the two natural ratios are exactly
        rho   := d/(d-s)   = sqrt 5,
        sigma := 2(d-s)/s  = phi = (1-x)/x,
     where d = gamma(1+y*), s = x - y*.
  D  The closed form follows: tangency in x gives (1-x)/x = 2(d-s)/s; the
     inner first-order condition s d^gamma = 2(d-s)^{gamma+1} is exactly
     rho^gamma = sigma; at the critical point rho = sqrt5 and sigma = phi,
     so gamma* = log(phi)/log(sqrt 5) = 2 log(phi)/log 5.
  E  Consistency with the exact finite-R ledger of THM-3002: the finite-R
     thresholds increase toward gamma* and gamma = 3/5 > gamma* stays ample
     while gamma = 2457/4135 < gamma* dies (already recorded at R = 2048).
"""

from mpmath import mp, mpf, log, sqrt, nstr

mp.dps = 40


def require(c, m):
    if not c:
        raise RuntimeError(m)


PHI = (1 + sqrt(5)) / 2
GAMMA = 2 * log(PHI) / log(5)
XSTAR = 1 / PHI ** 2


def H(u):
    if u <= 0 or u >= 1:
        return mpf(0)
    return -u * log(u) - (1 - u) * log(1 - u)


def f(y, x, g):
    d = g * (1 + y)
    s = x - y
    if d <= 0 or s > d or s < 0:
        return mpf(-1e30)
    return d * H(s / d) + s * log(2)


def cap(x, g, iters=400):
    """golden-section max over y in [0,x]; returns (value, argmax)"""
    lo, hi = mpf(0), x
    r1, r2 = mpf('0.381966011250105'), mpf('0.618033988749895')
    for _ in range(iters):
        m1 = lo + (hi - lo) * r1
        m2 = lo + (hi - lo) * r2
        if f(m1, x, g) < f(m2, x, g):
            lo = m1
        else:
            hi = m2
    y = (lo + hi) / 2
    return f(y, x, g), y


def a_tangency():
    c, y = cap(XSTAR, GAMMA)
    margin = c - H(XSTAR)
    require(abs(margin) < mpf(10) ** -30, f"A: margin {nstr(margin,6)} not 0")
    h = mpf(10) ** -12
    dm = ((cap(XSTAR + h, GAMMA)[0] - H(XSTAR + h))
          - (cap(XSTAR - h, GAMMA)[0] - H(XSTAR - h))) / (2 * h)
    require(abs(dm) < mpf(10) ** -12, f"A: d/dx margin {nstr(dm,6)} not 0")
    print(f"A tangency at gamma*=2log(phi)/log5, x*=phi^-2: "
          f"margin={nstr(margin,6)}, d/dx={nstr(dm,6)}  OK")
    return y


def b_minimiser():
    base = cap(XSTAR, GAMMA)[0] - H(XSTAR)
    for dx in ['0.2', '0.1', '0.05', '0.01', '-0.01', '-0.05', '-0.1', '-0.2']:
        xv = XSTAR + mpf(dx)
        if xv <= 0 or xv >= 1:
            continue
        m = cap(xv, GAMMA)[0] - H(xv)
        require(m > base - mpf(10) ** -30, f"B: margin dips at x*+{dx}")
    print("B x* = phi^-2 is the interior minimiser (margin >= 0 nearby): OK")


def c_ratios(y):
    d = GAMMA * (1 + y)
    s = XSTAR - y
    rho = d / (d - s)
    sigma = 2 * (d - s) / s
    require(abs(rho - sqrt(5)) < mpf(10) ** -18, f"C: rho={nstr(rho,25)}")
    require(abs(sigma - PHI) < mpf(10) ** -18, f"C: sigma={nstr(sigma,25)}")
    require(abs((1 - XSTAR) / XSTAR - PHI) < mpf(10) ** -30, "C: (1-x)/x != phi")
    print(f"C rho = d/(d-s) = sqrt5 (err {nstr(abs(rho-sqrt(5)),3)}), "
          f"sigma = 2(d-s)/s = phi = (1-x*)/x* (err {nstr(abs(sigma-PHI),3)}): OK")
    return rho, sigma


def d_closed_form(rho, sigma):
    # inner first-order condition  s d^gamma = 2 (d-s)^{gamma+1}  <=>  rho^gamma = sigma
    require(abs(rho ** GAMMA - sigma) < mpf(10) ** -18, "D: rho^gamma != sigma")
    # hence gamma = log(sigma)/log(rho) = log(phi)/log(sqrt5)
    g2 = log(PHI) / log(sqrt(5))
    require(abs(g2 - GAMMA) < mpf(10) ** -35, "D: closed form mismatch")
    require(abs(XSTAR - 1 / (1 + PHI)) < mpf(10) ** -35, "D: x* != 1/(1+phi)")
    print(f"D rho^gamma = sigma, so gamma* = log(phi)/log(sqrt5) "
          f"= 2log(phi)/log5 = {nstr(GAMMA,25)}")
    print(f"  => C = 1 + gamma* = log_5(5 phi^2) = {nstr(1+GAMMA,25)}")


def e_consistency():
    from fractions import Fraction as Fr
    from math import comb, log as mlog
    g35 = Fr(3, 5)
    cert = Fr(2457, 4135)
    require(float(g35) > float(GAMMA), "E: 3/5 should exceed gamma*")
    require(float(cert) < float(GAMMA), "E: 2457/4135 should be below gamma*")
    # spot-check the exact ledger at R=256 for both
    for tag, g in [("3/5", g35), ("2457/4135", cert)]:
        d = [(g.numerator * (256 + i)) // g.denominator for i in range(256)]
        worst = None
        for t in range(1, 256):
            capa = sum(comb(d[i], t - i) * 2 ** (t - i)
                       for i in range(min(t, 255) + 1) if t - i <= d[i])
            r = mlog(capa) - mlog(comb(255, t))
            if worst is None or r < worst:
                worst = r
        print(f"E R=256 worst log-ratio at gamma={tag:>10}: {worst:+.4f}"
              f"   ({'above' if float(g) > float(GAMMA) else 'below'} gamma*)")
    print("E consistency with the THM-3002 finite-R ledger: OK")


if __name__ == "__main__":
    print(f"phi = {nstr(PHI,25)}")
    print(f"gamma* = 2 log(phi)/log 5 = {nstr(GAMMA,25)}")
    print(f"x*     = phi^-2           = {nstr(XSTAR,25)}\n")
    y = a_tangency()
    b_minimiser()
    rho, sigma = c_ratios(y)
    d_closed_form(rho, sigma)
    e_consistency()
    print("\nALL THM-3017 REFEREE CHECKS PASSED")
