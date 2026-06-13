#!/usr/bin/env python3
"""
THM-438 ADDENDUM-13 follow-up (monad-explorer 15th session)
===========================================================
NUMERICAL CAUCHY-TRANSFORM INVERSION of the free factorial law mu_free
(= free compound Poisson of nu=e^{-x}dx; free cumulants kappa_n=n!; moments A088368).

ADD-12 characterized mu_free only by its EDGE (1/sqrt x at 0) and TAIL (~e e^{-x}).
The Bercovici-Pata frame (ADD-13) gives the explicit R-transform, so we can now
compute the actual density and VALIDATE it by recovering the A088368 moments.

  R(w) = \int_0^inf t e^{-t}/(1-w t) dt        (Borel sum of sum n! w^{n-1})
  K(w) = 1/w + R(w)                            (inverse Cauchy / K-transform)
  Cauchy transform G solves K(G(z)) = z.
  density  rho(x) = -(1/pi) Im G(x + i0).

We solve K(G) = x + i*eta for small eta>0 (physical sheet, Im G<0), marching x
downward with continuation. Then validate:
  total mass  \int rho dx        -> 1
  moments     \int x^k rho dx    -> A088368(k) = 1,1,3,13,69,421,...
  edge        rho(x) sqrt(x)     -> 1/pi  as x->0+   (critical free Poisson)
  tail        rho(x) e^{x}       -> e      as x->+inf (Levy-measure tail)
"""
import mpmath as mp

mp.mp.dps = 25
A088368 = [1,1,3,13,69,421,2867,21477,175769,1567273]

def R(w):
    # \int_0^inf t e^{-t}/(1 - w t) dt ; complex w off the positive real axis -> honest integral
    f = lambda t: t*mp.e**(-t)/(1 - w*t)
    return mp.quad(f, [0, mp.inf])

def K(w):
    return 1/w + R(w)

def solve_G(x, eta, guess):
    z = mp.mpf(x) + 1j*mp.mpf(eta)
    F = lambda G: K(G) - z
    return mp.findroot(F, guess)

def main():
    print("="*74)
    print("Free factorial law mu_free: Cauchy-transform inversion + moment validation")
    print("="*74)
    eta = mp.mpf('1e-6')

    # ---- build a grid in u with x=u^2 (tames the 1/sqrt x edge for integration) ----
    # x from ~0.0009 to ~49
    Nu = 240
    u_lo, u_hi = mp.mpf('0.03'), mp.mpf('7.0')
    us = [u_lo + (u_hi-u_lo)*mp.mpf(i)/(Nu) for i in range(Nu+1)]
    xs = [u*u for u in us]

    # ---- continuation from the largest x downward ----
    print(f"\nInverting K(G)=x+i{mp.nstr(eta,1)} on {Nu+1} nodes x in [{mp.nstr(xs[0],3)}, {mp.nstr(xs[-1],3)}] ...")
    # initial guess at large x: G ~ 1/x with small negative imaginary part
    guess = 1/(xs[-1]) - 1j*mp.mpf('1e-3')
    rho = [None]*(Nu+1)
    Gprev = guess
    fails = 0
    for i in range(Nu, -1, -1):
        x = xs[i]
        try:
            G = solve_G(x, eta, Gprev)
            im = mp.im(G)
            if im > 0:   # wrong sheet; reflect
                G = mp.conj(G); im = mp.im(G)
            rho[i] = max(mp.mpf(0), -im/mp.pi)
            Gprev = G
        except Exception:
            fails += 1
            rho[i] = rho[i+1] if i+1 <= Nu and rho[i+1] is not None else mp.mpf(0)
    print(f"  root-finder fallbacks: {fails}")

    # ---- integrate moments over x using the u-substitution: \int x^k rho dx = \int (u^2)^k rho * 2u du ----
    def moment(k):
        vals = [ (us[i]**(2*k)) * rho[i] * 2*us[i] for i in range(Nu+1) ]
        # trapezoid in u
        s = mp.mpf(0)
        for i in range(Nu):
            s += (us[i+1]-us[i])*(vals[i]+vals[i+1])/2
        return s

    print(f"\nMoment recovery  \\int x^k rho(x) dx  vs  A088368(k):")
    print(f"  {'k':>2} {'numeric':>16} {'A088368':>10} {'rel.err':>10}")
    for k in range(0, 6):
        mk = moment(k)
        ref = A088368[k]
        rel = abs(mk-ref)/ref
        print(f"  {k:>2} {mp.nstr(mk,10):>16} {ref:>10} {mp.nstr(rel,3):>10}")

    # ---- edge and tail diagnostics ----
    print(f"\nEdge  rho(x)*sqrt(x) near 0  (expect -> 1/pi = {mp.nstr(1/mp.pi,6)}):")
    for i in [0,1,2,4,8]:
        print(f"   x={mp.nstr(xs[i],4):>8}  rho*sqrt(x)={mp.nstr(rho[i]*mp.sqrt(xs[i]),6)}")
    print(f"\nTail  rho(x)*e^{{x}} for large x  (expect -> e = {mp.nstr(mp.e,6)}):")
    for i in [Nu, Nu-10, Nu-25, Nu-50]:
        print(f"   x={mp.nstr(xs[i],4):>8}  rho*e^x={mp.nstr(rho[i]*mp.e**(xs[i]),6)}")

    # ---- a few density samples for the record ----
    print(f"\nDensity samples rho(x):")
    for xt in ['0.1','0.5','1.0','2.0','4.0','8.0','16.0']:
        x = mp.mpf(xt)
        # find nearest node
        j = min(range(Nu+1), key=lambda i: abs(xs[i]-x))
        print(f"   x~{mp.nstr(xs[j],4):>8}  rho={mp.nstr(rho[j],6)}")

    print("\nDONE.")

if __name__ == "__main__":
    main()
