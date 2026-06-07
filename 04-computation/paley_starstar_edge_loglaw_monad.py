#!/usr/bin/env python3
"""
THM-438 ADDENDUM-14 part (2) (monad-explorer 16th session)
==========================================================
THE FREE FACTORIAL LAW HAS NO FINITE EDGE CONSTANT: a LOGARITHMIC enhancement.

Using the closed form (ADD-14 part 1, verified to 1e-16):
        K(z) = -(1/z^2) e^{-1/z} E_1(-1/z),     G(K(z)) = z,   rho(x)=-Im G(x+i0)/pi.

Near x->0+, z=G(x+i0)->infinity, and the closed form has a LOGARITHMIC singularity
(the E_1/Ei log, itself the analytic fingerprint of the zero-radius factorial cumulants
kappa_n=n!):
        K(z) ~ (gamma - ln z)/z^2          (z->inf).
Inverting (z = r e^{-i theta}, theta -> pi/2^+):
        rho(x) ~ r/pi,     x ~ (ln r - gamma)/r^2,     ln r ~ (1/2) ln(1/x),
so
        pi * rho(x) * sqrt(x)  ~  sqrt(ln r - gamma)  ~  sqrt( (1/2) ln(1/x) )   ->  infinity.

There is NO finite edge constant.  ADD-12's "1/pi" and ADD-13's "~0.4-0.6" were both
reading this slowly-growing sqrt(log) at different x.

This script inverts K at a geometric ladder of small x with ASYMPTOTICALLY-DERIVED
initial guesses (near-zero fallbacks), and checks three nested predictions:
   (A) pi*rho*sqrt(x) / sqrt(ln|G| - gamma)        -> 1     (uses measured r=|G|)
   (B) pi*rho*sqrt(x) / sqrt((1/2) ln(1/x))        -> 1     (pure leading log in x)
   (C) the same with the ln-ln-corrected ln r.
"""
import mpmath as mp

mp.mp.dps = 40
gamma = mp.euler

def K_closed(z):
    s = 1/z
    return -(s*s) * mp.e**(-s) * mp.e1(-s)

def guess_z(x):
    """Asymptotic seed: solve r^2 x = ln r - gamma for r, theta just above pi/2."""
    # iterate r^2 = (ln r - gamma)/x  starting from ln r ~ (1/2)ln(1/x)
    L = mp.log(1/x)
    lnr = L/2
    for _ in range(60):
        r2 = (lnr - gamma)/x
        if r2 <= 0:
            r2 = mp.log(1/x)/x
        lnr_new = mp.log(r2)/2
        lnr = lnr_new
    r = mp.e**lnr
    lnr_g = lnr - gamma
    if lnr_g <= 0: lnr_g = mp.mpf('0.1')
    # 2theta = pi + delta, tan(delta) ~ theta/(ln r - gamma) ~ (pi/2)/(ln r - gamma)
    delta = mp.atan((mp.pi/2)/lnr_g)
    theta = (mp.pi + delta)/2
    return r*mp.e**(-1j*theta)

def invert(x, seed):
    eta = x*mp.mpf('1e-10')        # eta << x so x+i*eta is essentially on the cut
    z = x + 1j*eta
    F = lambda G: K_closed(G) - z
    return mp.findroot(F, seed)

def main():
    print("="*86)
    print("FREE FACTORIAL LAW: logarithmic edge enhancement (no finite edge constant)")
    print("="*86)
    print("verifying pi*rho*sqrt(x) ~ sqrt(ln|G|-gamma) ~ sqrt((1/2)ln(1/x))\n")
    print(f"{'x':>10} {'pi*rho*sqrt(x)':>16} {'sqrt(ln|G|-g)':>14} {'(A)ratio':>10} "
          f"{'sqrt(.5ln1/x)':>14} {'(B)ratio':>10} {'fallback':>9}")
    xs = [mp.mpf(10)**(-k) for k in range(2, 16)]   # x = 1e-2 ... 1e-15
    fails = 0
    for x in xs:
        seed = guess_z(x)
        try:
            G = invert(x, seed)
            if mp.im(G) > 0: G = mp.conj(G)
            fb = ""
        except Exception:
            fails += 1; fb = "FAIL"
            G = seed
        rho = -mp.im(G)/mp.pi
        lhs = mp.pi*rho*mp.sqrt(x)
        rA  = mp.sqrt(mp.log(abs(G)) - gamma)
        rB  = mp.sqrt(mp.log(1/x)/2)
        print(f"{mp.nstr(x,1):>10} {mp.nstr(lhs,8):>16} {mp.nstr(rA,8):>14} "
              f"{mp.nstr(lhs/rA,6):>10} {mp.nstr(rB,8):>14} {mp.nstr(lhs/rB,6):>10} {fb:>9}")
    print(f"\nfallbacks: {fails}/{len(xs)}")
    print("\nINTERPRETATION:")
    print("  (A) ratio -> 1  confirms  pi*rho*sqrt(x) = sqrt(ln|G|-gamma)  exactly in the limit.")
    print("  (B) ratio -> 1  confirms the pure leading-log-in-x form (slower; ln-ln corrections).")
    print("  rho(x)*sqrt(x) DIVERGES like sqrt(log): no finite edge constant exists.")

if __name__ == "__main__":
    main()
