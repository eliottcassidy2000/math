#!/usr/bin/env python3
"""
THM-438 ADDENDUM-14 (monad-explorer 16th session)
=================================================
CLOSED-FORM K-TRANSFORM (inverse Cauchy transform) of the free factorial law
mu_free  (= free compound Poisson of nu=e^{-x}dx; free cumulants kappa_n=n!;
moments A088368 = 1,1,3,13,69,421,...).

The 15th session (ADD-13) inverted the Cauchy transform NUMERICALLY using the
convergent Borel integral
        R(z) = int_0^inf t e^{-t}/(1 - z t) dt           (Borel sum of  sum_{n>=1} n! z^{n-1})
        K(z) = 1/z + R(z),     G(K(z)) = z,    rho(x) = -(1/pi) Im G(x+i0).

CLAIM (this session): that Borel integral has a CLOSED FORM via the exponential
integral Ei.  Reducing the integrand,
        t/(1 - z t) = -1/z + 1/(z (1 - z t)),
so
        R(z) = -1/z + (1/z) int_0^inf e^{-t}/(1 - z t) dt
             = -1/z - (1/z^2) e^{-1/z} Ei(1/z),
and therefore
        K(z) = 1/z + R(z) = -(1/z^2) e^{-1/z} Ei(1/z).

The asymptotic series of  e^{-s}Ei(s) ~ sum_{n>=0} n!/s^{n+1}  (s=1/z -> +inf)
reproduces  K(z) = 1/z + 1 + 2z + 6z^2 + ... = 1/z + sum_{n>=1} n! z^{n-1}, so the
SIGN/branch that matches the cumulants is
        K(z) = (1/z^2) e^{-1/z} Ei(1/z).                 (*)

This is Euler's divergent factorial series resummed: the R-transform of mu_free
is an exponential-integral function.  The classical-side companion is the Gompertz
constant  delta = -e Ei(-1) = int_0^inf e^{-t}/(1+t) dt = 0.5963...

This script:
 (1) VERIFIES (*) against the convergent Borel integral 1/z+R(z) at many complex z.
 (2) Inverts (*) to the density rho and recovers the A088368 moments (closed-form,
     fast: Ei is a builtin, no per-point quadrature).
 (3) Pins the x->0 edge constant  C = lim rho(x) sqrt(x)  and tests a menu of
     candidate constants (refines ADD-12's loose 1/pi and ADD-13's "~0.4-0.6").
 (4) Pins the x->inf tail constant  lim rho(x) e^{x}  (ADD-12 predicted e).
"""
import mpmath as mp

mp.mp.dps = 30
A088368 = [1,1,3,13,69,421,2867,21477,175769,1567273]

# ---------- the two forms of K ----------
def R_integral(z):
    """Convergent Borel integral, valid for z off the positive real axis."""
    f = lambda t: t*mp.e**(-t)/(1 - z*t)
    return mp.quad(f, [0, mp.inf])

def K_integral(z):
    return 1/z + R_integral(z)

def K_closed(z):
    """(*)  K(z) = -(1/z^2) e^{-1/z} E_1(-1/z),  E_1 = principal exponential integral.

    Equivalent to (1/z^2) e^{-1/z} Ei(1/z) on the negative real axis, but E_1 is the
    correct analytic continuation off [0,inf): E_1 is analytic on C \ (-inf,0], and
    -1/z lands off that cut precisely when z is off the support [0,inf).
    """
    s = 1/z
    return -(s*s) * mp.e**(-s) * mp.e1(-s)

# ============================================================
def verify_closed_form():
    print("="*74)
    print("(1) VERIFY  K_closed(z) = (1/z^2) e^{-1/z} Ei(1/z)  ==  1/z + R_integral(z)")
    print("="*74)
    tests = [mp.mpf('-0.5'), mp.mpf('-0.1'), mp.mpf('-2.0'),
             mp.mpf('0.3')+0.2j, mp.mpf('0.1')+0.05j, mp.mpf('0.5')+1j,
             mp.mpf('-0.3')+0.4j, 0.2-0.3j, 1.0+0.7j]
    worst = mp.mpf(0)
    for z in tests:
        kc = K_closed(z)
        ki = K_integral(z)
        d = abs(kc-ki)
        worst = max(worst, d)
        print(f"  z={mp.nstr(z,4):>20}   K_closed={mp.nstr(kc,8):>26}   |closed-integral|={mp.nstr(d,3)}")
    print(f"\n  WORST |closed - integral| = {mp.nstr(worst,4)}   ->  {'MATCH' if worst<mp.mpf('1e-12') else 'MISMATCH'}")

    # asymptotic-series check: K(z) - 1/z should match 1+2z+6z^2+24z^3 (Borel-summed) for small |z|, z complex
    print("\n  Asymptotic-coefficient check (K(z)-1/z vs sum_{n>=1} n! z^{n-1}, z=0.05i):")
    z = mp.mpf('0.05')*1j
    approx = sum(mp.factorial(n)*z**(n-1) for n in range(1,9))   # truncated divergent series ~ ok at tiny z
    print(f"     K_closed-1/z = {mp.nstr(K_closed(z)-1/z,10)}")
    print(f"     trunc series = {mp.nstr(approx,10)}   |diff|={mp.nstr(abs(K_closed(z)-1/z-approx),3)}")
    return worst

# ============================================================
def invert_density():
    print("\n"+"="*74)
    print("(2) DENSITY via closed-form K + moment recovery (A088368)")
    print("="*74)
    eta = mp.mpf('1e-8')

    def solve_G(x, guess):
        z = mp.mpf(x) + 1j*eta
        F = lambda G: K_closed(G) - z
        return mp.findroot(F, guess)

    # grid x=u^2 to tame the 1/sqrt(x) edge; cap x_max~36 (density ~e^{-x}, beyond is below precision)
    Nu = 420
    u_lo, u_hi = mp.mpf('0.012'), mp.mpf('6.0')
    us = [u_lo + (u_hi-u_lo)*mp.mpf(i)/Nu for i in range(Nu+1)]
    xs = [u*u for u in us]

    rho = [mp.mpf(0)]*(Nu+1)
    Gs  = [None]*(Nu+1)
    Gprev = 1/xs[-1] - 1j*mp.mpf('1e-4')
    fails = 0
    for i in range(Nu, -1, -1):
        try:
            G = solve_G(xs[i], Gprev)
            if mp.im(G) > 0:
                G = mp.conj(G)
            Gs[i] = G
            rho[i] = max(mp.mpf(0), -mp.im(G)/mp.pi)
            Gprev = G
        except Exception:
            fails += 1
            Gs[i] = Gprev
            rho[i] = rho[i+1] if i+1 <= Nu else mp.mpf(0)
    print(f"  root-finder fallbacks: {fails} / {Nu+1}   x in [{mp.nstr(xs[0],3)},{mp.nstr(xs[-1],3)}]")

    def moment(k):
        # int x^k rho dx = int u^{2k} rho * 2u du  (trapezoid in u)
        g = [ us[i]**(2*k) * rho[i] * 2*us[i] for i in range(Nu+1) ]
        s = mp.mpf(0)
        for i in range(Nu):
            s += (g[i]+g[i+1])/2 * (us[i+1]-us[i])
        return s

    print("\n   k    integral      A088368     rel.err")
    for k in range(0,5):
        m = moment(k)
        tgt = A088368[k] if k < len(A088368) else None
        if tgt:
            print(f"   {k}   {mp.nstr(m,8):>12}   {tgt:>8}   {mp.nstr(abs(m-tgt)/tgt,3)}")
        else:
            print(f"   {k}   {mp.nstr(m,8):>12}")
    return us, xs, rho, Gs

# ============================================================
def edge_constant(us, xs, rho):
    print("\n"+"="*74)
    print("(3) EDGE constant  C = lim_{x->0+} rho(x) sqrt(x)")
    print("="*74)
    print("   x            rho*sqrt(x)")
    for i in range(1, 16):
        print(f"   {mp.nstr(xs[i],4):>12}   {mp.nstr(rho[i]*mp.sqrt(xs[i]),8)}")
    # extrapolate using the smallest few nodes
    C_est = rho[2]*mp.sqrt(xs[2])
    cands = {
        "1/pi":            1/mp.pi,
        "1/sqrt(2pi)":     1/mp.sqrt(2*mp.pi),
        "1/sqrt(pi e)":    1/mp.sqrt(mp.pi*mp.e),
        "Gompertz/pi":     (-mp.e*mp.ei(-1))/mp.pi,
        "Gompertz":        -mp.e*mp.ei(-1),
        "1/(pi) * e^{1/2}":mp.e**mp.mpf('0.5')/mp.pi,
        "e^{-1/2}":        mp.e**mp.mpf('-0.5'),
        "1/(2)":           mp.mpf('0.5'),
    }
    print(f"\n   C_est (smallest node) ~ {mp.nstr(C_est,8)}")
    print("   candidate            value        |C_est-cand|")
    for name,val in cands.items():
        print(f"   {name:<18} {mp.nstr(val,8):>10}   {mp.nstr(abs(C_est-val),4)}")

def tail_constant(us, xs, rho):
    print("\n"+"="*74)
    print("(4) TAIL constant  lim_{x->inf} rho(x) e^{x}   (ADD-12 predicted e=2.718)")
    print("="*74)
    print("   x            rho*e^{x}")
    N=len(xs)
    for i in range(N-1, N-12, -1):
        print(f"   {mp.nstr(xs[i],4):>12}   {mp.nstr(rho[i]*mp.e**xs[i],8)}")

if __name__ == "__main__":
    w = verify_closed_form()
    us,xs,rho,Gs = invert_density()
    edge_constant(us,xs,rho)
    tail_constant(us,xs,rho)
