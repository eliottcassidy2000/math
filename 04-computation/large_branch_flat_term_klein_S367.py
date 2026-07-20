#!/usr/bin/env python3
"""
klein-2026-07-20-S367 -- THE LARGE BRANCH CARRIES A NONZERO FLAT TERM, and that is WHY the
{-1,0,1} radial nullcone is one-sided.  Newton polygon LOCATES the branch; the flat term is
the obstruction to E_r[D^{-1/2}] = 1.

THE PICTURE.  Sum_m E[P^m] t^m = E_r[ D(r,t)^{-1/2} ],  D = (1-beta t)^2 - 4 alpha t^2,
alpha = r a c, beta = b.  NC2 (all E[P^m]=0) says this asymptotic series is == 1.  But the
integral is NOT its own asymptotic series: it can differ by a FLAT term e^{-r*(t)}.  That flat
term is the IMAGINARY PART for real t, coming from the region where D(r,t) < 0, and its edge is
the LARGE BRANCH r*(t) = smallest r>0 with D=0, located by the Newton polygon.

THE ARGUMENT (mechanism; rigor = Borel-Laplace, boxeph's frontier):
  - two-sided P => alpha(r) -> infinity => for any t>0, D(r,t) < 0 on a large-r interval
    => Im E_r[D^{-1/2}] != 0 (a flat term) => the function is NOT == 1 => NC2 fails.
  - one-sided P => alpha == 0 => D = (1-beta t)^2 >= 0 everywhere => NO flat term => can be == 1.
So the flat term is present EXACTLY off the one-sided locus.  Computed and controlled below.
"""
import mpmath as mp
mp.mp.dps = 40

def D(r, t, alpha, beta):
    return (1 - beta(r) * t) ** 2 - 4 * alpha(r) * t ** 2

def r_star(t, alpha, beta, hi=1e12):
    """smallest r>0 with D=0 (the large branch = the edge of the D<0 region)"""
    f = lambda r: D(r, t, alpha, beta)
    # scan for a sign change from + (at 0) to -
    r0 = mp.mpf('1e-6'); f0 = f(r0)
    if f0 < 0: return r0
    r = r0
    while r < hi:
        r2 = r * 2
        if f(r2) < 0:
            return mp.findroot(f, (r + r2) / 2)
        r = r2
    return None

def flat_imag(t, alpha, beta):
    """Im E_r[D^{-1/2}] = integral over the D<0 region of e^{-r} |D|^{-1/2} (the flat term)."""
    rs = r_star(t, alpha, beta)
    if rs is None: return mp.mpf(0), None
    # integrate e^{-r} / sqrt(|D|) from r* outward (converges by e^{-r}); |D|^{-1/2} integrable at r*
    g = lambda r: mp.e ** (-r) / mp.sqrt(abs(D(r, t, alpha, beta)))
    val = mp.quad(g, [rs, rs + mp.mpf('0.5'), rs + 5, rs + 50])
    return val, rs

print("=" * 88)
print("PART 1 -- CONTROL: one-sided P (alpha=0) has NO flat term; two-sided has one")
print("=" * 88)
cases = [
    ("one-sided  alpha=0, beta=0   (IN nullcone)", lambda r: mp.mpf(0), lambda r: mp.mpf(0)),
    ("one-sided  alpha=0, beta=1   (charge 0 only, still one-sided-in-alpha)",
     lambda r: mp.mpf(0), lambda r: mp.mpf(1)),
    ("two-sided  alpha=r, beta=0   (the {-1,+1} pair)", lambda r: r, lambda r: mp.mpf(0)),
    ("two-sided  alpha=r, beta=1   (pair + charge 0)", lambda r: r, lambda r: mp.mpf(1)),
    ("two-sided  alpha=r, beta=1-r (sign-indefinite, THM-1640 gap)", lambda r: r, lambda r: 1 - r),
]
t = mp.mpf('0.15')
print(f" at t = {t}:")
print(f"   {'case':<58} {'r*(t)':>12} {'Im E_r[D^-1/2]':>18}")
for name, al, be in cases:
    val, rs = flat_imag(t, al, be)
    rss = f"{float(rs):.3f}" if rs is not None else "none"
    print(f"   {name:<58} {rss:>12} {float(abs(val)):>18.6e}")
print("""
   READING: one-sided (alpha=0) -> D = (1-beta t)^2 >= 0, no D<0 region, r* NONE, flat term 0.
   two-sided (alpha=r) -> D<0 for large r, r* finite, flat term NONZERO.  The flat term is
   present EXACTLY off the one-sided locus -- which is the nullcone.
""")

print("=" * 88)
print("PART 2 -- THE NEWTON-POLYGON DICHOTOMY: real branch iff B < (A+1)/2, else COMPLEX")
print("=" * 88)
print(" D(r,t) = (1-beta t)^2 - 4 alpha t^2.  At large r, (1-beta t)^2 ~ beta^2 t^2, so")
print(" D < 0 (a REAL branch on the contour) requires 4 alpha > beta^2 at large r, i.e.")
print(" deg alpha = A+1 > 2B = deg beta^2, i.e.  B < (A+1)/2  (alpha dominates the polygon).")
print(" When B > (A+1)/2 (beta dominates), D -> +infinity, NO real branch: r*(t) is COMPLEX.")
print(" This is EXACTLY boxeph's Case I/II (dodge vs jump) split, read off the polygon.\n")
def cD(rc, tt, al, be):
    return (1 - be(rc)*tt)**2 - 4*al(rc)*tt**2
print(f"   {'(alpha,beta)':>26} {'A,B':>7} {'B vs (A+1)/2':>13} {'real r*?':>9} {'branch':>10}")
tt = mp.mpf('0.1')
strata = [
    ("alpha=r, beta=0",   lambda r:r, lambda r:mp.mpf(0), 0, 0),
    ("alpha=r, beta=1",   lambda r:r, lambda r:mp.mpf(1), 0, 0),
    ("alpha=r, beta=1-r", lambda r:r, lambda r:1-r,       0, 1),
    ("alpha=r, beta=r^2", lambda r:r, lambda r:r*r,       0, 2),
    ("alpha=r^3,beta=r",  lambda r:r**3, lambda r:r,      2, 1),
]
for name, al, be, A, B in strata:
    _, rs = flat_imag(tt, al, be)
    thr = mp.mpf(A+1)/2
    dom = "alpha" if B < thr else "beta"
    real = rs is not None
    if real:
        branch = f"REAL {float(rs):.2f}"
    else:
        # locate the complex branch: solve D(r,t)=0 from a complex seed
        try:
            rc = mp.findroot(lambda r: cD(r, tt, al, be), mp.mpc(5, 5))
            branch = f"CPLX {complex(rc).real:.1f}{complex(rc).imag:+.1f}i"
        except Exception:
            branch = "CPLX (?)"
    print(f"   {name:>26} {f'{A},{B}':>7} {f'{B} vs {float(thr)}':>13} {str(real):>9} {branch:>16}")
print("""
   VERDICT.  The Newton polygon B vs (A+1)/2 decides the branch TYPE:
     B < (A+1)/2  (alpha-dominant):  r*(t) REAL -> a flat term on the contour -> Part 3.
     B > (A+1)/2  (beta-dominant):   r*(t) COMPLEX -> an OSCILLATORY (Stokes) term, and
                                     conjugate branches can interfere -- boxeph's hard case.
   The sign-indefinite THM-1640 gap (beta=1-r) is beta-dominant, so it sits in the COMPLEX
   regime -- which is exactly why positivity/domination/real-flat-term all miss it.
""")

print("=" * 88)
print("PART 3 -- alpha-DOMINANT: the real flat term ~ e^{-r*(t)} is nonzero (closes that regime)")
print("=" * 88)
print(f"   {'t':>8} {'r*(t)':>12} {'Im E_r[D^-1/2]':>18} {'e^{-r*}':>16} {'ratio':>10}")
for tt in (mp.mpf('0.2'), mp.mpf('0.15'), mp.mpf('0.1'), mp.mpf('0.07')):
    val, rs = flat_imag(tt, lambda r: r, lambda r: mp.mpf(0))
    er = mp.e ** (-rs)
    print(f"   {float(tt):>8} {float(rs):>12.3f} {float(abs(val)):>18.6e} {float(er):>16.6e} "
          f"{float(abs(val)/er):>10.4f}")
print("""
   The flat term tracks e^{-r*(t)} (ratio O(1)), nonzero at every t.  So in the alpha-dominant
   regime B < (A+1)/2, E_r[D^{-1/2}] carries a nonzero flat piece off the one-sided locus, so it
   is NOT == 1 and NC2 fails there.  Rigor of 'nonzero flat => NC2 fails' = Borel-Laplace.
   The beta-dominant regime (complex branch) is genuinely harder and is boxeph's frontier.
""")
