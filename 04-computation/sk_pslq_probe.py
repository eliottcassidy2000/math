"""Lane A2 stage 2: verify known closed forms, then probe k=4,5 for elementary forms."""
from mpmath import mp, mpf, sqrt, pi, log, atan, asin, asinh, quad, hyp2f1, cos, sin, catalan, mpc, im, re, pslq, nstr
import sys

def Sval(k, dps):
    old = mp.dps
    mp.dps = dps + 15
    kk = mpf(k)
    g = lambda th: hyp2f1(mpf(1)/2, 2/kk, 1+2/kk, cos(2*th))
    v = (2/pi) * quad(g, [0, pi/4, pi/2])
    mp.dps = old
    return +v

def pslq_report(name, vec, names, dps, maxcoeff=10**12, maxsteps=10**6):
    """Run pslq at dps, then re-run at dps+30 to see if the relation survives."""
    out = []
    for D in (dps, dps+30):
        mp.dps = D
        v = [x(D) for x in vec]
        r = pslq(v, maxcoeff=maxcoeff, maxsteps=maxsteps, tol=mpf(10)**(-(D-12)))
        out.append(r)
    mp.dps = dps
    print("== %s ==" % name)
    print("   basis: %s" % ", ".join(names))
    if out[0] is None:
        print("   NO RELATION at dps=%d (coeff bound %s)" % (dps, maxcoeff))
        return None
    mp.dps = dps + 15
    v = [x(dps+15) for x in vec]
    resid = sum(mpf(c)*x for c, x in zip(out[0], v))
    print("   relation @dps=%d : %s" % (dps, out[0]))
    print("   residual          : %s" % nstr(resid, 8))
    print("   survives +30 dps? : %s   (%s)" % (out[1] == out[0] or (out[1] is not None and _prop(out[0], out[1])), out[1]))
    return out[0]

def _prop(a, b):
    # proportional check
    if len(a) != len(b):
        return False
    ratio = None
    for x, y in zip(a, b):
        if x == 0 and y == 0:
            continue
        if x == 0 or y == 0:
            return False
        r = mpf(y)/mpf(x)
        if ratio is None:
            ratio = r
        elif abs(r-ratio) > mpf(10)**-20:
            return False
    return ratio is not None
