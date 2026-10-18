"""K-directions with poles above 19 (mac-mini-2026-10-18-sixcolour): for each prime q of K = Q(zeta15, sqrt(-11)) above 19
(residue field F_361 = F_19(t), t^2 = -11), gamma_q = j - c0 - c1 zeta with t = c0 + c1 zeta_bar in F_361 vanishes at q;
delta_q = conj(gamma_q)/gamma_q is a unit vector of K with a pole at the three primes of F above q (one of type h = 7)
unless gamma_q also vanishes at conj(q). Returns the list of (label, delta) for one prime of each conjugate pair."""
from fractions import Fraction as Q
P = 19; NR = (-11) % P
def fm(x, y): return ((x[0] * y[0] + NR * x[1] * y[1]) % P, (x[0] * y[1] + x[1] * y[0]) % P)
def fpow(x, k):
    r = (1, 0)
    for _ in range(k): r = fm(r, x)
    return r
def kdirections(geo):
    f = geo.f
    els = [(a, b) for a in range(P) for b in range(P)]
    z15 = [x for x in els if fpow(x, 15) == (1, 0) and fpow(x, 5) != (1, 0) and fpow(x, 3) != (1, 0)]
    def hval(x, z, j):                          # evaluate a K-element (w-exponent 0) at zeta -> z, j -> j
        s = (0, 0)
        for c, e in zip(x, f.basis):
            if c:
                assert e[2] == 0 and c.denominator % P
                cc = c.numerator * pow(c.denominator, -1, P) % P
                v = fm(fpow(z, e[0]), fpow(j, e[1])); s = ((s[0] + cc * v[0]) % P, (s[1] + cc * v[1]) % P)
        return s
    zeta = f.gen(0); out = []; seen = set()
    for z0 in z15:
        for j0 in [(0, 1), (0, P - 1)]:
            key = frozenset([(z0, j0), (fpow(z0, P), (0, (-j0[1]) % P))])
            ckey = frozenset([(fpow(z0, 14), (0, (-j0[1]) % P)), (fpow(fpow(z0, 14), P), j0)])   # complex conjugate prime
            if key in seen or ckey in seen: continue
            seen.add(key)
            # solve j0 = c0 + c1 z0 in F_361 (basis 1, t)
            c1 = (j0[1] * pow(z0[1], -1, P)) % P; c0 = (j0[0] - c1 * z0[0]) % P
            lift = lambda c: c - P if c > P // 2 else c
            gamma = f.sub(geo.j, f.add(f.scalar(lift(c0)), f.scale(zeta, lift(c1))))
            assert hval(gamma, z0, j0) == (0, 0)
            zc, jc = fpow(z0, 14), (0, (-j0[1]) % P)
            if hval(gamma, zc, jc) == (0, 0): continue
            delta = f.mul(geo.bar(gamma), f.inv(gamma))
            assert f.mul(delta, geo.bar(delta)) == f.one
            out.append(('q(z=%s,j=%s) c0=%d c1=%d' % (z0, j0, lift(c0), lift(c1)), delta))
    return out
