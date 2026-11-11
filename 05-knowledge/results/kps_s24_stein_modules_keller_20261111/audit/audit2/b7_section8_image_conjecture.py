"""b7_section8_image_conjecture.py - second audit, README section 8, first bullet.

 A. W = x^2 y: 1 = 2 (d/dy - W_y)(y) - (d/dx - W_x)(x); the functional x^k y^l -> l! [k = 2l + 1] kills
    M = Im(d/dx - W_x) + Im(d/dy - W_y) and is 1 on x.  So M is proper and contains 1: not a Mathieu-Zhao subspace.
 B. W = x^p y^q, p, q >= 1: (d/dx - W_x)(alpha x) + (d/dy - W_y)(beta y) = alpha + beta - (p alpha + q beta) W.
    All linear functionals killing M: m(i, j) supported on q (i + 1) = p (j + 1), with gcd(p, q) free parameters
    (elementary; no comparison theorem needed).  So codim M = gcd(p, q) >= 1, and 1 in M iff p != q.
    Machine: dimension of C[x, y]/M in a box, class of 1, for 1 <= p, q <= 4.
 C. The same mechanism in n variables: W = z^a (a monomial with all a_i >= 1, not all equal).
"""
import itertools
from math import gcd
import sympy as sp
import eng2
from eng2 import RingL, add, scal, Echelon

x, y = sp.symbols("x y")


def partA():
    W = x ** 2 * y
    Dx = lambda h: sp.expand(sp.diff(h, x) - sp.diff(W, x) * h)
    Dy = lambda h: sp.expand(sp.diff(h, y) - sp.diff(W, y) * h)
    print("A. W = x^2 y: 2 (d/dy - W_y)(y) - (d/dx - W_x)(x) = %s;  [d/dx - W_x, d/dy - W_y] on x^3 y^2: %s"
          % (sp.expand(2 * Dy(y) - Dx(x)), sp.expand(Dx(Dy(x ** 3 * y ** 2)) - Dy(Dx(x ** 3 * y ** 2)))))

    def lam(h):
        tot = 0
        for (k, l), c in sp.Poly(h, x, y).terms():
            if k == 2 * l + 1:
                tot += c * sp.factorial(l)
        return tot
    bad = 0
    for k, l in itertools.product(range(0, 12), range(0, 8)):
        if lam(Dx(x ** k * y ** l)) != 0 or lam(Dy(x ** k * y ** l)) != 0:
            bad += 1
    print("   functional x^k y^l -> l! [k = 2l + 1] on the generators, k < 12, l < 8: %d failures; value on x: %s; value on 1: %s" % (bad, lam(x), lam(sp.Integer(1))))


def partB():
    al, be, p, q = sp.symbols("alpha beta p q")
    W = x ** p * y ** q
    e = sp.simplify(sp.diff(al * x, x) - sp.diff(W, x) * al * x + sp.diff(be * y, y) - sp.diff(W, y) * be * y - (al + be - (p * al + q * be) * W))
    print("B. (d/dx - W_x)(alpha x) + (d/dy - W_y)(beta y) - [alpha + beta - (p alpha + q beta) W] = %s" % e)
    R = RingL()
    N, Np = 14, 30
    for pp in range(1, 5):
        row = []
        for qq in range(1, 5):
            Wd = R.mono(pp, qq)
            Wx, Wy = R.dx(Wd), R.dy(Wd)
            E = Echelon(lambda m: (m[0] + m[1], m))
            for i in range(Np + 1):
                for j in range(Np + 1 - i):
                    m = R.mono(i, j)
                    E.insert(add(R.dx(m), R.mul(Wx, m), -1))
                    E.insert(add(R.dy(m), R.mul(Wy, m), -1))
            dimV = (N + 1) * (N + 2) // 2
            cod = dimV - sum(1 for m in E.piv if m[0] + m[1] <= N)
            one_in = E.contains(R.one())
            ok = (cod == gcd(pp, qq)) and (one_in == (pp != qq))
            row.append("(%d,%d): codim %d, 1 in M: %s%s" % (pp, qq, cod, "yes" if one_in else "no", "" if ok else "  <== MISMATCH"))
        print("   " + ";  ".join(row))
    print("   expected: codim = gcd(p, q); 1 in M iff p != q.")


def partC():
    z = sp.symbols("z1 z2 z3")
    for a in ((2, 1, 1), (1, 1, 1), (3, 2, 1)):
        W = z[0] ** a[0] * z[1] ** a[1] * z[2] ** a[2]
        D = [lambda h, i=i: sp.expand(sp.diff(h, z[i]) - sp.diff(W, z[i]) * h) for i in range(3)]
        # two weight vectors with sum a_i w_i = 1
        w1 = [sp.Rational(1, a[0]), 0, 0]
        w2 = [0, sp.Rational(1, a[1]), 0]
        e1 = sum(D[i](w1[i] * z[i]) for i in range(3))
        e2 = sum(D[i](w2[i] * z[i]) for i in range(3))
        print("C. W = z^%s: sum D_i(w_i z_i) for two weights: %s and %s; difference (in M): %s" % (a, e1, e2, sp.expand(e1 - e2)))


if __name__ == "__main__":
    partA()
    partB()
    partC()
