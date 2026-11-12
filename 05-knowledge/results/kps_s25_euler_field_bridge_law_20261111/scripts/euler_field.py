"""euler_field.py - the Euler-field lemma for twisted de Rham classes (kind-pasteur-2026-11-11-S25).

Lemma.  W a polynomial on C^2, xi a polynomial vector field with xi(W) = W and constant divergence delta.  Then
W h = (xi + delta) h modulo M_W = Im(d/dx - W_x) + Im(d/dy - W_y), for every polynomial h.  So for an eigen-polynomial
(xi h = e h):  [W^m h] = (e + delta)(e + delta + 1)...(e + delta + m - 1) [h].

Checks (exact arithmetic modulo the prime of the S24 engine; membership is relative to the truncation):
  T1  W = x + x^2 y (Broughton), xi = x d/dx - y d/dy, delta = 0: [W^m] = 0 and [x W^m] = m! [x] != 0, m <= 4:
      the kernel is not a Mathieu-Zhao subspace although it does not contain 1 (HYP-9668 (E') is false);
  T2  W = x + x^3 y, xi = x d/dx - 2 y d/dy, delta = -1: [W] = -[1], [W^2] = 0, [x^2 W^n] = n! [x^2];
  T3  W = x y, xi = (x d/dx + y d/dy)/2, delta = 1;
  T4  a Keller automorphism pair (f, g) = (x + y^2, y), centre (2, 3): [U V h] = [(xi + 1) h], xi = (U d_u + V d_v)/2;
  T5  the slice map of THM-1300 (S24 note, section 6.2), centre (2, 3): the same identity.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..", "kps_s24_stein_modules_keller_20261111", "scripts")))
from stein import Ring, twisted, codim, P          # noqa: E402
from stein_more import slice_data                   # noqa: E402

INV2 = pow(2, P - 2, P)


def fact(n):
    r = 1
    for i in range(2, n + 1):
        r *= i
    return r


def zero(E, el, R, wmax):
    """True if el reduces to zero; None if its weight exceeds wmax."""
    if not el:
        return True
    if max(Ring.weight(m) for m in el) > wmax:
        return None
    return not E.normal_form(el)


def t1():
    R = Ring("A")
    x, y, one = R.mono(1, 0), R.mono(0, 1), R.mono(0, 0)
    W = R.add(x, R.mul(R.mul(x, x), y))
    E = twisted(R, W, 40)
    print("T1  W = x + x^2 y: codimension c(N) =", [codim(R, E, N) for N in (6, 12, 18, 24)],
          "; class of 1:", "nonzero" if E.normal_form(one) else "zero", "; class of x:", "nonzero" if E.normal_form(x) else "zero")
    bad = tot = 0
    for k in range(0, 5):
        for l in range(0, 4):
            h = R.mono(k, l)
            z = zero(E, R.add(R.mul(W, h), R.scal(h, k - l), -1), R, 30)
            if z is None:
                continue
            tot += 1
            bad += (not z)
    print("    W h = (k - l) h modulo M for monomials h = x^k y^l: %d checked, %d failures" % (tot, bad))
    res1, res2, res3 = [], [], []
    Wm = one
    for m in range(1, 5):
        Wm = R.mul(Wm, W)
        res1.append(zero(E, Wm, R, 30))
        xWm = R.mul(x, Wm)
        res2.append(zero(E, R.add(xWm, R.scal(x, fact(m)), -1), R, 30))
        res3.append(zero(E, xWm, R, 30))
    print("    W^m in M, m = 1..4:", res1, "; x W^m - m! x in M:", res2, "; x W^m in M:", res3)


def t2():
    R = Ring("A")
    x, y, one = R.mono(1, 0), R.mono(0, 1), R.mono(0, 0)
    W = R.add(x, R.mul(R.power(x, 3), y))
    E = twisted(R, W, 44)
    x2 = R.mul(x, x)
    print("T2  W = x + x^3 y: codimension c(N) =", [codim(R, E, N) for N in (6, 12, 18, 24)],
          "; class of 1:", "nonzero" if E.normal_form(one) else "zero", "; class of x^2:", "nonzero" if E.normal_form(x2) else "zero")
    bad = tot = 0
    for k in range(0, 5):
        for l in range(0, 3):
            h = R.mono(k, l)
            z = zero(E, R.add(R.mul(W, h), R.scal(h, k - 2 * l - 1), -1), R, 32)
            if z is None:
                continue
            tot += 1
            bad += (not z)
    print("    W h = (k - 2l - 1) h modulo M for monomials: %d checked, %d failures" % (tot, bad))
    W2 = R.mul(W, W)
    print("    W + 1 in M:", zero(E, R.add(W, one), R, 32), "; W^2 in M:", zero(E, W2, R, 32), "; W^4 in M:", zero(E, R.mul(W2, W2), R, 32))
    res = []
    Wn = one
    for n in range(1, 5):
        Wn = R.mul(Wn, W)
        res.append(zero(E, R.add(R.mul(x2, Wn), R.scal(x2, fact(n)), -1), R, 32))
    print("    x^2 W^n - n! x^2 in M, n = 1..4:", res)


def t3():
    R = Ring("A")
    x, y = R.mono(1, 0), R.mono(0, 1)
    W = R.mul(x, y)
    E = twisted(R, W, 30)
    bad = tot = 0
    for k in range(0, 5):
        for l in range(0, 5):
            h = R.mono(k, l)
            # 2 W h = (k + l + 2) h
            z = zero(E, R.add(R.scal(R.mul(W, h), 2), R.scal(h, k + l + 2), -1), R, 22)
            if z is None:
                continue
            tot += 1
            bad += (not z)
    print("T3  W = x y: 2 W h = (k + l + 2) h modulo M for monomials: %d checked, %d failures" % (tot, bad))


def t4():
    R = Ring("A")
    x, y, one = R.mono(1, 0), R.mono(0, 1), R.mono(0, 0)
    f = R.add(x, R.mul(y, y))
    g = y
    du = lambda h: R.dx(h)                                           # d/df at fixed g
    dv = lambda h: R.add(R.dy(h), R.scal(R.mul(y, R.dx(h)), 2), -1)   # d/dg at fixed f
    assert du(f) == one and dv(g) == one and not du(g) and not dv(f)
    U = R.add(f, R.scal(one, 2), -1)
    V = R.add(g, R.scal(one, 3), -1)
    Wp = R.mul(U, V)
    E = twisted(R, Wp, 30)
    bad = tot = 0
    for h in (x, y, R.mul(x, y), R.mul(x, x), R.mul(y, y), R.mul(R.mul(x, x), y), one):
        lhs = R.scal(R.mul(Wp, h), 2)
        rhs = R.add(R.add(R.mul(U, du(h)), R.mul(V, dv(h))), R.scal(h, 2))
        z = zero(E, R.add(lhs, rhs, -1), R, 22)
        if z is None:
            continue
        tot += 1
        bad += (not z)
    print("T4  Keller pair (x + y^2, y), centre (2, 3): 2 U V h = (U d_u + V d_v + 2) h modulo M_z: %d checked, %d failures" % (tot, bad))


def t5():
    R, x, y, w, one, u, v, J = slice_data()
    Jinv = R.scal(R.power(R.add(one, R.mul(x, y)), 3), -INV2)
    du = lambda h: R.mul(R.add(R.mul(R.dx(h), R.dy(v)), R.mul(R.dy(h), R.dx(v)), -1), Jinv)
    dv = lambda h: R.mul(R.add(R.mul(R.dx(u), R.dy(h)), R.mul(R.dy(u), R.dx(h)), -1), Jinv)
    U = R.add(u, R.scal(one, 2), -1)
    V = R.add(v, R.scal(one, 3), -1)
    Wp = R.mul(U, V)
    E = twisted(R, Wp, 34)
    bad = tot = 0
    for h in (one, w, x, y, R.mul(w, w), R.mul(x, w)):
        lhs = R.scal(R.mul(Wp, h), 2)
        rhs = R.add(R.add(R.mul(U, du(h)), R.mul(V, dv(h))), R.scal(h, 2))
        z = zero(E, R.mul(R.add(lhs, rhs, -1), J), R, 26)
        if z is None:
            continue
        tot += 1
        bad += (not z)
    print("T5  slice map, centre (2, 3): 2 U V h = (U d_u + V d_v + 2) h modulo M_z: %d checked, %d failures" % (tot, bad))


if __name__ == "__main__":
    t1()
    t2()
    t3()
    t4()
    t5()
