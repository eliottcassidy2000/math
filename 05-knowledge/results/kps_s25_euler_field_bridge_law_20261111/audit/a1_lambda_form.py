"""a1_lambda_form.py - audit of Theorem A2 (1): the linear form lambda on C[x, y] for W = x + x^(l+1) y.

Independent of the author's engine: exact rational arithmetic (fractions), polynomials as dicts {(k, j): coeff}.

Checks
  L1  lambda((d/dx - W_x) x^k y^j) = 0 and lambda((d/dy - W_y) x^k y^j) = 0 for l = 1..4, k < 14, j < 6
  L2  lambda(1) = (-1)^l / l!, lambda(x^l) = 1
  L3  lambda(W h) = (e + delta) lambda(h) for monomials h = x^k y^j, e = k - l j, delta = 1 - l   (Lemma A1, Euler)
  L4  lambda(W^m) = 0 exactly for m >= l, nonzero for m < l;  lambda(x^l W^n) = n!
  L5  the formula for lambda is the residue of the exponential integral: sympy residue of
      j! x^(k - (l+1)(j+1)) exp(-x) at x = 0 (the inner integral int_0^oo y^j exp(-x^(l+1) y) dy = j!/x^((l+1)(j+1)))
  L6  symbolic identity behind L1 for general l (hand computation, checked here with sympy gamma functions):
      k/N! + 1/(N-1)! - (l+1)(j+1)/N! = 0  with N = (l+1)(j+1) - k   [after dividing by (-1)^N j!]
"""
from fractions import Fraction
from math import factorial

import sympy as sp


def lam(k, j, l):
    n = (l + 1) * (j + 1) - 1 - k
    if n < 0:
        return Fraction(0)
    return Fraction((-1) ** n * factorial(j), factorial(n))


def lam_poly(P, l):
    return sum(c * lam(k, j, l) for (k, j), c in P.items())


def add(P, Q, s=1):
    R = dict(P)
    for m, c in Q.items():
        R[m] = R.get(m, 0) + s * c
        if R[m] == 0:
            del R[m]
    return R


def mul(P, Q):
    R = {}
    for (a, b), c in P.items():
        for (a2, b2), c2 in Q.items():
            m = (a + a2, b + b2)
            R[m] = R.get(m, 0) + c * c2
    return {m: c for m, c in R.items() if c != 0}


def dx(P):
    return {(k - 1, j): c * k for (k, j), c in P.items() if k > 0}


def dy(P):
    return {(k, j - 1): c * j for (k, j), c in P.items() if j > 0}


def Wpoly(l):
    return {(1, 0): Fraction(1), (l + 1, 1): Fraction(1)}


def main():
    # L1
    tot = bad = 0
    for l in range(1, 5):
        W = Wpoly(l)
        Wx, Wy = dx(W), dy(W)
        for k in range(14):
            for j in range(6):
                h = {(k, j): Fraction(1)}
                gx = add(dx(h), mul(Wx, h), -1)
                gy = add(dy(h), mul(Wy, h), -1)
                for gen in (gx, gy):
                    tot += 1
                    bad += (lam_poly(gen, l) != 0)
    print("L1  lambda kills (d/dx - W_x) x^k y^j and (d/dy - W_y) x^k y^j, l = 1..4, k < 14, j < 6: %d generators, %d failures" % (tot, bad))

    # L2
    ok = all(lam(0, 0, l) == Fraction((-1) ** l, factorial(l)) and lam(l, 0, l) == 1 for l in range(1, 9))
    print("L2  lambda(1) = (-1)^l/l!, lambda(x^l) = 1 for l = 1..8:", ok,
          "; values of lambda(1):", [str(lam(0, 0, l)) for l in range(1, 6)])

    # L3
    tot = bad = 0
    for l in range(1, 5):
        W = Wpoly(l)
        for k in range(14):
            for j in range(6):
                h = {(k, j): Fraction(1)}
                lhs = lam_poly(mul(W, h), l)
                rhs = (k - l * j + 1 - l) * lam(k, j, l)
                tot += 1
                bad += (lhs != rhs)
    print("L3  lambda(W h) = (k - l j + 1 - l) lambda(h) on monomials: %d checked, %d failures" % (tot, bad))

    # L4
    for l in range(1, 5):
        W = Wpoly(l)
        Wm = {(0, 0): Fraction(1)}
        vals = []
        xl = {(l, 0): Fraction(1)}
        okx = True
        for m in range(1, 9):
            Wm = mul(Wm, W)
            vals.append(lam_poly(Wm, l))
            okx = okx and (lam_poly(mul(xl, Wm), l) == factorial(m))
        pred = []
        for m in range(1, 9):
            v = Fraction((-1) ** l, factorial(l))
            for t in range(m):
                v *= (1 - l + t)
            pred.append(v)
        print("L4  l = %d: lambda(W^m), m = 1..8: %s ; equals (1-l)...(m-l) lambda(1): %s ; lambda(x^l W^n) = n! for n = 1..8: %s"
              % (l, [str(v) for v in vals], vals == pred, okx))

    # L5
    x = sp.symbols("x")
    tot = bad = 0
    for l in range(1, 4):
        for k in range(0, 9):
            for j in range(0, 4):
                r = sp.residue(sp.factorial(j) * x ** (k - (l + 1) * (j + 1)) * sp.exp(-x), x, 0)
                tot += 1
                bad += (sp.nsimplify(r) != sp.Rational(lam(k, j, l).numerator, lam(k, j, l).denominator))
    print("L5  lambda(x^k y^j) = Res_(x=0) j! x^(k-(l+1)(j+1)) e^(-x), l = 1..3, k < 9, j < 4: %d checked, %d failures" % (tot, bad))
    yy, t = sp.symbols("y t", positive=True)
    inner = [sp.integrate(yy ** j * sp.exp(-t * yy), (yy, 0, sp.oo)) for j in range(4)]
    print("    inner integral int_0^oo y^j e^(-t y) dy, j = 0..3:", inner)

    # L6
    l_, j_, k_ = sp.symbols("l j k")
    N = (l_ + 1) * (j_ + 1) - k_
    expr = k_ / sp.gamma(N + 1) + 1 / sp.gamma(N) - (l_ + 1) * (j_ + 1) / sp.gamma(N + 1)
    print("L6  k/N! + 1/(N-1)! - (l+1)(j+1)/N! with N = (l+1)(j+1) - k simplifies to:", sp.simplify(sp.gammasimp(expr)))


if __name__ == "__main__":
    main()
