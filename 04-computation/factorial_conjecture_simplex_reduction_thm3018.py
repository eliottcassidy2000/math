#!/usr/bin/env python3
"""Referee for THM-3018: the Factorial Conjecture is a simplex moment problem.

FC(n) (van den Essen-Wright-Zhao): L(x^alpha) = alpha! = prod_i alpha_i!;
for f HOMOGENEOUS in n variables, L(f^m) = 0 for all m >= 1 implies f = 0.

Checks:
  R1  Integral form of L and the polar reduction:
      L(f^m) = (dm+n-1)! * int_{Delta_{n-1}} g^m dA,  g = f restricted to the
      simplex, verified symbolically for n = 2 (d = 1,2,3; m = 1,2,3) and
      n = 3 (d = 1,2; m = 1,2).
  R2  f |-> f|_Delta is a linear BIJECTION from degree-d forms in n variables
      onto polynomials of degree <= d in n-1 variables (Bernstein basis), so
      the reduction loses nothing.
  R3  FC(2) for low degree: the moment system int_0^1 g^m du = 0 (m up to
      d+3) has ONLY g = 0, by exact Groebner elimination, for d = 1,2,3.
  R4  The S_3 mechanism for FC(3): the 3-cycle sigma on barycentric
      coordinates is an area-preserving self-map of the triangle; if
      g o sigma = omega g (omega = e^{2 pi i/3}) then int_T g^m dA = 0
      automatically for every m not divisible by 3.  Verified exactly, in
      rational arithmetic over Q(omega), on the degree-1 eigenvector
      g = u + omega v + omega^2 w, together with the fact that the surviving
      moments m = 3,6,9 are NONZERO for that g.
"""

from fractions import Fraction as Fr
from itertools import combinations
from math import factorial

import sympy as sp

u, v, x, y, z = sp.symbols('u v x y z')


def require(c, m):
    if not c:
        raise RuntimeError(m)


def L_poly(poly, gens):
    P = sp.Poly(sp.expand(poly), *gens)
    tot = 0
    for mon, c in P.terms():
        t = c
        for a in mon:
            t *= factorial(a)
        tot += t
    return sp.expand(tot)


def r1():
    for d, cs in [(1, [2, -3]), (2, [1, sp.Rational(5, 2), -1]), (3, [1, -1, 2, 1])]:
        f = sum(cs[i] * x ** i * y ** (d - i) for i in range(d + 1))
        g = sp.expand(f.subs({x: u, y: 1 - u}))
        for m in (1, 2, 3):
            lhs = L_poly(sp.expand(f ** m), (x, y))
            rhs = factorial(d * m + 1) * sp.integrate(sp.expand(g ** m), (u, 0, 1))
            require(sp.simplify(lhs - rhs) == 0, f"R1 n=2 d={d} m={m}")
    for d, cs in [(1, {(1, 0, 0): 2, (0, 1, 0): -1, (0, 0, 1): 3}),
                  (2, {(2, 0, 0): 1, (1, 1, 0): -2, (0, 0, 2): 1, (0, 1, 1): 3})]:
        f = sum(c * x ** a * y ** b * z ** cc for (a, b, cc), c in cs.items())
        g = sp.expand(f.subs({x: u, y: v, z: 1 - u - v}))
        for m in (1, 2):
            lhs = L_poly(sp.expand(f ** m), (x, y, z))
            inner = sp.integrate(sp.expand(g ** m), (v, 0, 1 - u))
            rhs = factorial(d * m + 2) * sp.integrate(sp.expand(inner), (u, 0, 1))
            require(sp.simplify(lhs - rhs) == 0, f"R1 n=3 d={d} m={m}")
    print("R1 L(f^m) = (dm+n-1)! * int_Delta g^m  (n=2 d<=3, n=3 d<=2): OK")


def r2():
    for n_, d_ in [(2, 3), (2, 5), (3, 2), (3, 3)]:
        if n_ == 2:
            forms = [x ** i * y ** (d_ - i) for i in range(d_ + 1)]
            restr = [sp.expand(F.subs({x: u, y: 1 - u})) for F in forms]
            M = sp.Matrix([[sp.Poly(r, u).coeff_monomial(u ** k)
                            for k in range(d_ + 1)] for r in restr])
        else:
            forms = [x ** a * y ** b * z ** (d_ - a - b)
                     for a in range(d_ + 1) for b in range(d_ + 1 - a)]
            restr = [sp.expand(F.subs({x: u, y: v, z: 1 - u - v})) for F in forms]
            mons = [(a, b) for a in range(d_ + 1) for b in range(d_ + 1 - a)]
            M = sp.Matrix([[sp.Poly(r, u, v).coeff_monomial(u ** a * v ** b)
                            for (a, b) in mons] for r in restr])
        require(M.rank() == M.shape[0] == M.shape[1], f"R2 n={n_} d={d_}")
    print("R2 restriction to the simplex is a bijection onto polys of deg <= d: OK")


def r3():
    for d in (1, 2, 3):
        cs = sp.symbols(f'c0:{d+1}')
        g = sum(cs[i] * u ** i for i in range(d + 1))
        eqs, gm = [], sp.Integer(1)
        for m in range(1, d + 4):
            gm = sp.expand(gm * g)
            eqs.append(sp.expand(sp.integrate(gm, (u, 0, 1))))
        sols = sp.solve(eqs, cs, dict=True)
        require(sols and all(all(sp.simplify(s.get(c, 0)) == 0 for c in cs)
                             for s in sols), f"R3 d={d}: nonzero solution!")
        print(f"R3 FC(2) degree {d}: moment system forces g = 0 (exact Groebner): OK")


def r4():
    w = sp.Rational(-1, 2) + sp.sqrt(3) * sp.I / 2          # primitive cube root
    require(sp.simplify(w ** 3 - 1) == 0 and sp.simplify(w - 1) != 0, "R4 omega")

    def bint(a, b, c):
        return Fr(factorial(a) * factorial(b) * factorial(c),
                  factorial(a + b + c + 2))

    def mult(p, q):
        r = {}
        for k1, c1 in p.items():
            for k2, c2 in q.items():
                k = (k1[0] + k2[0], k1[1] + k2[1], k1[2] + k2[2])
                r[k] = sp.expand(r.get(k, 0) + c1 * c2)
        return r

    def integ(p):
        return sp.expand(sum(c * sp.Rational(bint(*k)) for k, c in p.items()))

    g = {(1, 0, 0): sp.Integer(1), (0, 1, 0): w, (0, 0, 1): w ** 2}
    # sigma: (a,b,c) -> cycle; check g o sigma = omega^2 g  (so moments with
    # 3 not dividing m must vanish)
    gs = {(k[2], k[0], k[1]): c for k, c in g.items()}
    for k in set(list(g) + list(gs)):
        require(sp.simplify(gs.get(k, 0) - w ** 2 * g.get(k, 0)) == 0,
                "R4 g o sigma != omega^2 g")
    cur = {(0, 0, 0): sp.Integer(1)}
    vals = []
    for m in range(1, 10):
        cur = mult(cur, g)
        vals.append(sp.simplify(integ(cur)))
    for m, val in enumerate(vals, 1):
        if m % 3:
            require(sp.simplify(val) == 0, f"R4 moment m={m} should vanish")
        else:
            require(sp.simplify(val) != 0, f"R4 moment m={m} should NOT vanish")
    print("R4 S_3 mechanism: g o sigma = omega^2 g forces int_T g^m = 0 for "
          "3 not| m, exactly (m<=9); m = 3,6,9 survive nonzero: OK")
    print(f"   surviving values: m=3 -> {vals[2]}, m=6 -> {vals[5]}, m=9 -> {vals[8]}")


if __name__ == "__main__":
    r1(); r2(); r3(); r4()
    print("ALL THM-3018 REFEREE CHECKS PASSED")
