#!/usr/bin/env python3
"""Exact virtual/actual doubling controls; universal rootwise transport is OPEN.

The companion note proves the virtual sign and coefficientwise path inclusion.
SymPy supplies exact polynomial operations and real-root isolation. Independent
path dynamic programming and Fraction quotient recurrences check the two new
maps. No repository mathematical implementation is imported.
"""

from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from math import comb, factorial, gcd
from random import Random

import sympy as sp


t, parameter = sp.symbols("t g")
CHECKS = 0
TRACE = sha256()


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def record(*row):
    TRACE.update((repr(row)+"\n").encode())


def binomial(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def polynomial(values):
    return sp.Poly.from_list(list(reversed(values)), t, domain=sp.QQ)


def rows(A, B, h, r, z, x):
    delta, y = B-A, B*h+r
    mass, L = x+y+z, x//delta

    def alpha(j):
        return binomial(mass, z+A*j) if j >= 0 else 0

    def beta(j):
        u, v = x+delta*j, y-B*j
        return binomial(u+v, u) if min(u, v) >= 0 else 0

    P = polynomial([alpha(j)*beta(j) for j in range(h+1)])
    V = polynomial([sum(alpha(i)*alpha(j-i) for i in range(j+1)) *
                    sum(beta(i)*beta(j-i) for i in range(max(-L, j-h), min(h, j+L)+1))
                    for j in range(2*h+1)])
    actual = {}
    for j in range(-1, 2*h+2):
        u, v, w = 2*x+delta*j, 2*y-B*j, 2*z+A*j
        if min(u, v, w) >= 0:
            actual[j] = factorial(2*mass)//(factorial(u)*factorial(v)*factorial(w))
    epsilon = 2*z//A
    Q = polynomial([actual.get(j-epsilon, 0) for j in range(2*h+2+epsilon)])
    return mass, L, P, V, Q, epsilon, actual, alpha, beta


def interval_product(left, right):
    values = [a*b for a in left for b in right]
    return min(values), max(values)


def interval_value(p, a, b):
    value = (F(0), F(0))
    for c in p.all_coeffs():
        value = interval_product(value, (F(a), F(b)))
        value = value[0]+F(c), value[1]+F(c)
    return value


def signs_at_roots(P, target):
    require(P.gcd(target).degree() == 0, "sign certificate has no shared root")
    remainder = target.rem(P)
    signs, certificates = [], []
    for (a, b), multiplicity in P.intervals():
        require(multiplicity == 1 and b <= 0 and P.eval(0) != 0, "simple negative first root")
        for _ in range(240):
            lo, hi = interval_value(remainder, a, b)
            if lo > 0 or hi < 0:
                break
            a, b = P.refine_root(a, b, eps=(b-a)/4)
        else:
            raise RuntimeError("root-value refinement exhausted")
        signs.append(1 if lo > 0 else -1)
        certificates.append(tuple(map(str, (a, b, lo, hi))))
    require(len(signs) == P.degree(), "every first root certified")
    return signs, certificates


def path_count(U, V, selected):
    if min(U, V) < 0:
        return 0, 0
    table = {}
    for i in range(U+1):
        for j in range(V+1):
            if i == j == 0:
                pair = (1, 0)
            else:
                left = table.get((i-1, j), (0, 0))
                below = table.get((i, j-1), (0, 0))
                pair = left[0]+below[0], left[1]+below[1]
            if selected(i, j):
                pair = 0, pair[0]+pair[1]
            table[i, j] = pair
    no_hit, hit = table[U, V]
    return no_hit+hit, hit


def path_controls():
    cases = endpoints = 0
    for B in range(2, 5):
        for A in range(1, B):
            delta = B-A
            for h in (1, 2):
                for x in (0, 1, 2):
                    for r in (0, B-1):
                        for z in (0, A-1):
                            mass, L, P, V, Q, epsilon, actual, alpha, beta = rows(A, B, h, r, z, x)
                            y = B*h+r
                            for j in range(-1, 2*h+2):
                                ac = sum(alpha(i)*alpha(j-i) for i in range(max(0, j+1)))
                                bc = sum(beta(i)*beta(j-i) for i in range(max(-L, j-h), min(h, j+L)+1))
                                U, W = 2*(mass-z)-A*j, 2*z+A*j
                                at, ah = path_count(U, W, lambda i, q: i+q == mass and (q-z) % A == 0)
                                require(ah == ac <= at, "literal alpha midpoint path injection")
                                U, W = 2*y-B*j, 2*x+delta*j
                                level = delta*y+B*x
                                def beta_selected(i, q):
                                    if delta*i+B*q != level or (y-i) % B:
                                        return False
                                    ell = (y-i)//B
                                    return q == x+delta*ell
                                bt, bh = path_count(U, W, beta_selected)
                                require(bh == bc <= bt, "literal beta linear-level path injection")
                                require(actual.get(j, 0) == at*bt, "two full paths equal actual multinomial")
                                virtual = V.nth(j) if j >= 0 else 0
                                require(ac*bc == virtual, "selected path pair equals virtual coefficient")
                                require(actual.get(j, 0) >= virtual, "coefficientwise actual dominates virtual")
                                endpoints += 1
                            cases += 1
    print(f"midpoint DP: {cases} indexed parameter cases, {endpoints} endpoints; noncoprime and x=0 controls retained")


def eligible(data):
    A, B, h, r, z, x = data
    mass, a = x+B*h+r+z, A*(B*h+r)+B*z
    return gcd(A, B) == 1 and A*x-(B-A)*z > 0 and gcd(a, mass) == 1


def sign_bank():
    head = [(A, B, h, r, z, x) for B in range(2, 7) for A in range(1, B)
            for h in (1, 2, 3) for r in range(B) for z in range(A) for x in (1, 2, 3)]
    wide = [(1, 2, h, r, 0, x) for h in (2, 3, 4, 5, 6, 8, 10, 12)
            for r in (0, 1) for x in (5, 7, 11, 101, 503, 997)]
    rng = Random(4440)
    for _ in range(300):
        B = rng.randrange(3, 13)
        A = rng.randrange(1, B)
        h, r, z = rng.randrange(1, 9), rng.randrange(B), rng.randrange(A)
        x = rng.choice((1, 2, 3, 5, 7, 11, 23, 101, 997))
        wide.append((A, B, h, r, z, x))
    distinct = set()
    for name, bank in (("head", head), ("wide", wide)):
        count = roots = 0
        for data in bank:
            if not eligible(data):
                continue
            mass, L, P, V, Q, epsilon, actual, _, _ = rows(*data)
            D = sp.Poly(Q.as_expr()-t**epsilon*V.as_expr(), t)
            ds, dc = signs_at_roots(P, D)
            vs, vc = signs_at_roots(P, V)
            require(all(s == -1 for s in vs), "virtual Hadamard square has negative root values")
            require(all((-1)**epsilon*s == -1 for s in ds), "FINITE-EXACT actual-minus-virtual root signs")
            for j, value in actual.items():
                require(value >= (V.nth(j) if j >= 0 else 0), "full signed bank coefficientwise inclusion")
            record("root signs", data, epsilon, dc, vc)
            count += 1
            roots += P.degree()
            distinct.add(data)
        print(f"{name} sign bank: {count} rows, {roots} roots; every raw actual-minus-virtual and virtual value negative")
    print(f"distinct parameter tuples across sign banks: {len(distinct)}")


def numeric_difference_remainder(g):
    _, _, _, V, Q, _, _, _, _ = rows(2, 3, 2, 0, 1, g-7)
    b, c = F(20*(g-5), 6), F((g-5)*(g-6), 6)
    powers = {0: (F(1), F(0)), 1: (F(0), F(1)), -1: (-b/c, -1/c)}
    for j in range(2, 6):
        u, v = powers[j-1]
        powers[j] = -c*v, u-b*v
    C = sum(F(Q.nth(j))*powers[j-1][0] for j in range(6))-sum(F(V.nth(j))*powers[j][0] for j in range(5))
    E = sum(F(Q.nth(j))*powers[j-1][1] for j in range(6))-sum(F(V.nth(j))*powers[j][1] for j in range(5))
    return C, E, 2*C-b*E, C*C-b*C*E+c*E*E


def width15_symbolic():
    g = parameter
    field = sp.QQ.frac_field(g)
    def bc(v, k):
        return sp.prod(v-j for j in range(k))*sp.Rational(1, factorial(k)) if k >= 0 else sp.Integer(0)
    alpha = lambda j: bc(g, 1+2*j)
    beta = lambda j: bc(g-1-2*j, 6-3*j)
    P = sp.Poly(6*t*t+20*(g-5)*t+(g-5)*(g-6), t, domain=field)
    V = sp.Poly(sum(sp.expand(sum(alpha(i)*alpha(j-i) for i in range(j+1)) *
                             sum(beta(i)*beta(j-i) for i in range(j-2, 3)))*t**j for j in range(5)), t, domain=field)
    Q = sp.Poly(sum(sp.prod(2*g-i for i in range(15-j))*
                    sp.Rational(1, factorial(15-3*j)*factorial(2*j))*t**j for j in range(6)), t, domain=field)
    require(not V.as_expr().atoms(sp.Float) and not Q.as_expr().atoms(sp.Float), "symbolic producer has no floating coefficients")
    inverse = -(6*t+20*(g-5))/((g-5)*(g-6))
    remainder = sp.Poly(inverse*(Q.as_expr()-t*V.as_expr()), t, domain=field).rem(P)
    C, E = [sp.cancel(remainder.nth(i)) for i in (0, 1)]
    trace = sp.cancel(2*C-sp.Rational(10, 3)*(g-5)*E)
    norm = sp.cancel(C*C-sp.Rational(10, 3)*(g-5)*C*E+(g-5)*(g-6)*E*E/6)
    require(sp.denom(C) == sp.denom(C).subs(g, 0) and sp.degree(C, g) <= 14, "constant remainder polynomial degree bound")
    require(sp.denom(E) == sp.denom(E).subs(g, 0) and sp.degree(E, g) <= 13, "linear remainder polynomial degree bound")
    require(sp.degree(trace, g) <= 14 and sp.degree(norm, g) <= 28, "trace/norm exact interpolation degree bounds")
    K = g*(g-1)*(g-2)*(g-3)*(g-4)
    J = sp.Poly(sp.cancel(-trace*294226732800/(K*(g-5))), g, domain=sp.ZZ)
    H = sp.Poly(sp.cancel(norm*961881892157362176000000/(K*K*(g-5))), g, domain=sp.ZZ)
    require(J.degree() == 8 and H.degree() == 17, "normalized transport trace/norm degrees")
    Jshift = sp.Poly(J.as_expr().subs(g, g+8), g)
    Hshift = sp.Poly(H.as_expr().subs(g, g+8), g)
    require(all(c > 0 for c in Jshift.all_coeffs()), "all nine shifted trace coefficients positive")
    require(all(c > 0 for c in Hshift.all_coeffs()), "all eighteen shifted norm coefficients positive")
    for value in range(8, 37):
        numeric = numeric_difference_remainder(value)
        symbolic = tuple(F(expr.subs(g, value)) for expr in (C, E, trace, norm))
        require(numeric == symbolic, "29 independent Fraction quotient evaluations certify degree<=28 identities")
    print("width15 transport: exact QQ(g) remainder; negative trace degree14 and positive norm degree28")
    print("normalized J8(g+8) coefficients ascending="+str(list(reversed(Jshift.all_coeffs()))))
    print("normalized H17(g+8) coefficients ascending="+str(list(reversed(Hshift.all_coeffs()))))
    print("width15 independent replay: 29 exact Fraction points g=8..36; algebraic nonprimitive points are not first-return claims")
    record("width15", J.all_coeffs(), H.all_coeffs())


def hostiles():
    _, _, P, V, Q, epsilon, _, _, _ = rows(1, 3, 1, 0, 0, 1)
    require(P == polynomial([4, 4]) and V == polynomial([16, 64, 28]), "virtual row exact hostile")
    require(Q == polynomial([28, 280, 28]) and epsilon == 0, "actual row differs from virtual, no carry needed")
    require(V.eval(-1) == -20 and Q.eval(-1) == -224, "negative rows are not proportional or identical")
    require(polynomial([1, 1]).eval(-1) == 0 and polynomial([1]).eval(-1) > 0, "nonnegative coefficients do not force negative value at a first root")
    _, _, P, V, Q, epsilon, _, _, _ = rows(1, 2, 2, 0, 0, 1)
    require(Q-sp.Poly(t**epsilon*V.as_expr(), t) == polynomial([8, 220, 1080, 1080]), "three-channel missed-midpoint defect")
    print("hostiles: actual!=virtual at(-3,1,9); positive coefficient ordering has no negative-root implication")


def main():
    print("HADAMARD VIRTUAL DOUBLING AND ACTUAL MIDPOINT DEFECT")
    print("scope=PROVED virtual sign and coefficientwise injection; universal rootwise defect sign OPEN")
    path_controls()
    sign_bank()
    width15_symbolic()
    hostiles()
    print("trace_sha256="+TRACE.hexdigest())
    print(f"PASS explicit_gates={CHECKS}")


if __name__ == "__main__":
    main()
