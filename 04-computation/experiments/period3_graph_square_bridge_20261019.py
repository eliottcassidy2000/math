#!/usr/bin/env python3
"""Exact controls for period-three graph actions and square-sum AP embeddings.

Requires SymPy. No floating-point arithmetic and no truth-bearing asserts.
The square-sum minimum is a complete finite search with a proved root bound.
"""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, permutations
from math import isqrt, prod

import sympy as S


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def orbit_type(p):
    seen, lengths = set(), []
    for start in range(len(p)):
        if start in seen:
            continue
        v, length = start, 0
        while v not in seen:
            seen.add(v)
            length += 1
            v = p[v]
        lengths.append(length)
    return tuple(sorted(lengths))


def exact_order_three(p):
    return p != tuple(range(len(p))) and all(p[p[p[v]]] == v for v in range(len(p)))


def graph_automorphisms(n, edges):
    """Backtrack directly from adjacency, independent of Kneser coordinates."""
    adj = [{j for j in range(n) if frozenset((i, j)) in edges} for i in range(n)]
    mapping, used = [], set()

    def visit():
        i = len(mapping)
        if i == n:
            yield tuple(mapping)
            return
        for image in range(n):
            if image in used or len(adj[i]) != len(adj[image]):
                continue
            if any((j in adj[i]) != (mapping[j] in adj[image]) for j in range(i)):
                continue
            mapping.append(image)
            used.add(image)
            yield from visit()
            used.remove(image)
            mapping.pop()

    return list(visit())


def dynatomic_checks():
    x, c = S.symbols("x c")
    f = x*x + c
    f2 = S.expand(f*f+c)
    f3 = S.expand(f2*f2+c)
    phi, rem = S.div(f3-x, f-x, x)
    require(rem == 0 and S.degree(phi, x) == 6, "dynatomic degree")
    P = x**3 + x**2/2 - S.Rational(9, 4)*x - S.Rational(1, 8)
    Q = x**3 + x**2/4 - S.Rational(41, 16)*x + S.Rational(23, 64)
    R = (x+S.Rational(7, 4))*(x-S.Rational(5, 4))*(x+S.Rational(1, 4))
    require(S.expand(phi.subs(c, -S.Rational(7, 4))-P**2) == 0, "parabolic factorization")
    require(S.expand(phi.subs(c, -S.Rational(29, 16))-Q*R) == 0, "rational parameter factorization")
    for label, parameter, factors, count in [
        ("parabolic", -S.Rational(7, 4), [P], 5),
        ("rational", -S.Rational(29, 16), [Q, R], 8),
    ]:
        squarefree = S.Poly((f3-x).subs(c, parameter), x).sqf_part()
        require(squarefree.degree() == count, label+" distinct root count")
        require(squarefree.count_roots(-S.oo, S.oo) == count, label+" real root count")
        require(S.Poly((f-x).subs(c, parameter), x).is_irreducible,
                label+" fixed points are not rational")
        for factor in factors:
            require(S.gcd(factor, (f-x).subs(c, parameter)) == 1, "cycle/fixed separation")
            require(S.rem(factor.subs(x, x*x+parameter), factor, x) == 0, "cycle preservation")
        print(f"{label}: c={parameter}; distinct real roots of f^3-id={count}; fixed=2; exact-3={count-2}")
    require(S.Poly(P, x).is_irreducible, "parabolic cubic irreducible")
    require(S.Poly(Q, x).is_irreducible, "other rational-parameter cubic irreducible")
    require(S.discriminant(P, x) == 49, "parabolic cubic discriminant")
    orbit = [F(-7, 4), F(5, 4), F(-1, 4)]
    require(all(orbit[i]**2-F(29, 16) == orbit[(i+1)%3] for i in range(3)), "rational orbit")
    require(prod(orbit[i]+orbit[(i+1)%3] for i in range(3)) == 1, "secant product")
    require(prod(2*v for v in orbit) == F(35, 8), "multiplier differs from secant product")
    print("rational audit: parabolic 3-cycle irreducible; other -29/16 cycle irreducible; secant product=1, multiplier=35/8")


def graph_checks():
    vertices = list(combinations(range(5), 2))
    index = {v: i for i, v in enumerate(vertices)}
    E5 = {frozenset(e) for e in combinations(range(5), 2)}
    E33 = {frozenset((a, b)) for a in range(3) for b in range(3, 6)}
    EP = {frozenset((i, j)) for i, j in combinations(range(10), 2)
          if set(vertices[i]).isdisjoint(vertices[j])}
    expected = {"K5": (120, Counter({(1, 1, 3): 20})),
                "K3,3": (72, Counter({(1, 1, 1, 3): 4, (3, 3): 4})),
                "Petersen": (120, Counter({(1, 3, 3, 3): 20}))}
    for name, n, edges in [("K5", 5, E5), ("K3,3", 6, E33), ("Petersen", 10, EP)]:
        autos = graph_automorphisms(n, edges)
        census = Counter(orbit_type(p) for p in autos if exact_order_three(p))
        require((len(autos), census) == expected[name], name+" automorphism census")
        print(f"{name}: automorphisms={len(autos)}; order-three orbit census={dict(sorted(census.items()))}")
        if name == "Petersen":
            induced = {tuple(index[tuple(sorted((p[a], p[b])))] for a, b in vertices)
                       for p in permutations(range(5))}
            require(set(autos) == induced, "independent Kneser/full-adjacency agreement")
            independent_four = [v for v in combinations(range(10), 4)
                                if all(frozenset(e) not in EP for e in combinations(v, 2))]
            require(len(independent_four) == 5, "five maximum independent stars")
            require(all(set.intersection(*(set(vertices[i]) for i in v)) for v in independent_four),
                    "every independent four is a star")
    print("Petersen: direct adjacency automorphisms equal the 120 induced S5 actions; five independent 4-stars")


def triangle_minimum():
    N = 6242
    hits = []
    # For a<b<c<=N, all pair sums are <=2N-1. The root bound is complete.
    for p in range(1, isqrt(2*N-1)+1):
        for q in range(p+1, isqrt(2*N-1)+1):
            r2 = 2*q*q-p*p
            r = isqrt(r2)
            if r*r != r2:
                continue
            numerators = (2*p*p-q*q, q*q, 3*q*q-2*p*p)
            if any(v % 2 for v in numerators):
                continue
            a, b, c = (v//2 for v in numerators)
            if 0 < a < b < c <= N:
                hits.append(((a, b, c), (p, q, r)))
    require(hits == [((482, 3362, 6242), (62, 82, 98))], "minimal AP square-sum triangle")
    # Independent enumeration by the middle label: 2b must be an even square.
    independent = []
    for q in range(2, isqrt(2*N-1)+1, 2):
        b = q*q//2
        if b >= N:
            continue
        for p in range(1, q):
            a = p*p-b
            c = 2*b-a
            if 0 < a < b < c <= N and isqrt(b+c)**2 == b+c:
                independent.append((a, b, c))
    require(independent == [(482, 3362, 6242)], "independent minimum search")
    orbit = [F(-7, 4), F(5, 4), F(-1, 4)]
    image = [1920*v+3842 for v in orbit]
    require(image == [482, 6242, 3362], "affine orbit image")
    require(all(isqrt(int(a+b))**2 == a+b for a, b in combinations(image, 2)), "all square sums")
    print("AP square-sum minimum: unique triple with max<=6242 is (482,3362,6242); edge roots=(62,82,98)")
    print("affine orbit map x -> 1920*x+3842: 482 -> 6242 -> 3362 -> 482")


def symbolic_bridge_checks():
    x, t, k, j, trace = S.symbols("x t k j trace")
    par = -trace*trace-trace-2
    cubic = x**3-trace*x*x+(par-trace-1)*x-par*trace-par-1
    require(S.factor(S.discriminant(cubic, x)) == (4*trace*trace+6*trace+9)**2,
            "trace cubic discriminant")
    require(S.factor(S.resultant(cubic, x*x+par-x, x)) == 4*trace*trace+6*trace+9,
            "trace cubic no real fixed root")
    require(S.rem(cubic.subs(x, x*x+par), cubic, x) == 0, "trace cubic preservation")
    f1 = x*x+par
    f2 = S.expand(f1*f1+par)
    f3 = S.expand(f2*f2+par)
    require(S.expand((f1-x)*cubic*cubic.subs(trace, -1-trace)-(f3-x)) == 0,
            "trace factorization all parameters")
    print("general real classification: trace-cubic discriminant=(4*s^2+6*s+9)^2; f^3-id factorization verified")
    q = S.interpolate([(0, 1), (1, t), (t, 0)], x)
    A, B, C = [S.factor(S.expand(q).coeff(x, power)) for power in (2, 1, 0)]
    parameter = S.factor(A*C+B/2-B*B/4)
    defect = -(t**3-2*t**2-t+1)**2/(4*t*t*(t-1)**2)
    require(S.factor(parameter+S.Rational(7, 4)-defect) == 0, "interpolation discriminant bridge")
    require(S.factor(A+(t*t-t+1)/(t*(t-1))) == 0, "quadratic coefficient")
    require(S.Poly(t**3-2*t**2-t+1, t).is_irreducible, "rational parabolic-shape obstruction")
    p, r, s = 2*k*k-1, 2*k*k+2*k+1, 2*k*k+4*k+1
    require(S.expand(p*p+s*s-2*r*r) == 0, "infinite AP square family")
    a, b, c = 2*(2*p*p-r*r), 2*r*r, 2*(3*r*r-2*p*p)
    require(S.expand(a+c-2*b) == 0, "AP vertex family")
    for value, root in [(a+b, 2*p), (a+c, 2*r), (b+c, 2*s)]:
        require(S.expand(value-root*root) == 0, "family square edge")
    shifted = S.Poly(S.expand((a/2).subs(k, j+4)), j)
    require(all(v > 0 for v in shifted.all_coeffs()), "family positivity for k>=4")
    require([v.subs(k, 4) for v in (a, b, c)] == [482, 3362, 6242], "first family triple")
    print("symbolic interpolation: c(t)+7/4=-(t^3-2t^2-t+1)^2/[4*t^2*(t-1)^2]")
    print("infinite square-sum AP family: verified identically, positive for k>=4, k=4 gives the minimum")
    # A square-sum path is not automatically graceful, even after its sign gauge.
    path = [8, 1, 15, 10, 6, 3, 13, 12, 4, 5, 11, 14, 2, 7, 9]
    require(sorted(path) == list(range(1, 16)), "Q15 spanning path")
    edge_sums = [a+b for a, b in zip(path, path[1:])]
    require(set(edge_sums) == {9, 16, 25}, "Q15 squares")
    gauged = [(-1)**i*v for i, v in enumerate(path)]
    require([abs(a-b) for a, b in zip(gauged, gauged[1:])] == edge_sums, "sum/difference gauge")
    require(len(set(edge_sums)) < len(edge_sums), "graceful collision hostile")
    print("Q15 hostile: sign gauge preserves absolute edge magnitudes {9,16,25}, with repetitions; no graceful labeling follows")
    print("PASS: all exact symbolic identities and complete declared finite universes")


if __name__ == "__main__":
    dynatomic_checks()
    graph_checks()
    triangle_minimum()
    symbolic_bridge_checks()
