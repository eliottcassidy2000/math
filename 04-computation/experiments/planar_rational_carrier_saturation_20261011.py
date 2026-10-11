#!/usr/bin/env python3
"""Exact controls for rational-carrier saturation; no planar JC claim."""
import itertools
import json
import sympy as s

p, y, x, t, lam = s.symbols('p y x t lam')
D = p**3-y**2
u = y**2/D
g0 = -u/2
checks = 0


def check(ok, label):
    global checks
    checks += 1
    if not ok:
        raise RuntimeError(label)


def bracket(f, g):
    return s.cancel(-D/p*(s.diff(f, p)*s.diff(g, y)-s.diff(f, y)*s.diff(g, p)))


def poly(f):
    return s.cancel(f).is_polynomial(p, y) is True


def main():
    denominators = [p, y, D, p-y, p*p+y+1]
    pole_cases = 0
    for h in denominators:
        for order in range(1, 5):
            for numerator in [s.Integer(1), 1+2*p+3*y]:
                S = numerator/h**order
                check(not poly(S), 'genuine hostile denominator')
                check(not (poly(bracket(p, S)) and poly(bracket(y, S))),
                      'at least one coordinate pole survives')
                pole_cases += 1

    # This bounded universe is independent of the source-test classifier.
    bases = [p, y, p*p, p*D, p*p*D]
    polynomial_cases = 0
    for coefficients in itertools.product([0, 1], repeat=len(bases)):
        S = 7+sum(c*b for c, b in zip(coefficients, bases))
        coordinate_ok = poly(bracket(p, S)) and poly(bracket(y, S))
        all_tests = all(poly(bracket(g0+H, S)) for H in [0, p, y])
        check(coordinate_ok == poly((S-7)/p**2), 'coordinate criterion iff')
        check(all_tests == poly((S-7)/(p**2*D)), 'three-test criterion iff')
        polynomial_cases += 1

    positive_cases = 0
    for a in range(4):
        for b in range(4):
            S = 11+p*p*D*(1+p**a*y**b)
            for H in [0, p, y, p*p*y+y**3+1]:
                check(poly(bracket(g0+H, S)), 'positive carrier response')
            positive_cases += 1

    check(poly(bracket(g0, D)) and poly(bracket(g0+y, D))
          and not poly(bracket(g0+p, D)), 'missing p response hostile')
    check(poly(bracket(g0, p*D)) and poly(bracket(g0+p, p*D))
          and not poly(bracket(g0+y, p*D)), 'missing y response hostile')
    check(poly(bracket(p, p*p)) and poly(bracket(y, p*p))
          and not poly(bracket(g0, p*p)), 'missing baseline-source hostile')
    check(s.cancel(bracket(g0, p*p)+2*p**3*y/D) == 0,
          'baseline pole coefficient')

    # With repeated Poisson numerator p^2, S=1/p is a counterexample.
    repeated_y_response = -p**2*s.diff(1/p, p)
    check(s.cancel(repeated_y_response-1) == 0, 'squarefree hypothesis hostile')
    check(bracket(g0, u) == 0, 'fixed-source rational survivor')
    check(not poly(bracket(p, u)) and not poly(bracket(y, u)),
          'fixed-source survivor leaves carrier')
    check(s.cancel(bracket(p, u)+2*y*p**2/D) == 0, 'u p-response')
    check(s.cancel(bracket(y, u)+3*p*y**2/D) == 0, 'u y-response')

    # Independent literal source derivatives, using the birational chart.
    px = t*(1+x*x*t)
    yx = x*t*px
    def literal(f, g):
        return s.cancel(s.diff(f, x)*s.diff(g, t)-s.diff(f, t)*s.diff(g, x))
    chart = {p:px, y:yx}
    for S in [p*p*D, p*D, D, 1/p, 1/D, u]:
        sx = s.cancel(S.subs(chart, simultaneous=True))
        for F in [p, y, g0]:
            fx = s.cancel(F.subs(chart, simultaneous=True))
            expected = s.cancel(bracket(F, S).subs(chart, simultaneous=True))
            check(s.cancel(literal(fx, sx)-expected) == 0, 'literal source bracket')

    X = x/(1-lam*x)
    T = t*(1-lam*x)**2
    check(s.cancel(X*X*T-x*x*t) == 0, 'rational flow fixes u')
    check(s.cancel(literal(X, T)-1) == 0, 'rational flow symplectic')
    check(s.cancel(X.subs(x, x/(1+lam*x))-x) == 0, 'rational x inverse')
    check(s.cancel(T.subs({x:x/(1+lam*x), t:t*(1+lam*x)**2},
                          simultaneous=True)-t) == 0, 'rational t inverse')
    print(json.dumps({
        'status':'FINITE-EXACT interfaces; universal rational-time obstruction inherited; JC OPEN',
        'checks':checks,
        'pole_cases':pole_cases,
        'polynomial_cases':polynomial_cases,
        'positive_carriers':positive_cases,
        'independent_source_chart_cases':18,
        'squarefree_numerator_hostile':'{p,y}=p^2, S=1/p, {y,S}=1'
    }, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
