#!/usr/bin/env python3
"""Independent symbolic controls of the charge and source-domain boundaries."""
import json
from itertools import product
import sympy as s

x, p = s.symbols('x p')


def require(ok, msg):
    if not ok:
        raise RuntimeError(msg)


def bracket(f, g):
    return s.expand(s.diff(f, x) * s.diff(g, p) - s.diff(f, p) * s.diff(g, x))


def sector(f, q, m):
    return sum(c * x**a * p**b for (a, b), c in s.Poly(f, x, p).terms() if (a - b) % m == q % m)


checks = 0
for m in range(2, 10):
    u, v, w = x**m, p**m, x*p
    for actual, wanted in ((bracket(u, v), m*m*w**(m-1)),
                           (bracket(w, u), -m*u), (bracket(w, v), m*v)):
        require(s.expand(actual-wanted) == 0, 'invariant bracket')
        checks += 1
    F, G = x**m, p/(m*x**(m-1))
    require(s.simplify(bracket(F, G)-1) == 0, 'punctured cotangent Jacobian')
    checks += 1
    require(sector(bracket(x, p), 0, m) == 1 and bracket(sector(x, 0, m), sector(p, 0, m)) == 0,
            'opposite-charge contribution')
    checks += 1
    require(bracket(x*(1+p), p) == 1+p and sector(1+p, 0, m) == 1, 'averaged Keller hostile')
    checks += 1
    # Verify bilinear charge contraction against direct symbolic derivatives.
    for a, b, c, d in product(range(4), repeat=4):
        f, g = x**a*p**b, x**c*p**d
        actual = sector(bracket(f, g), 0, m)
        wanted = sum(bracket(sector(f, q, m), sector(g, -q, m)) for q in range(m))
        require(s.expand(actual-wanted) == 0, 'charge-pair contraction')
        checks += 1
require(bracket(x*p, x) == -x and bracket(x*p, p) == p, 'left-bracket sign')
require(bracket(x, x*p) == x and bracket(p, x*p) == -p, 'right-bracket sign')
require(bracket(x, p) == 1, 'zero total weight is not cancellation')
checks += 3
print(json.dumps(dict(status='FINITE-EXACT PASS', checks=checks, sympy=s.__version__,
                      cyclic_orders=list(range(2, 10)), monomial_exponents=list(range(4)),
                      averaged_bracket_hostile='P=x(1+p), Q=p: det=1+p, every C_m average=1'),
                 indent=2, sort_keys=True))
