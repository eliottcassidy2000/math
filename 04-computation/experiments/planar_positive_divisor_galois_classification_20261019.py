"""Small exact identity witnesses for packet 135; no Galois-group census.

The universal theorem is proved in the companion note. These specified
witnesses check only formulas used there, including one square-class example.
"""
import json
import math
import sympy as sp

v, e, y = sp.symbols('v E Y')
checks = []


def check(name, value):
    if not bool(value):
        raise RuntimeError(name)
    checks.append(name)


for a, r in [(1, 2), (2, 1), (2, 2), (3, 3)]:
    m = a*r
    n = m+a+1
    tag = f'a{a}_r{r}'
    p = -sum(sp.binomial(r, j)*(-1)**(r-j)*e**j*v**(a*(j+1)+1)
             / sp.Integer(a*(j+1)+1) for j in range(r+1))
    kappa = sp.Rational((-1)**(r+1)*a**r*math.factorial(r),
                        math.prod(a*j+1 for j in range(1, r+2)))
    check(tag+'_primitive_derivative', sp.expand(sp.diff(p, v)+v**a*(e*v**a-1)**r) == 0)
    check(tag+'_degree_and_leading', sp.Poly(p, v).degree() == n and sp.Poly(p, v).LC() == -e**r/n)
    check(tag+'_critical_constant', sp.simplify(p.subs({e: 1, v: 1})-kappa) == 0)
    check(tag+'_all_critical_values', sp.rem(p.subs(e, 1)-kappa*v, v**a-1, v) == 0)
    u = e*v**a-1
    chi = 1/(v*u**r)
    theta = v**2*u**(2*r+1)
    check(tag+'_birational_u', sp.cancel(chi**2*theta-u) == 0)
    check(tag+'_birational_energy', sp.cancel(chi**a*u**m*(1+u)-e) == 0)
    check(tag+'_coordinate_degrees', sp.degree(sp.denom(chi), v) == m+1
          and sp.degree(theta, v) == 2*n-a and (m+1) % n != 0 and (2*n-a) % n != 0)
    if (a, r) == (2, 2):
        f = sp.expand(p/(-e**r/n)-y/(-e**r/n))
        h = v**(a//2)*(e*v**a-1)**(r//2)
        res = sp.resultant(f, h, v)
        disc = sp.discriminant(f, v)
        exact_disc = (-1)**(n*(n-1)//2)*(sp.Integer(n)/e**r)**n*res**2
        check('a2_r2_discriminant_resultant', sp.cancel(disc-exact_disc) == 0)
        square_part = sp.Integer(n)**((n-1)//2)/e**(r*n//2)*res
        check('a2_r2_square_class_minus7', sp.cancel(disc-(-7)*square_part**2) == 0)

check('a8_r2_arithmetic_A25_constant', (8*(2+1)+1, (-1)**(8//2)*(8*(2+1)+1)) == (25, 25))
print(json.dumps({'packet': 135, 'status': 'PASS', 'checks': len(checks),
                  'witnesses': [[1, 2], [2, 1], [2, 2], [3, 3]], 'verified': checks,
                  'scope': 'Exact identity witnesses only; no computation of or census over Galois groups.'},
                 ensure_ascii=False, indent=2, sort_keys=True))
