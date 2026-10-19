"""Small exact identities for packet 137; no Galois-group computation."""
import json
import math
import sympy as sp

v, e, y, z, lam, q, zeta = sp.symbols('v E Y z lambda q zeta')
checks = []


def check(name, value):
    if not bool(value):
        raise RuntimeError(name)
    checks.append(name)


witnesses = [(2, 2), (5, 2), (2, 3), (2, 4), (3, 3)]
for b, r in witnesses:
    m, n = b*r-1, b*(r-1)-1
    tag = f'b{b}_r{r}'
    p = -sum(sp.binomial(r-2, j)*(-1)**(r-2-j)*e**j*v**(b*(j+1)-1)
             / sp.Integer(b*(j+1)-1) for j in range(r-1))
    kappa = sp.Rational((-1)**(r-1)*b**(r-2)*math.factorial(r-2),
                        math.prod(b*j-1 for j in range(1, r)))
    check(tag+'_primitive_derivative', sp.expand(sp.diff(p, v)+v**(b-2)*(e*v**b-1)**(r-2)) == 0)
    check(tag+'_degree_and_leading', sp.degree(p, v) == n and sp.Poly(p, v).LC() == -e**(r-2)/n)
    check(tag+'_critical_constant', sp.simplify(p.subs({e: 1, v: 1})-kappa) == 0)
    check(tag+'_critical_scaling_remainder', sp.rem(v*p.subs(e, 1)-kappa, v**b-1, v) == 0)
    u = 1/(e*v**b-1)
    chi = v/(e*v**b-1)**r
    theta = (e*v**b-1)**(2*r-1)/v**2
    check(tag+'_birational_u', sp.cancel(chi**2*theta-u) == 0)
    check(tag+'_birational_energy', sp.cancel(chi**(-b)*u**m*(1+u)-e) == 0)
    check(tag+'_parameter_recovery', sp.cancel(chi/u**r-v) == 0)
    if r == 2:
        check(tag+'_kummer_equation', sp.expand(-n*(p-p.subs(v, z)-lam)-(v**n-z**n+n*lam)) == 0)
    if (b, r) == (2, 4):
        f = sp.expand((p-y)/(-e**(r-2)/n))
        h = v**((b-2)//2)*(e*v**b-1)**((r-2)//2)
        res = sp.resultant(f, h, v)
        disc = sp.discriminant(f, v)
        square_part = sp.Integer(n)**((n-1)//2)/e**((r-2)*n//2)*res
        check('b2_r4_discriminant_class5', sp.cancel(disc-5*square_part**2) == 0)

check('kummer_X_trivial_stabilizer_coefficient', sp.Poly(zeta*(q-1)**2-(zeta*q-1)**2, q).coeff_monomial(1) == zeta-1)
check('kummer_T_trivial_stabilizer_coefficient', sp.Poly((zeta*q-1)**3-zeta**2*(q-1)**3, q).LC() == zeta**3-zeta**2)
check('b2_r6_arithmetic_A9_constant', (2*(6-1)-1, (-1)**(2//2-1)*(2*(6-1)-1)) == (9, 9))
print(json.dumps({'packet': 137, 'status': 'PASS', 'checks': len(checks),
                  'witnesses': [list(pair) for pair in witnesses], 'verified': checks,
                  'scope': 'Small exact primitive, coordinate, Kummer, and discriminant identities only.'},
                 ensure_ascii=False, indent=2, sort_keys=True))
