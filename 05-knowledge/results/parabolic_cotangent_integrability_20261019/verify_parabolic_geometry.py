#!/usr/bin/env python3
"""Exact, always-active controls; not a finite test for analytic divergence."""
from pathlib import Path
import json
import sympy as s

q, p, r, z = s.symbols('q p r z')
checks = 0


def check(value, label):
    global checks
    if value is not True:
        raise RuntimeError(label)
    checks += 1


def eq(a, b, label):
    check(s.cancel(a-b) == 0, label)


def trunc(poly, n):
    poly = s.Poly(s.expand(poly), q)
    return s.Add(*(v*q**i[0] for i, v in poly.terms() if i[0] <= n))


def compose(poly, f, n):
    result = 0
    for (k,), a in s.Poly(poly, q).terms():
        result += a*trunc(f**k, n)
    return trunc(result, n)


def exp_vector(v, n):
    value, term = q, q
    for j in range(1, n):
        term = trunc(v*s.diff(term, q), n)
        if term == 0:
            break
        value += term/s.factorial(j)
    return trunc(value, n)


f = q+q*q
N = 14
v, term = 0, q
for j in range(1, N):
    term = trunc(compose(term, f, N)-term, N)
    v += s.Rational((-1)**(j+1), j)*term
v = trunc(v, N)

# Independent triangular inversion of the vector-field exponential.
v_recursive = 0
for n in range(2, N+1):
    candidate = exp_vector(v_recursive, n)
    cn = s.expand(f-candidate).coeff(q, n)
    v_recursive += cn*q**n
    eq(exp_vector(v_recursive, n), trunc(f, n), f'exp recursion order {n}')
eq(v, v_recursive, 'two independent iterative-log constructions')
eq(trunc(compose(v, f, N)-s.diff(f, q)*v, N), 0, 'Julia equation')

expected = [1, -1, s.Rational(3, 2), s.Rational(-8, 3), s.Rational(31, 6),
            s.Rational(-157, 15), s.Rational(649, 30), s.Rational(-9427, 210)]
for k, a in enumerate(expected, start=2):
    eq(s.expand(v).coeff(q, k), a, f'cusp log coefficient {k}')

# Independent p-coordinate logarithm, computed directly from pullback C_f.
inv_fp = sum((-2*q)**k for k in range(N))
term, log_p_coefficient = s.Integer(1), s.Integer(0)
for j in range(1, N):
    term = trunc(compose(term, f, N-1)*inv_fp-term, N-1)
    log_p_coefficient += s.Rational((-1)**(j+1), j)*term
eq(trunc(log_p_coefficient, N-1), -s.diff(v, q), 'full cotangent logarithm')

time_form = s.series(1/v, q, 0, 7).removeO()
eq(time_form.coeff(q, -1), 1, 'cusp iterative residue')
for exponent, coefficient in [(0, s.Rational(-1, 2)), (1, s.Rational(2, 3)),
                              (2, s.Rational(-13, 12))]:
    eq(time_form.coeff(q, exponent), coefficient, f'time form term {exponent}')

# Hamiltonian first-integral powers and canonical form.
Q, P = f, p/s.diff(f, q)
eq(s.det(s.Matrix([Q, P]).jacobian([q, p])), 1, 'cotangent determinant')
eq(P*s.diff(Q, q), p, 'canonical one-form')
for k in range(1, 5):
    # Julia equation to finite order gives the powers; evaluate independently.
    defect = s.series((compose(v, f, N)/s.diff(f, q))**k-v**k, q, 0, N+1).removeO()
    eq(defect, 0, f'formal invariant power {k}')
eq(Q.subs(q, -1-q), Q, 'quadratic deck base')
eq(P.subs({q: -1-q, p: -p}, simultaneous=True), P, 'quadratic deck momentum')

# Parabolic Möbius positive control and arbitrary-order polynomial hostiles.
M, MP = q/(1-q), p*(1-q)**2
eq(M**2*MP, p*q*q, 'source-flow invariant')
jet_receipts = []
for n in range(2, 25):
    fn = sum(q**k for k in range(1, n+1))
    fp = s.diff(fn, q)
    eq(fn-M, -q**(n+1)/(1-q), f'hostile base difference {n}')
    eq(1/fp-(1-q)**2, q**n*((n+1)-n*q)/fp, f'hostile momentum difference {n}')
    eq(s.det(s.Matrix([fn, p/fp]).jacobian([q, p])), 1, f'hostile Jacobian {n}')
    eq(fp.subs(q, 0), 1, f'hostile denominator unit {n}')
    jet_receipts.append({'N': n, 'first_difference_q': f'-q^{n+1}',
                         'first_difference_p': f'{n+1}*p*q^{n}'})

# Exact rational eigenobservables must not be mistaken for invariants.
for d in range(2, 9):
    fd = q**d
    eq(fd*(p/s.diff(fd, q)), p*q/d, f'monomial eigenobservable {d}')
    td = s.chebyshevt(d, q)
    eq((td**2-1)*(p/s.diff(td, q))**2,
       (q*q-1)*p*p/d**2, f'Chebyshev eigenobservable {d}')

# Cubic quotient: direct derivatives versus truncated iterative composition.
c = -s.Rational(7, 4)
base = q*q+c
F = base
for _ in range(2):
    F = s.expand(F.subs(q, base))
Pr = r**3+r*r/2-s.Rational(9, 4)*r-s.Rational(1, 8)
Pq = Pr.subs(r, q)
eq(F-q, (base-q)*Pq**2, 'period-three factorization')
eq(s.discriminant(Pr, r), 49, 'cycle discriminant')
check(s.Poly(Pr, r, domain=s.QQ).is_irreducible, 'cycle cubic irreducible')


def red(a):
    return s.rem(s.expand(a), Pr, r)


a = red(s.diff(F, q, 2).subs(q, r)/2)
b = red(s.diff(F, q, 3).subs(q, r)/6)
eq(a, 14*r*r-28*r-s.Rational(21, 2), 'return quadratic coefficient')
eq(b, -56*r*r+56*r, 'return cubic coefficient')
eq(red(a*a), 1372*r*(r-1), 'square of return coefficient')
eq(red(49*b+2*a*a), 0, 'ratio minus two forty-ninths')
eq(red(1-b*s.invert(a*a, Pr, r)), s.Rational(51, 49), 'iterative residue 51/49')
eq(red(F.subs(q, r)-r), 0, 'return fixed point')
eq(red(s.diff(F, q).subs(q, r)), 1, 'return multiplier')
check(Pr.subs(r, 0) != 0 and Pr.subs(r, 1) != 0, 'nonzero quadratic coefficient')

u = s.symbols('u')
local = r+u
for _ in range(3):
    local = s.Poly(s.expand(local*local+c), u)
    local = s.Add(*(red(coeff)*u**degree[0]
                    for degree, coeff in local.terms() if degree[0] <= 3))
for k, expected_coefficient in [(0, r), (1, 1), (2, a), (3, b)]:
    eq(red(s.expand(local).coeff(u, k)-expected_coefficient), 0,
       f'independent composed return coefficient {k}')

# General simple-parabolic residue, without importing a normal-form sign.
aa, bb = s.symbols('a b', nonzero=True)
vv = aa*q*q+(bb-aa*aa)*q**3
eq(exp_vector(vv, 3), q+aa*q*q+bb*q**3, 'general cubic log')
eq(s.series(1/vv, q, 0, 0).removeO().coeff(q, -1),
   1-bb/aa**2, 'general iterative residue formula')

# Coalescing poles and the exact local discrete linearization period.
delta = s.symbols('delta', positive=True)
omega = (1+q)/(q*q+delta*delta)
res_plus = s.cancel((1+q)/(q+s.I*delta)).subs(q, s.I*delta)
res_minus = s.cancel((1+q)/(q-s.I*delta)).subs(q, -s.I*delta)
eq(2*s.pi*s.I*res_plus, s.pi/delta+s.pi*s.I, 'upper model period')
eq(2*s.pi*s.I*res_minus, -s.pi/delta+s.pi*s.I, 'lower model period')
eq(res_plus+res_minus, 1, 'coalescing total residue')
discrete_period = s.series(2*s.pi*s.I/s.log(1+2*s.I*delta), delta, 0, 5).removeO()
expected_period = (s.pi/delta+s.I*s.pi+s.pi*delta/3-s.I*s.pi*delta**2/3
                   -19*s.pi*delta**3/45+3*s.I*s.pi*delta**4/5)
eq(discrete_period, expected_period, 'exact local discrete period series')

# Catalan inverse: an exact algebraic branch, not a rational inverse.
inverse = sum((-1)**(n-1)*s.catalan(n-1)*q**n for n in range(1, 21))
eq(trunc(inverse+inverse*inverse, 20), q, 'Catalan right inverse')
eq(compose(inverse, f, 20), q, 'Catalan left inverse')
for n in range(1, 21):
    check(inverse.coeff(q, n) != 0, f'nonzero inverse coefficient {n}')

receipt = {
    'status': 'FINITE-EXACT; general and analytic claims require the proofs/citations',
    'checks': checks,
    'iterative_log_order': N,
    'iterative_log_coefficients': {str(k): str(v.coeff(q, k)) for k in range(2, N+1)},
    'time_form_through_order_6': str(time_form),
    'period_three': {'P': str(Pr), 'a': str(a), 'b': str(b), 'residue': '51/49'},
    'jet_hostiles': jet_receipts,
    'model_periods': ['pi/delta+pi*i', '-pi/delta+pi*i'],
    'exact_local_discrete_period_through_order_4': str(discrete_period),
    'catalan_inverse_order': 20,
    'no_computational_claim': ['divergence', 'all-degree rational invariant exclusion',
                               'analytic first-integral nonexistence'],
}
output = Path(__file__).with_name('exact-results.json')
output.write_text(json.dumps(receipt, indent=2)+'\n', encoding='utf-8')
print(f'{checks} exact checks passed; iterative logarithm through degree {N}; period-three residue 51/49.')
