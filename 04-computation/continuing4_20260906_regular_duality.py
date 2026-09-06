"""Exact symbolic certificate for reflected regular trinomial rows.

Universe: H=1,...,6, all polynomial parameter z; arithmetic in QQ[z].
Run normally or with -O. No assertion is used as a validity gate.
"""
from pathlib import Path
from math import factorial, gcd, comb
import json
import sys
import sympy as S

sys.stdout.reconfigure(newline='\n')
z = S.Symbol('z')
Z = S.Poly(0, z, domain=S.QQ)
O = S.Poly(1, z, domain=S.QQ)
gates = 0

def require(ok, label):
    global gates
    gates += 1
    if not ok:
        raise ArithmeticError(label)

def rising(a, n):
    p = O
    for i in range(n):
        p *= S.Poly(a+i, z, domain=S.QQ)
    return p

def rows(h):
    p = [rising(z+2*j+1, 2*(h-j)).mul_ground(
        S.Rational(factorial(h), factorial(j)*factorial(3*h-3*j)))
         for j in range(h+1)]
    q = [rising(2*z+2*j+1, 4*h-2*j).mul_ground(
        S.Rational(factorial(2*h), factorial(j)*factorial(6*h-3*j)))
         for j in range(2*h+1)]
    return p, q

def rem(a, p):
    a = list(a)
    h = len(p)-1
    while len(a) >= len(p):
        c = a.pop()
        for j in range(h):
            a[len(a)-h+j] -= c*p[j]
    return a+[Z]*(h-len(a))

def mm(a, b):
    h = len(a)
    return [[sum((a[i][j]*b[j][k] for j in range(h)), Z)
             for k in range(h)] for i in range(h)]

def divisor(h, k):
    p = O
    for ell in range(1, h+1):
        p *= S.Poly((z+2*ell-1)*(z+2*ell), z,
                     domain=S.QQ)**max(0, k-h+ell)
    return p

def literal_counts(h, zz, mass):
    a, b, c = -3*h, zz, 6*h+3*zz
    return {(i, j, mass-i-j) for i in range(mass+1)
            for j in range(mass-i+1)
            if a*i+b*j+c*(mass-i-j) == 0}

certificate = {'status': 'FINITE-EXACT symbolic identity certificate',
               'universe': 'H=1..6; parameter z indeterminate over QQ',
               'coefficients': 'ascending powers of z, exact rational strings',
               'rows': []}
count = 0
for h in range(1, 7):
    p, q = rows(h)
    require(p[-1] == O and q[-1] == O, 'monic rows')
    r = rem(q, p)
    columns = [rem([Z]*j+r, p) for j in range(h)]
    m = [[columns[j][i] for j in range(h)] for i in range(h)]
    b = [[O if i == j else Z for j in range(h)] for i in range(h)]
    cs = []
    for k in range(1, h+1):
        b = mm(m, b)
        c = sum((b[i][i] for i in range(h)), Z).mul_ground(S.Rational(-1, k))
        cs.append(c)
        for i in range(h):
            b[i][i] += c
    for row in b:
        for entry in row:
            require(entry.is_zero, 'complete Cayley-Hamilton remainder')
    residuals = []
    for k, c in enumerate(cs, 1):
        require(c.degree() <= 4*h*k, 'weighted degree')
        d = divisor(h, k)
        require(d.degree() == k*(k+1), 'paired divisor degree')
        quotient, remainder = S.div(c, d)
        require(remainder.is_zero, 'all-height forced divisor in exact bank')
        require(quotient.degree() <= 4*h*k-k*(k+1), 'residual degree')
        coefficients = list(reversed(quotient.all_coeffs()))
        for value in coefficients:
            require(value > 0, 'strict residual coefficient positivity')
        residuals.append({'k': k, 'degree': quotient.degree(),
                          'coefficients': [str(v) for v in coefficients]})
        count += len(coefficients)
    certificate['rows'].append({'H': h, 'residuals': residuals})
    print('H', h, 'residual_degrees', [v['degree'] for v in residuals],
          'positive_coefficients', sum(len(v['coefficients']) for v in residuals))
    # Literal complete-fibre and normalization controls, including both returns.
    for zz in (1, 2, 5):
        g = 3*h+zz
        if gcd(zz, 3*h) != 1:
            continue
        for multiplier in (1, 2):
            mass = multiplier*g
            triples = [(multiplier*zz+2*j, multiplier*3*h-3*j, j)
                       for j in range(multiplier*h+1)]
            require(literal_counts(h, zz, mass) == set(triples), 'complete fibre')
            row = p if multiplier == 1 else q
            for j, counts in enumerate(triples):
                literal = S.Rational(factorial(mass),
                                     S.prod(factorial(v) for v in counts))
                require(literal == comb(mass, multiplier*h)*row[j].eval(zz),
                        'literal multinomial normalization')
        require(all(not literal_counts(h, zz, mass) for mass in range(1, g)),
                'first-return positive control')
    require(literal_counts(h, 3*h, 1) == set(), 'gcd hostile mass one')
    require((1, 1, 0) in literal_counts(h, 3*h, 2), 'gcd hostile first mass two')
    # Boundary reflection identity, including the H=0 complement.
    for rr in range(1, h+1):
        hh, zz = h-rr, 2*rr+1
        pp, qq = rows(hh)
        for j in range(h+1):
            old = S.Rational(factorial(2*h+1),
                             factorial(3*h-3*j)*factorial(1+2*j))
            old *= S.ff(h-rr, h-j)
            new = 0 if j < rr else pp[j-rr].eval(zz)
            require(old == new, 'first boundary reflection')
        for e in range(-1, 2*h+1):
            old = S.ff(2*h-2*rr, 2*h-e)/(
                factorial(6*h-3*e)*factorial(2+2*e))
            new = 0 if e < 2*rr else qq[e-2*rr].eval(zz)/factorial(4*h+2)
            require(old == new, 'doubled boundary reflection')
certificate['positive_coefficient_count'] = count
require(count == 833, 'certificate coefficient inventory')
destination = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).with_name(
    'continuing4_20260906_regular_duality_certificate.json')
destination.write_text(json.dumps(certificate, indent=2)+'\n', encoding='utf-8', newline='\n')
print('positive_coefficients', count)
print('always_active_gates', gates)
print('PASS')
