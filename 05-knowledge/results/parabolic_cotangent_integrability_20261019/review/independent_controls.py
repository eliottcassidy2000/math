#!/usr/bin/env python3
"""Independent rational-arithmetic controls, with no SymPy/import of author code.

These finite identities do not prove analytic divergence or universal invariant
exclusion. Those claims require the accompanying proof and cited theorem.
"""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import json

checks = 0

def require(test, label):
    global checks
    if not test:
        raise RuntimeError(label)
    checks += 1

def trim(p):
    p = list(map(F, p))
    while len(p) > 1 and not p[-1]:
        p.pop()
    return p or [F(0)]

def add(a, b):
    return trim([(a[i] if i < len(a) else 0)+(b[i] if i < len(b) else 0)
                 for i in range(max(len(a), len(b)))])

def scale(a, k):
    return trim([k*x for x in a])

def sub(a, b):
    return add(a, scale(b, -1))

def mul(a, b):
    c = [F(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i+j] += x*y
    return trim(c)

def deriv(a):
    return trim([i*a[i] for i in range(1, len(a))])

def divrem(a, b):
    a, b = trim(a), trim(b)
    quotient = [F(0)]*max(1, len(a)-len(b)+1)
    while a != [0] and len(a) >= len(b):
        k, c = len(a)-len(b), a[-1]/b[-1]
        quotient[k] += c
        a = sub(a, [F(0)]*k+scale(b, c))
    return trim(quotient), a

def rem(a, b):
    return divrem(a, b)[1]

def gcd(a, b):
    while trim(b) != [0]:
        a, b = b, rem(a, b)
    return scale(a, 1/a[-1])

def reciprocal(a, n):
    require(a[0] != 0, 'reciprocal has nonzero constant')
    result = [1/a[0]]
    for k in range(1, n+1):
        result.append(-sum(a[j]*result[k-j] for j in range(1, min(k, len(a)-1)+1))/a[0])
    return result

# Solve the coefficient of q^(n+1) in v(q+q^2)=(1+2q)v(q).
# This differs from both the author's operator logarithm and flow inversion.
v = [F(0), F(0), F(1)]
for n in range(3, 21):
    v.append(-sum(v[k]*comb(k, n+1-k) for k in range(2, n)
                  if 0 <= n+1-k <= k)/F(n-2))
source = Path(__file__).resolve().parent.parent
author = json.loads((source/'exact-results.json').read_text())
for n in range(2, 15):
    require(v[n] == F(author['iterative_log_coefficients'][str(n)]),
            f'independent Julia recurrence coefficient {n}')
time_unit = reciprocal(v[2:], 8)
expected_time = {0:F(1), 1:F(1), 2:F(-1,2), 3:F(2,3), 4:F(-13,12),
                 5:F(113,60), 6:F(-1187,360), 7:F(1754,315), 8:F(-14569,1680)}
for k, expected in expected_time.items():
    require(time_unit[k] == expected, f'independent time-form coefficient {k-2}')

# Recover return coefficients from the inherited square factorization, rather
# than differentiating or locally iterating the cubic return.
x = [F(0), F(1)]
P = [F(-1,8), F(-9,4), F(1,2), F(1)]
f = [F(-7,4), F(0), F(1)]
D = sub(f, x)
F3 = x
for unused in range(3):
    F3 = add(mul(F3, F3), [F(-7,4)])
require(sub(F3, x) == mul(D, mul(P, P)), 'independent period-three factorization')
require(gcd(P, D) == [1], 'period-three roots are not fixed')
require(gcd(P, deriv(P)) == [1], 'period-three roots are distinct')
bb, cc, dd = P[2], P[1], P[0]
require(bb*bb*cc*cc-4*cc**3-4*bb**3*dd-27*dd**2+18*bb*cc*dd == 49,
        'independent discriminant')
a = rem(mul(D, mul(deriv(P), deriv(P))), P)
b = rem(add(mul(deriv(D), mul(deriv(P), deriv(P))),
            mul(D, mul(deriv(P), deriv(deriv(P))))), P)
require(a == [F(-21,2), F(-28), F(14)], 'factor-derived quadratic coefficient')
require(b == [F(0), F(56), F(-56)], 'factor-derived cubic coefficient')
require(rem(mul(a,a), P) == [F(0), F(-1372), F(1372)], 'factor-derived a squared')
require(gcd(P, a) == [1], 'quadratic coefficient unit in cubic field')
require(rem(add(scale(mul(a,a), 2), scale(b,49)), P) == [0], 'ratio b/a^2=-2/49')
require(F(1)-F(-2,49) == F(51,49), 'residue 51/49')

# The polynomial identity implies the exact rational momentum defect and its
# first term, for a larger independent finite universe N=2,...,40.
for n in range(2, 41):
    fn = [F(0)]+[F(1)]*n
    require(mul(fn, [F(1),F(-1)]) == [F(0),F(1)]+[F(0)]*(n-1)+[F(-1)],
            f'geometric base identity {n}')
    expected = [F(1)]+[F(0)]*(n-1)+[F(-(n+1)),F(n)]
    require(mul(deriv(fn), [F(1),F(-2),F(1)]) == expected,
            f'exact momentum numerator {n}')

# Invert log(1+t)/t with rational arithmetic. Multiplication by 2*pi*i and
# substitution t=2*i*delta then recover the discrete period coefficients.
log_unit = [F((-1)**k,k+1) for k in range(10)]
inv_log_unit = reciprocal(log_unit, 9)
require(inv_log_unit[:6] == [F(1),F(1,2),F(-1,12),F(1,24),F(-19,720),F(3,160)],
        'reciprocal logarithm coefficients')

def imul(a,b):
    return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])

expected_period = { -1:(F(1),F(0)), 0:(F(0),F(1)), 1:(F(1,3),F(0)),
                    2:(F(0),F(-1,3)), 3:(F(-19,45),F(0)), 4:(F(0),F(3,5)) }
power = (F(0),F(-1,2)) # (2i)^(-1)
for j in range(-1,5):
    value = imul((F(0),F(2)), power)
    value = (value[0]*inv_log_unit[j+1],value[1]*inv_log_unit[j+1])
    require(value == expected_period[j], f'discrete period delta exponent {j}')
    power = imul(power,(F(0),F(2)))

# The unfolding singular coefficients follow from the odd/even reciprocal
# logarithm parts. This checks A through delta^8 and B through delta^6.
expected_A = [F(1),F(1,3),F(-19,45),F(863,945),F(-33953,14175)]
expected_B = [F(1),F(-1,3),F(3,5),F(-275,189)]
for k, expected in enumerate(expected_A):
    require(F((-4)**k)*inv_log_unit[2*k] == expected,
            f'unfolding A delta exponent {2*k}')
for k, expected in enumerate(expected_B):
    require(2*F((-4)**k)*inv_log_unit[2*k+1] == expected,
            f'unfolding B delta exponent {2*k}')

# The shift from a local return to one original step would divide time by 3;
# no such change is made in the return residue 51/49.
require(F(51,49) != 1, 'period-three return and cusp have different residue')

receipt = {'status':'PASS; FINITE-EXACT, not a universal or analytic proof',
           'checks':checks,
           'independent_method':'fractions only; Julia recurrence, factor-derived cubic quotient, geometric identities, reciprocal logarithm',
           'iterative_log_constructed_through':20,
           'iterative_log_compared_through':14,
           'time_form_checked_through':6,
           'jet_family_checked_N':[2,40],
           'period_three_residue':'51/49',
           'unfolding_A_checked_through':8,'unfolding_B_checked_through':6,
           'no_computational_claim':['divergence','all-degree rational invariant exclusion','meromorphic first-integral nonexistence']}
Path(__file__).with_name('independent-results.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(f'{checks} independent exact controls passed.')
