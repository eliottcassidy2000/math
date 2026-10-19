#!/usr/bin/env python3
"""Exact finite-jet algebra and explicitly numerical orbit controls.

Reproduce: python3 verify_finite_jet_clocks.py
Also valid under python3 -O. No mathematical claim depends on float controls.
Dependencies: SymPy and mpmath. All outputs are written beside this script.
"""
from __future__ import annotations

from decimal import Decimal, localcontext
from fractions import Fraction as F
import json
import math
from pathlib import Path

import mpmath as mp
import sympy as sp

ROOT = Path(__file__).resolve().parent
checks = 0

def check(condition, label):
    global checks
    if not condition:
        raise RuntimeError(label)
    checks += 1


def mul(a, b, n):
    ans = [F(0)] * (n + 1)
    for i, ai in enumerate(a[:n+1]):
        for j, bj in enumerate(b[:n+1-i]):
            ans[i+j] += ai * bj
    return ans


def power(a, exponent, n):
    """Integer nonnegative power, independently of SymPy."""
    out = [F(1)] + [F(0)] * n
    for _ in range(exponent):
        out = mul(out, a, n)
    return out


def compose(a, b, n):
    out = [F(0)] * (n + 1)
    for k, ak in enumerate(a[:n+1]):
        p = power(b, k, n)
        for j in range(n + 1):
            out[j] += ak * p[j]
    return out


def binomial_series(unit, exponent, n):
    """Rational coefficient generalized binomial of a unit with constant 1."""
    check(unit[0] == 1, 'unit constant')
    z = list(unit) + [F(0)] * (n + 1 - len(unit))
    z[0] -= 1
    out = [F(1)] + [F(0)] * n
    term = out[:]
    binom = F(1)
    for k in range(1, n + 1):
        term = mul(term, z, n)
        binom *= (exponent - k + 1) / k
        for j in range(n + 1):
            out[j] += binom * term[j]
    return out


# A second construction of the inverse: solve U(h(u))=u recursively.
# The Lagrange construction below is compared coefficient by coefficient.
def inverse_recursive(u_series, n):
    h = [F(0)] * (n + 1)
    h[1] = 1 / u_series[1]
    for k in range(2, n + 1):
        residual = compose(u_series, h, k)[k]
        h[k] = -residual / u_series[1]
    check(compose(u_series, h, n) == [F(0), F(1)] + [F(0)] * (n - 1),
          'recursive inverse composition')
    return h


y = sp.Symbol('y')
a = sp.Symbol('a', positive=True)
b5, b6 = sp.symbols('b5 b6')
q1, q2, q3, q4 = sp.symbols('q1 q2 q3 q4')

quartic = sp.expand(sp.series((a + b5*y + b6*y*y)**(-sp.Rational(3,4)),
                             y, 0, 3).removeO()).coeff(y, 2)
quartic_expected = a**(-sp.Rational(3,4)) * (
    sp.Rational(21,32)*b5*b5/(a*a) - sp.Rational(3,4)*b6/a)
check(sp.simplify(quartic-quartic_expected) == 0, 'asymmetric quartic coefficient')
check(sp.simplify(quartic.subs({a:1,b5:0,b6:1})*sp.pi/sp.sqrt(2)
                  +3*sp.pi*sp.sqrt(2)/8) == 0, 'lambda quartic coefficient')
check(quartic.subs({a:1,b5:1,b6:sp.Rational(7,8)}) == 0,
      'different jets, equal divergent quartic clocks')

unit = 1+q1*y+q2*y*y+q3*y**3+q4*y**4
sextic2 = sp.expand(sp.series(unit**(-sp.Rational(1,2)), y, 0, 3).removeO()).coeff(y,2)
check(sp.expand(sextic2-(sp.Rational(3,8)*q1*q1-sp.Rational(1,2)*q2)) == 0,
      'sextic c2')
sextic4 = sp.expand(sp.series(unit**(-sp.Rational(5,6)), y, 0, 5).removeO()).coeff(y,4)
sextic_expected = (sp.Rational(21505,31104)*q1**4-sp.Rational(935,432)*q1*q1*q2
                  +sp.Rational(55,36)*q1*q3+sp.Rational(55,72)*q2*q2
                  -sp.Rational(5,6)*q4)
check(sp.expand(sextic4-sextic_expected) == 0, 'sextic c4')

inverse_results = []
for r in [2, 4, 6, 8, 10]:
    # Q=2^r times a rational unit gives U'(0)=2 exactly.
    n = r-1
    rational_unit = [F(1)] + [F((-1)**(k+1), k+2) for k in range(1,n+1)]
    phi = [2*x for x in binomial_series(rational_unit, F(1,r), n)]
    u_series = [F(0)] + phi[:n]
    h = inverse_recursive(u_series, n)
    row = {'r':r, 'h_prime_even':[]}
    for k in range(n):
        lhs = (k+1)*h[k+1]
        rhs = binomial_series(rational_unit, F(-(k+1),r), k)[k] / 2**(k+1)
        check(lhs == rhs, f'Lagrange formula r={r}, k={k}')
        if k%2 == 0:
            row['h_prime_even'].append(str(lhs))
    inverse_results.append(row)

    # A change to g in degree 2r-1 changes Q only in degree r-1.
    perturbed = rational_unit[:]
    perturbed[r-1] += F(11,7)
    for k in range(0,r-1,2):
        c = binomial_series(rational_unit, F(-(k+1),r), k)[k]
        d = binomial_series(perturbed, F(-(k+1),r), k)[k]
        check(c == d, f'jet independence r={r}, k={k}')

    # Embedded flow and polynomial increment agree through degree 2r-2.
    emb = y*(1-(r-1)*a*y**(r-1))**(-sp.Rational(1,r-1)) - y
    expansion = sp.series(emb,y,0,2*r).removeO().expand()
    check(sp.expand(expansion-a*y**r-sp.Rational(r,2)*a*a*y**(2*r-1)) == 0,
          f'embedded flow coefficient r={r}')

for r in [4,6,8,10]:
    lam = sp.Symbol('lambda')
    last = sp.series((a+lam*y**(r-2))**(-sp.Rational(r-1,r)), y,0,r-1).removeO()
    last = sp.expand(last).coeff(y,r-2)
    expected = -sp.Rational(r-1,r)*lam*a**(-sp.Rational(2*r-1,r))
    check(sp.simplify(last-expected)==0,f'last jet sharpness r={r}')

# Borderline term is exactly a bounded endpoint log ratio.
eps,u = sp.symbols('epsilon u', positive=True)
for r in [2,4,6,8]:
    check(sp.simplify(sp.diff(sp.log(eps+u**r)/r,u)-u**(r-1)/(eps+u**r)) == 0,
          f'borderline primitive r={r}')

# Beta integrals: independent high-precision quadrature is diagnostic only.
mp.mp.dps=60
integral_controls=[]
for r in [2,4,6,8]:
    for j in range(r//2):
        numeric=2*mp.quad(lambda t:t**(2*j)/(1+t**r),[0,1,mp.inf])
        closed=2*mp.pi/r/mp.sin((2*j+1)*mp.pi/r)
        check(abs(numeric-closed)<mp.mpf('1e-45'),f'beta numerical control {r},{j}')
        integral_controls.append({'r':r,'j':j,'absolute_error':mp.nstr(abs(numeric-closed),8)})


def numerical_count(r, qs, n, left, right):
    epsilon=2.0**(-r*n)
    x=float(left)
    target=float(right)
    coeff=[float(q) for q in qs]
    count=0
    while x<target:
        poly=0.0
        for q in reversed(coeff):
            poly=poly*x+q
        x += epsilon+(x**r)*poly
        count+=1
        if count>10_000_000:
            raise RuntimeError('numerical orbit control unexpectedly large')
    return count


def decimal_count(r,qs,n,left,right):
    with localcontext() as ctx:
        ctx.prec=70
        def dec(q):
            q=F(q)
            return Decimal(q.numerator)/Decimal(q.denominator)
        epsilon=Decimal(2)**(-r*n)
        x,target=dec(left),dec(right)
        coeff=[dec(q) for q in qs]
        count=0
        while x<target:
            poly=Decimal(0)
            for q in reversed(coeff):
                poly=poly*x+q
            x += epsilon+(x**r)*poly
            count+=1
        return count


def clock_coefficients(r,qs):
    check(qs[0] == 1, 'numerical examples leading coefficient one')
    out=[]
    for j in range(r//2):
        k=2*j
        ck=binomial_series(qs,F(-(k+1),r),k)[k]
        ckmp=mp.mpf(ck.numerator)/ck.denominator
        out.append(ckmp*2*mp.pi/r/mp.sin((k+1)*mp.pi/r))
    return out

examples=[
    ('quadratic_asymmetric',2,[F(1),F(1,2),F(1,4)],[7,9,11],F(-1,16),F(1,20)),
    ('quartic_pure',4,[F(1)],[3,4,5,6],F(-1,8),F(1,10)),
    ('quartic_equal_clock',4,[F(1),F(1),F(7,8)],[3,4,5,6],F(-1,8),F(1,10)),
    ('quartic_asymmetric',4,[F(1),F(1),F(1,2)],[3,4,5,6],F(-1,8),F(1,10)),
    ('quartic_lambda',4,[F(1),F(0),F(1)],[3,4,5,6],F(-1,8),F(1,10)),
    ('sextic_asymmetric',6,[F(1),F(1,2),F(1,4),F(-1,8),F(1,8)],[3,4],F(-1,8),F(1,10)),
]
orbit_controls=[]
for name,r,qs,ns,left,right in examples:
    coeff=clock_coefficients(r,qs)
    for index,n in enumerate(ns):
        count=numerical_count(r,qs,n,left,right)
        main=sum(A*mp.mpf(2)**((r-1-2*j)*n) for j,A in enumerate(coeff))
        row={'case':name,'r':r,'n':n,'count_float':count,
             'full_clock_residual':float(count-main),
             'leading_only_residual':float(count-coeff[0]*mp.mpf(2)**((r-1)*n))}
        if index==0:
            dec_count=decimal_count(r,qs,n,left,right)
            check(count==dec_count,f'float/70-digit Decimal orbit agreement {name}')
            row['count_decimal70']=dec_count
        orbit_controls.append(row)

# The no-log assertion fails for a one-sided pure-quadratic passage.
one_sided=[]
for n in [6,8,10,12]:
    count=numerical_count(2,[F(1)],n,F(0),F(1,16))
    main=mp.pi/2*mp.mpf(2)**n+n*mp.log(2)
    one_sided.append({'n':n,'count_float':count,
                      'with_log_residual':float(count-main),
                      'without_log_residual':float(count-mp.pi/2*mp.mpf(2)**n)})

result={
    'status':'PASS; exact algebra plus explicitly numerical diagnostics',
    'checks':checks,
    'scope':'Fixed-endpoint constant-additive families; coefficients are not asserted to be conjugacy invariants.',
    'quartic_c2':str(quartic_expected),
    'sextic_c4_normalized':str(sextic_expected),
    'independent_rational_series_inverse_checks':inverse_results,
    'beta_numerical_controls':integral_controls,
    'orbit_controls':orbit_controls,
    'one_sided_log_controls':one_sided,
    'limitations':['Finite calculations do not prove the uniform analytic theorem.',
                  'Floating-point and Decimal counts are supporting numerical controls, not certified exact integer counts.',
                  'No empirical residual bound is promoted to an all-parameter bound.'],
}
(ROOT/'verification-results.json').write_text(json.dumps(result,indent=2)+'\n')
print(f'PASS: {checks} explicit checks (exact algebra and labelled numerical controls).')
print('Independent rational-series inverse checks: r=2,4,6,8,10.')
print('Quartic c2:',quartic_expected)
print('Sextic normalized c4:',sextic_expected)
print('Numerical orbit diagnostics:')
for row in orbit_controls:
    print(f"  {row['case']:24} n={row['n']} N={row['count_float']:8} "
          f"full residual={row['full_clock_residual']: .9f} "
          f"leading residual={row['leading_only_residual']: .9f}")
print('One-sided controls (log retained / log omitted):')
for row in one_sided:
    print(f"  n={row['n']}: {row['with_log_residual']:.9f} / {row['without_log_residual']:.9f}")
print('Detailed output: verification-results.json')
