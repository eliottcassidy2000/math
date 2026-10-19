#!/usr/bin/env python3
"""Weighted formal controls for q -> q+q^2+delta^2, delta fixed."""
from pathlib import Path
import json
import sympy as s

q, d = s.symbols('q delta')
checks = 0


def eq(a, b, label):
    global checks
    if s.cancel(a-b) != 0:
        raise RuntimeError(label + ': ' + str(s.cancel(a-b)))
    checks += 1


def trunc(a, degree):
    poly = s.Poly(s.expand(a), q, d)
    return s.Add(*(coef*q**i*d**j for (i, j), coef in poly.terms() if i+j <= degree))


N = 10
h = q*q+d*d
f = q+h
term, v = q, s.Integer(0)
for j in range(1, N):
    term = trunc(term.subs(q, f)-term, N)
    v += s.Rational((-1)**(j+1), j)*term
v = trunc(v, N)
W, remainder = s.div(v, h, q)
eq(remainder, 0, 'formal logarithm divisible by h')
eq(W.subs({q: 0, d: 0}), 1, 'unit W')
eq(W.subs(d, -d), W, 'even parameter symmetry')
eq(trunc(v, 4), h*(1-q+s.Rational(3, 2)*q*q+d*d/6), 'weighted low-order log')

# Independently reconstruct the original map by the formal flow.
flow, term = q, q
for j in range(1, N):
    term = trunc(v*s.diff(term, q), N)
    flow += term/s.factorial(j)
eq(trunc(flow, N), f, 'weighted exp-log inversion')

reciprocal_unit = s.Integer(0)
term = s.Integer(1)
for j in range(N-1):
    reciprocal_unit += term
    term = trunc(-term*(W-1), N-2)
reciprocal_unit = trunc(reciprocal_unit, N-2)
eq(trunc(W*reciprocal_unit, N-2), 1, 'inverse formal unit')
U, remainder = s.div(reciprocal_unit, h, q)
A = s.expand(remainder).coeff(q, 0)
B = s.expand(remainder).coeff(q, 1)
eq(remainder, A+B*q, 'Weierstrass remainder degree less than two')
eq(reciprocal_unit, A+B*q+h*U, 'full formal division identity')

lp, lm = s.log(1+2*s.I*d), s.log(1-2*s.I*d)
Ae = s.series(s.I*d*(1/lp-1/lm), d, 0, N-1).removeO()
Be = s.series(1/lp+1/lm, d, 0, N-2).removeO()
eq(A, Ae, 'explicit A from logarithmic multipliers')
eq(B, Be, 'explicit B from logarithmic multipliers')
for sign in [-1, 1]:
    actual = s.expand(s.diff(v, q).subs(q, sign*s.I*d))
    expected = s.series(s.log(1+sign*2*s.I*d), d, 0, N).removeO()
    eq(actual, expected, f'fixed-section logarithmic derivative {sign}')

expected_A = 1+d*d/3-19*d**4/45+s.Rational(863, 945)*d**6-s.Rational(33953, 14175)*d**8
expected_B = 1-d*d/3+3*d**4/5-s.Rational(275, 189)*d**6
eq(A, expected_A, 'A coefficients through degree eight')
eq(B, expected_B, 'B coefficients through degree six')
eq(U.subs(d, 0).coeff(q, 0), -s.Rational(1, 2), 'regular time form begins minus one-half')

receipt = {'status': 'FINITE-EXACT controls for a formal theorem, not an analytic transit bound',
           'checks': checks, 'joint_log_degree': N, 'v': str(v), 'W': str(W),
           'A_through_delta_8': str(A), 'B_through_delta_6': str(B),
           'U_truncated': str(U),
           'convergence_claim': 'A and B converge by explicit formulas; no convergence claim for U'}
Path(__file__).with_name('formal-unfolding-exact-results.json').write_text(
    json.dumps(receipt, indent=2)+'\n', encoding='utf-8')
print(f'{checks} exact checks passed; joint logarithm through degree {N}; singular coefficients A and B verified.')
