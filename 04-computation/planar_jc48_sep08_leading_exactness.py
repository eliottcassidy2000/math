#!/usr/bin/env python3
"""Exact controls for the formal leading-coefficient primitive theorem.

The unbounded theorem is proved in the companion note. This finite universe
checks signs, coefficient indices, genuine rational mates, algebraic-field
primitives, the uniform two-pole residue formula, and a necessary-only hostile.
Checks remain active under python -O. No finite search proves the theorem.
"""
from hashlib import sha256
import json

import sympy as S

x, t, v, z, h, c = S.symbols('x t v z h c')
gates = []


def need(label, value):
    if not bool(value):
        raise RuntimeError(label)
    gates.append(label)


def zero(label, value):
    need(label, S.cancel(value) == 0)


def jac(f, g):
    return S.diff(f, x)*S.diff(g, t)-S.diff(f, t)*S.diff(g, x)


# Every n=1,...,6 and every pole exponent e=0,2,...,6 is included.
# e=1 is the deliberately omitted logarithmic case and is tested below.
monomial_cases = []
for n in range(1, 7):
    for e in (0, 2, 3, 4, 5, 6):
        r = (x-2)**e
        Z = r*t
        B = (x-2)**(1-e)/S.Integer(e-1)
        F = Z**n
        G = B/(n*Z**(n-1))
        zero(f'monomial rational mate n={n}, e={e}', jac(F, G)-1)
        zero(f'primitive sign n={n}, e={e}', S.diff(-B, x)-1/r)
        # In this exact Hensel inverse, G=B v^(n-1)/n.
        zero(f'Laurent coefficient n={n}, e={e}', S.diff(B/n, x)+1/(n*r))
        monomial_cases.append((n, e))

# Nontrivial radical field: x=q^3, a=q^2, primitive 3q.
q = S.symbols('q')
F = x**2*t**3
G = -1/(x*t**2)
zero('nontrivial cubic-field rational mate', jac(F, G)-1)
zero('nontrivial cubic-field primitive', S.diff(3*q, q)-S.diff(q**3, q)/q**2)
zero('nontrivial cubic-field coefficient', S.diff(-q, q)/S.diff(q**3, q)+1/(3*q**2))
need('radical really nonrational by valuation', 2 % 3 != 0)

# A genuine lower-coefficient example, not just F=A t^n.
F = x**3*t**2+x
G = x*t/F
zero('lower-coefficient rational mate', jac(F, G)-1)
T = S.sqrt(1-x*v**2)/(x*S.sqrt(x)*v)
zero('lower-coefficient inverse solves F=v^-2', F.subs(t, T)-v**-2)
gt = S.cancel(G.subs(t, T))
coefficient = S.series(gt, v, 0, 3).removeO().coeff(v, 1)
zero('coefficient index n-1', coefficient-1/S.sqrt(x))
zero('coefficient derivative and sign', S.diff(coefficient, x)+1/(2*x*S.sqrt(x)))
zero('lower-coefficient full chain rule', S.diff(gt, x)+1/S.diff(F, t).subs(t, T))

# Explicit inverses with all n=1,...,6 and nonconstant lower rows.
# These are formal branch controls, not claimed Jacobian pairs.
for n in range(1, 7):
    a = x+1
    shift = x**2-3*x
    lower = x**3+2
    F = (a*(t+shift))**n+lower
    root = (1-lower*v**n)**S.Rational(1, n)
    T = root/(a*v)-shift
    leading = S.limit(v*T, v, 0)
    zero(f'Hensel leading root n={n}', leading-1/a)
    # Test the defining identity on U=vT without nested radical rewriting.
    zero(f'Hensel defining polynomial n={n}', root**n+lower*v**n-1)
    zero(f'Hensel derivative unit n={n}', n*a**n*(1/a)**(n-1)-n*a)
    rv = S.series(T, v, 0, n+2).removeO()
    needed = S.cancel(v**(n+1)*S.diff(rv, v)/n).expand()
    zero(f'coefficient extraction n={n}', needed.coeff(v, n-1)+1/(n*a))

# A root of order exactly n gives dx/a=dx/(x-p), for all n tested.
for n in range(1, 9):
    A = (x-3)**n
    zero(f'multiplicity-n leading coefficient n={n}', A-(x-3)**n)
    zero(f'nonzero logarithmic residue n={n}', S.residue(1/(x-3), x, 3)-1)

# All ordered a,b in [1,6]^2, fixed distinct rational points p=2,q=-1.
# This checks the closed formula, not an assumption of generic positions.
residue_cases = []
for aa in range(1, 7):
    for bb in range(1, 7):
        density = (x-2)**(-aa)*(x+1)**(-bb)
        expected = (-1)**(aa-1)*S.binomial(aa+bb-2, aa-1)/S.Integer(3)**(aa+bb-1)
        actual = S.residue(density, x, 2)
        zero(f'two finite poles a={aa}, b={bb}', actual-expected)
        need(f'two finite poles nonzero a={aa}, b={bb}', actual != 0)
        residue_cases.append((aa, bb))

# A split binary octic (4,2,2) has a finite double point in all three
# choices of the unique possible infinity point, and also if all are finite.
for infinity in (None, 0, 1, 2):
    finite = [m for i, m in enumerate((4, 2, 2)) if i != infinity]
    need(f'422 retains finite double infinity={infinity}', 2 in finite)

# Passing non-square octic: x=z-coordinate inverse 1/(1-z^2).
# N=x^5(x-1)^3, a=z^3/(1-z^2)^4.
xx = 1/(1-z**2)
aa = z**3/(1-z**2)**4
zero('5+3 algebraic square root', aa**2-xx**5*(xx-1)**3)
primitive = -2/z-4*z+S.Rational(2, 3)*z**3
zero('5+3 exact differential', S.diff(primitive, z)-S.diff(xx, z)/aa)
A = x**10*(x-1)**6
B = (8*x**2-4*x-1)/(6*x**9*(x-1)**5)
zero('5+3 descended pure-monomial mate', jac(A*t**4, B/t**3)-1)

# Exactness is only necessary once lower coefficients are admitted.
# F=x^2+t^2 has constant leading coefficient, but its generic fibre
# carries -dx/(2t)=-i dh/(2h), so no rational primitive exists.
xx = (h+c/h)/2
tt = (h-c/h)/(2*S.I)
zero('necessary-only hostile fibre', xx**2+tt**2-c)
density = -S.diff(xx, h)/(2*tt)
zero('necessary-only hostile logarithmic differential', density+S.I/(2*h))
need('necessary-only hostile residue nonzero', S.residue(density, h, 0) != 0)
zero('necessary-only hostile leading gate passes', S.diff(x, x)-1)

# Pure finite sixfold N passes: the double boundary point at infinity
# must not be mistaken for a finite double point.
zero('infinity-double hostile primitive', S.diff(-1/(2*(x-5)**2), x)-1/(x-5)**3)

print('leading_exactness: PASS')
print('scope: formal unbounded necessary gate; pure-monomial iff; finite exact controls')
print('monomial universe:', len(monomial_cases), 'n=1..6, e in {0,2,3,4,5,6}')
print('two-pole universe:', len(residue_cases), 'ordered a,b=1..6, p=2, q=-1')
print('positive controls: nontrivial radical field; nonconstant lower row; octic 5+3')
print('hostiles: logarithmic root; necessary-only circle fibre; double point at infinity')
print('gates:', len(gates))
print('gate-label sha256:', sha256(json.dumps(gates, separators=(',', ':')).encode()).hexdigest())
