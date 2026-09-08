#!/usr/bin/env python3
"""Exact controls for dx/sqrt(N), N in C[x], degree at most eight.

The theorem is analytic; this enumerates all 67 multiplicity partitions and
checks its uniform formulas. No source Jacobian existence is inferred.
"""
from collections import Counter
from hashlib import sha256
import json

import sympy as S

gates = []


def need(label, predicate):
    if not bool(predicate):
        raise RuntimeError(label)
    gates.append(label)


def zero(label, expression):
    need(label, S.cancel(expression) == 0)


def partitions(n, upper=None):
    if n == 0:
        yield ()
        return
    if upper is None:
        upper = n
    for first in range(min(n, upper), 0, -1):
        for rest in partitions(n - first, first):
            yield (first,) + rest


def classify(part):
    n = sum(part)
    if n == 0:
        return 'constant'
    if 2 in part:
        return 'reject:finite-simple-pole'
    square = all(m % 2 == 0 for m in part)
    if square:
        return 'pure' if len(part) == 1 else 'reject:split-rational-residue'
    if n == 1:
        return 'pure'
    if n == 2:
        return 'reject:infinity-simple-pole'
    capacity = sum(max(m - 2, 0) for m in part)
    required = n - 2 if n % 2 else n // 2 - 1
    if capacity < required:
        return 'reject:primitive-degree'
    if len(part) == 1:
        return 'pure'
    if len(part) == 2 and all(m % 2 for m in part):
        return 'odd-pair'
    return {
        (4, 1, 1): '411:midpoint',
        (6, 1, 1): '611:quadratic-residue',
        (4, 3, 1): '431:weighted-midpoint',
        (5, 1, 1, 1): '5111:elliptic-two-coefficients',
    }[part]


expected = {
    0: {(): 'constant'},
    1: {(1,): 'pure'},
    2: {},
    3: {(3,): 'pure'},
    4: {(4,): 'pure', (3, 1): 'odd-pair'},
    5: {(5,): 'pure'},
    6: {(6,): 'pure', (5, 1): 'odd-pair',
        (4, 1, 1): '411:midpoint', (3, 3): 'odd-pair'},
    7: {(7,): 'pure'},
    8: {(8,): 'pure', (7, 1): 'odd-pair',
        (6, 1, 1): '611:quadratic-residue', (5, 3): 'odd-pair',
        (5, 1, 1, 1): '5111:elliptic-two-coefficients',
        (4, 3, 1): '431:weighted-midpoint'},
}
bank = []
for n in range(9):
    pp = list(partitions(n))
    need(f'partition count degree {n}', len(pp) == [1, 1, 2, 3, 5, 7, 11, 15, 22][n])
    got = {}
    for part in pp:
        tag = classify(part)
        need(f'partition sums correctly {part}', sum(part) == n)
        bank.append((part, tag))
        if not tag.startswith('reject:'):
            got[part] = tag
            if part and not all(m % 2 == 0 for m in part):
                branch_count = sum(m % 2 for m in part) + n % 2
                genus = (branch_count - 2) // 2
                need(f'surviving genus {part}', genus == (1 if part == (5, 1, 1, 1) else 0))
    need(f'complete survivor list degree {n}', got == expected[n])
need('entire partition universe', len(bank) == 67)

x, z, d = S.symbols('x z d', nonzero=True)
d1, d2, ss, tt = S.symbols('d1 d2 ss tt', nonzero=True)
a, e = S.symbols('a e', nonzero=True)
b, c = S.symbols('b c')

# Finite and infinite local orders are derived from the literal covering maps.
for m in range(1, 9):
    if m % 2:
        order = 1 - m  # x-p=z^2, sqrt(N)=unit*z^m.
        need(f'odd local primitive capacity m={m}', max(-order - 1, 0) == max(m - 2, 0))
    else:
        order = -m // 2  # Two unramified points, each with x-p as parameter.
        need(f'even local primitive capacity m={m}', 2 * max(-order - 1, 0) == max(m - 2, 0))
        need(f'even logarithmic threshold m={m}', (order == -1) == (m == 2))
for n in range(1, 9):
    order = n - 3 if n % 2 else n // 2 - 2
    if n >= 3:
        need(f'infinity local degree n={n}', order + 1 == (n - 2 if n % 2 else n // 2 - 1))
    else:
        need(f'infinity low-degree exception n={n}', order == (-2 if n == 1 else -1))

# If y^2=N and B=R*y, dB=(N R'+N'R/2) dx/y.
def radical_derivative(N, R):
    return S.diff(R, x) * N + R * S.diff(N, x) / 2


for m in [0, 1, 3, 4, 5, 6, 7, 8]:
    N = x ** m
    R = S.Rational(2, 2 - m) * x ** (1 - m)
    zero(f'pure power primitive degree {m}', radical_derivative(N, R) - 1)
need('pure degree two has nonzero logarithmic residue', S.residue(1 / x, x, 0) == 1)
zero('split 4+4 residue at zero', S.residue(1 / (x**2 * (x - d)**2), x, 0) - 2 / d**3)
zero('split 4+4 opposite residue', S.residue(1 / (x**2 * (x - d)**2), x, d) + 2 / d**3)

# All orientations of every two-odd-root survivor; no logarithmic monomial.
for k in range(1, 4):
    for right in range(k + 1):
        left = k - right
        eta = -2 * (z**2 - 1)**(k - 1) * z**(-2 * right) / d**k
        primitive = sum(
            -2 * S.binomial(k - 1, j) * (-1)**(k - 1 - j)
            * z**(2*j - 2*right + 1) / ((2*j - 2*right + 1) * d**k)
            for j in range(k)
        )
        zero(f'odd-pair primitive {(2*left+1, 2*right+1)}', S.diff(primitive, z) - eta)
        need(f'odd-pair no log {(left, right)}', all(2*j - 2*right + 1 != 0 for j in range(k)))
        # Independent substitution x-p=d/(z^2-1), sqrt(N)=x^(k+1)z^(2right+1).
        u = d / (z**2 - 1)
        zero(f'odd-pair change of variable {(left, right)}', S.diff(u, z) / (u**(k+1) * z**(2*right+1)) - eta)

D = x*x + ss*x + tt
# Normalize the residue by removing the nonzero constant square root.
zero('411 residue coefficient', S.diff(D**S.Rational(-1, 2), x).subs(x, 0) + ss / (2 * tt**S.Rational(3, 2)))
zero('611 residue coefficient', S.diff(D**S.Rational(-1, 2), x, 2).subs(x, 0) / 2 - (3*ss**2 - 4*tt) / (8 * tt**S.Rational(5, 2)))
DD = (x + d1) * (x + d2)
# Algebraic logarithmic derivative avoids symbolic square-root branch choices.
zero('431 normalized residue', (-1 / (x+d1) - S.diff(DD, x) / (2*DD)).subs(x, 0) + (d1 + 3*d2) / (2*d1*d2))

R411 = -1 / (tt*x)
zero('411 all-parameter primitive', (D*S.diff(R411, x) + S.diff(D, x)*R411/2 - 1/x**2).subs(ss, 0))
R611 = -1/(2*tt*x**2) + 3*ss/(4*tt**2*x)
res611 = S.factor(D*S.diff(R611, x) + S.diff(D, x)*R611/2 - 1/x**3)
zero('611 exact obstruction factor', res611 - (4*tt - 3*ss**2)/(8*tt**2*x))
AA, BB = -2/(d1*d2*(d1-d2)), -1/(d1*d2)
R431 = (AA*x + BB) / (x*(x+d1))
zero('431 all-parameter primitive', (DD*S.diff(R431, x) + S.diff(DD, x)*R431/2 - 1/(x**2*(x+d1))).subs(d2, -d1/3))

# Rational-coefficient positive controls in each exceptional genus-zero family.
for label, poly, RR, denominator in [
    ('411', x*x + 1, -1/x, x*x),
    ('611', x*x + 2*x + 3, (x-1)/(6*x*x), x**3),
    ('431', (x+3)*(x-1), (x+2)/(6*x*(x+3)), x*x*(x+3)),
]:
    need(f'{label} control squarefree', S.discriminant(poly, x) != 0)
    need(f'{label} control avoids even root', poly.subs(x, 0) != 0)
    zero(f'{label} literal positive primitive', poly*S.diff(RR, x) + S.diff(poly, x)*RR/2 - 1/denominator)
need('411 named hostile residue nonzero', (-ss/(2*tt)).subs({ss: -3, tt: 2}) != 0)
need('611 named hostile residue nonzero', (3*ss**2 - 4*tt).subs({ss: 0, tt: 1}) != 0)
need('431 named hostile residue nonzero', (d1 + 3*d2).subs({d1: 1, d2: 2}) != 0)

# Elliptic stratum: N=x^5(a x^3+b x^2+c x+e), X=1/x,
# Y=sqrt(N)/x^4. The same symbol x now denotes X in the cubic model.
cubic = a + b*x + c*x*x + e*x**3
old_N = z**5 * (a*z**3 + b*z*z + c*z + e)
zero('elliptic inversion curve', old_N.subs(z, 1/x)*x**8 - cubic)
zero('elliptic inversion differential', S.diff(1/x, x)*x**4 + x*x)
zero('elliptic explicit primitive all allowed parameters', (S.diff(cubic, x) * (-S.Rational(2, 3)/e) / 2 + x*x).subs({b: 0, c: 0}))

# Complete pole-space argument: after finite regularity, A(X)+B(X)Y;
# pole orders 2i and 2j+3 have opposite parity, so no top cancellation.
need('elliptic primitive basis at pole budget three',
     [i for i in range(5) if 2*i <= 3] == [0, 1]
     and [j for j in range(5) if 2*j + 3 <= 3] == [0])
V = S.symbols('V')
eq = S.Poly(V * S.diff(cubic, x) / 2 + x*x, x)
zero('elliptic degree-two coefficient fixes V', eq.coeff_monomial(x*x).subs(V, -S.Rational(2, 3)/e))
zero('elliptic degree-one obstruction', eq.coeff_monomial(x).subs(V, -S.Rational(2, 3)/e) + 2*c/(3*e))
zero('elliptic constant obstruction', eq.coeff_monomial(1).subs(V, -S.Rational(2, 3)/e) + b/(3*e))
for label, aa, bb, cc, ee, exact in [
    ('positive', 1, 0, 0, 1, True),
    ('b-hostile', 1, 1, 0, 1, False),
    ('c-hostile', 1, 0, 1, 1, False),
]:
    pol = cubic.subs({a: aa, b: bb, c: cc, e: ee})
    need(f'elliptic {label} actual genus-one control', S.discriminant(pol, x) != 0 and aa*ee != 0)
    need(f'elliptic {label} exactness predicate', (bb == cc == 0) == exact)
zero('elliptic primitive directly in original x', radical_derivative(x**5*(a*x**3+e), -2/(3*e*x**4)) - 1)

# Type and scope controls: the radical extension matters, and a passed leading
# condition does not supply a full mate. The latter elliptic fibre is recorded
# analytically in the note rather than replaced by a bounded mate search.
t = S.symbols('t', nonzero=True)
F, G = x**3*t*t, 1/(x*x*t)
zero('actual rational mate requiring radical leading field', S.diff(F,x)*S.diff(G,t)-S.diff(F,t)*S.diff(G,x)-1)
zero('radical primitive in sqrt(x)', S.diff(-2/z, z) - S.diff(z*z, z)/z**3)
for degree in [1, 2, 4]:
    F, G = t**degree, -x/(degree*t**(degree-1))
    zero(f'constant leading coefficient positive mate degree {degree}', S.diff(F,x)*S.diff(G,t)-S.diff(F,t)*S.diff(G,x)-1)
need('simple all-root elliptic differential is nonzero holomorphic', classify((1,1,1,1)) == 'reject:primitive-degree')

rows = [{'partition': list(p), 'decision': tag} for p, tag in bank]
semantic = sha256(json.dumps({'rows': rows, 'gates': gates}, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
counts = dict(sorted(Counter(tag for _, tag in bank).items()))
print('Boundary differential dx/sqrt(N), C[x], 0 <= degree <= 8')
print('Universe: 67 multiplicity partitions; distinct labelled roots; nonzero scalar')
print('Decisions: ' + json.dumps(counts, sort_keys=True, separators=(',', ':')))
print('Exceptional conditions and explicit primitives: symbolic PASS')
print('Elliptic primitive space: 1, X, Y; exact iff b=c=0')
print(f'Gates: {len(gates)}')
print('Semantic SHA256: ' + semantic)
print('RESULT: PASS')
