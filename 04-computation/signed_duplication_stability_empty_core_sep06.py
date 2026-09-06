#!/usr/bin/env python3
"""Tiny exact controls for the quartic zero-coefficient stability lemma.

Universe: first three reciprocal roots are multisets from {-3,-2,-1,1,2,3};
tune the fourth rational root to e_2=0, reject undefined/zero fourth roots,
and deduplicate the resulting sorted four-root tuples. No census claim.
"""

from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement
from hashlib import sha256


checks = 0
trace = []


def gate(condition, label):
    global checks
    checks += 1
    if not condition:
        raise RuntimeError(label)


def mul(a, b):
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def elementary(roots):
    out = [F(1)]
    for r in roots:
        out = mul(out, [F(1), F(r)])
    return out


def audit(roots, label):
    a = elementary(roots)
    gate(a[2] == 0, label + ': cancellation')
    e = sum(x*x*y*y for x, y in combinations(roots, 2))
    d = -mul(a, a)[4]
    gate(e > 0, label + ': denominator')
    gate(e == -2*a[1]*a[3] + 2*a[0]*a[4], label + ': denominator identity')
    gate(d == -2*a[1]*a[3] - 2*a[0]*a[4], label + ': square identity')
    gate(d >= F(7, 9)*e, label + ': sharp bound')
    pos = [r for r in roots if r > 0]
    neg = [-r for r in roots if r < 0]
    if len(pos) == len(neg) == 2:
        p, q = pos
        u, v = neg
        pp, uu = p*q, u*v
        alpha, beta = (p-q)**2/pp, (u-v)**2/uu
        z = 3*(alpha + beta) + alpha*beta
        correction = F(4, 9)*(3*uu*(p-q)**2 + 3*pp*(u-v)**2 + (p-q)**2*(u-v)**2)
        gate(d-F(7, 9)*e == correction, label + ': exact stability')
        gate(d/e == (7+z)/(9+z), label + ': normalized stability')
        eps = d/e-F(7, 9)
        gate(z == 81*eps/(2-9*eps), label + ': inverse stability')
        gate(alpha+beta <= 27*eps/(2-9*eps), label + ': pair imbalance')
        # t+1/t >= 4 without introducing a square root.
        gate((pp+uu)**2 >= 16*pp*uu, label + ': root-mean separation')
    else:
        gate(a[4] < 0 and d/e > 1, label + ': sign-count hostile')
    # Three independent real scalar/variable gauges on the coefficient form.
    for scalar, variable in [(F(2), F(-3)), (F(-5, 2), F(2, 3)), (F(1, 7), F(5, 4))]:
        b = [scalar*variable**i*x for i, x in enumerate(a)]
        dd = -mul(b, b)[4]
        ee = b[0]**2 * sum((variable*x)**2*(variable*y)**2 for x, y in combinations(roots, 2))
        gate(dd/ee == d/e, label + ': real gauge')
    # Literal monomial shifts, with the extracted coefficient shifted too.
    for shift in [1, 3]:
        b = [F(0)]*shift + a
        gate(-mul(b, b)[2*(2+shift)] == d, label + ': monomial shift')
    trace.append((label, tuple(str(r) for r in roots), str(d/e)))
    return d/e


rows = set()
rejected = 0
for triple in combinations_with_replacement(map(F, [-3, -2, -1, 1, 2, 3]), 3):
    a = elementary(triple)
    if a[1] == 0 or a[2] == 0:
        rejected += 1
        continue
    fourth = -a[2]/a[1]
    rows.add(tuple(sorted(triple + (fourth,))))
ratios = [audit(row, 'tuned:'+','.join(map(str, row))) for row in sorted(rows)]

audit(tuple(map(F, [1, 1, 1, -1])), 'three-one sign hostile')
audit((F(1), F(1, 2), F(-1, 10), F(-1, 4)), 'inherited positive opposite-coefficient product')

# Exact arithmetic in Q(sqrt(2)); pairs encode a+b*sqrt(2).
def qa(a, b):
    return (a[0]+b[0], a[1]+b[1])


def qm(a, b):
    return (a[0]*b[0]+2*a[1]*b[1], a[0]*b[1]+a[1]*b[0])


def qpoly(a, b):
    out = [(F(0), F(0))]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] = qa(out[i+j], qm(x, y))
    return out


one, zero, root2 = (F(1), F(0)), (F(0), F(0)), (F(0), F(1))
q = [one, root2, (-F(1), F(0))]
h = qpoly(q, q)
hh = qpoly(h, h)
gate(h == [one, (F(0), F(2)), zero, (F(0), -F(2)), one], 'sharp quadratic-field carrier')
gate(hh[4] == (-F(14), F(0)), 'sharp doubled coefficient')
gate(qa(qm((-F(2), F(0)), qm(h[1], h[3])), (F(2), F(0))) == (F(18), F(0)), 'sharp energy')
gate(F(14, 18) == F(7, 9), 'sharp exact ratio')

# Degree five defeats transporting the improved quartic constant unchanged.
five = (F(1), F(1), -F(1, 6), -F(1, 6), -F(13, 60))
five_a = elementary(five)
five_e = sum(x*x*y*y for x, y in combinations(five, 2))
five_ratio = -mul(five_a, five_a)[4]/five_e
gate(five_a[2] == 0, 'degree-five cancellation hostile')
gate(five_ratio == F(18501, 26101) < F(7, 9), 'degree-five margin hostile')

# Existing Gaussian/Hermite k=2 limit at x=1: square coefficient -5/6,
# limiting reciprocal-root product energy 1/2. This is a different degree limit.
gate(-(F(2, 3)-2+F(1, 2))/F(1, 2) == F(5, 3), 'Hermite comparison')

# Boundaries: root-realness and interior condition cannot be dropped.
gate(mul([F(1), F(0), F(1)], [F(1), F(0), F(1)])[2] == 2, 'nonreal-core hostile')
gate(mul([F(0), F(0), F(1)], [F(0), F(0), F(1)])[2] == 0, 'noninterior hostile')

print('PASS signed duplication quartic stability')
print('Tuned universe: 56 root triples; rejected', rejected, '; unique quartic rows', len(rows))
print('Tuned normalized ratio range:', min(ratios), max(ratios))
print('Sharp Q(sqrt(2)) carrier: H=(1+sqrt(2)s-s^2)^2; E=18; doubled=-14; ratio=7/9')
print('Degree-five hostile: (1,1,-1/6,-1/6,-13/60); ratio=18501/26101<7/9')
print('Inherited Hermite k=2 limit at first-coefficient zero: normalized ratio=5/3')
print('Additional controls: inherited mixed pair signs, 3/1 root signs, gauges, shifts, nonreal and noninterior boundaries')
print('Exact gates:', checks)
print('Trace SHA256:', sha256(repr(trace).encode()).hexdigest())
