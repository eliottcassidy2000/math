#!/usr/bin/env python3
"""Exact controls for the quadratic-anchor compactness and negative-tail proof.

Finite universe: one rational strict-interlacing shape; all four of its
original row roots; the excluded boundary (0,0,3,5,5); two exact evaluations
of an identity proved linear in its symbolic parameter. No producer imports,
floating point, disabled assertions, or claim of a full anchored census.
"""
from fractions import Fraction as F
from itertools import combinations
from math import comb, gcd, lcm
from functools import reduce
import hashlib
import json

GATES = 0


def require(ok, message):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(message)


def trim(p):
    p = list(map(F, p))
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def at(p, s):
    out = F(0)
    for c in reversed(p):
        out = out*s+c
    return out


def coeff(p, j):
    return p[j] if 0 <= j < len(p) else F(0)


def conv(p, q):
    out = [F(0)]*(len(p)+len(q)-1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i+j] += a*b
    return trim(out)


def plus(p, q):
    return trim([coeff(p, j)+coeff(q, j) for j in range(max(len(p), len(q)))])


def scale(p, z):
    return [z*a for a in p]


def interval(p, ab):
    lo = hi = F(0)
    for c in reversed(p):
        ends = [lo*ab[0], lo*ab[1], hi*ab[0], hi*ab[1]]
        lo, hi = min(ends)+c, max(ends)+c
    return lo, hi


def elementary(roots):
    out = []
    for k in range(len(roots)+1):
        total = F(0)
        for subset in combinations(roots, k):
            product = F(1)
            for a in subset:
                product *= a
            total += product
        out.append(total)
    return out


def kernels(e):
    b = [(-1)**(5-j)*e[5-j] for j in range(6)]
    c = [F(3*(j+1), 7+2*j)*b[j+1] for j in range(5)]
    d = [F(2+3*j, 6+2*j)*c[j] for j in range(5)]
    return b, c, d


def response(e):
    e3, e4, e5 = e[3:]
    p = [F(182), F(-20020), 2002*e3, -3432*e4, 2002*e5]
    GB = e
    GC = [F(1), F(12), F(45), F(2, 3)*e3, F(3, 7)*e4]
    GD = [F(1), F(11), F(36), F(5, 12)*e3, F(1, 7)*e4]
    GG, CD = conv(GB, GB), conv(GC, GD)
    q = [F(comb(28, 2*j+2))*(coeff(GG, j+2)+2*coeff(CD, j+1))
         for j in range(-1, 9)]
    # q[0] is raw coefficient at exponent -1. Keep it before clearing.
    qbar = [(-1 if j % 2 else 1)*q[j+1] for j in range(-1, 9)]
    require(q[0] == 28, 'negative-support carry')
    require(q[-1] == comb(28, 18)*e5*e5, 'complete top coefficient')
    require(q[-2] == comb(28, 16)*(2*e4*e5+F(6, 49)*e4*e4), 'complete crossing-dominated coefficient')
    relation = [F(-1, 11), F(10), -e3, F(12, 7)*e4, -e5]
    require(scale(relation, -2002) == p, 'exact original-root reciprocal identity')
    return p, qbar


def carrier(p, s):
    n = len(p)-1
    even = [F(0)]*(2*n+1)
    for j, a in enumerate(p):
        even[2*j] = a*s**(n-j)
    return conv(even, [F(comb(14, j)) for j in range(15)])


def primitive(p):
    denominator = lcm(*(a.denominator for a in p))
    integer = [int(a*denominator) for a in p]
    factor = reduce(gcd, map(abs, integer))
    return [a//factor for a in integer]


roots = [F(a, 5) for a in (1, 3, 9, 22, 30)]
e = elementary(roots)
require(e[1:3] == [13, 55] and sum(a*a for a in roots) == 59, 'two genuine anchors')
require(e[3:] == [F(2127, 25), F(27144, 625), F(3564, 625)], 'named exact shape')
b, c, d = kernels(e)
db = [j*b[j] for j in range(1, 6)]
residues = []
for a in roots:
    require(at(b, a) == 0, 'literal B root')
    pair = [at(p, a)/at(db, a) for p in (c, d)]
    require(all(z > 0 for z in pair), 'both strict interlacers')
    residues.append(list(map(str, pair)))
require(c == [F(3, 7)*e[4], -F(2, 3)*e[3], F(45), F(-12), F(1)], 'induced C anchor')
require(d == [F(1, 7)*e[4], -F(5, 12)*e[3], F(36), F(-11), F(1)], 'induced D anchor')
p, qbar = response(e)
brackets = [(F(9921, 10**6), F(9922, 10**6)),
            (F(121555, 10**6), F(121556, 10**6)),
            (F(1119825, 10**6), F(1119826, 10**6)),
            (F(11804974, 10**6), F(11804975, 10**6))]
controls = []
for lo, hi in brackets:
    require(at(p, lo)*at(p, hi) < 0, 'original row sign bracket')
    bounds = interval(qbar, (lo, hi))
    require(bounds[1] < 0, 'negative full response on entire bracket')
    s = (lo+hi)/2
    HB, HC, HD = [carrier(poly, s) for poly in (b, c, d)]
    require(coeff(HB, 9) == -s*at(p, s), 'independent ordinary first-row extraction')
    require(s*at(qbar, s) == coeff(conv(HB, HB), 18)-2*s*coeff(conv(HC, HD), 16),
            'independent full ordinary response with carry')
    controls.append({'phase': list(map(str, (lo, hi))), 'sQ_bounds': list(map(str, bounds))})
require(len(p)-1 == len(brackets) == 4 and
        all(brackets[j][1] < brackets[j+1][0] for j in range(3)), 'degree exhaustion of all original roots')

# The proposed two-zero boundary preserves C interlacing and fails D.
boundary = list(map(F, (0, 0, 3, 5, 5)))
eb = elementary(boundary)
BB, CC, DD = kernels(eb)
require(eb[1:3] == [13, 55] and eb[3:] == [75, 0, 0], 'boundary anchors and root count')
require(CC == conv(conv([0, 1], [-2, 1]), conv([-5, 1], [-5, 1])), 'C weak boundary interlacing')
require(at(DD, 5) == F(-25, 4), 'D excludes the repeated-root boundary')
require(13**2 > 2*59, 'at least three positive boundary roots')
require(3+5 == 8 and 3*5 == 15, 'forced boundary pair')
# Both sides below have coefficient degree <=1 in E; two evaluations prove
# the identity identically, in addition to its displayed direct factorization.
target = scale(conv([0, 0, 1], conv([-5, 1], [-5, 1])), F(1, 3))
for E in (F(0), F(75)):
    C = [F(0), -F(2, 3)*E, F(45), F(-12), F(1)]
    cubic = [-E, F(55), F(-13), F(1)]
    require(plus(C, scale(target, -1)) == scale([0]+cubic, F(2, 3)), 'generic boundary identity, degree-one certificate')

O = [F(comb(14, 2*j+1)) for j in range(7)]
E = [F(comb(14, 2*j)) for j in range(8)]
OO, EE = conv(O, O), conv(E, E)
for j in range(-1, 14):
    require(coeff(OO, j)+coeff(EE, j+1) == comb(28, 2*j+2), 'full alpha completion identity')

bank = {'roots': list(map(str, roots)), 'elementary_coefficients': list(map(str, e)),
        'interlacer_residues': residues, 'P_minus_s_primitive': primitive(p),
        'all_four_original_root_controls': controls, 'excluded_boundary': [0, 0, 3, 5, 5],
        'D_at_5': '-25/4', 'negative_tail_factor': str(-F(6, 49)*comb(28, 16))}
print('STATUS: exact controls for qualitative uniform negative tail with two anchors; all-root model negativity OPEN')
print(json.dumps(bank, sort_keys=True, separators=(',', ':')))
print('gates', GATES)
print('semantic_sha256', hashlib.sha256(json.dumps(bank, sort_keys=True, separators=(',', ':')).encode()).hexdigest())
