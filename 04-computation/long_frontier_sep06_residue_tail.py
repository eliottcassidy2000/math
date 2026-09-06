#!/usr/bin/env python3
"""Exact all-phase certificate on one genuine four-anchor coefficient fibre.

Universe: formal D-residue determinants, complete raw carried response,
four uniform phase intervals, three named fibre controls, and f=6 hostile.
Only standard-library rational arithmetic; no numerical roots or imports
from other producers. Every gate remains active under Python -O.
"""
from fractions import Fraction as F
from itertools import permutations
from math import comb, factorial
import hashlib
import json

GATES = 0


def require(ok, label):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(label)


# Sparse polynomials in A=e3,b=e4,f=e5,s. Negative phase exponents are
# retained during the same-original-root coefficient elimination.
def mon(k=(0, 0, 0, 0), c=1):
    return {k: F(c)} if c else {}


def add(*ps):
    out = {}
    for p in ps:
        for k, c in p.items():
            out[k] = out.get(k, F(0))+c
    return {k: c for k, c in out.items() if c}


def scale(p, c):
    return {k: a*c for k, a in p.items() if a*c}


def mul(p, q):
    out = {}
    for i, a in p.items():
        for j, b in q.items():
            k = tuple(x+y for x, y in zip(i, j))
            out[k] = out.get(k, F(0))+a*b
    return {k: c for k, c in out.items() if c}


def power(p, n):
    out = mon()
    for _ in range(n):
        out = mul(out, p)
    return out


def substitute(p, index, value):
    out = {}
    for k, c in p.items():
        n = k[index]
        require(n >= 0, 'nonnegative eliminated-variable exponent')
        rest = list(k)
        rest[index] = 0
        out = add(out, mul(mon(tuple(rest), c), power(value, n)))
    return out


def determinant(matrix):
    out = {}
    for order in permutations(range(len(matrix))):
        inv = sum(order[i] > order[j] for i in range(len(order)) for j in range(i+1, len(order)))
        p = mon(c=(-1)**inv)
        for i, j in enumerate(order):
            p = mul(p, matrix[i][j])
        out = add(out, p)
    return out


def phase_poly(cs):
    return add(*(mon((0, 0, 0, j), c) for j, c in enumerate(cs)))


def moments(den, num, count):
    out = []
    for n in range(count):
        out.append(add(num[n] if n < len(num) else {},
                       *(scale(mul(den[j], out[n-j]), -1) for j in range(1, min(n+1, len(den))))))
    return out


A, B, Z = [mon(tuple(int(j == i) for j in range(4))) for i in range(3)]
den = [mon(), mon(c=-13), mon(c=55), scale(A, -1), B, scale(Z, -1)]
num = [mon(), mon(c=-11), mon(c=36), scale(A, F(-5, 12)), scale(B, F(1, 7))]
nu = moments(den, num, 5)
require(nu == [mon(), mon(c=2), mon(c=7), add(scale(A, F(7, 12)), mon(c=-19)),
               add(scale(A, F(115, 12)), scale(B, F(-6, 7)), mon(c=-632))], 'general D residue moments')
require(add(mul(nu[1], nu[3]), scale(power(nu[2], 2), -1)) ==
        add(scale(A, F(7, 6)), mon(c=-87)), 'shifted two-moment constraint')
require(determinant([[nu[i+j] for j in range(3)] for i in range(3)]) ==
        scale(add(scale(power(A, 2), 343), scale(A, -67788), scale(B, 2592), mon(c=3157056)), F(-1, 1008)),
        'general three-moment constraint')
require(F(-343*75**2+67788*75-3157056, 2592) == F(-259, 288), 'excluded C-only boundary')

fixed_den = [mon(), mon(c=-13), mon(c=55), mon(c=-84), mon(c=35), scale(Z, -1)]
fixed_num = list(map(lambda c: mon(c=c), [1, -11, 36, -35, 5]))
fixed_nu = moments(fixed_den, fixed_num, 9)
det5 = determinant([[fixed_nu[i+j] for j in range(5)] for i in range(5)])
factor1 = add(power(Z, 2), scale(Z, -15), mon(c=-225))
factor2 = add(power(Z, 2), scale(Z, 32), mon(c=-136))
require(det5 == mul(factor1, factor2), 'factored five-by-five residue determinant')
require(7**5 < 12**4, 'AM-GM fourth-coefficient bound places f below twelve')
require(F(19, 5)**2+32*F(19, 5)-136 > 0, 'algebraic endpoint below nineteen fifths')
require((6**2-15*6-225)*(6**2+32*6-136) == -25668, 'hostile violates D moment positivity')


def convolution(p, q):
    out = [F(0)]*(len(p)+len(q)-1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i+j] += a*b
    return out


def coefficient(p, j):
    return p[j] if 0 <= j < len(p) else F(0)


def row_product(p, q):
    out = [{} for _ in range(len(p)+len(q)-1)]
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i+j] = add(out[i+j], mul(a, b))
    return out


GB = [mon(c=c) for c in [1, 13, 55, 84, 35]]+[Z]
GC = [mon(c=c) for c in [1, 12, 45, 56, 15]]
GD = [mon(c=c) for c in [1, 11, 36, 35, 5]]
O = [F(comb(14, 2*j+1)) for j in range(7)]
E = [F(comb(14, 2*j)) for j in range(8)]
OO, EE = convolution(O, O), convolution(E, E)
GG, CD = row_product(GB, GB), row_product(GC, GD)
raw = {}
Qbar = {}
for j in range(-1, 9):
    alpha = coefficient(OO, j)+coefficient(EE, j+1)
    require(alpha == comb(28, 2*j+2), 'complete alpha binomial identity')
    beta = add(GG[j+2], scale(CD[j+1], 2) if j+1 < len(CD) else {})
    raw[j] = scale(beta, alpha)
    Qbar = add(Qbar, mul(scale(raw[j], -1 if j % 2 else 1), mon((0, 0, 0, j+1))))
require(raw[-1] == mon(c=28), 'retained negative-support carry')
P = add(*(mul(scale(GB[j+1], O[j]*(-1)**j), mon((0, 0, 0, j))) for j in range(5)))
p = add(phase_poly([F(1, 11), -10, 84, -60]), mul(Z, mon((0, 0, 0, 4))))
require(P == scale(p, 2002), 'original first equation')
replace_f = add(mon((0, 0, 0, -1), 60), mon((0, 0, 0, -2), -84),
                mon((0, 0, 0, -3), 10), mon((0, 0, 0, -4), F(-1, 11)))
require(substitute(p, 2, replace_f) == {}, 'elimination at original zero')
h = list(map(F, [-443993, 73031400, -3278853435, 46232902140, -234760993560,
                343030019640, -83518139925, -26087589000, 3585421125]))
require(substitute(Qbar, 2, replace_f) == scale(phase_poly(h), F(-14, 11)), 'full same-zero response independent of f')


def value(poly, s):
    out = F(0)
    for c in reversed(poly):
        out = out*s+c
    return out


def f_value(poly, f):
    require(all(k[0] == k[1] == k[3] == 0 for k in poly), 'univariate coefficient extraction')
    return sum(c*f**k[2] for k, c in poly.items())


def specialize(poly, f):
    require(all(k[0] == k[1] == 0 and k[3] >= 0 for k in poly), 'ordinary phase polynomial')
    out = [F(0)]*(1+max(k[3] for k in poly))
    for k, c in poly.items():
        out[k[3]] += c*f**k[2]
    return out


intervals = [(F(99, 10000), F(1, 100)), (F(1, 9), F(13, 100)), (F(1), F(3, 2))]
endpoint_bank = []
for i, (a, b) in enumerate(intervals):
    signs = (1, -1) if i % 2 == 0 else (-1, 1)
    corners = []
    for f in [F(0), F(5)]:
        pp = specialize(p, f)
        corners.append([value(pp, a), value(pp, b)])
        require(all(x*sign > 0 for x, sign in zip(corners[-1], signs)), 'uniform finite phase signs')
    endpoint_bank.append([[str(x) for x in row] for row in corners])
require(all(value(specialize(p, f), F(10)) < 0 for f in [F(0), F(5)]), 'uniform fourth phase starts after ten')


def transformed(poly, a, b):
    n = len(poly)-1
    out = [F(0)]*(n+1)
    for k, c in enumerate(poly):
        first = [F(comb(k, j))*a**(k-j)*b**j for j in range(k+1)]
        second = [F(comb(n-k, j)) for j in range(n-k+1)]
        for j, x in enumerate(convolution(first, second)):
            out[j] += c*x
    return out


transforms = []
for a, b in intervals:
    row = transformed(h, a, b)
    require(all(x > 0 for x in row), 'positive full finite-interval coefficient transform')
    require(row[0] == value(h, a) and row[-1] == value(h, b), 'transform endpoints independently agree')
    transforms.append(row)
tail = [sum(h[k]*comb(k, j)*10**(k-j) for k in range(j, len(h))) for j in range(len(h))]
require(all(x > 0 for x in tail), 'positive h ten-plus-phase shift')

# The genuine f=1 point is independently reconstructed from all literal
# factorial channels. Other coefficient deformations are not actual rows.
for j in range(5):
    literal = F(factorial(14), factorial(1+j)*factorial(12-3*j)*factorial(1+2*j))
    require(O[j]*f_value(GB[j+1], F(1)) == literal, 'genuine full first charge row')
for j in range(-1, 9):
    literal = F(factorial(28), factorial(2+j)*factorial(24-3*j)*factorial(2+2*j))
    require(f_value(raw[j], F(1)) == literal, 'genuine full doubled charge row')


def enclosure(poly, bounds):
    lo = hi = F(0)
    for c in reversed(poly):
        corners = [lo*bounds[0], lo*bounds[1], hi*bounds[0], hi*bounds[1]]
        lo, hi = min(corners)+c, max(corners)+c
    return lo, hi


def carrier(poly, phase):
    n = len(poly)-1
    out = [F(0)]*(2*n+1)
    for j, c in enumerate(poly):
        out[2*j] = c*phase**(n-j)
    return convolution(out, [F(comb(14, j)) for j in range(15)])


controls = []
for f in [F(0), F(1), F(5), F(6)]:
    pp, qq = specialize(p, f), specialize(Qbar, f)
    bounds = intervals.copy()
    if f == 6:
        bounds[2] = (F(3, 2), F(8, 5))
        bounds.append((F(83465, 10000), F(83466, 10000)))
    elif f:
        high = F(20)
        for _ in range(8):
            if value(pp, high) > 0:
                break
            high *= 2
        require(value(pp, high) > 0, 'named fourth-root finite bracket')
        bounds.append((F(10), high))
    require(len(bounds) == (4 if f else 3), 'degree-drop control root count')
    rows = []
    for index, (lo, hi) in enumerate(bounds):
        require(value(pp, lo)*value(pp, hi) < 0, 'control disjoint sign bracket')
        positive_left = value(pp, lo) > 0
        for _ in range(40):
            mid = (lo+hi)/2
            if (value(pp, mid) > 0) == positive_left:
                lo = mid
            else:
                hi = mid
        qbound = enclosure(qq, (lo, hi))
        if f == 6 and index == 3:
            require(qbound[0] > 0, 'outside-class largest-root sign hostile')
        else:
            require(qbound[1] < 0, 'named full response negative')
        phase = (lo+hi)/2
        HB, HC, HD = [carrier(row, phase) for row in [[-f, 35, -84, 55, -13, 1], [15, -56, 45, -12, 1], [5, -35, 36, -11, 1]]]
        require(coefficient(HB, 9) == -2002*phase*value(pp, phase), 'independent ordinary first-row carrier')
        ordinary = coefficient(convolution(HB, HB), 18)-2*phase*coefficient(convolution(HC, HD), 16)
        require(ordinary == phase*value(qq, phase), 'independent ordinary complete response carrier')
        rows.append({'phase': list(map(str, [lo, hi])), 'response_sign': 1 if qbound[0] > 0 else -1})
    controls.append({'f': str(f), 'roots': rows})

bank = {'universe': 'four-anchor fibre f in [0,5], exact endpoint and coefficient certificates',
        'moment_det': 'f^4+17f^3-841f^2-5160f+30600',
        'finite_intervals': [[str(a), str(b)] for a, b in intervals], 'tail_start': '10',
        'endpoint_corners': endpoint_bank, 'finite_transforms': [list(map(str, row)) for row in transforms],
        'tail_transform': list(map(str, tail)), 'controls': controls}
expanded_det = [f_value(det5, F(z)) for z in range(5)]
require(expanded_det == [F(z**4+17*z**3-841*z**2-5160*z+30600) for z in range(5)], 'reported determinant expansion')
print('STATUS: all original roots have Q(-s)<0 throughout f in [0,5]; general two-anchor class OPEN')
print(json.dumps(bank, sort_keys=True, separators=(',', ':')))
print('gates', GATES)
print('semantic_sha256', hashlib.sha256(json.dumps(bank, sort_keys=True, separators=(',', ':')).encode()).hexdigest())
