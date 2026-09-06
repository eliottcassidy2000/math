#!/usr/bin/env python3
"""Exact residue-moment and original-first-root certificate.

Universe: formal four-variable carry/elimination identities, displayed
rational inequalities, two named root shapes and the excluded old shape.
No producer imports, floating point, disabled assertions or sampled theorem.
"""
from fractions import Fraction as F
from itertools import combinations, permutations
from math import comb
import hashlib
import json

GATES = 0


def require(ok, label):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(label)


# Sparse polynomials in (e3,e4,e5,s); the s exponent may be negative during
# substitution. This makes the original Laurent normalization explicit.
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


def derivative(p, index):
    out = {}
    for k, c in p.items():
        if k[index]:
            new = list(k)
            new[index] -= 1
            out[tuple(new)] = c*k[index]
    return out


def phase_poly(coefficients):
    return add(*(mon((0, 0, 0, j), c) for j, c in enumerate(coefficients)))


def product_rows(p, q):
    out = [{} for _ in range(len(p)+len(q)-1)]
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i+j] = add(out[i+j], mul(a, b))
    return out


def determinant(matrix):
    out = {}
    for order in permutations(range(len(matrix))):
        inversions = sum(order[i] > order[j] for i in range(len(order)) for j in range(i+1, len(order)))
        p = mon(c=(-1)**inversions)
        for i, j in enumerate(order):
            p = mul(p, matrix[i][j])
        out = add(out, p)
    return out


A, B, Z = [mon(tuple(int(j == i) for j in range(4))) for i in range(3)]
GB = [mon(c=1), mon(c=13), mon(c=55), A, B, Z]
GC = [mon(c=1), mon(c=12), mon(c=45), scale(A, F(2, 3)), scale(B, F(3, 7))]
GD = [mon(c=1), mon(c=11), mon(c=36), scale(A, F(5, 12)), scale(B, F(1, 7))]

# Residues at infinity are compiled from B*(C/B)=C, not assumed moments.
den = [mon(), mon(c=-13), mon(c=55), scale(A, -1), B, scale(Z, -1)]
num = [mon(), mon(c=-12), mon(c=45), scale(A, F(-2, 3)), scale(B, F(3, 7))]
mu = []
for n in range(5):
    mu.append(add(num[n], *(scale(mul(den[j], mu[n-j]), -1) for j in range(1, n+1))))
expected_mu = [mon(), mon(), mon(c=3), add(scale(A, F(1, 3)), mon(c=-16)),
               add(scale(A, F(16, 3)), mon(c=-373), scale(B, F(-4, 7)))]
require(mu == expected_mu, 'exact residue moments')
require(add(mul(mu[1], mu[3]), scale(mul(mu[2], mu[2]), -1)) ==
        add(scale(A, F(1, 3)), mon(c=-25)), 'nonnegative-axis Cauchy determinant')
envelope = add(scale(mul(add(A, mon(c=-75)), add(mon(c=135), scale(A, -1))), F(1, 9)), scale(B, F(-8, 7)))
require(determinant([[mu[i+j] for j in range(3)] for i in range(3)]) == envelope,
        'three-by-three moment determinant')
require(F(7, 72)*900 == F(175, 2), 'coefficient envelope maximum')
require(F(13, 5)**5 < 119, 'fifth coefficient AM-GM bound')

# Compile the full alpha completion from its two parity rows.
O = [F(comb(14, 2*j+1)) for j in range(7)]
E = [F(comb(14, 2*j)) for j in range(8)]


def convolution(p, q):
    out = [F(0)]*(len(p)+len(q)-1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i+j] += a*b
    return out


def coefficient(p, j):
    return p[j] if 0 <= j < len(p) else F(0)


OO, EE = convolution(O, O), convolution(E, E)
GG, CD = product_rows(GB, GB), product_rows(GC, GD)
Qbar = {}
for j in range(-1, 9):
    alpha = coefficient(OO, j)+coefficient(EE, j+1)
    require(alpha == comb(28, 2*j+2), 'full alpha parity identity')
    beta = add(GG[j+2], scale(CD[j+1], 2) if j+1 < len(CD) else {})
    if j == -1:
        require(scale(beta, alpha) == mon(c=28), 'raw lower carry at exponent -1')
    Qbar = add(Qbar, mul(scale(beta, alpha*(-1 if j % 2 else 1)), mon((0, 0, 0, j+1))))
P = add(*(mul(scale(GB[j+1], O[j]*(-1)**j), mon((0, 0, 0, j))) for j in range(5)))
original = add(phase_poly([182, -20020]), mul(scale(A, 2002), mon((0, 0, 0, 2))),
               mul(scale(B, -3432), mon((0, 0, 0, 3))), mul(scale(Z, 2002), mon((0, 0, 0, 4))))
require(P == original, 'original first-row polynomial')
replace_A = add(mon((0, 0, 0, -1), 10), mon((0, 0, 0, -2), F(-1, 11)),
                mon((0, 1, 0, 1), F(12, 7)), mon((0, 0, 1, 2), -1))
require(substitute(P, 0, replace_A) == {}, 'same-zero third-coefficient elimination')
reduced = substitute(Qbar, 0, replace_A)
H = [-63866, 19755450, -1420269695, 19125419885, 9335990950]
expected = add(
    mul(power(B, 2), phase_poly([0]*7+[F(-39330*31920, 49), F(-39330*19601, 49)])),
    mul(mul(B, Z), phase_poly([0]*8+[F(26220*24708, 7), F(26220*9605, 7)])),
    mul(B, phase_poly([0]*3+[F(30*z, 77) for z in [1753752, -166142340, 1898076564, 600168371]])),
    mul(power(Z, 2), phase_poly([0]*9+[-2185*24708, -2185*7735])),
    mul(Z, phase_poly([0]*4+[F(-5*z, 11) for z in [13585572, -716499174, 5260283632, 1724814091]])),
    scale(phase_poly(H), F(-7, 121)))
require(reduced == expected, 'complete exact six-term response identity')
require(all(k[3] >= 0 for k in reduced), 'all Laurent poles cancel after first-row elimination')

# The discarded e5-elimination route really is indefinite for every s>0.
replace_Z = add(mon((0, 1, 0, -1), F(12, 7)), mon((1, 0, 0, -2), -1),
                mon((0, 0, 0, -3), 10), mon((0, 0, 0, -4), F(-1, 11)))
other = substitute(Qbar, 2, replace_Z)
hessian = add(mul(derivative(derivative(other, 0), 0), derivative(derivative(other, 1), 1)),
              scale(power(derivative(derivative(other, 0), 1), 2), -1))
require(hessian == phase_poly([0]*12+[F(-1031232600*z, 49) for z in [144097056, 72692884, 10966105]]),
        'indefinite alternate elimination')


def value(p, s):
    out = F(0)
    for c in reversed(p):
        out = out*s+c
    return out


left, right = F(1, 110), F(1, 90)
require(75-F(12, 7)*F(175, 2)*left == F(810, 11), 'positive first-root left remainder')
right_bound = 182-20020*right+2002*135*right**2+2002*119*right**4
derivative_bound = -20020+4004*135*right+8008*119*right**3
require(right_bound < 0 and derivative_bound < 0, 'first-root existence and uniqueness bounds')
Hprime = [j*H[j] for j in range(1, 5)]
require(12*H[4]*right**2+6*H[3]*right+2*H[2] < 0, 'quartic second derivative negative')
require(value(Hprime, left) < 0, 'quartic derivative starts negative')
require(value(H, right) == F(4379140411, 656100) and value(H, right) > 6000, 'quartic lower bound')
require(13585572-716499174*right > 0, 'negative linear fifth-coefficient term')
linear_bound = F(30, 77)*right**3*(600168371*right**3+1898076564*right**2+1753752)
cross_bound = F(26220, 7)*right**8*(9605*right+24708)
require(linear_bound < 2 and cross_bound < F(1, 1000), 'positive response terms bounded')
margin = -F(42000, 121)+175+F(175, 2)*119/1000
require(margin < -160, 'uniform original-root full-response margin')


def specialize(p, e):
    require(all(k[3] >= 0 for k in p), 'ordinary specialization')
    out = [F(0)]*(max(k[3] for k in p)+1)
    for k, c in p.items():
        out[k[3]] += c*e[3]**k[0]*e[4]**k[1]*e[5]**k[2]
    return out


def elementary(roots):
    out = []
    for j in range(6):
        total = F(0)
        for subset in combinations(roots, j):
            term = F(1)
            for a in subset:
                term *= a
            total += term
        out.append(total)
    return out


def enclosure(p, ab):
    lo = hi = F(0)
    for c in reversed(p):
        ends = [lo*ab[0], lo*ab[1], hi*ab[0], hi*ab[1]]
        lo, hi = min(ends)+c, max(ends)+c
    return lo, hi


def carrier(p, s):
    n = len(p)-1
    even = [F(0)]*(2*n+1)
    for j, c in enumerate(p):
        even[2*j] = c*s**(n-j)
    return convolution(even, [F(comb(14, j)) for j in range(15)])


controls = []
for name, roots in [('strict', [F(a, 5) for a in [1, 3, 9, 22, 30]]),
                    ('C_only_boundary', list(map(F, [0, 0, 3, 5, 5])))]:
    e = elementary(roots)
    require(e[1:3] == [13, 55], 'control anchors')
    b = [(-1)**(5-j)*e[5-j] for j in range(6)]
    c = [F(3*(j+1), 7+2*j)*b[j+1] for j in range(5)]
    d = [F(2+3*j, 6+2*j)*c[j] for j in range(5)]
    if name == 'strict':
        db = [j*b[j] for j in range(1, 6)]
        for a in roots:
            require(value(b, a) == 0 and all(value(poly, a)*value(db, a) > 0 for poly in [c, d]),
                    'strict control interlacers')
    else:
        require(c == convolution([0, 1], convolution([-2, 1], convolution([-5, 1], [-5, 1]))), 'weak C factorization')
        require(value(d, 5) == F(-25, 4), 'D is deliberately not required')
    pp, qq = specialize(P, e), specialize(Qbar, e)
    lo, hi = left, right
    require(value(pp, lo) > 0 and value(pp, hi) < 0, 'independent control phase bracket')
    for _ in range(48):
        mid = (lo+hi)/2
        if value(pp, mid) > 0:
            lo = mid
        else:
            hi = mid
    bounds = enclosure(qq, (lo, hi))
    require(bounds[1] < -160, 'literal full carried response at control first root')
    phase = (lo+hi)/2
    HB, HC, HD = [carrier(poly, phase) for poly in [b, c, d]]
    require(coefficient(HB, 9) == -phase*value(pp, phase), 'ordinary original-row extraction')
    require(phase*value(qq, phase) == coefficient(convolution(HB, HB), 18)-2*phase*coefficient(convolution(HC, HD), 16),
            'independent ordinary full response')
    controls.append({'name': name, 'roots': list(map(str, roots)), 'phase_bracket': list(map(str, [lo, hi])),
                     'sQ_bounds': list(map(str, bounds))})

old_roots = [F(13*512**i, sum(512**j for j in range(5))) for i in range(5)]
require(elementary(old_roots)[1] == 13 and elementary(old_roots)[2] != 55, 'old hostile excluded by second anchor')
bank = {'phase_interval': ['1/110', '1/90'], 'e3_interval': [75, 135], 'e4_upper': '175/2',
        'final_margin_upper': str(margin), 'P_right_upper': str(right_bound), 'P_derivative_upper': str(derivative_bound),
        'H_right': str(value(H, right)), 'positive_linear_upper': str(linear_bound),
        'positive_cross_upper': str(cross_bound), 'controls': controls}
print('STATUS: exact first-phase closure; sQ(-s)<-160 at the unique smallest original root; other phases OPEN')
print(json.dumps(bank, sort_keys=True, separators=(',', ':')))
print('gates', GATES)
print('semantic_sha256', hashlib.sha256(json.dumps(bank, sort_keys=True, separators=(',', ':')).encode()).hexdigest())
