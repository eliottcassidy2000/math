#!/usr/bin/env python3
"""Exact anchored Euler/pencil model obstruction; actual path law stays open.

No repository producer, floating arithmetic, or assertion statements are used.
The finite universe is geometric shapes R=256,512 at q=5,h=4 and an actual
q=2,h=1 control. Exact polynomial interval certificates additionally prove
admissibility for all real R>=256 and positive full response for all R>=384.
"""
from fractions import Fraction as F
from itertools import combinations
from math import comb, gcd, lcm
from functools import reduce
import hashlib
import json

GATES = 0


def require(truth, message):
    global GATES
    GATES += 1
    if not truth:
        raise RuntimeError(message)


def trim(p):
    p = list(map(F, p))
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def conv(p, q):
    out = [F(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i+j] += a*b
    return trim(out)


def coeff(p, j):
    return p[j] if 0 <= j < len(p) else F(0)


def plus(*ps):
    return trim([sum(coeff(p, j) for p in ps) for j in range(max(map(len, ps)))])


def times(p, c):
    return trim([c*a for a in p])


def value(p, x):
    out = F(0)
    for a in reversed(p):
        out = out*x+a
    return out


def interval(p, ab):
    lo = hi = F(0)
    a, b = ab
    for c in reversed(p):
        ends = [lo*a, lo*b, hi*a, hi*b]
        lo, hi = min(ends)+c, max(ends)+c
    return lo, hi


def sign(p, ab, wanted, label):
    lo, hi = interval(p, ab)
    require(lo > 0 if wanted > 0 else hi < 0, label)
    # Exact rational bounds are rounded outwards only for compact printing.
    unit = 10**6
    low = (lo*unit).numerator//(lo*unit).denominator
    high = -((-hi*unit).numerator//(-hi*unit).denominator)
    require(low > 0 if wanted > 0 else high < 0, label+' outward bound')
    return [str(F(low, unit)), str(F(high, unit))]


def from_roots(roots):
    # Independent elementary-symmetric subset construction of monic B.
    q = len(roots)
    b = []
    for i in range(q+1):
        total = F(0)
        for subset in combinations(roots, q-i):
            product = F(1)
            for a in subset:
                product *= a
            total += product
        b.append((-1)**(q-i)*total)
    c = [F(3*(i+1), q+2+2*i)*b[i+1] for i in range(q)]
    d = [F(2+3*i, q+1+2*i)*c[i] for i in range(q)]
    return b, c, d


def residue_certificate(roots, b, c, d):
    derivative = [i*b[i] for i in range(1, len(b))]
    out = []
    for a in roots:
        require(value(b, a) == 0, 'literal B root')
        pair = [value(p, a)/value(derivative, a) for p in (c, d)]
        require(all(v > 0 for v in pair), 'strict C,D/B interlacing')
        out.append(list(map(str, pair)))
    return out


def reverse_sign(p):
    return [(-1)**j*a for j, a in enumerate(reversed(p))]


def laurent_product(p, q):
    out = {}
    for i, a in p.items():
        for j, b in q.items():
            out[i+j] = out.get(i+j, F(0))+a*b
    return {i: a for i, a in out.items() if a}


def laurent_sum(*ps):
    keys = set().union(*(p.keys() for p in ps))
    return {j: sum(p.get(j, F(0)) for p in ps) for j in keys}


def shifted(p, j, factor=1):
    return {i+j: factor*a for i, a in p.items()}


def hadamard(p, q):
    return {j: p[j]*q[j] for j in p.keys() & q.keys()}


def clear_negative_phase(p):
    # Every model response has lower exponent at least -1, including the carry.
    require(min(p) >= -1, 'full raw response lower support')
    return trim([(-1 if j % 2 else 1)*p.get(j, F(0)) for j in range(-1, max(p)+1)])


def make_model(b, c, d, h):
    q = len(b)-1
    x, g, k = q-h, q+2*h+1, 2*h+1
    GB, GC, GD = map(reverse_sign, (b, c, d))
    require(all(a > 0 for p in (GB, GC, GD) for a in p), 'full positive auxiliary coefficients')
    O = {j: F(comb(g, 2*j+1)) for j in range((g-1)//2+1)}
    E = {j: F(comb(g, 2*j)) for j in range(g//2+1)}
    beta, cr, dr = [{j-x: a for j, a in enumerate(p)} for p in (GB, GC, GD)]
    OO, EE = laurent_product(O, O), shifted(laurent_product(E, E), -1)
    BB = laurent_product(beta, beta)
    cross = shifted(laurent_product(cr, dr), 1, 2)
    raw = {'P': hadamard(O, beta), 'V': hadamard(OO, BB),
           'G1': hadamard(EE, BB), 'G2': hadamard(OO, cross), 'G3': hadamard(EE, cross)}
    raw['W'] = laurent_sum(raw['V'], raw['G1'])
    raw['skip'] = laurent_sum(raw['G2'], raw['G3'])
    raw['Q'] = laurent_sum(raw['W'], raw['skip'])
    raw['Q-V'] = laurent_sum(raw['G1'], raw['skip'])
    raw['R0'] = laurent_sum(raw['G1'], raw['G3'])
    raw['K'] = {j: a/2 for j, a in raw['G2'].items()}
    p = [(-1)**j*raw['P'].get(j, F(0)) for j in range(max(raw['P'])+1)]
    rows = {name: clear_negative_phase(v) for name, v in raw.items() if name != 'P'}
    return p, rows, (x, g, k)


def ordinary(p, s, g):
    n = len(p)-1
    even = [F(0)]*(2*n+1)
    for j, a in enumerate(p):
        even[2*j] = a*s**(n-j)
    return conv(even, [F(comb(g, j)) for j in range(g+1)])


def primitive(p):
    denominator = lcm(*(a.denominator for a in p))
    ints = [int(a*denominator) for a in p]
    factor = reduce(gcd, map(abs, ints))
    return [a//factor for a in ints]


def poly_remainder(p, q):
    p = trim(p)
    while len(p) >= len(q) and any(p):
        offset, factor = len(p)-len(q), p[-1]/q[-1]
        for j, a in enumerate(q):
            p[offset+j] -= factor*a
        p = trim(p)
    return p


def gcd_degree(p, q):
    while any(q):
        p, q = q, poly_remainder(p, q)
    return len(trim(p))-1


def symbolic_kernels():
    # Each entry is a polynomial in R, not a specialization or interpolation.
    bs = []
    for i in range(6):
        p = [F(0)]*11
        for subset in combinations(range(5), 5-i):
            p[sum(subset)] += (-1)**(5-i)
        bs.append(trim(p))
    cs = [times(bs[i+1], F(3*(i+1), 7+2*i)) for i in range(5)]
    ds = [times(cs[i], F(2+3*i, 6+2*i)) for i in range(5)]
    return bs, cs, ds


def poly_array_product(p, q):
    out = [[F(0)] for _ in range(len(p)+len(q)-1)]
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i+j] = plus(out[i+j], conv(a, b))
    return out


def reciprocal_shift(p, power):
    require(len(p)-1 <= power, 'reciprocal normalization has no negative power')
    out = [F(0)]*(power+1)
    for j, a in enumerate(p):
        out[power-j] = a
    return trim(out)


def box(p_by_u, t_interval, u_interval):
    lo = hi = F(0)
    for p in reversed(p_by_u):
        low, high = interval(p, t_interval)
        ends = [lo*u_interval[0], lo*u_interval[1], hi*u_interval[0], hi*u_interval[1]]
        lo, hi = min(ends)+low, max(ends)+high
    return lo, hi


def substitute_u(p_by_u, u_polynomial):
    out = [F(0)]
    for p in reversed(p_by_u):
        out = plus(conv(out, u_polynomial), p)
    return out


def bivariate_value(p_by_u, t, u):
    return value([value(p, t) for p in p_by_u], u)


def uniform_certificate():
    bs, cs, ds = symbolic_kernels()
    It = (F(0), F(1, 256))
    interlacer_bounds = []
    for polynomials in (cs, ds):
        for i in range(5):
            signed = times(plus(*([F(0)]*(i*j)+p for j, p in enumerate(polynomials))), (-1)**i)
            normalized = list(reversed(signed))
            lo, hi = interval(normalized, It)
            require(lo > 0, 'uniform strict interlacing for all R>=256')
            require(lo > F(1, 14), 'compact uniform interlacer lower bound')
            interlacer_bounds.append([str(lo), str(hi)])

    # f(t,u)=t^4 p_base(u t^3), with p_base from the ordinary first extraction.
    f = [times(reciprocal_shift(bs[4-j], 4+3*j), comb(14, 2*j+1)) for j in range(5)]
    require(trim([coeff(p, 0) for p in f]) == [F(-14), F(364)], 'exact limiting first root')
    lower = substitute_u(f, [F(1, 26)])
    upper = substitute_u(f, [F(1, 26), F(1, 100)])
    require(lower[0] == upper[0] == 0, 'tube endpoint cancellation at t=0')
    lower_bounds, upper_bounds = interval(lower[1:], It), interval(upper[1:], It)
    require(lower_bounds[1] < 0 and upper_bounds[0] > 0, 'uniform phase tube 1/26<u<1/26+t/100')
    require(lower_bounds[0] > -3 and lower_bounds[1] < F(-29, 10), 'compact lower endpoint bounds')
    require(upper_bounds[0] > F(2, 3) and upper_bounds[1] < F(3, 4), 'compact upper endpoint bounds')
    derivative = [times(f[j], j) for j in range(1, len(f))]
    derivative_bounds = box(derivative, It, (F(0), F(1, 26)+F(1, 25600)))
    require(derivative_bounds[0] > 0, 'uniform monotonicity and smallest-root continuation')
    require(derivative_bounds[0] > 363 and derivative_bounds[1] < 366, 'compact derivative bounds')

    GB = [times(bs[5-j], (-1)**j) for j in range(6)]
    GC = [times(cs[4-j], (-1)**j) for j in range(5)]
    GD = [times(ds[4-j], (-1)**j) for j in range(5)]
    GG, CD = poly_array_product(GB, GB), poly_array_product(GC, GD)
    Sbar = [F(1)]*5
    T = []
    for j in range(-1, 9):
        alpha = (-1 if j % 2 else 1)*comb(28, 2*j+2)
        hit = times(reciprocal_shift(GG[j+2], 6+3*(j+1)), 13*alpha)
        mixed = CD[j+1] if j+1 < len(CD) else [F(0)]
        skip = times(conv(Sbar, reciprocal_shift(mixed, 2+3*(j+1))), 2*alpha)
        T.append(plus(hit, skip))
    limiting_T = trim([coeff(p, 0) for p in T])
    require(limiting_T == [F(0), F(0), F(-1512*3025, 143),
                          F(1512*92391, 143), F(-1512*314847, 143)], 'exact limiting full-response polynomial')
    require(value(limiting_T, F(1, 26)) == F(47439, 48334), 'positive limiting full response')

    # Blow up the known first root: u=1/26+t*z, 0<z<1/100.
    tube_T = [[F(0)] for _ in T]
    for i, p in enumerate(T):
        for j in range(i+1):
            tube_T[j] = plus(tube_T[j], [F(0)]*j+times(p, comb(i, j)*F(1, 26)**(i-j)))
    full_bounds = box(tube_T, (F(0), F(1, 384)), (F(0), F(1, 100)))
    require(full_bounds[0] > 0, 'all real R>=384 have positive anchored full response')
    require(full_bounds[0] > F(1, 40) and full_bounds[1] < F(33, 32), 'compact uniform full-response margin')

    # Verify the symbolic transport against the separate raw Laurent engine.
    # At normalized phase sigma, u=13*t*sigma/Sbar and T=t^2*Sbar*sigma*Q(-sigma).
    for R, phases in [(256, [F(760984, 10**6), F(760985, 10**6)]),
                      (512, [F(1518379, 10**6), F(1518380, 10**6)])]:
        t = F(1, R)
        roots = [F(13*R**i, sum(R**j for j in range(5))) for i in range(5)]
        b, c, d = from_roots(roots)
        p, rows, _ = make_model(b, c, d, 4)
        Sb = value(Sbar, t)
        for sigma in phases:
            u = 13*t*sigma/Sb
            require(value(p, sigma) == -F(13, 1)/Sb*bivariate_value(f, t, u), 'symbolic-to-literal original row')
            require(bivariate_value(T, t, u) == t*t*Sb*value(rows['Q'], sigma), 'symbolic-to-literal full carried response')
            z = (u-F(1, 26))/t
            require(bivariate_value(T, t, u) == bivariate_value(tube_T, t, z), 'exact blown-up response identity')
        if R == 256:
            require(gcd_degree(p, rows['Q']) == 0, 'one coprime specialization proves resultant nonzero')

    return {'interlacer_bounds': interlacer_bounds,
            'f_by_u_ascending': [list(map(str, p)) for p in f],
            'tube_lower_divided_by_t_bounds': list(map(str, lower_bounds)),
            'tube_upper_divided_by_t_bounds': list(map(str, upper_bounds)),
            'first_derivative_bounds': list(map(str, derivative_bounds)),
            'T_limit': list(map(str, limiting_T)), 'T_limit_at_root': '47439/48334',
            'Q_over_R_limit': str(F(338)*value(limiting_T, F(1, 26))),
            'T_tube_positive_bounds_R_ge_384': list(map(str, full_bounds)),
            'T_degree_t': max(map(len, T))-1, 'T_degree_u': len(T)-1,
            'T_terms': sum(bool(a) for p in T for a in p)}


def case(R, bracket, full_sign):
    roots = [F(13*R**i, sum(R**j for j in range(5))) for i in range(5)]
    b, c, d = from_roots(roots)
    residue = residue_certificate(roots, b, c, d)
    require([-b[-2], -c[-2], -d[-2]] == [13, 12, 11], 'genuine three linear anchors')
    p, rows, (x, g, k) = make_model(b, c, d, 4)
    # p is literally P(-s); its sign is the reverse of the ordinary extraction.
    require(value(p, bracket[0])*value(p, bracket[1]) < 0, 'original first-row root bracket')
    derivative = [i*p[i] for i in range(1, len(p))]
    require(value(p, 0) > 0 and interval(derivative, (F(0), bracket[1]))[1] < 0,
            'unique smallest positive phase')
    signs = {name: sign(rows[name], bracket, wanted, str(R)+' '+name)
             for name, wanted in [('V', -1), ('W', -1), ('skip', 1), ('K', 1), ('Q', full_sign)]}
    for j in range(5):
        require(3*(j+1)*b[j+1] == (7+2*j)*c[j], 'first coefficient Euler identity')
        require((6+2*j)*d[j] == (2+3*j)*c[j], 'second coefficient Euler identity')
    for s in bracket:
        HB, HC, HD = [ordinary(poly, s, g) for poly in (b, c, d)]
        require(coeff(HB, k) == (-s)**x*value(p, s), 'literal original common zero map')
        require(coeff(conv(HB, HB), 2*k) == s**(2*x-1)*value(rows['W'], s), 'literal hit extraction')
        require(-2*coeff(conv(HC, HD), 2*k-2) == s**(2*x-2)*value(rows['skip'], s), 'literal skip extraction and carry')
    return {'R': R, 'roots': list(map(str, roots)), 'phase': list(map(str, bracket)),
            'primitive_P_minus_s': primitive(p), 'C_D_residues': residue, 's_times_response_bounds': signs}


bank = [case(256, (F(760984, 10**6), F(760985, 10**6)), -1),
        case(512, (F(1518379, 10**6), F(1518380, 10**6)), 1)]
bank.append({'uniform': uniform_certificate()})
ap, ar, _ = make_model([F(1), F(-4), F(1)], [F(-3), F(1)], [F(-2), F(1)], 1)
require(value(ap, 2) == 0 and value(ar['Q'], 2) == -12610, 'actual (-9,1,6) negative control')

print('STATUS: EXACT anchored-model full Q positive for all R>=384; algebraic same-anchor double cancellation in (256,384); actual path transport OPEN')
for row in bank:
    print(json.dumps(row, sort_keys=True, separators=(',', ':')))
print('gates', GATES)
print('semantic_sha256', hashlib.sha256(json.dumps(bank, sort_keys=True, separators=(',', ':')).encode()).hexdigest())
