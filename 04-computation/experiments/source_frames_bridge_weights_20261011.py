#!/usr/bin/env python3
"""Exact controls for source frames, branch collisions, and bridge weights.

No external preprint theorem, universal Collatz claim, or JC claim is tested.
The inequalities for all parameters are proved in the companion notes.
"""
import itertools
import json
from fractions import Fraction as Q
from math import gcd, prod

import sympy as s

p, y = s.symbols('p y')
D = p**3-y**2
u = y**2/D
G0 = -u/2
checks = 0


def check(ok, label):
    global checks
    checks += 1
    if not ok:
        raise RuntimeError(label)


def bracket(f, g):
    return s.cancel(-D/p*(s.diff(f, p)*s.diff(g, y)-s.diff(f, y)*s.diff(g, p)))


def polynomial(f):
    return s.cancel(f).is_polynomial(p, y) is True


def genus(a, b, cs):
    """Geometric connected component; caller supplies distinct cusp factors."""
    d = gcd(a, b)
    for c in cs:
        d = gcd(d, c)
    a, b, cs = a//d, b//d, [c//d for c in cs]
    weight = 2*a+3*b+6*sum(cs)
    k = a+b+3*sum(cs)
    orders = [k]+[-c for c in cs]+[sum(cs)-k]
    check(gcd(weight, gcd(*orders)) == 1, 'connected Kummer cover')
    # Independently enumerate the cycles of each sheet permutation.
    cycles = []
    for e in orders:
        unseen = set(range(weight))
        number = 0
        while unseen:
            v = next(iter(unseen))
            while v in unseen:
                unseen.remove(v)
                v = (v+e) % weight
            number += 1
        check(number == gcd(weight, e), 'branch permutation cycle count')
        cycles.append(number)
    euler = 2*weight-sum(weight-n for n in cycles)
    check(euler % 2 == 0, 'even curve Euler characteristic')
    return {'weight': weight, 'orders': orders, 'cycles': cycles,
            'genus': (2-euler)//2}


def nu(r, m):
    return Q(1, (m+1)**r)-Q(1, (m+2)**r)


def source(j, t):
    k = 1 << (18*j)
    return (2264*k-23)//27+96*k*t


def replay(n, word):
    for a in word:
        v = 3*n+1
        actual = (v & -v).bit_length()-1
        if n == 1 or actual != a:
            raise ValueError('non-native or ROOT-padded word')
        n = v >> a
    return n


def main():
    frames = [(p, y), (p+y**2, y), (p, y+p**3),
              (p+y**2, y+(p+y**2)**2), (2*p+y, p+y)]
    candidates = [s.Integer(7), p**2*D, 3+p**2*D*(p+y),
                  p**2, D, p*D, u, 1/p, 1/(p-y)]
    for h1, h2 in frames:
        jac = s.expand(s.diff(h1, p)*s.diff(h2, y)
                       - s.diff(h1, y)*s.diff(h2, p))
        check(jac != 0 and not jac.free_symbols, 'Keller frame')
        for S in candidates:
            b1, b2 = bracket(h1, S), bracket(h2, S)
            rp = s.cancel((s.diff(h2, y)*b1-s.diff(h1, y)*b2)/jac)
            ry = s.cancel((-s.diff(h2, p)*b1+s.diff(h1, p)*b2)/jac)
            check(rp == bracket(p, S) and ry == bracket(y, S),
                  'adjugate recovers coordinate responses')
            direct = all(polynomial(bracket(G0+h, S)) for h in [0, p, y])
            changed = all(polynomial(bracket(G0+h, S)) for h in [0, h1, h2])
            check(direct == changed, 'unit-Jacobian three-test equivalence')

    # Algebraically independent but ramified frames can hide both defects.
    for h1, h2, S in [(p**3, y**2, D), (D, p*D, u)]:
        check(all(polynomial(bracket(G0+h, S)) for h in [0, h1, h2]),
              'ramified frame passes')
        check(not all(polynomial(bracket(G0+h, S)) for h in [0, p, y]),
              'ramified frame loses full source condition')
    check(bracket(D, u) == -6*p*y, 'rational pole hidden by D frame')
    check(bracket(p*D, u) == -8*p**2*y, 'second hidden-pole response')
    check(s.expand(s.diff(D, p)*s.diff(p*D, y)
                   - s.diff(D, y)*s.diff(p*D, p)-2*y*D) == 0,
          'ramification divisor of the rational hostile frame')
    check(s.cancel(bracket(p**2*D, p**2*y*D)
                   + p**2*D**2*(5*p**3-2*y**2)) == 0, 'leading commutator')

    curves = {
        'p2_y_D': genus(2, 1, [1]),
        'p2_D_E': genus(2, 0, [1, 1]),
        'p2_y2_D': genus(2, 2, [1]),
        'p2_y_D2': genus(2, 1, [2]),
        'p2_D2_E_generic': genus(2, 0, [2, 1]),
        'beta_1_p2_D3': genus(2, 0, [3]),
        'beta_0_p5_D2': genus(5, 0, [2]),
    }
    check([v['genus'] for v in curves.values()] == [6, 11, 7, 9, 15, 5, 5],
          'same exponent boxes need not have same genus; collision lowers genus')
    v, z, beta = s.symbols('v z beta')
    transformed = (p**2*D**2*(p**3-beta*y**2)).subs(
        {p: v**2/z, y: v**3/z}, simultaneous=True)
    check(s.cancel(transformed-v**22*(1-z)**2*(1-beta*z)/z**11) == 0,
          'independent rational chart gives the branch orders')
    check(15 == 5*6//2 == 4**2-1 == 3*5 and 3*5+1 == 16,
          'triangular arithmetic boundary')
    profiles = {}
    for exps in [(1,), (3,), (2, 1, 1), (2, 2, 1)]:
        divisors = list(itertools.product(*(range(a+1) for a in exps)))
        proper = [d for d in divisors if any(d) and d != exps]
        squarefree = sum(all(e <= 1 for e in d) for d in proper)
        primes = sum(sum(d) == 1 for d in proper)
        formula = prod(a+1 for a in exps)-2
        check(len(proper) == formula, 'divisor box cardinality')
        check((formula == squarefree+primes) == (exps != (2, 2, 1)),
              'inherited divisor-balance profiles and hostile')
        profiles[str(exps)] = [len(proper), squarefree, primes]

    # Exact polygon controls: holomorphic moments vanish, area current survives.
    triangle = [s.Integer(0), s.Integer(1), s.I]
    tau = s.symbols('tau', real=True)
    area_current = 0
    for a, b in zip(triangle, triangle[1:]+triangle[:1]):
        z = a+(b-a)*tau
        area_current += s.integrate(s.conjugate(z)*(b-a), (tau, 0, 1))
    check(s.simplify(area_current) == s.I, 'triangle current equals twice i area')
    for m in range(25):
        total = sum((b**(m+1)-a**(m+1))/Q(m+1)
                    for a, b in zip(triangle, triangle[1:]+triangle[:1]))
        check(s.simplify(total) == 0, 'closed triangle invisible to holomorphic moments')
        opened = sum((b**(m+1)-a**(m+1))/Q(m+1)
                     for a, b in zip(triangle[:-1], triangle[1:]))
        check(s.simplify(opened-s.I**(m+1)/Q(m+1)) == 0,
              'cut-open tree retains endpoint moment')

    theta = Q(9056, 9153)
    sampled_ratios = []
    for j in range(1, 7):
        k = 1 << (18*j)
        scale = Q(3*k, 8)
        shift = Q(2*(k-1), 27)
        check(scale.denominator == shift.denominator == 1, 'integer affine index map')
        for t in list(range(21))+[127, 1024, 10**9]:
            n = source(j, t)
            h = 223+256*t
            m, b = (n-1)//2, (h-1)//2
            check(Q(m+1) == scale*(b+1)-shift, 'exact source-index transport')
            check(replay(n, (1, 2, 18*j, 1))
                  == replay(h, (1, 1, 1, 1, 3)) == 425+486*t,
                  'independent native common-future replay')
            v = 27*n+23
            check((v & -v).bit_length()-1 == 18*j+3, 'disjoint bridge address')
            check(n % 256 == 75 and h % 256 == 223, 'entry and child guards differ')
            for r in [1, 2, 6, 10]:
                ratio = nu(r, m)/nu(r, b)*scale**(r+1)
                check(1 <= ratio <= theta**(-r-1), 'pointwise distortion bound')
                if r == 6:
                    sampled_ratios.append(ratio)
        # Hostile: deliberately changing one native valuation is rejected.
        try:
            replay(source(j, 0), (1, 2, 18*j+1, 1))
        except ValueError:
            check(True, 'wrong bridge exponent rejected')
        else:
            check(False, 'wrong bridge exponent accepted')

    # Every subset of a small ray box; normalization uses its actual total mass.
    for j in [1, 2]:
        for r in [1, 2, 6]:
            child_weights = [nu(r, 111+128*t) for t in range(8)]
            source_weights = [nu(r, (source(j, t)-1)//2) for t in range(8)]
            cw, sw = sum(child_weights), sum(source_weights)
            for mask in range(256):
                pc = sum(v for t, v in enumerate(child_weights) if mask >> t & 1)/cw
                ps = sum(v for t, v in enumerate(source_weights) if mask >> t & 1)/sw
                check(theta**(r+1)*pc <= ps <= theta**(-r-1)*pc,
                      'every-subset conditional distortion in finite control')

    c6 = Q(8, 3)**7/((1 << 126)-1)
    tail6 = Q(8, 3)**7*Q(1, 1 << (126*6))/((1 << 126)-1)
    partial6 = sum(Q(3*(1 << (18*j)), 8)**-7 for j in range(1, 7))
    check(partial6+tail6 == c6, 'infinite branch sum with exact geometric tail')
    check(theta**-7 < Q(1078, 1000), 'less than 7.8 percent conditional distortion')
    check(0 < c6 < theta**-7*c6 < Q(1, 10**34), 'tiny import mass is not expansion')

    print(json.dumps({'checks': checks, 'curves': curves, 'divisor_profiles': profiles,
                      'bridge_distortion_theta': str(theta),
                      'r6_conditional_distortion_upper': str(theta**-7),
                      'r6_all_branch_scale': str(c6),
                      'r6_sample_scaled_ratio_min': str(min(sampled_ratios)),
                      'r6_sample_scaled_ratio_max': str(max(sampled_ratios)),
                      'scope': 'Finite controls; general claims use companion proofs; Collatz and JC open.'},
                     indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
