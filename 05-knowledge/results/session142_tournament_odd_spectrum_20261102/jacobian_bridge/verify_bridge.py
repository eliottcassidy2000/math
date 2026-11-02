#!/usr/bin/env python3
"""Exact bounded checks; standard library only. No repository writes or network."""
from collections import defaultdict
from fractions import Fraction as Q
from itertools import permutations
from math import gcd
import json


def check(ok, message):
    if not ok:
        raise RuntimeError(message)


# Independent sparse coefficient arithmetic in Z[x,y,z].
ZERO = (0, 0, 0)


def const(a):
    return {ZERO: a} if a else {}


def add(*ps):
    out = defaultdict(int)
    for p in ps:
        for mon, coef in p.items():
            out[mon] += coef
    return {m: c for m, c in out.items() if c}


def scale(p, a):
    return {m: a*c for m, c in p.items() if a*c}


def mul(p, q):
    out = defaultdict(int)
    for a, ca in p.items():
        for b, cb in q.items():
            out[tuple(x+y for x, y in zip(a, b))] += ca*cb
    return {m: c for m, c in out.items() if c}


def power(p, n):
    ans = const(1)
    for _ in range(n):
        ans = mul(ans, p)
    return ans


def deriv(p, i):
    out = {}
    for m, c in p.items():
        if m[i]:
            mm = list(m)
            mm[i] -= 1
            out[tuple(mm)] = c*m[i]
    return out


def evaluate(p, xyz):
    return sum((c * xyz[0]**m[0] * xyz[1]**m[1] * xyz[2]**m[2]
                for m, c in p.items()), Q(0))


def determinant(matrix):
    n = len(matrix)
    ans = {}
    for sigma in permutations(range(n)):
        inv = sum(sigma[i] > sigma[j] for i in range(n) for j in range(i+1, n))
        term = const((-1)**inv)
        for i in range(n):
            term = mul(term, matrix[i][sigma[i]])
        ans = add(ans, term)
    return ans


def keller_certificate():
    x, y, z = ({(1, 0, 0): 1}, {(0, 1, 0): 1}, {(0, 0, 1): 1})
    xy = mul(x, y)
    u = add(const(1), xy)
    v = add(const(4), scale(xy, 3))
    f = [add(mul(power(u, 3), z), mul(mul(power(y, 2), u), v)),
         add(y, scale(mul(mul(x, power(u, 2)), z), 3),
             scale(mul(mul(x, power(y, 2)), v), 3)),
         add(scale(x, 2), scale(mul(power(x, 2), y), -3),
             scale(mul(power(x, 3), z), -1))]
    j = [[deriv(p, i) for i in range(3)] for p in f]
    det = determinant(j)
    check(det == const(-2), 'Full polynomial Jacobian determinant failed')
    pts = [(Q(0), Q(0), Q(-1, 4)), (Q(1), Q(-3, 2), Q(13, 2)),
           (Q(-1), Q(3, 2), Q(13, 2))]
    images = [tuple(evaluate(p, t) for p in f) for t in pts]
    check(all(im == (Q(-1, 4), Q(0), Q(0)) for im in images), 'Collision failed')
    # Independent compact 2-variable identity, embedded in the same ring.
    # Here x=w, y=c. Original coordinates have det d(x,w,c)/d(x,y,z)=-x^3.
    w, c = x, y
    a = add(power(w, 2), w, scale(mul(power(w, 3), c), -1))
    b = add(scale(w, 4), const(2), scale(mul(power(w, 2), c), -3))
    d = determinant([[scale(a, -2), deriv(a, 0), deriv(a, 1)],
                     [scale(b, -1), deriv(b, 0), deriv(b, 1)],
                     [c, {}, const(1)]])
    check(d == const(2), 'Weighted-coordinate hand identity failed')
    return {'determinant': -2, 'weighted_coordinate_determinant': 2,
            'F_term_counts': [len(p) for p in f],
            'F_coefficients': [[[list(m), v] for m, v in sorted(p.items())] for p in f],
            'points': [[str(x) for x in p] for p in pts],
            'images': [[str(x) for x in p] for p in images]}


def edges(n):
    return [(i, j) for i in range(n) for j in range(i+1, n)]


def matrix(n, mask, varied=None, t=0):
    a = [[0]*n for _ in range(n)]
    for k, (i, j) in enumerate(edges(n)):
        v = t if k == varied else (mask >> k) & 1
        a[i][j], a[j][i] = v, 1-v
    return a


def ham_dp(a):
    n = len(a)
    dp = [[0]*n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for s in range(1, 1 << n):
        for i in range(n):
            if dp[s][i]:
                for j in range(n):
                    if not (s >> j) & 1 and a[i][j]:
                        dp[s | (1 << j)][j] += dp[s][i]*a[i][j]
    return sum(dp[-1])


def ham_perm(a):
    n = len(a)
    ans = 0
    for p in permutations(range(n)):
        term = 1
        for i in range(n-1):
            term *= a[p[i]][p[i+1]]
        ans += term
    return ans


def tournament_certificate():
    controls = []
    for n, lo, hi, expect in [(5, 4, 260, 5), (6, 20, 28, 19)]:
        bit = (lo ^ hi).bit_length()-1
        check(lo ^ hi == 1 << bit, 'Not a one-arc deformation')
        check(not ((lo >> bit) & 1), 'Required orientation convention failed')
        values = []
        for t in [Q(0), Q(1, 4), Q(1, 2), Q(3, 4), Q(1)]:
            a = matrix(n, lo, bit, t)
            v = ham_dp(a)
            check(v == ham_perm(a) == expect + 4*t, 'Independent weighted path counters differ')
            values.append([str(t), str(v)])
        controls.append({'n': n, 'masks': [lo, hi], 'flipped_edge': edges(n)[bit],
                         'affine_response': [expect, 4], 'values': values,
                         'midpoint': expect+2})
    unimodular_controls = []
    for n, start, end, h0 in [(5, 4, 5, 5), (6, 20, 16, 19)]:
        bit = (start ^ end).bit_length()-1
        check(start ^ end == 1 << bit, 'Unimodular control is not a one-edge flip')
        initial_bit = (start >> bit) & 1
        values = []
        for t in [Q(-1), Q(0), Q(1, 2), Q(1)]:
            p = (1-t) if initial_bit else t
            a = matrix(n, start, bit, p)
            h = ham_dp(a)
            check(h == ham_perm(a) == h0 - 2*t, 'Unimodular guard control failed')
            values.append([str(t), str(h), str((h-1)/2)])
        unimodular_controls.append({'n': n, 'masks': [start, end],
            'flipped_edge': edges(n)[bit], 'normalized_affine_response': [(h0-1)//2, -1],
            'forbidden_target': h0+2, 'unique_integral_inverse': -1,
            'legal_parameter_set': [0, 1], 'values_t_H_Q': values})
    check(ham_dp(matrix(4, 0)) == 1, 'Transitive positive control failed')
    check(ham_dp(matrix(3, 2)) == 3, 'Directed triangle positive control failed')
    finite = {}
    for n in range(1, 7):
        e = n*(n-1)//2
        table = [ham_dp(matrix(n, mask)) for mask in range(1 << e)]
        hs = sorted(set(table))
        check(all(h % 2 for h in hs), 'Parity control failed')
        check(7 not in hs and 21 not in hs, 'Finite forbidden-value check failed')
        coeffs = list(table)
        for i in range(e):
            for s in range(1 << e):
                if (s >> i) & 1:
                    coeffs[s] -= coeffs[s ^ (1 << i)]
        check(coeffs[0] == 1 and all(c % 2 == 0 for c in coeffs[1:]),
              'Coefficientwise parity normalization failed')
        cgcd = 0
        for c in coeffs[1:]:
            cgcd = gcd(cgcd, c)
        finite[str(n)] = {'labeled_tournaments': 1 << e, 'spectrum': hs,
                         'nonconstant_coefficient_gcd': cgcd,
                         'nonzero_coefficients': sum(c != 0 for c in coeffs),
                         'polynomial_degree': max(s.bit_count() for s,c in enumerate(coeffs) if c)}
    return {'arc_interpolation_controls': controls,
            'unimodular_source_guard_controls': unimodular_controls,
            'all_labeled_n_le_6': finite,
            'scope': 'The n<=6 census is bounded, not a proof of all-order gaps or completeness.'}


def arithmetic_controls():
    check(Q(7-5, 4) == Q(21-19, 4) == Q(1, 2), 'Rational preimages failed')
    # Unimodular matrix sends (u,v) to (u+v,u+2v).
    # Target (1,3) has unique inverse (-1,2), violating nonnegative source guard.
    check((-1)+2 == 1 and (-1)+2*2 == 3, 'Unimodular guard hostile failed')
    x, y = {(1, 0, 0): 1}, {(0, 1, 0): 1}
    f = [add(x, power(y, 2)), y]
    check(determinant([[deriv(p, i) for i in range(2)] for p in f]) == const(1),
          'Positive-graded shear determinant failed')
    check(all(all(2*m[0]+m[1] == wt for m in p) for p, wt in zip(f, [2, 1])),
          'Positive-graded shear weighting failed')
    check([evaluate(p, (-1, 1, 0)) for p in f] == [0, 1],
          'Positive-graded shear source-guard witness failed')
    attainable = {2*u+3*v for u in range(20) for v in range(20)}
    check(1 not in attainable and all(v in attainable for v in range(2, 40)),
          'Primitive-lattice nonsaturation control failed')
    return {'unimodular_guard_hostile': {'matrix': [[1, 1], [1, 2]],
            'target': [1, 3], 'integral_inverse': [-1, 2]},
            'positive_weight_unimodular_guard_hostile': {
                'map': ['x+y^2', 'y'], 'weights': [2, 1], 'determinant': 1,
                'target': [0, 1], 'integral_inverse': [-1, 1]},
            'primitive_lattice_hole': {'matrix': [[2, 3]], 'hole': 1,
            'integer_solution': [-1, 1], 'nonnegative_solution': None}}


def main():
    print(json.dumps({'status': 'PASS', 'keller': keller_certificate(),
                     'tournaments': tournament_certificate(),
                     'arithmetic': arithmetic_controls()}, sort_keys=True, indent=2))


if __name__ == '__main__':
    main()
