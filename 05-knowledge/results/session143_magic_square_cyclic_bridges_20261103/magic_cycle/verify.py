#!/usr/bin/env python3
"""Portable exact controls for the signed-AP obstruction and nine-state cycle.

No external packages, source checkout, floating-point arithmetic, or assertions.
The all-real AP theorem is proved in THEOREM.md; finite controls are not its proof.
"""
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from itertools import permutations, product
from math import isqrt
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def line_totals(grid):
    require(len(grid) == 3 and all(len(row) == 3 for row in grid), "grid shape")
    return ([sum(row) for row in grid]
            + [sum(grid[i][j] for i in range(3)) for j in range(3)]
            + [sum(grid[i][i] for i in range(3)),
               sum(grid[i][2-i] for i in range(3))])


def is_magic(grid):
    return len(set(line_totals(grid))) == 1


def compose(p, q):
    return tuple(p[q[i]] for i in range(len(p)))


def order(p):
    identity = tuple(range(len(p)))
    q = identity
    for k in range(1, 1000):
        q = compose(p, q)
        if q == identity:
            return k
    raise RuntimeError("permutation order exceeded declared finite bound")


def source_pins():
    pins = json.loads((HERE / 'SOURCE_PINS.json').read_text())['pins']
    for pin in pins:
        raw = (HERE / pin['snapshot_path']).read_bytes()
        require(len(raw) == pin['bytes'], 'source byte length: ' + pin['snapshot_path'])
        require(sha256(raw).hexdigest() == pin['sha256'], 'source hash: ' + pin['snapshot_path'])
    return len(pins)


def main():
    checked_sources = source_pins()
    # A polynomial identity, not sampling coefficients: expand each of the
    # nine formal cubes ((k^2-20/3)A+kB)^3 and compare all monomial coefficients.
    moment_coefficients = { (3, 0): Q(0), (2, 1): Q(0),
                            (1, 2): Q(0), (0, 3): Q(0) }
    for k in range(-4, 5):
        x = Q(k*k)-Q(20, 3)
        moment_coefficients[3, 0] += x**3
        moment_coefficients[2, 1] += 3*x*x*k
        moment_coefficients[1, 2] += 3*x*k*k
        moment_coefficients[0, 3] += k**3
    require(moment_coefficients == {(3, 0): Q(2860, 3), (2, 1): 0,
                                     (1, 2): 924, (0, 3): 0},
            'formal cubic centered-moment polynomial identity')
    power_sums = {r: sum(k**r for k in range(-4, 5)) for r in (2, 4, 6)}
    require(power_sums == {2: 60, 4: 708, 6: 9780}, 'independent power sums')
    # Concrete hostile A<0 / nonintegral coefficient and linear positive control.
    A, B, C = Q(-7, 5), Q(11, 3), Q(-19, 2)
    vals = [A*k*k+B*k+C for k in range(-4, 5)]
    mean = sum(vals)/9
    direct = sum((x-mean)**3 for x in vals)
    require(direct == Q(44, 3)*A*(65*A*A+63*B*B) and direct < 0,
            'negative-quadratic rational hostile control')
    centered_lo_shu = [[3, -4, 1], [-2, 0, 2], [-1, 4, -3]]
    require(is_magic([[B*k+C for k in row] for row in centered_lo_shu]),
            'linear rational assignment positive control')

    # General packet moment controls; the all-m formula is proved by power sums.
    for m in range(1, 31):
        q = m*(m+1)
        n = 2*m+1
        av = Q(q, 3)
        a3 = sum((Q(k*k)-av)**3 for k in range(-m, m+1))
        ab2 = 3*sum((Q(k*k)-av)*k*k for k in range(-m, m+1))
        factor = Q(n*q*(4*q-3), 945)
        require(a3 == factor*(4*q-15) and ab2 == 63*factor,
                'general packet moment control m='+str(m))
    require(sorted(3*k*k+k for k in (-1, 0, 1)) == [0, 2, 4], 'three-point sharp boundary')

    # COMPLETE real-parameter odd-cubic classification. The normalized parameter
    # t=b/a solves two linear equations in each of all 384 signed assignments.
    cubic_survivors = []
    cubic_rejections = Counter()
    cubic_branches = 0
    for inds in permutations(range(1, 5)):
        for signs in product((-1, 1), repeat=4):
            cubic_branches += 1
            i, j, h, l = inds
            si, sj, sh, sl = signs
            constants = (si*i**3+sj*j**3-sh*h**3, si*i**3-sj*j**3-sl*l**3)
            slopes = (si*i+sj*j-sh*h, si*i-sj*j-sl*l)
            if any(c != 0 and s == 0 for c, s in zip(constants, slopes)):
                cubic_rejections['nonzero_constant_equation'] += 1
                continue
            candidates = [-Q(c, s) for c, s in zip(constants, slopes) if s != 0]
            require(candidates, 'unexpected identically soluble odd-cubic branch')
            if len(set(candidates)) > 1:
                cubic_rejections['inconsistent_unique_parameters'] += 1
                continue
            t = candidates[0]
            require(all(c+s*t == 0 for c, s in zip(constants, slopes)), 'cubic simultaneous equations')
            vals = [Q(k)**3+t*k for k in range(-4, 5)]
            require(t == -15 and len(set(vals)) == 9, 'only nondegenerate normalized cubic parameter')
            u, v = si*(i**3+t*i), sj*(j**3+t*j)
            grid = [[u, -u-v, v], [-u+v, 0, u-v], [-v, u+v, -u]]
            require(is_magic(grid) and sorted(x for row in grid for x in row) == sorted(vals),
                    'surviving cubic consequence object')
            cubic_survivors.append({'indices': list(inds), 'signs': list(signs), 'b_over_a': str(t)})
    require(cubic_branches == 384 and len(cubic_survivors) == 8
            and sum(cubic_rejections.values()) == 376, 'complete odd-cubic branch partition')
    cubic_control = [[18, -22, 4], [-14, 0, 14], [-4, 22, -18]]
    require(line_totals(cubic_control) == [0]*8, 'root cubic positive witness')
    require(sorted(x for row in cubic_control for x in row)
            == sorted(k**3-15*k for k in range(-4, 5)), 'root cubic witness multiset')
    cycle = [4, 8, 3, 7, 2, 6, 10, 5, 9, 4]
    require(len(set(cycle[:-1])) == 9 and set(cycle[:-1]) == set(range(2, 11)), 'nine states')
    successor = lambda x: 2 + ((x - 2 + 4) % 9)
    require(all(successor(a) == b for a, b in zip(cycle, cycle[1:])), 'cycle transition')
    jumps = [b-a for a, b in zip(cycle, cycle[1:])]
    require(Counter(jumps) == Counter({4: 5, -5: 4}), 'integer wrap invoice')

    points = list(product(range(3), repeat=2))
    index = {p: i for i, p in enumerate(points)}
    T = lambda u, v: ((u+1) % 3, (v+1+u+2*u*u) % 3)
    perm = tuple(index[T(*p)] for p in points)
    for u, v in points:
        require((T(u, v)[0] + 3*T(u, v)[1]) % 9 == (u+3*v+4) % 9, 'digit coordinate map')
        require(T(*T(*T(u, v))) == (u, (v+1) % 3), 'third iterate')
    require(order(perm) == 9, 'cycle order')

    # Exhaust ALL 48 invertible linear parts and all nine translations, no filter.
    affine_orders = Counter()
    linear_parts = 0
    for a, b, c, d in product(range(3), repeat=4):
        if (a*d-b*c) % 3 == 0:
            continue
        linear_parts += 1
        for x, y in points:
            p = tuple(index[((a*u+b*v+x) % 3, (c*u+d*v+y) % 3)] for u, v in points)
            require(len(set(p)) == 9, 'affine map invertibility')
            affine_orders[order(p)] += 1
    expected_orders = {1: 1, 2: 45, 3: 80, 4: 54, 6: 144, 8: 108}
    require(linear_parts == 48 and sum(affine_orders.values()) == 432, 'affine universe')
    require(dict(affine_orders) == expected_orders, 'complete affine order distribution')
    require(9 not in affine_orders, 'no affine nine-cycle')
    suppressed = tuple(index[((u+1) % 3, (v+1) % 3)] for u, v in points)
    require(order(suppressed) == 3, 'carry-suppressed hostile control')

    # Independent exact center-equation classification of all possible signed APs.
    # Normalize d=1. j=0 would say 20/3=0. Each other j gives this unique a/d.
    ap_cases = []
    for j in range(-4, 5):
        if j == 0:
            continue
        a = Q(20-3*j*j, 6*j)
        squares = [(a+k)**2 for k in range(-4, 5)]
        mean = sum(squares) / 9
        require(mean == a*a+Q(20, 3) == squares[j+4], 'AP center equation')
        require(max(squares) > 2*mean, 'AP maximum/opposite obstruction')
        pairs = [(u, v) for u in range(-4, 5) for v in range(u+1, 5)
                 if u != j and v != j and squares[u+4]+squares[v+4] == 2*mean]
        require(not pairs, 'AP center candidates have no noncentral opposite pair')
        ap_cases.append({'center_index': j, 'a_over_d': str(a), 'center_over_d2': str(mean),
                         'maximum_minus_twice_center_over_d2': str(max(squares)-2*mean),
                         'noncentral_opposite_pairs': pairs})
    require(sum(range(2, 11)) == 54 and sum(k*k for k in range(2, 11)) == 384,
            'literal cycle moments')
    require(Q(384, 9) not in {k*k for k in range(2, 11)}, 'literal squared-cycle center absent')

    # Non-vacuity and semantic controls: an ordinary magic square, a constant
    # squared AP at d=0, the inherited square-line-total witness, and a mutation.
    lo_shu = [[9, 2, 7], [4, 6, 8], [5, 10, 3]]
    require(is_magic(lo_shu) and line_totals(lo_shu) == [18]*8, 'shifted Lo Shu positive control')
    squared_lo_shu = [[x*x for x in row] for row in lo_shu]
    require(not is_magic(squared_lo_shu), 'entrywise squaring does not preserve magic')
    require(is_magic([[49]*3 for _ in range(3)]), 'd=0 exception positive control')
    inherited_roots = [[2655, 2976, 320], [2624, 2520, 1665], [1440, 895, 3624]]
    inherited_grid = [[x*x for x in row] for row in inherited_roots]
    flat = [x for row in inherited_grid for x in row]
    totals = line_totals(inherited_grid)
    require(len(set(flat)) == 9 and min(flat) > 0, 'square-total witness distinct positivity')
    require(totals == [4001**2]*6 + [5151**2, 2920**2], 'square-total witness exact totals')
    require(all(isqrt(t)**2 == t for t in totals) and not is_magic(inherited_grid),
            'square-total versus equal-total distinction')
    damaged = [row[:] for row in lo_shu]
    damaged[0][0] += 1
    require(not is_magic(damaged), 'mutated magic witness rejected')

    result = {
        'status': 'PASS',
        'scope': 'Self-contained all-real proof plus independent finite-exact controls; general MSOS remains OPEN.',
        'source_snapshots_verified': checked_sources,
        'quadratic_centered_cubic_coefficients': {
            'A^3': str(moment_coefficients[3, 0]), 'A^2 B': str(moment_coefficients[2, 1]),
            'A B^2': str(moment_coefficients[1, 2]), 'B^3': str(moment_coefficients[0, 3])},
        'index_power_sums': power_sums,
        'negative_quadratic_control_centered_cubic': str(direct),
        'general_packet_moment_control_m_range': [1, 30],
        'odd_cubic_branch_count': cubic_branches,
        'odd_cubic_branch_rejections': dict(sorted(cubic_rejections.items())),
        'odd_cubic_survivors': cubic_survivors,
        'odd_cubic_positive_witness': cubic_control,
        'cycle': cycle, 'cycle_order': 9, 'cycle_jump_counts': dict(sorted(Counter(jumps).items())),
        'cycle_sum': 54, 'cycle_square_sum': 384, 'cycle_square_mean': '128/3',
        'affine_group_universe': 432,
        'affine_element_orders': dict(sorted(affine_orders.items())),
        'carry_suppressed_order': 3,
        'signed_ap_center_candidates': ap_cases,
        'shifted_lo_shu_totals': line_totals(lo_shu),
        'squared_shifted_lo_shu_totals': line_totals(squared_lo_shu),
        'square_line_total_witness_totals': totals,
        'controls': ['ordinary Lo Shu passes', 'd=0 constant-square boundary passes',
                     'mutated Lo Shu rejected', 'carry suppression drops order to3',
                     'square line totals do not imply equality'],
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
