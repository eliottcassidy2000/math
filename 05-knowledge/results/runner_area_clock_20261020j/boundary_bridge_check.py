"""Closed polygon endpoint guard and a literal failed colouring equality."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import lcm
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
CHECKS = 0


def need(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(label)


def dist(t):
    z = t % 1
    return min(z, 1-z)


def adjacent_sides(n, t):
    z = n*(t % 1)
    right = z.numerator // z.denominator
    left = (right-1) % n if z.denominator == 1 else right
    return left, right


def closed_area(n, t):
    return any(s not in (0, n-1) for s in adjacent_sides(n, t))


def colour(z):
    r = z % 8
    return 2*(r//4)+r % 2


def run():
    for n in range(3, 31):
        for j in range(24*n):
            t = F(j, 24*n)
            need(closed_area(n, t) == (dist(t) >= F(1, n)), 'closed endpoint predicate')
    t = F(1, 3)
    need(all(closed_area(3, d*t) for d in (1, 2)), 'AP tight endpoint retained')
    need(not all(adjacent_sides(3, d*t)[1] not in (0, 2) for d in (1, 2)), 'right-rate-only hostile')
    D = (1, 3, 4, 7)
    for z in range(-8, 16):
        for d in D:
            need(colour(z) != colour(z+d), 'exact periodic four-colouring')
    need(all(abs(x-y) in D for x, y in combinations((0, 3, 4, 7), 2)), 'literal four-clique')
    intervals = ((F(0), F(1, 5), 1), (F(1, 5), F(3, 10), 4),
                 (F(3, 10), F(2, 5), 3), (F(2, 5), F(9, 20), 7),
                 (F(9, 20), F(1, 2), 4))
    for lo, hi, d in intervals:
        points = {lo, hi} | {F(j, 2*d) for j in range(d+1) if lo <= F(j, 2*d) <= hi}
        need(max(dist(d*x) for x in points) <= F(1, 5), 'whole interval obstruction')
    candidates = {F(0), F(1)}
    for d in D:
        candidates.update(F(j, 2*d) for j in range(2*d+1))
    for a, b in combinations(D, 2):
        for den in (a+b, abs(a-b)):
            candidates.update(F(j, den) for j in range(den+1))
    values = {t: min(dist(d*t) for d in D) for t in candidates}
    optimum = max(values.values())
    need(optimum == F(1, 5), 'independent affine-envelope maximum')
    need(min(dist(d*F(1, 5)) for d in D) == optimum, 'attainment')
    need(colour(2) != (2*colour(1)) % 4, 'additive equivariance lost')
    for k in range(2, 9):
        for L in range(3, 13):
            m = lcm(*range(2, L+1))
            speeds = (1, *(j*m for j in range(1, k)))
            witness = F(1, 2)+F(1, k*m)
            need(len(set(speeds)) == k, 'distinct LCM speeds')
            need(min(dist(d*witness) for d in speeds) >= F(1, k+1),
                 'finite lonely witness above the denominator obstruction')
            for den in range(1, L+1):
                need(m % den == 0, 'every denominator at most L collides at speed m')
    return {'status': 'FINITE-EXACT controls; general proofs in boundary_and_colouring.md',
            'checks': CHECKS, 'kappa': str(optimum), 'rotation_circumference': 5,
            'fractional_and_ordinary_chromatic': 4,
            'finite_unbounded_denominator_controls': {'k': [2, 8], 'L': [3, 12]},
            'maximizing_candidate_times': sorted(str(t) for t, v in values.items() if v == optimum),
            'checker_sha256_lf': sha256(Path(__file__).read_text(encoding='utf-8').encode()).hexdigest()}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = run()
    path = HERE/'boundary_bridge_receipt.json'
    if args.verify:
        need(result == json.loads(path.read_text(encoding='utf-8')), 'saved boundary receipt')
    else:
        path.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, indent=2))
