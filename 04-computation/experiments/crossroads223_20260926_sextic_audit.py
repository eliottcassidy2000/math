"""Independent exhaustive projective-chart audit; no geometry-lane imports."""
from math import isqrt
from pathlib import Path
import hashlib


def check(ok, message):
    if not ok:
        raise AssertionError(message)


def main():
    primes = [p for p in range(2, 434)
              if all(p % d for d in range(2, isqrt(p) + 1))]
    empty, empty_torus = [], []
    chart_cases = 0
    controls = {}
    for p in primes:
        sixth = [pow(a, 6, p) for a in range(p)]
        affine = torus = 0
        for x in range(p):
            for y in range(p):
                hit = (sixth[x] + sixth[y] + 1) % p == 0
                affine += hit
                torus += hit and x != 0 and y != 0
                chart_cases += 1
        # At z=0, y cannot be zero; normalize y=1.
        infinity = sum((v + 1) % p == 0 for v in sixth)
        total = affine + infinity
        if not total:
            empty.append(p)
        if not torus:
            empty_torus.append(p)
        if p in (2, 3, 223, 277, 397, 401, 433):
            controls[p] = (total, torus)
    check(empty == [7, 31, 67, 79, 139, 223], 'projective-empty census')
    check(empty_torus == [2, 5, 7, 13, 31, 61, 67, 79, 97, 139, 157, 223, 277],
          'torus-empty census')
    check(402 ** 2 > 400 * 401 and (439 - 17) ** 2 > 400 * 439,
          'exact Hasse-Weil cutoff comparisons')
    print('source_sha256', hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    print('FINITE-EXACT:', len(primes), 'primes <=433;', chart_cases, 'affine-chart cases')
    print('projective_empty', empty)
    print('torus_empty', empty_torus)
    print('controls (projective,torus)', controls)
    print('cutoffs verified using integer squares; extension uses CITED genus/Hasse-Weil')
    print('PASS')


if __name__ == '__main__':
    main()
