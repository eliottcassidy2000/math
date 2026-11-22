"""Invert perfect-number decimal suffixes and decode a finite Volchkov reading."""
import argparse
from collections import Counter
from fractions import Fraction
from functools import cache
from itertools import combinations, permutations
from math import factorial, prod
from pathlib import Path

# GIMPS snapshot 2026-11-22, https://www.mersenne.org/primes/
MERSENNE = [
    2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203,
    2281, 3217, 4253, 4423, 9689, 9941, 11213, 19937, 21701, 23209, 44497,
    86243, 110503, 132049, 216091, 756839, 859433, 1257787, 1398269, 2976221,
    3021377, 6972593, 13466917, 20996011, 24036583, 25964951, 30402457,
    32582657, 37156667, 42643801, 43112609, 57885161, 74207281, 77232917,
    82589933, 136279841,
]


def suffix(p, k):
    modulus = 10**k
    x = pow(2, p - 1, modulus)
    return x * (2*x - 1) % modulus


def parity(r, k):
    return ''.join('EO'[int(c) % 2] for c in str(r).zfill(k))


def unit_roots(d, k):
    roots = [i for i in range(1, 5) if (i*i - d) % 5 == 0]
    mod = 5
    for _ in range(1, k):
        new = 5 * mod
        roots = [x + t*mod for x in roots for t in range(5) if ((x + t*mod)**2 - d) % new == 0]
        mod = new
    return roots


def discrete_log_two(unit, k):
    order = 4
    value = next(x for x in range(order) if pow(2, x, 5) == unit % 5)
    for j in range(2, k + 1):
        mod = 5**j
        value = next(value + t*order for t in range(5) if pow(2, value + t*order, mod) == unit % mod)
        order *= 5
    return value, order


def inverse(r, k):
    if r % (1 << k):
        return []
    d = 8*r + 1
    v = 0
    while d % 5 == 0 and v < k:
        v += 1
        d //= 5
    if v == k:
        period = 4 * 5**((k-1)//2)
        return [(period-1, period)]
    if v % 2:
        return []
    a = v // 2
    result = []
    for root in unit_roots(d, k - 2*a):
        w = (5**a) * root
        if w % 5 in (0, 3):
            x, period = discrete_log_two(1 + w, k-a)
            result.append(((x-1) % period, period))
    return sorted(result)


def packing_slots(rows):
    n = len(rows)

    @cache
    def coeff(vertices):
        if not vertices:
            return Counter({0: 1})
        root, *rest = vertices
        result = Counter({m+1: a for m, a in coeff(tuple(rest)).items()})
        for size in range(2, len(rest) + 1, 2):
            for chosen in combinations(rest, size):
                cycles = 0
                for tail in permutations(chosen):
                    cycle = (root, *tail, root)
                    cycles += all(rows[i] >> j & 1 for i, j in zip(cycle, cycle[1:]))
                if cycles:
                    remaining = tuple(i for i in rest if i not in chosen)
                    for m, a in coeff(remaining).items():
                        result[m+1] += 2 * cycles * a
        return result

    return coeff(tuple(range(n)))


def moment(poly, b):
    return sum((Fraction(a) / (m+b+1) for m, a in poly.items()), Fraction(0))


def decoded(n, reading, signed=False):
    L = factorial(n)
    levels = range(n % 2 or 2, n+1, 2)
    denoms = {m: 1+(m+1)*L for m in levels}
    product = prod(denoms.values())
    w = reading * product / L
    assert w.denominator == 1
    result = Counter()
    for m,d in denoms.items():
        a = int(w) * pow(product//d, -1, d) % d
        if signed and a > d//2:
            a -= d
        if a:
            result[m] = a
    return result


def run():
    for k in range(1, 6):
        period = 4*5**(k-1)
        observed = {}
        for p in range(1, period+1, 2):
            r = suffix(p + (k+1)*period, k)
            observed.setdefault(r, set()).add(p)
        for r in range(0, 10**k, 1 << k):
            predicted = {p for a, d in inverse(r, k) for p in range(a, period, d)}
            assert predicted == observed.get(r, set())
    phases = Counter(parity(suffix(p, 4), 4) for p in MERSENNE)
    assert phases == dict(EOEE=15, EOOE=9, OEOE=9, EEOE=7, OOOE=7, OOEE=3, EEEE=2)
    assert sum(suffix(p, 1) == 6 for p in MERSENNE) == 33
    assert inverse(2176, 4) == [(61, 500)] and inverse(1328, 4) == [(19, 20)]
    target = (1 << 60) * ((1 << 61)-1) - 1
    for (c,h,k0,k1), expected in [
        ((81,25,43,68), [[22],[22]]),
        ((411,75,367,392), [[],[]]),
        ((1021,125,915,940), [[60],[35]]),
    ]:
        assert [[y for y in range(h) if (c*y+base-target) % h == 0] for base in (k0,k1)] == expected
    assert [y for y in range(4000) if y % 125 == 35 and y % 32 == 31] == [1535]
    table = Path(__file__).resolve().parents[3] / 'session212/cycle-levels'
    for n in range(1, 6):
        rows = [(0,)] if n == 1 else [tuple(map(int, s.split())) for s in
                                       (table / f'classes{n}.rows').read_text().splitlines()]
        for graph in rows:
            p = packing_slots(graph)
            direct = sum(all(graph[i] >> j & 1 for i,j in zip(w, w[1:]))
                         for w in permutations(range(n)))
            assert sum(p.values()) == direct
            assert decoded(n, moment(p, Fraction(1, factorial(n)))) == p
    a = (6,12,24,33,11,23)
    b = list(a)
    b[0] ^= 1 << 3
    b[3] ^= 1
    pa, pb = packing_slots(a), packing_slots(b)
    diff = Counter({m: pa[m]-pb[m] for m in pa.keys() | pb.keys()})
    assert diff == Counter({4: 2, 2: -2}) and moment(diff, 0) == Fraction(-4,15)
    assert decoded(6, moment(diff, Fraction(1,720)), signed=True) == diff
    print('Decimal inverse exact for every suffix through five digits; known four-digit parity counts:', dict(phases))
    print('Receiver congruences, finite moment decodings and marked six-vertex difference passed.')


if __name__ == '__main__':
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--suffix', help='decimal suffix, with leading zeros retained')
    cli.add_argument('--table', action='store_true', help='show the 52 recorded perfect-number suffixes')
    args = cli.parse_args()
    if args.suffix is not None:
        if not args.suffix.isascii() or not args.suffix.isdigit():
            cli.error('supply one or more decimal digits')
        print(inverse(int(args.suffix), len(args.suffix)))
    elif args.table:
        for i, p in enumerate(MERSENNE, 1):
            r = suffix(p, 4)
            print(i, p, f'{r:04d}', parity(r, 4))
    else:
        run()
