"""Reproduce the orbit bounds and the modularly prime protected source."""
from itertools import combinations
from math import factorial, gcd, isqrt, prod
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'chanpak_gap_star_ports_20261122'))
from replay import subpaths


def orbit_floor(q):
    triangles = q * (q * q - 1) // 24
    through_vertex = (q * q - 1) // 8
    ordered = total = 1
    for j in range((q + 8) // 9):
        ordered *= triangles - 3 * j * through_vertex
        denominator = factorial(j + 1)
        total += ((ordered << (j + 1)) + denominator - 1) // denominator
    return total


def small_prime_factors(value, bound):
    return [p for p in range(2, bound + 1) if value % p == 0
            and all(p % d for d in range(2, isqrt(p) + 1))]


def perfect_shape(p):
    return (1 << (p - 1)) * ((1 << p) - 1)


def carousel(size):
    return tuple(sum(1 << ((v + j) % size) for j in range(1, size // 2 + 1))
                 for v in range(size))


def protected_source():
    rows = [0] * 18
    for a, b in combinations(range(18), 2):
        i, u = divmod(a, 3)
        j, v = divmod(b, 3)
        if i == j:
            forward = (v - u) % 3 == 1
        elif i // 2 == j // 2:
            forward = u == v
        else:
            forward = (i, j) in {(0, 2), (0, 4)}
        winner, loser = (a, b) if forward else (b, a)
        rows[winner] |= 1 << loser
    return rows


def pair_closure(rows, a, b):
    members = (1 << a) | (1 << b)
    while True:
        enlarged = members
        for v, row in enumerate(rows):
            if not members >> v & 1 and row & members not in (0, members):
                enlarged |= 1 << v
        if enlarged == members:
            return members
        members = enlarged


def check():
    for p, first_excluded_prime, candidates in ((31, 53, [29]), (61, 89, [])):
        h = perfect_shape(p) + 1
        assert orbit_floor(first_excluded_prime) > h
        assert small_prime_factors(h, first_excluded_prime - 1) == candidates
        print(p, 'upper count', h, 'orbit floor at cutoff prime', orbit_floor(first_excluded_prime),
              'possible symmetry primes', candidates)
    assert (perfect_shape(31) + 1) % (29 * 29) != 0
    assert orbit_floor(29) ** 2 > perfect_shape(31) + 1
    factor = prod([11, 13, 31, 41, 61, 151, 331, 1321])
    assert 9 * 175 * factor == (1 << 60) - 1
    assert [subpaths(carousel(q))[-1] for q in (3, 5, 7)] == [3, 15, 175]
    for p in range(3, 200, 2):
        for i in range(8):
            f = (1 << (1 << i)) + 1
            assert (perfect_shape(p) + 1) % f
            assert ((perfect_shape(p) - 1) % f == 0) == ((p - 1) % (2 << i) == 0)
    lucas_trace = 1
    q = 13631489
    for i in range(19):
        two_to_n = 1 << (1 << i)
        assert gcd(two_to_n + 1, lucas_trace) == (q if i == 18 else 1)
        lucas_trace = lucas_trace * lucas_trace - 2 * two_to_n
    assert all(q % d for d in range(2, isqrt(q) + 1))
    period = 524288
    assert (pow(2, period, q * q) - 1) // q == 725967
    for exponent, lift in ((96230, 3578626406374), (428059, 3568199698459)):
        root = pow(2, exponent, q)
        lifted_root = pow(2, lift, q * q)
        assert (root * root - root + 2) % q == 0 and lift % period == exponent
        assert (lifted_root * lifted_root - lifted_root + 2) % (q * q) == 0
    prime_exponent = 9865243
    assert prime_exponent % period == 428059
    assert all(prime_exponent % d for d in range(2, isqrt(prime_exponent) + 1))
    print('first Fermat overlap', 18, q, 'phases', 96230, 428059)
    rows = protected_source()
    for i in range(3):
        perm = [v - v % 3 + (v + 1) % 3 if v // 6 == i else v for v in range(18)]
        assert all((rows[u] >> v & 1) == (rows[perm[u]] >> perm[v] & 1)
                   for u, v in combinations(range(18), 2))
    assert all(pair_closure(rows, a, b) == (1 << 18) - 1 for a, b in combinations(range(18), 2))
    h = subpaths(rows)[-1]
    assert h == 480496617 and (h + 2) % 3
    print('18-vertex modularly prime source', h, 'upper scalar', h + 2)


if __name__ == '__main__':
    check()
