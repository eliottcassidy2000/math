"""Finite sanity checks for circle_antiperiodicity.md (not a formal proof).

Uses only the Python standard library. Exhaustively enumerates signings for
the claimed pinned-angle optimum at small denominators, and selected full
maximum-cut problems. No network or repository writes.
"""
from fractions import Fraction
from math import gcd, lcm

def rational_pair(p, q):
    """Return the finite cyclic group and steps for pi/3 and p*pi/q."""
    if not (q > 0 and gcd(abs(p), q) == 1):
        raise RuntimeError('Exact verification gate failed: circle/check_finite_cases.py:15')
    order_gamma = 2 * q // gcd(abs(p), 2 * q)
    size = lcm(6, order_gamma)
    return (size, size // 6, size * p // (2 * q) % size)

def alpha_orbits(size, a):
    remaining = set(range(size))
    result = []
    while remaining:
        k = min(remaining)
        orbit = []
        while k not in orbit:
            orbit.append(k)
            k = (k + a) % size
        if not len(orbit) == 6:
            raise RuntimeError('Exact verification gate failed: circle/check_finite_cases.py:30')
        remaining.difference_update(orbit)
        result.append(orbit)
    return result

def pinned_optimum(p, q):
    size, a, b = rational_pair(p, q)
    orbits = alpha_orbits(size, a)
    best = size
    for mask in range(1 << len(orbits)):
        sigma = [0] * size
        for i, orbit in enumerate(orbits):
            base = 1 if mask & 1 << i else -1
            for j, k in enumerate(orbit):
                sigma[k] = base if j % 2 == 0 else -base
        if not all((sigma[(k + a) % size] != sigma[k] for k in range(size))):
            raise RuntimeError('Exact verification gate failed: circle/check_finite_cases.py:46')
        bad = sum((sigma[(k + b) % size] == sigma[k] for k in range(size)))
        best = min(best, bad)
    return Fraction(best, size)

def unweighted_optimum(p, q):
    size, a, b = rational_pair(p, q)
    best = 2 * size
    best_pair = None
    for mask in range(1 << size - 1):
        sigma = [1] + [1 if mask & 1 << k else -1 for k in range(size - 1)]
        da = sum((sigma[(k + a) % size] == sigma[k] for k in range(size)))
        db = sum((sigma[(k + b) % size] == sigma[k] for k in range(size)))
        if da + db < best:
            best = da + db
            best_pair = (Fraction(da, size), Fraction(db, size))
    return (size, Fraction(best, size), best_pair)

def main():
    checked = 0
    for q in range(1, 13):
        for p in range(0, 2 * q):
            if gcd(p, q) != 1:
                continue
            expected = Fraction(0) if p % 2 == q % 2 == 1 else Fraction(gcd(3, q), q)
            actual = pinned_optimum(p, q)
            if not actual == expected:
                raise RuntimeError((p, q, actual, expected))
            size, a, b = rational_pair(p, q)
            admissible = [n for n in range(size) if n * a % size == size // 2 and n * b % size == size // 2]
            if not bool(admissible) == (p % 2 == q % 2 == 1):
                raise RuntimeError((p, q, admissible))
            checked += 1
    print(f'PASS: {checked} reduced rational pairs, denominators 1 through 12.')
    print('Pinned pi/3 defect optima and common-frequency compatibility agree.')
    for p, q in [(2, 3), (1, 2), (1, 3), (1, 6)]:
        size, defect, pair = unweighted_optimum(p, q)
        print(f'gamma/pi={p}/{q}; |H|={size}; min(D_pi/3 + D_gamma)={defect}; one minimizing pair={pair}')
if __name__ == '__main__':
    main()
