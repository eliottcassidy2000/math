#!/usr/bin/env python3
"""Find the first physical owner-distinct defect +/-3 beyond m+h<=13."""

from fractions import Fraction
from itertools import product
from math import gcd


R = Fraction(3, 14)


def nint(x):
    return (2 * x.numerator + x.denominator) // (2 * x.denominator)


def cells(speeds):
    walls = {Fraction(0), Fraction(1)}
    for w in speeds:
        for n in range(-1, w + 2):
            for sign in (-1, 1):
                x = (Fraction(n) + sign * R) / w
                if 0 <= x <= 1:
                    walls.add(x)
    walls = sorted(walls)
    for left, right in zip(walls, walls[1:]):
        if left < right:
            yield (left + right) / 2


def witness(p, b, q, h, m, s, t):
    for y in cells((p, b, q)):
        ns = tuple(nint(w * y) for w in (p, b, q))
        errors = tuple(w * y - n for w, n in zip((p, b, q), ns))
        if not all(abs(e) < R for e in errors):
            continue
        owners = tuple((-pow(w, -1, 3) * n) % 3
                       for w, n in zip((p, b, q), ns))
        if set(owners) != {0, 1, 2}:
            continue
        delta = m * ns[1] - s * h * ns[0] - t * ns[2]
        if delta:
            return y, ns, errors, owners, delta
    return None


def main():
    pairs = [(2, 13), (4, 11), (5, 10), (7, 8)]
    for h, m in pairs:
        found = None
        for max_speed in range(3, 200, 2):
            pool = [w for w in range(1, max_speed + 1, 2) if w % 3]
            for p in pool:
                for q in pool:
                    for s, t in product((-1, 1), repeat=2):
                        numerator = s * h * p + t * q
                        if numerator <= 0 or numerator % m:
                            continue
                        b = numerator // m
                        if b not in pool or len({p, b, q}) != 3:
                            continue
                        if gcd(gcd(p, b), q) != 1:
                            continue
                        hit = witness(p, b, q, h, m, s, t)
                        if hit:
                            found = (max_speed, p, b, q, s, t, hit)
                            break
                    if found:
                        break
                if found:
                    break
            if found:
                break
        print((h, m), found)
    print("PASS")


if __name__ == "__main__":
    main()
