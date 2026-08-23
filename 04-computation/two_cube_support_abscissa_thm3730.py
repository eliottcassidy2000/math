#!/usr/bin/env python3
"""Exact hostile probe for the positive distinct two-cube support.

This is session scratch.  It compares the indexed pair fibre with the
deduplicated integer support and independently rechecks every represented
integer in the bounded universe with THM-463's good-divisor criterion, plus
an unrepresented hostile.  Floating Gamma and harmonic diagnostics are
controls; the fibre and divisor gates are exact.
"""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from math import gamma, isqrt, fsum


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def icbrt(n: int) -> int:
    """Floor cube root, using a floating seed followed by exact repair."""
    if n < 0:
        raise ValueError("icbrt expects n>=0")
    x = int(round(n ** (1.0 / 3.0))) if n else 0
    while (x + 1) ** 3 <= n:
        x += 1
    while x ** 3 > n:
        x -= 1
    return x


def indexed_count(cutoff: int) -> int:
    """#{1<=x<y:x^3+y^3<=cutoff}, without materializing the fibre."""
    total = 0
    for x in range(1, icbrt(cutoff) + 1):
        rem = cutoff - x**3
        if rem <= 0:
            break
        ymax = icbrt(rem)
        if ymax > x:
            total += ymax - x
    return total


def pair_fibre(cutoff: int) -> Counter[int]:
    fibre: Counter[int] = Counter()
    for x in range(1, icbrt(cutoff) + 1):
        rem = cutoff - x**3
        if rem <= 0:
            break
        ymax = icbrt(rem)
        for y in range(x + 1, ymax + 1):
            fibre[x**3 + y**3] += 1
    return fibre


def divisors(n: int) -> list[int]:
    lo: list[int] = []
    hi: list[int] = []
    for d in range(1, isqrt(n) + 1):
        if n % d == 0:
            lo.append(d)
            if d * d != n:
                hi.append(n // d)
    return lo + hi[::-1]


def positive_distinct_good_divisors(n: int) -> list[tuple[int, int, int]]:
    """Return (d,x,y), x<y positive, from the exact THM-463 criterion."""
    reps: list[tuple[int, int, int]] = []
    for d in divisors(n):
        numerator = d**3 - n
        denominator = 3 * d
        if numerator % denominator:
            continue
        s = numerator // denominator
        delta = d * d - 4 * s
        if delta <= 0:
            continue
        e = isqrt(delta)
        if e * e != delta or (d - e) % 2:
            continue
        x = (d - e) // 2
        y = (d + e) // 2
        if 0 < x < y and x**3 + y**3 == n:
            reps.append((d, x, y))
    return reps


def tau(n: int) -> int:
    return len(divisors(n))


def main() -> None:
    print("== positive distinct two-cube support / indexed-fibre probe ==")

    kappa = gamma(4.0 / 3.0) ** 2 / (2.0 * gamma(5.0 / 3.0))
    residue = (2.0 / 3.0) * kappa
    print(f"area_constant={kappa:.15f}")
    print(f"critical_residue={(residue):.15f}")

    lattice_rows: list[tuple[int, int, str, str]] = []
    for exponent in range(3, 16):
        cutoff = 10**exponent
        count = indexed_count(cutoff)
        scale = cutoff ** (2.0 / 3.0)
        boundary = cutoff ** (1.0 / 3.0)
        ratio = count / scale
        error = (count - kappa * scale) / boundary
        lattice_rows.append((exponent, count, f"{ratio:.12f}", f"{error:.9f}"))
        print(
            f"X=1e{exponent:02d} indexed={count} "
            f"ratio={ratio:.12f} boundary_error={error:.9f}"
        )

    cutoff = 2_000_000
    fibre = pair_fibre(cutoff)
    indexed = sum(fibre.values())
    support = len(fibre)
    collision_tax = sum(v - 1 for v in fibre.values())
    require(indexed == indexed_count(cutoff), "pair/count paths disagree")
    require(indexed == support + collision_tax, "collision-tax identity")

    divisor_rep_total = 0
    max_mult = 0
    for n, multiplicity in fibre.items():
        good = positive_distinct_good_divisors(n)
        require(len(good) == multiplicity, f"good-divisor mismatch at {n}")
        require(multiplicity <= tau(n), f"divisor multiplicity bound at {n}")
        divisor_rep_total += len(good)
        max_mult = max(max_mult, multiplicity)
    require(divisor_rep_total == indexed, "divisor and pair fibres disagree")

    collisions = sorted((n, m) for n, m in fibre.items() if m >= 2)
    require(collisions and collisions[0][0] == 1729, "1729 positive control")
    require(dict(collisions).get(4104) == 2, "4104 collision control")
    require(2**3 + 2**3 not in fibre, "diagonal hostile leaked into distinct fibre")
    require(1**3 + 12**3 in fibre, "ordinary positive control missing")

    multi_mass = fsum(m / n for n, m in fibre.items())
    support_mass = fsum(1.0 / n for n in fibre)
    collision_mass = fsum((m - 1) / n for n, m in fibre.items())
    require(abs(multi_mass - support_mass - collision_mass) < 2e-14, "weighted tax")

    print(
        f"bounded_universe=1..{cutoff} indexed={indexed} support={support} "
        f"collision_tax={collision_tax} collisions={len(collisions)} max_mult={max_mult}"
    )
    print(f"first_collisions={collisions[:12]}")
    print(
        f"harmonic_prefix_multiset={multi_mass:.15f} "
        f"support={support_mass:.15f} collision={collision_mass:.15f}"
    )

    hostile_n = 1728
    require(not positive_distinct_good_divisors(hostile_n), "1728 hostile")
    good_1729 = positive_distinct_good_divisors(1729)
    require(good_1729 == [(13, 1, 12), (19, 9, 10)], "1729 divisor ledger")
    print(f"good_divisors_1729={good_1729}")

    semantic = repr(
        {
            "lattice": lattice_rows,
            "bounded": (cutoff, indexed, support, collision_tax, len(collisions), max_mult),
            "first_collisions": collisions[:12],
            "good_1729": good_1729,
        }
    ).encode()
    print(f"semantic_sha256={sha256(semantic).hexdigest()}")
    print("ALL GATES PASS")


if __name__ == "__main__":
    main()
