#!/usr/bin/env python3
"""Independent exact probes for THM-2100's recursive guard-ratio box.

All geometry is rebuilt from rational arrangement endpoints. No repository
helper is imported, and every check remains active under ``python -O``.
"""
from fractions import Fraction as F
from math import gcd
import random


def dist0(x: F) -> F:
    x %= 1
    return min(x, 1 - x)


def endpoints(speed: int, denominator: int):
    """Boundaries ``||speed*t|| = 1/denominator`` on the unit circle."""
    return {
        (F(n, speed) + sign * F(1, denominator * speed)) % 1
        for n in range(speed)
        for sign in (-1, 1)
    }


def carrier_atoms(h: int, known: tuple[int, ...]):
    ep = {F(0), F(1)} | endpoints(h, 7)
    for q in known:
        ep |= endpoints(q, 14)
    xs = sorted(ep)
    atoms = []
    for a, b in zip(xs, xs[1:]):
        if a == b:
            continue
        t = (a + b) / 2
        live = dist0(h * t) >= F(1, 7) and all(
            dist0(q * t) >= F(1, 14) for q in known
        )
        atoms.append((a, b, live))
    measure = sum((b - a for a, b, live in atoms if live), F(0))
    flags = [live for _, _, live in atoms]
    if not any(flags):
        components = 0
    elif all(flags):
        components = 1
    else:
        components = sum(flags[i] and not flags[i - 1] for i in range(len(flags)))
    return atoms, measure, components


def intersect_danger(atoms, q: int) -> F:
    qep = sorted(endpoints(q, 14) | {F(0), F(1)})
    total = F(0)
    for a, b, live in atoms:
        if not live:
            continue
        cuts = [a] + [x for x in qep if a < x < b] + [b]
        for u, v in zip(cuts, cuts[1:]):
            if dist0(q * (u + v) / 2) < F(1, 14):
                total += v - u
    return total


def overlap(h: int, q: int) -> F:
    ep = sorted({F(0), F(1)} | endpoints(h, 7) | endpoints(q, 14))
    answer = F(0)
    for a, b in zip(ep, ep[1:]):
        t = (a + b) / 2
        if dist0(h * t) < F(1, 7) and dist0(q * t) < F(1, 14):
            answer += b - a
    return answer


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def main():
    overlap_checks = 0
    for h in range(1, 32, 2):
        for q in range(1, 81):
            value = overlap(h, q)
            if q == 6 * h:
                require(value == F(1, 42), (h, q, value, "6h"))
            elif 11 * q == h:
                require(value == F(2, 77), (h, q, value, "h/11"))
            else:
                require(value >= F(1, 35), (h, q, value, "generic"))
            overlap_checks += 1

    cases = []
    for h in range(1, 16, 2):
        pool = [q for q in range(1, 22) if q != h]
        for k in range(1, min(6, len(pool)) + 1):
            cases.append((h, tuple(pool[:k])))
            cases.append((h, tuple(pool[-k:])))
    rng = random.Random(2097)
    for _ in range(220):
        h = rng.randrange(1, 30, 2)
        k = rng.randrange(1, 7)
        pool = [q for q in range(1, 45) if q != h]
        cases.append((h, tuple(sorted(rng.sample(pool, k)))))

    carrier_checks = 0
    discrepancy_checks = 0
    min_slack = None
    for h, known in cases:
        atoms, measure, components = carrier_atoms(h, known)
        k = len(known)
        if k == 1:
            floor = F(5 - k, 7) + F(1, 42)
        else:
            floor = F(5 - k, 7) + F(1, 42) + F(2, 77) + (k - 2) * F(1, 35)
        require(measure >= floor, (h, known, measure, floor, "mass floor"))
        require(1 <= components <= h + sum(known), (h, known, components, "components"))
        min_slack = measure - floor if min_slack is None else min(min_slack, measure - floor)
        carrier_checks += 1
        for q in (1, 2, 3, 5, 7, 11, 17, 31, 47, 73):
            got = intersect_danger(atoms, q)
            rhs = measure / F(7) + F(6 * components, 49 * q)
            require(got <= rhs, (h, known, q, got, rhs, "discrepancy"))
            discrepancy_checks += 1

    divisor = 3**4 * 5**2 * 11 * 17 * 23 * 29
    guard_bound = 57 * divisor
    selected_sum = guard_bound
    expected = [
        248777478576,
        1243418355401,
        4774568778091,
        16202606721059,
        56928264210988,
        534912686676376,
    ]
    rows = []
    for k in range(1, 7):
        remaining = 7 - k
        if k == 1:
            delta = F(5 - k, 7) + F(1, 42)
        else:
            delta = F(5 - k, 7) + F(1, 42) + F(2, 77) + (k - 2) * F(1, 35)
        coefficient = F(6 * remaining, 7 * k) / delta
        component_bound = guard_bound + selected_sum
        next_bound = coefficient.numerator * component_bound // coefficient.denominator
        require(next_bound == expected[k - 1], (k, next_bound, expected[k - 1]))
        rows.append((k, delta, coefficient, component_bound, next_bound))
        selected_sum += next_bound

    terminal_bound = max([guard_bound] + expected)
    full_bound = 128 * terminal_bound // 3
    require(terminal_bound == 534912686676376, "terminal box")
    require(full_bound == 22822941298192042, "full box")

    print(f"overlap_checks={overlap_checks}")
    print(f"carrier_checks={carrier_checks}")
    print(f"discrepancy_checks={discrepancy_checks}")
    print(f"minimum_mass_floor_slack={min_slack}")
    for k, delta, coefficient, component_bound, next_bound in rows:
        print(
            f"k={k} delta={delta} coefficient={coefficient} "
            f"components_bound={component_bound} next_speed_bound={next_bound}"
        )
    print(f"terminal_max_bound={terminal_bound}")
    print(f"full_row_max_bound={full_bound}")
    print("PASS")


if __name__ == "__main__":
    main()
