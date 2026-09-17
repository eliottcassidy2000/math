"""Exact controls for the distinct-summand/Collatz graph connection.

No trajectory census is used to assert Collatz convergence.  All universes
and finite cutoffs are printed. Run from repository root with Python 3.
"""

from fractions import Fraction
from itertools import combinations
import json


def valuation(n, p):
    assert n > 0 and p >= 2
    a = 0
    while n % p == 0:
        n //= p
        a += 1
    return a


def odd_core(n):
    return n // 2 ** valuation(n, 2)


def shortcut(n, p=3, sign=1):
    return (p * n + sign) // 2 if n % 2 else n // 2


def odd_step(n, p=3, sign=1):
    return odd_core(p * n + sign)


def closure_step(values):
    return values | {a + b for a, b in combinations(sorted(values), 2)}


def check():
    bound = 256
    # Build the strict shadow from genuine unordered parent pairs, not
    # from the formula which is being tested.
    strict_arcs = set()
    for z in range(1, bound + 1):
        for x in range(1, z):
            y = z - x
            if x < y:
                strict_arcs.update(((x, z), (y, z)))
    ambient = {(x, z) for z in range(2, bound + 1) for x in range(1, z)}
    missing = ambient - strict_arcs
    assert missing == {(x, 2 * x) for x in range(1, bound // 2 + 1)}
    assert (1, 2) not in strict_arcs and (1, 3) in strict_arcs

    multiplicative_weak = set()
    multiplicative_strict = set()
    for x, z in ambient:
        if z % x == 0:
            multiplicative_weak.add((x, z))
            if z // x != x:
                multiplicative_strict.add((x, z))
    multiplicative_missing = multiplicative_weak - multiplicative_strict
    assert multiplicative_missing == {
        (x, x * x) for x in range(2, bound + 1) if x * x <= bound
    }

    # Independent synchronous closure proves finite controls for the
    # all-height interval argument recorded in the accompanying note.
    depths = {1: 0, 2: 0}
    live = {1, 2}
    frontiers = []
    for depth in range(1, 9):
        enlarged = closure_step(live)
        for n in enlarged - live:
            depths[n] = depth
        live = enlarged
        assert live == set(range(1, 2 ** depth + 2))
        frontiers.append(max(live))
    assert depths[34] == depths[35] == depths[55] == 6
    assert depths[27] == depths[28] == 5

    # Countermodels have the same colour rule as shortcut Collatz:
    # even nodes halve, and odd nodes >1 use strict summand arcs.
    def divergent(n):
        return 2 if n == 1 else n + 2 if n % 2 else n // 2

    def extra_cycle(n):
        return 6 if n == 5 else divergent(n)

    for fn in (divergent, extra_cycle):
        assert fn(1) == 2 and fn(2) == 1
        for n in range(3, 201, 2):
            z = fn(n)
            assert z > n and z - n != n
    assert [extra_cycle(n) for n in (3, 5, 6)] == [5, 6, 3]
    assert [divergent(n) for n in range(3, 20, 2)] == list(range(5, 22, 2))

    # Shortcut and odd-only versions of the explicit signed/5n+1 cycles.
    minus_cycles = ((1,), (5, 7), (17, 25, 37, 55, 41, 61, 91))
    for cycle in minus_cycles:
        assert tuple(odd_step(n, sign=-1) for n in cycle) == cycle[1:] + cycle[:1]
    five_cycle = (13, 33, 83, 208, 104, 52, 26)
    assert tuple(shortcut(n, p=5) for n in five_cycle) == five_cycle[1:] + five_cycle[:1]

    # Independent braid audit using iteration (not the closed formula).
    braid_records = []
    for p, q in ((3, 4), (5, 16), (7, 8), (17, 256)):
        c = (q - 1) // p
        assert valuation(q - 1, p) == 1
        for s in range(1, 5):
            modulus = 2 * p ** s
            visited = set()
            n = 1
            while n not in visited:
                visited.add(n)
                n = (q * n + c) % modulus
            assert n == 1
            assert visited == set(range(1, modulus, 2))
            braid_records.append([p, s, len(visited)])
        for n in range(1, 60, 2):
            a, b, z = n, ((p - 2) * n + 1) // 2, (p * n + 1) // 2
            transformed = (q * a + c, q * b - c, q * z)
            nn = q * n + c
            assert transformed == (nn, ((p - 2) * nn + 1) // 2, (p * nn + 1) // 2)
            assert transformed[0] + transformed[1] == transformed[2]
            assert odd_core(z) == odd_core(transformed[2])

    # Height sidecar exactly recovers the full shortcut graph.
    for n in range(1, 1001):
        u, k = odd_core(n), valuation(n, 2)
        if k:
            reconstructed = 2 ** (k - 1) * u
        else:
            a = valuation(3 * u + 1, 2)
            reconstructed = 2 ** (a - 1) * odd_core(3 * u + 1)
        assert reconstructed == shortcut(n)

    # If the half-step affine map is applied even at even inputs, each
    # further step adds a denominator factor two and never returns to Z.
    rational_records = []
    for n in (2, 4, 6, 8):
        x = Fraction(n)
        row = []
        for t in range(1, 7):
            x = (3 * x + 1) / 2
            assert x.denominator == 2 ** t
            row.append(str(x))
        rational_records.append([n, row])

    return {
        "status": "all exact controls passed; no global Collatz claim",
        "strict_summand_universe": [1, bound],
        "ascending_arcs": len(ambient),
        "strict_summand_arcs": len(strict_arcs),
        "missing_doublings": len(missing),
        "missing_multiplicative_squarings": len(multiplicative_missing),
        "seed_1_2_synchronous_frontiers_depths_1_to_8": frontiers,
        "depths": {n: depths[n] for n in (27, 28, 34, 35, 55)},
        "minus_odd_cycles": minus_cycles,
        "five_shortcut_cycle": five_cycle,
        "braid_records_p_s_period": braid_records,
        "shortcut_reconstruction_checked_through": 1000,
        "half_affine_even_input_examples": rational_records,
    }


if __name__ == "__main__":
    print(json.dumps(check(), indent=2))
