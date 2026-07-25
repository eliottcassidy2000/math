#!/usr/bin/env python3
"""Exact half-time hidden-state Bellman referee for THM-2227.

The state retains the three current unit-core bits at every multiplication by
13.  Every transition optimizes over all joint current-bit couplings with the
three exact one-core conditional marginals.  Exact primal and dual vertices
are independently enumerated.
"""

from fractions import Fraction
from functools import lru_cache
from itertools import combinations, product


CLAUSES = 3
ALL = 7
BIT_VECTORS = tuple(product((0, 1), repeat=3))


def solve_square(matrix, rhs):
    """Solve a square rational system, returning None when singular."""
    n = len(rhs)
    aug = [list(row) + [rhs[i]] for i, row in enumerate(matrix)]
    for col in range(n):
        pivot = next((r for r in range(col, n) if aug[r][col]), None)
        if pivot is None:
            return None
        aug[col], aug[pivot] = aug[pivot], aug[col]
        scale = aug[col][col]
        aug[col] = [x / scale for x in aug[col]]
        for r in range(n):
            if r == col:
                continue
            scale = aug[r][col]
            if scale:
                aug[r] = [
                    aug[r][j] - scale * aug[col][j] for j in range(n + 1)
                ]
    return tuple(aug[i][-1] for i in range(n))


@lru_cache(maxsize=None)
def maximize_joint_bits(values, marginals):
    """Exact three-Bernoulli coupling LP with primal/dual parity."""
    columns = tuple(
        (Fraction(1),) + tuple(Fraction(bit) for bit in bits)
        for bits in BIT_VECTORS
    )
    rhs = (Fraction(1),) + tuple(marginals)
    rank = 4
    primal = None
    for basis in combinations(range(8), rank):
        matrix = tuple(
            tuple(columns[basis[col]][row] for col in range(rank))
            for row in range(rank)
        )
        solution = solve_square(matrix, rhs)
        if solution is None or any(x < 0 for x in solution):
            continue
        candidate = sum(values[basis[i]] * solution[i] for i in range(rank))
        if primal is None or candidate > primal:
            primal = candidate
    if primal is None:
        raise RuntimeError("joint-bit primal has no basic feasible solution")

    dual = None
    for active in combinations(range(8), rank):
        matrix = tuple(columns[i] for i in active)
        solution = solve_square(matrix, tuple(values[i] for i in active))
        if solution is None:
            continue
        if any(
            sum(columns[z][i] * solution[i] for i in range(rank)) < values[z]
            for z in range(8)
        ):
            continue
        candidate = sum(rhs[i] * solution[i] for i in range(rank))
        if dual is None or candidate < dual:
            dual = candidate
    if dual is None:
        raise RuntimeError("joint-bit dual has no enumerated vertex")
    if primal != dual:
        raise RuntimeError(f"joint-bit primal/dual drift: {primal} != {dual}")
    return primal


def eliminate_half_time(value, n, valuations):
    """Eliminate one multiplication-by-13 root digit."""
    active = []
    for j, valuation in enumerate(valuations):
        delta = n - valuation
        if delta >= 0 and delta % 2 == 0 and delta // 2 < CLAUSES:
            active.append((j, delta // 2))
    out = [[Fraction(0) for _ in range(8)] for _ in range(8)]
    for later in range(8):
        for next_mask, next_bits in enumerate(BIT_VECTORS):
            marginals = tuple(Fraction(2 - bit, 13) for bit in next_bits)
            costs = []
            for current_mask, current_bits in enumerate(BIT_VECTORS):
                newly = 0
                for j, clause in active:
                    if current_bits[j]:
                        newly |= 1 << clause
                costs.append(value[later | newly][current_mask])
            out[later][next_mask] = maximize_joint_bits(
                tuple(costs), marginals
            )
    return tuple(tuple(row) for row in out)


def hidden_bound(valuations):
    """Robust hidden-state upper bound for one relative valuation triple."""
    value = tuple(
        tuple(Fraction(int(mask == ALL)) for _ in range(8))
        for mask in range(8)
    )
    for n in range(max(valuations) + 5):
        value = eliminate_half_time(value, n, valuations)
    return maximize_joint_bits(value[0], (Fraction(1, 7),) * 3)


def power_tower_actual(valuations):
    """Exact same-core Markov control for a 13-power coefficient tower."""
    indices = tuple(
        tuple(2 * r + valuations[j] for j in range(3))
        for r in range(3)
    )
    length = max(max(row) for row in indices) + 1
    transition = (
        (Fraction(11, 13), Fraction(2, 13)),
        (Fraction(12, 13), Fraction(1, 13)),
    )
    stationary = (Fraction(6, 7), Fraction(1, 7))
    answer = Fraction(0)
    for word in product((0, 1), repeat=length):
        if not all(any(word[t] for t in row) for row in indices):
            continue
        mass = stationary[word[0]]
        for left, right in zip(word, word[1:]):
            mass *= transition[left][right]
        answer += mass
    return answer


def main():
    target = Fraction(961, 6930)
    expected_bounds = {
        (0, 2, 4): Fraction(1_086_371_907, 5_710_115_047),
        (0, 2, 5): Fraction(541_866_697, 5_710_115_047),
        (0, 3, 4): Fraction(464_433_205, 5_710_115_047),
        (0, 3, 5): Fraction(1_005_912_252, 10_604_499_373),
    }
    expected_actual = {
        (0, 2, 4): Fraction(916_159, 4_826_809),
        (0, 2, 5): Fraction(895_649_983, 10_604_499_373),
        (0, 3, 4): Fraction(57_121_111, 815_730_721),
        (0, 3, 5): Fraction(895_649_983, 10_604_499_373),
    }
    expected_margins = {
        (0, 2, 4): Fraction(291_590_965_049, 5_653_013_896_530),
        (0, 2, 5): Fraction(247_469_192_851, 5_653_013_896_530),
        (0, 3, 4): Fraction(324_128_349_931, 5_653_013_896_530),
        (0, 3, 5): Fraction(3_219_951_991_093, 73_489_180_654_890),
    }
    bounds = {}
    for valuations in expected_bounds:
        bound = hidden_bound(valuations)
        actual = power_tower_actual(valuations)
        if bound != expected_bounds[valuations]:
            raise RuntimeError(f"hidden bound drift for {valuations}: {bound}")
        if actual != expected_actual[valuations]:
            raise RuntimeError(f"power-tower drift for {valuations}: {actual}")
        margin = bound - target if bound >= target else target - bound
        if margin != expected_margins[valuations]:
            raise RuntimeError(f"target margin drift for {valuations}")
        bounds[valuations] = bound
        print(
            f"valuations={valuations} bound={bound} "
            f"decimal={float(bound):.15f} "
            f"tower={actual} tower_decimal={float(actual):.15f} "
            f"verdict={'PASS' if bound < target else 'FAIL'}"
        )

    failing = tuple(v for v, bound in bounds.items() if bound >= target)
    if failing != ((0, 2, 4),):
        raise RuntimeError(f"hidden-parity failure set drift: {failing}")

    eight_profiles = (
        (4, 6, 8),
        (4, 6, 9),
        (4, 7, 8),
        (4, 7, 9),
        (5, 7, 9),
        (5, 7, 10),
        (5, 8, 9),
        (5, 8, 10),
    )
    closed = tuple(
        profile
        for profile in eight_profiles
        if bounds[
            (0, profile[1] - profile[0], profile[2] - profile[0])
        ]
        < target
    )
    remaining = tuple(profile for profile in eight_profiles if profile not in closed)
    expected_closed = (
        (4, 6, 9),
        (4, 7, 8),
        (4, 7, 9),
        (5, 7, 10),
        (5, 8, 9),
        (5, 8, 10),
    )
    if closed != expected_closed:
        raise RuntimeError(f"closed profile drift: {closed}")
    if remaining != ((4, 6, 8), (5, 7, 9)):
        raise RuntimeError(f"remaining profile drift: {remaining}")
    if 458 - len(closed) != 452:
        raise RuntimeError("current scalar ledger drift")

    print(f"target={target}")
    print(f"closed_profiles={closed}")
    print(f"remaining_profiles={remaining}")
    print("standalone_thm2219_thm2226_scalar_ledger=452")
    print("exact_primal_dual_vertex_parity=PASS")
    print("status=THM2227_HALF_TIME_HIDDEN_PARITY_BELLMAN")


if __name__ == "__main__":
    main()
