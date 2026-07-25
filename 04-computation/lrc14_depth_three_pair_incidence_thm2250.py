#!/usr/bin/env python3
"""Exact Lagrange-Bellman certificates for non-all-equal (3,4,5) cores."""

from fractions import Fraction
from functools import lru_cache
from itertools import combinations, product


Q = Fraction
P = 13
TARGET = Q(961, 6930)
PROFILE = (3, 4, 5)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def invert_square(matrix):
    n = len(matrix)
    augmented = [
        list(row) + [Q(int(i == j)) for j in range(n)]
        for i, row in enumerate(matrix)
    ]
    for column in range(n):
        pivot = next(
            (row for row in range(column, n) if augmented[row][column]),
            None,
        )
        if pivot is None:
            return None
        augmented[column], augmented[pivot] = (
            augmented[pivot],
            augmented[column],
        )
        scale = augmented[column][column]
        augmented[column] = [value / scale for value in augmented[column]]
        for row in range(n):
            if row == column:
                continue
            scale = augmented[row][column]
            if scale:
                augmented[row] = [
                    augmented[row][j] - scale * augmented[column][j]
                    for j in range(2 * n)
                ]
    return tuple(tuple(augmented[i][n:]) for i in range(n))


def matrix_vector(matrix, vector):
    return tuple(
        sum((matrix[i][j] * vector[j] for j in range(len(vector))), Q())
        for i in range(len(matrix))
    )


def dot(left, right):
    return sum((x * y for x, y in zip(left, right)), Q())


class CouplingOracle:
    def __init__(self, coordinate_count):
        self.bits = tuple(product((0, 1), repeat=coordinate_count))
        self.rank = coordinate_count + 1
        self.columns = tuple(
            (Q(1),) + tuple(Q(bit) for bit in vector)
            for vector in self.bits
        )
        inverses = {}
        for basis in combinations(range(len(self.bits)), self.rank):
            matrix = tuple(
                tuple(
                    self.columns[basis[column]][row]
                    for column in range(self.rank)
                )
                for row in range(self.rank)
            )
            inverse = invert_square(matrix)
            if inverse is not None:
                inverses[basis] = inverse
        self.inverses = inverses
        self.local_calls = 0
        self.dual_checks = 0

    def rhs(self, next_key):
        if next_key == len(self.bits):
            return (
                Q(1),
                Q(5, 7),
                *((Q(1, 7) for _ in range(len(self.bits[0]) - 1))),
            )
        next_bits = self.bits[next_key]
        return (
            Q(1),
            Q(10 - next_bits[0], P),
            *((Q(2 - bit, P) for bit in next_bits[1:])),
        )

    @lru_cache(maxsize=None)
    def vertices(self, next_key):
        rhs = self.rhs(next_key)
        by_distribution = {}
        for basis, inverse in self.inverses.items():
            weights = matrix_vector(inverse, rhs)
            if any(weight < 0 for weight in weights):
                continue
            distribution = [Q() for _ in self.bits]
            for key, weight in zip(basis, weights):
                distribution[key] = weight
            distribution = tuple(distribution)
            by_distribution.setdefault(distribution, []).append(basis)
        require(by_distribution, f"empty coupling polytope at {next_key}")
        return tuple(
            (
                distribution,
                tuple(bases),
                tuple(
                    (index, weight)
                    for index, weight in enumerate(distribution)
                    if weight
                ),
            )
            for distribution, bases in by_distribution.items()
        )

    @lru_cache(maxsize=None)
    def maximize_cached(self, values, next_key):
        rhs = self.rhs(next_key)
        optimum = None
        optimizers = []
        for distribution, bases, support in self.vertices(next_key):
            value = sum(
                (weight * values[index] for index, weight in support), Q()
            )
            if optimum is None or value > optimum:
                optimum = value
                optimizers = [(distribution, bases)]
            elif value == optimum:
                optimizers.append((distribution, bases))
        require(optimum is not None, "no primal optimum")
        for distribution, bases in optimizers:
            for basis in bases:
                inverse = self.inverses[basis]
                basic_values = tuple(values[index] for index in basis)
                dual = tuple(
                    sum(
                        (
                            inverse[row][column] * basic_values[row]
                            for row in range(self.rank)
                        ),
                        Q(),
                    )
                    for column in range(self.rank)
                )
                if any(
                    dot(column, dual) < value
                    for column, value in zip(self.columns, values)
                ):
                    continue
                require(dot(rhs, dual) == optimum, "dual value mismatch")
                require(
                    dot(distribution, values) == optimum,
                    "primal value mismatch",
                )
                self.dual_checks += 1
                return optimum
        raise RuntimeError("no dual-feasible optimal basis")

    def maximize(self, values, next_key):
        self.local_calls += 1
        return self.maximize_cached(tuple(values), next_key)


PARTITIONS = {
    "distinct": (0, 1, 2),
    "eq01": (0, 0, 1),
    "eq02": (0, 1, 0),
    "eq12": (0, 1, 1),
    "all_equal": (0, 0, 0),
}


def exact_bound(name, penalties):
    group_map = PARTITIONS[name]
    group_count = max(group_map) + 1
    oracle = CouplingOracle(1 + group_count)
    all_mask = 3

    def original_bits(hidden_bits):
        return tuple(hidden_bits[1 + group_map[j]] for j in range(3))

    def update(time, bits, pmask, nmask):
        for blocker, valuation in enumerate(PROFILE):
            if not bits[blocker]:
                continue
            if time == valuation:
                pmask |= 1
            if time == valuation - 2:
                pmask |= 2
            if time == valuation - 1:
                nmask |= 1
            if time == valuation - 3:
                nmask |= 2
        return pmask, nmask

    unique_pairs = tuple(combinations(range(group_count), 2))
    require(
        all(pair in unique_pairs for time, pair in penalties),
        "penalty on nonexistent/equal group pair",
    )

    def penalty(time, hidden_bits):
        return sum(
            coefficient * hidden_bits[1 + left] * hidden_bits[1 + right]
            for (penalty_time, (left, right)), coefficient in penalties.items()
            if penalty_time == time
        )

    @lru_cache(maxsize=None)
    def recurse(time, pmask, nmask, next_key):
        if time < 0:
            hidden_bits = oracle.bits[next_key]
            if not hidden_bits[0] or pmask != all_mask:
                return Q()
            return Q(bool(nmask & 1)) + Q(not bool(nmask & 2), 169)
        values = []
        for current_key, hidden_bits in enumerate(oracle.bits):
            bits = original_bits(hidden_bits)
            next_pmask, next_nmask = update(
                time, bits, pmask, nmask
            )
            values.append(
                recurse(
                    time - 1,
                    next_pmask,
                    next_nmask,
                    current_key,
                )
                - penalty(time, hidden_bits)
            )
        return oracle.maximize(values, next_key)

    terminal_values = tuple(
        recurse(5, 0, 0, terminal_key)
        for terminal_key in range(len(oracle.bits))
    )
    relaxed = oracle.maximize(terminal_values, len(oracle.bits))
    cap_charge = sum(penalties.values(), Q(0)) / 14
    return (
        relaxed + cap_charge,
        relaxed,
        cap_charge,
        oracle.local_calls,
        oracle.dual_checks,
        len(oracle.inverses),
        tuple(
            len(oracle.vertices(key))
            for key in range(len(oracle.bits) + 1)
        ),
    )


# Coarse rational penalties, chosen near the float occupation-LP dual.
CANDIDATES = {
    "distinct": {
        (5, (1, 2)): Q(1, 128),
        (4, (0, 1)): Q(1, 256),
        (4, (1, 2)): Q(1, 14),
        (3, (0, 1)): Q(7, 64),
        (3, (1, 2)): Q(1, 2),
    },
    "eq01": {
        (5, (0, 1)): Q(1, 128),
        (4, (0, 1)): Q(1, 8),
        (3, (0, 1)): Q(1, 2),
    },
    "eq02": {
        (5, (0, 1)): Q(1, 128),
        (4, (0, 1)): Q(3, 32),
        (3, (0, 1)): Q(2, 3),
    },
    "eq12": {
        (4, (0, 1)): Q(1, 256),
        (3, (0, 1)): Q(2, 3),
    },
}


if __name__ == "__main__":
    print("target", TARGET, float(TARGET))
    for name in PARTITIONS:
        uncontrolled = exact_bound(name, {})
        print(
            name,
            "uncontrolled",
            uncontrolled[0],
            float(uncontrolled[0]),
            "FAIL" if uncontrolled[0] >= TARGET else "PASS",
        )
        if name == "all_equal":
            require(
                uncontrolled[0] == Q(52940, 371293),
                "all-equal hostile control drift",
            )
            require(
                uncontrolled[0] > TARGET,
                "all-equal hostile control unexpectedly passed",
            )
            continue
        penalties = CANDIDATES[name]
        result = exact_bound(name, penalties)
        value = result[0]
        require(value < TARGET, f"{name} exact certificate failed")
        print(
            name,
            "penalties",
            tuple(sorted(penalties.items())),
            "certified",
            value,
            float(value),
            "PASS",
            "gap",
            TARGET - value,
            "relaxed",
            result[1],
            "cap",
            result[2],
            "calls",
            result[3],
            "duals",
            result[4],
            "bases",
            result[5],
            "vertex_census",
            result[6],
        )
