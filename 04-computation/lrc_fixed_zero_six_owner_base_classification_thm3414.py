#!/usr/bin/env python3
"""Exact global finite-profile atlas for THM-3414.

This standard-library companion proves an all-modulus statement without a
modulus cutoff.  It starts from THM-3408's globally exhaustive exceptional
order set B and constructs exact rational probability measures on sheet-order
strata q/d.  A finite-companion lemma reduces every arbitrary extra owner to
a finite set depending only on d and B; all remaining owners have inherited
load at most 7/44.  A deterministic Fraction simplex discovers the measures,
and independent exact inequalities verify every certificate before it is
accepted.  The program also checks the five positive base covers and the
q=21, q=22, and q=102 hostile controls.
"""

from __future__ import annotations

import ast
from collections import defaultdict
from fractions import Fraction
from functools import cache
from hashlib import sha256
from itertools import combinations, combinations_with_replacement
from math import gcd, lcm
from pathlib import Path


BASES = (15, 16, 18, 20, 24)
B_EXCEPTIONS = (
    21, 22, 26, 28, 33, 35, 39, 42, 44, 46, 50, 52, 56, 57, 63, 65,
    70, 74, 76, 77, 78, 84, 91, 99, 102, 104, 114, 117, 130, 143,
    156, 186, 190,
)
H_CORRIDOR = (22, 33, 44, 46, 50, 102)
RESIDUAL_MAXIMAL_CLIQUES = (
    (21, 26, 39, 42, 46, 57, 74, 78, 91, 102, 114, 186),
    (21, 33, 39, 57, 63, 77, 91, 99, 117, 143),
    (26, 35, 46, 50, 65, 70, 74, 91, 130, 190),
    (26, 28, 46, 52, 56, 74, 76, 91, 104),
    (26, 46, 57, 74, 76, 91, 102, 114, 186),
    (21, 33, 39, 46, 57, 91, 102),
    (22, 46, 102),
    (44, 46, 102),
)
NONREPEATABLE = (22, 44, 84, 156)
BAD_ANCHORS = {
    2: (),
    3: (),
    4: ((46, 46), (46, 102), (102, 102)),
    5: ((46, 46), (46, 74), (46, 91), (46, 102), (102, 102)),
}
POSITIVE_WITNESSES = {
    15: (1, 2, 3, 4, 5, 7),
    16: (1, 3, 5, 7, 8),
    18: (1, 5, 6, 7, 9),
    20: (1, 3, 4, 7, 9, 10),
    24: (1, 5, 7, 8, 11, 12),
}
Q102_SHEET_WEIGHTS = {
    6: Fraction(1, 25),
    34: Fraction(4, 25),
    51: Fraction(4, 25),
    102: Fraction(16, 25),
}
Q102_OWNER_WEIGHTS = {
    2: Fraction(7, 150),
    3: Fraction(7, 150),
    17: Fraction(4, 25),
    102: Fraction(56, 75),
}
EXPECTED_SEMANTIC_DIGEST = "681e1b1f8066b54970ed791f14867df708e9808c8a66dce034002f0404e66b3f"


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def lf_hash(path):
    return sha256(path.read_bytes().replace(bytes((13, 10)), bytes((10,)))).hexdigest()


class ExactDigest:
    def __init__(self):
        self._hash = sha256()

    def add(self, value):
        self._hash.update(repr(value).encode("ascii"))
        self._hash.update(bytes((10,)))

    def hexdigest(self):
        return self._hash.hexdigest()


@cache
def prime_factorization(number):
    require(number >= 1, number)
    answer = []
    divisor = 2
    while divisor * divisor <= number:
        if number % divisor == 0:
            exponent = 0
            while number % divisor == 0:
                number //= divisor
                exponent += 1
            answer.append((divisor, exponent))
        divisor += 1
    if number > 1:
        answer.append((number, 1))
    return tuple(answer)


@cache
def prime_divisors(number):
    return tuple(prime for prime, _ in prime_factorization(number))


@cache
def phi(number):
    answer = number
    for prime in prime_divisors(number):
        answer = answer // prime * (prime - 1)
    return answer


@cache
def divisors(number):
    answer = [1]
    for prime, exponent in prime_factorization(number):
        old = tuple(answer)
        power = 1
        for _ in range(exponent):
            power *= prime
            answer.extend(divisor * power for divisor in old)
    return tuple(sorted(answer))


@cache
def coprime_prefix_count(number, bound):
    primes = prime_divisors(number)
    answer = 0
    for mask in range(1 << len(primes)):
        divisor = 1
        sign = 1
        for index, prime in enumerate(primes):
            if mask & (1 << index):
                divisor *= prime
                sign = -sign
        answer += sign * (bound // divisor)
    return answer


@cache
def alpha(order):
    if order == 1:
        return Fraction(1)
    cutoff = (order - 1) // 14
    return Fraction(2 * coprime_prefix_count(order, cutoff), phi(order))


def base_free(number):
    return not any(number % base == 0 for base in BASES)


class ExactLinearProgram:
    """Deterministic exact two-phase simplex.

    Maximize c*x subject to A*x<=b and x>=0.  The caller uses the solver only
    to discover a rational certificate and then checks every relevant
    inequality directly with Fraction arithmetic.
    """

    def __init__(self, matrix, bounds, objective):
        self.row_count = len(bounds)
        self.variable_count = len(objective)
        rows = self.row_count
        columns = self.variable_count
        require(rows >= 1 and columns >= 1, (rows, columns))
        require(all(len(row) == columns for row in matrix), "ragged matrix")
        self.basic = [columns + row for row in range(rows)]
        self.nonbasic = list(range(columns)) + [-1]
        self.tableau = [
            [Fraction(0) for _ in range(columns + 2)]
            for _ in range(rows + 2)
        ]
        for row in range(rows):
            for column in range(columns):
                self.tableau[row][column] = Fraction(matrix[row][column])
            self.tableau[row][columns] = Fraction(-1)
            self.tableau[row][columns + 1] = Fraction(bounds[row])
        for column in range(columns):
            self.tableau[rows][column] = -Fraction(objective[column])
        self.tableau[rows + 1][columns] = Fraction(1)

    def pivot(self, leaving_row, entering_column):
        rows = self.row_count
        columns = self.variable_count
        tableau = self.tableau
        pivot_value = tableau[leaving_row][entering_column]
        require(pivot_value != 0, ("zero pivot", leaving_row, entering_column))
        for row in range(rows + 2):
            if row == leaving_row:
                continue
            multiplier = tableau[row][entering_column]
            if multiplier == 0:
                continue
            for column in range(columns + 2):
                if column != entering_column:
                    tableau[row][column] -= (
                        tableau[leaving_row][column] * multiplier / pivot_value
                    )
        for column in range(columns + 2):
            if column != entering_column:
                tableau[leaving_row][column] /= pivot_value
        for row in range(rows + 2):
            if row != leaving_row:
                tableau[row][entering_column] /= -pivot_value
        tableau[leaving_row][entering_column] = Fraction(1, 1) / pivot_value
        self.basic[leaving_row], self.nonbasic[entering_column] = (
            self.nonbasic[entering_column],
            self.basic[leaving_row],
        )

    def simplex(self, phase):
        rows = self.row_count
        columns = self.variable_count
        objective_row = rows + 1 if phase == 1 else rows
        while True:
            eligible = tuple(
                column
                for column in range(columns + 1)
                if not (phase == 2 and self.nonbasic[column] == -1)
            )
            entering = min(
                eligible,
                key=lambda column: (
                    self.tableau[objective_row][column],
                    self.nonbasic[column],
                ),
            )
            if self.tableau[objective_row][entering] >= 0:
                return True
            possible = tuple(
                row for row in range(rows) if self.tableau[row][entering] > 0
            )
            if not possible:
                return False
            leaving = min(
                possible,
                key=lambda row: (
                    self.tableau[row][columns + 1]
                    / self.tableau[row][entering],
                    self.basic[row],
                ),
            )
            self.pivot(leaving, entering)

    def solve(self):
        rows = self.row_count
        columns = self.variable_count
        first_leaving = min(
            range(rows), key=lambda row: self.tableau[row][columns + 1]
        )
        if self.tableau[first_leaving][columns + 1] < 0:
            self.pivot(first_leaving, columns)
            require(self.simplex(1), "infeasible phase-I ray")
            require(self.tableau[rows + 1][columns + 1] == 0, "infeasible LP")
            for row in range(rows):
                if self.basic[row] == -1:
                    candidates = tuple(
                        column
                        for column in range(columns + 1)
                        if self.nonbasic[column] != -1
                        and self.tableau[row][column] != 0
                    )
                    if candidates:
                        entering = max(
                            candidates,
                            key=lambda column: (
                                abs(self.tableau[row][column]),
                                -self.nonbasic[column],
                            ),
                        )
                        self.pivot(row, entering)
                    break
        require(self.simplex(2), "unbounded LP")
        solution = [Fraction(0) for _ in range(columns)]
        for row, label in enumerate(self.basic):
            if label < columns:
                solution[label] = self.tableau[row][columns + 1]
        return self.tableau[rows][columns + 1], tuple(solution)


def residual_pairs():
    answer = set()
    for clique in RESIDUAL_MAXIMAL_CLIQUES:
        answer.update(combinations_with_replacement(clique, 2))
    for vertex in NONREPEATABLE:
        answer.discard((vertex, vertex))
    return frozenset(answer)


def adjacency_from_pairs(pairs):
    answer = {vertex: set() for vertex in B_EXCEPTIONS}
    for left, right in pairs:
        answer[left].add(right)
        answer[right].add(left)
    return {vertex: frozenset(neighbours) for vertex, neighbours in answer.items()}


def finite_candidates(modulus, support_divisors, ordinary_only):
    """All companions with an exceptional projected order on the support."""
    answer = set()
    for divisor in support_divisors:
        for projected in (1,) + B_EXCEPTIONS:
            for owner in divisors(divisor * projected):
                if owner < 2:
                    continue
                if not base_free(lcm(modulus, owner)):
                    continue
                if owner // gcd(owner, divisor) != projected:
                    continue
                if ordinary_only and owner in B_EXCEPTIONS:
                    continue
                answer.add(owner)
    return tuple(sorted(answer))


def profile_load(weights, support_divisors, owner):
    return sum(
        weight * alpha(owner // gcd(owner, divisor))
        for weight, divisor in zip(weights, support_divisors)
    )


def solve_pair_pruning_group(modulus, fixed_pairs, digest):
    support = tuple(divisor for divisor in divisors(modulus) if divisor < modulus)
    candidates = finite_candidates(modulus, support, ordinary_only=False)
    dimension = len(support)
    matrix = []
    bounds = []
    for owner in candidates:
        matrix.append(
            [alpha(owner // gcd(owner, divisor)) for divisor in support]
            + [Fraction(-1), Fraction(0)]
        )
        bounds.append(Fraction(0))
    matrix.append(
        [Fraction(7, 44)] * dimension + [Fraction(-1), Fraction(0)]
    )
    bounds.append(Fraction(0))
    for fixed in fixed_pairs:
        matrix.append(
            [
                sum(alpha(owner // gcd(owner, divisor)) for owner in fixed)
                for divisor in support
            ]
            + [Fraction(4), Fraction(-1)]
        )
        bounds.append(Fraction(0))
    matrix.append([Fraction(1)] * dimension + [Fraction(0), Fraction(0)])
    bounds.append(Fraction(1))
    matrix.append([Fraction(-1)] * dimension + [Fraction(0), Fraction(0)])
    bounds.append(Fraction(-1))
    objective_value, solution = ExactLinearProgram(
        matrix,
        bounds,
        [Fraction(0)] * (dimension + 1) + [Fraction(-1)],
    ).solve()
    weights = solution[:dimension]
    companion_bound = solution[dimension]
    certificate_value = solution[-1]
    require(all(value >= 0 for value in solution), ("negative pair solution", modulus))
    require(sum(weights) == 1, ("pair normalization", modulus, sum(weights)))
    require(objective_value == -certificate_value, (modulus, objective_value))
    candidate_loads = tuple(profile_load(weights, support, owner) for owner in candidates)
    require(
        max(candidate_loads + (Fraction(7, 44),)) <= companion_bound,
        ("pair companion bound", modulus),
    )
    fixed_loads = tuple(
        sum(profile_load(weights, support, owner) for owner in fixed)
        + 4 * companion_bound
        for fixed in fixed_pairs
    )
    require(max(fixed_loads) == certificate_value, ("pair value", modulus))
    require(certificate_value < 1, ("pair certificate failure", modulus, certificate_value))
    digest.add((
        "pair",
        modulus,
        tuple(fixed_pairs),
        tuple(zip(support, weights)),
        companion_bound,
        certificate_value,
        candidates,
    ))
    return certificate_value, sum(weight != 0 for weight in weights), len(candidates)


def solve_anchor_group(modulus, fixed_pairs, b_count, adjacency, digest):
    support = tuple(divisor for divisor in divisors(modulus) if divisor < modulus)
    candidates = finite_candidates(modulus, support, ordinary_only=True)
    dimension = len(support)
    pair_count = len(fixed_pairs)
    matrix = []
    bounds = []
    for owner in candidates:
        matrix.append(
            [alpha(owner // gcd(owner, divisor)) for divisor in support]
            + [Fraction(-1)]
            + [Fraction(0)] * pair_count
            + [Fraction(0)]
        )
        bounds.append(Fraction(0))
    matrix.append(
        [Fraction(7, 44)] * dimension
        + [Fraction(-1)]
        + [Fraction(0)] * pair_count
        + [Fraction(0)]
    )
    bounds.append(Fraction(0))
    neighbour_sets = []
    for index, fixed in enumerate(fixed_pairs):
        neighbours = tuple(sorted(
            owner
            for owner in adjacency[fixed[0]] & adjacency[fixed[1]]
            if base_free(lcm(modulus, owner))
        ))
        neighbour_sets.append(neighbours)
        for owner in neighbours:
            row = (
                [alpha(owner // gcd(owner, divisor)) for divisor in support]
                + [Fraction(0)]
                + [Fraction(0)] * pair_count
                + [Fraction(0)]
            )
            row[dimension + 1 + index] = Fraction(-1)
            matrix.append(row)
            bounds.append(Fraction(0))
        row = (
            [
                sum(alpha(owner // gcd(owner, divisor)) for owner in fixed)
                for divisor in support
            ]
            + [Fraction(6 - b_count)]
            + [Fraction(0)] * pair_count
            + [Fraction(-1)]
        )
        row[dimension + 1 + index] = Fraction(b_count - 2)
        matrix.append(row)
        bounds.append(Fraction(0))
    matrix.append(
        [Fraction(1)] * dimension + [Fraction(0)] * (pair_count + 2)
    )
    bounds.append(Fraction(1))
    matrix.append(
        [Fraction(-1)] * dimension + [Fraction(0)] * (pair_count + 2)
    )
    bounds.append(Fraction(-1))
    objective_value, solution = ExactLinearProgram(
        matrix,
        bounds,
        [Fraction(0)] * (dimension + pair_count + 1) + [Fraction(-1)],
    ).solve()
    weights = solution[:dimension]
    ordinary_bound = solution[dimension]
    neighbour_bounds = solution[dimension + 1:dimension + 1 + pair_count]
    certificate_value = solution[-1]
    require(all(value >= 0 for value in solution), ("negative anchor solution", b_count, modulus))
    require(sum(weights) == 1, ("anchor normalization", b_count, modulus))
    require(objective_value == -certificate_value, (b_count, modulus, objective_value))
    ordinary_loads = tuple(profile_load(weights, support, owner) for owner in candidates)
    require(
        max(ordinary_loads + (Fraction(7, 44),)) <= ordinary_bound,
        ("ordinary bound", b_count, modulus),
    )
    fixed_totals = []
    for fixed, neighbours, neighbour_bound in zip(
        fixed_pairs, neighbour_sets, neighbour_bounds
    ):
        loads = tuple(profile_load(weights, support, owner) for owner in neighbours)
        require(
            max(loads, default=Fraction(0)) <= neighbour_bound,
            ("neighbour bound", b_count, modulus, fixed),
        )
        fixed_totals.append(
            sum(profile_load(weights, support, owner) for owner in fixed)
            + (b_count - 2) * neighbour_bound
            + (6 - b_count) * ordinary_bound
        )
    require(max(fixed_totals) == certificate_value, ("anchor value", b_count, modulus))
    require(certificate_value < 1, ("anchor failure", b_count, modulus, certificate_value))
    digest.add((
        "anchor",
        b_count,
        modulus,
        tuple(fixed_pairs),
        tuple(zip(support, weights)),
        ordinary_bound,
        tuple(zip(fixed_pairs, neighbour_sets, neighbour_bounds)),
        certificate_value,
        candidates,
    ))
    return (
        certificate_value,
        sum(weight != 0 for weight in weights),
        len(candidates),
        sum(len(neighbours) for neighbours in neighbour_sets),
    )


def exceptional_multisets(b_count, pairs):
    bad = frozenset(BAD_ANCHORS[b_count])
    answer = []
    for family in combinations_with_replacement(B_EXCEPTIONS, b_count):
        if not any(owner in H_CORRIDOR for owner in family):
            continue
        if not base_free(lcm(*family)):
            continue
        index_pairs = tuple(
            (family[left], family[right])
            for left, right in combinations(range(b_count), 2)
        )
        if not all(pair in pairs for pair in index_pairs):
            continue
        h_pairs = tuple(
            pair
            for pair in index_pairs
            if pair[0] in H_CORRIDOR or pair[1] in H_CORRIDOR
        )
        if h_pairs and all(pair in bad for pair in h_pairs):
            answer.append(family)
    return tuple(answer)


def solve_exception_group(modulus, fixed_families, b_count, digest):
    support = tuple(divisor for divisor in divisors(modulus) if divisor < modulus)
    candidates = finite_candidates(modulus, support, ordinary_only=True)
    dimension = len(support)
    matrix = []
    bounds = []
    for owner in candidates:
        matrix.append(
            [alpha(owner // gcd(owner, divisor)) for divisor in support]
            + [Fraction(-1), Fraction(0)]
        )
        bounds.append(Fraction(0))
    matrix.append(
        [Fraction(7, 44)] * dimension + [Fraction(-1), Fraction(0)]
    )
    bounds.append(Fraction(0))
    for family in fixed_families:
        matrix.append(
            [
                sum(alpha(owner // gcd(owner, divisor)) for owner in family)
                for divisor in support
            ]
            + [Fraction(6 - b_count), Fraction(-1)]
        )
        bounds.append(Fraction(0))
    matrix.append([Fraction(1)] * dimension + [Fraction(0), Fraction(0)])
    bounds.append(Fraction(1))
    matrix.append([Fraction(-1)] * dimension + [Fraction(0), Fraction(0)])
    bounds.append(Fraction(-1))
    objective_value, solution = ExactLinearProgram(
        matrix,
        bounds,
        [Fraction(0)] * (dimension + 1) + [Fraction(-1)],
    ).solve()
    weights = solution[:dimension]
    ordinary_bound = solution[dimension]
    certificate_value = solution[-1]
    require(all(value >= 0 for value in solution), ("negative exception solution", b_count, modulus))
    require(sum(weights) == 1, ("exception normalization", b_count, modulus))
    require(objective_value == -certificate_value, (b_count, modulus, objective_value))
    ordinary_loads = tuple(profile_load(weights, support, owner) for owner in candidates)
    require(
        max(ordinary_loads + (Fraction(7, 44),)) <= ordinary_bound,
        ("exception ordinary bound", b_count, modulus),
    )
    totals = tuple(
        sum(profile_load(weights, support, owner) for owner in family)
        + (6 - b_count) * ordinary_bound
        for family in fixed_families
    )
    require(max(totals) == certificate_value, ("exception value", b_count, modulus))
    require(certificate_value < 1, ("exception failure", b_count, modulus, certificate_value))
    digest.add((
        "exception",
        b_count,
        modulus,
        tuple(fixed_families),
        tuple(zip(support, weights)),
        ordinary_bound,
        certificate_value,
        candidates,
    ))
    return certificate_value, sum(weight != 0 for weight in weights), len(candidates)


def residual_six_families(pairs, digest):
    count = 0
    worst = (Fraction(0), None)
    witness_divisor_counts = defaultdict(int)
    for family in combinations_with_replacement(B_EXCEPTIONS, 6):
        if not any(owner in H_CORRIDOR for owner in family):
            continue
        modulus = lcm(*family)
        if not base_free(modulus):
            continue
        if not all(
            (family[left], family[right]) in pairs
            for left, right in combinations(range(6), 2)
        ):
            continue
        choices = tuple(
            (
                sum(alpha(owner // gcd(owner, divisor)) for owner in family),
                divisor,
            )
            for divisor in divisors(modulus)
            if divisor < modulus
        )
        best_value, best_divisor = min(choices)
        require(best_value < 1, ("six-family failure", family, modulus, best_value))
        digest.add(("six", family, modulus, best_divisor, best_value))
        witness_divisor_counts[best_divisor] += 1
        count += 1
        if best_value > worst[0]:
            worst = (best_value, (family, modulus, best_divisor))
    return count, worst, tuple(sorted(witness_divisor_counts.items()))


def cyclic_distance_numerator(residue, modulus):
    residue %= modulus
    return min(residue, modulus - residue)


def danger_set(modulus, speed):
    require(modulus >= 2 and speed % modulus != 0, (modulus, speed))
    return frozenset(
        sheet
        for sheet in range(modulus)
        if 14 * cyclic_distance_numerator(speed * sheet, modulus) < modulus
    )


def literal_coverage_types(modulus):
    by_mask = {}
    for speed in range(1, modulus):
        block = danger_set(modulus, speed)
        by_mask.setdefault(block, speed)
    return tuple(sorted((speed, block) for block, speed in by_mask.items()))


def minimum_literal_cover(modulus):
    coverage_types = literal_coverage_types(modulus)
    full = frozenset(range(modulus))
    tested = 0
    for rank in range(1, len(coverage_types) + 1):
        for chosen in combinations(coverage_types, rank):
            tested += 1
            if frozenset().union(*(block for _, block in chosen)) == full:
                return rank, tuple(speed for speed, _ in chosen), len(coverage_types), tested
    return None, (), len(coverage_types), tested


def image_order(modulus, owner_order, sheet_order):
    return owner_order // gcd(owner_order, modulus // sheet_order)


def main():
    source = Path(__file__)
    tree = ast.parse(source.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "assert node")
    require(
        not any(
            isinstance(node, ast.Constant) and isinstance(node.value, float)
            for node in ast.walk(tree)
        ),
        "floating literal",
    )

    semantic = ExactDigest()
    certificate_digest = ExactDigest()
    b_set = frozenset(B_EXCEPTIONS)
    h_set = frozenset(H_CORRIDOR)
    require(h_set < b_set, "H must be a proper subset of B")
    require(max(alpha(owner) for owner in B_EXCEPTIONS) == Fraction(1, 5), "B maximum")
    unit_b_one_bound = Fraction(1, 5) + 5 * Fraction(7, 44)
    require(unit_b_one_bound == Fraction(219, 220) < 1, unit_b_one_bound)

    pairs = residual_pairs()
    adjacency = adjacency_from_pairs(pairs)
    admissible_pairs = tuple(
        pair
        for pair in combinations_with_replacement(B_EXCEPTIONS, 2)
        if base_free(lcm(*pair))
    )
    require(len(admissible_pairs) == 368, len(admissible_pairs))
    require(len(pairs) == 213 and pairs <= frozenset(admissible_pairs), len(pairs))
    excluded_pairs = tuple(pair for pair in admissible_pairs if pair not in pairs)
    require(len(excluded_pairs) == 155, len(excluded_pairs))
    repeatable = tuple(owner for owner in B_EXCEPTIONS if (owner, owner) in pairs)
    require(len(repeatable) == 29, repeatable)

    # Verify that the displayed cliques are exactly the nontrivial maximal
    # simple cliques; 84 and 156 are the two isolated singleton cliques.
    simple_neighbours = {
        owner: frozenset(neighbour for neighbour in adjacency[owner] if neighbour != owner)
        for owner in B_EXCEPTIONS
    }
    found_maximal = []

    def bron_kerbosch(chosen, possible, excluded):
        if not possible and not excluded:
            found_maximal.append(tuple(sorted(chosen)))
            return
        for vertex in tuple(sorted(possible)):
            bron_kerbosch(
                chosen | {vertex},
                possible & set(simple_neighbours[vertex]),
                excluded & set(simple_neighbours[vertex]),
            )
            possible.remove(vertex)
            excluded.add(vertex)

    bron_kerbosch(set(), set(B_EXCEPTIONS), set())
    expected_maximal = tuple(sorted(RESIDUAL_MAXIMAL_CLIQUES + ((84,), (156,))))
    require(tuple(sorted(found_maximal)) == expected_maximal, tuple(sorted(found_maximal)))

    pair_groups = defaultdict(list)
    for pair in excluded_pairs:
        pair_groups[lcm(*pair)].append(pair)
    require(len(pair_groups) == 64, len(pair_groups))
    pair_worst = (Fraction(0), None)
    pair_support_total = 0
    pair_candidate_total = 0
    for modulus, fixed_pairs in sorted(pair_groups.items()):
        value, support_size, candidate_count = solve_pair_pruning_group(
            modulus, tuple(fixed_pairs), certificate_digest
        )
        pair_support_total += support_size
        pair_candidate_total += candidate_count
        if value > pair_worst[0]:
            pair_worst = (value, modulus)
    require(pair_worst == (Fraction(7903, 7920), 814), pair_worst)

    anchor_pairs = tuple(sorted(
        pair for pair in pairs if pair[0] in h_set or pair[1] in h_set
    ))
    require(len(anchor_pairs) == 60, len(anchor_pairs))
    anchor_rows = []
    anchor_support_total = 0
    anchor_candidate_total = 0
    anchor_neighbour_total = 0
    for b_count in range(2, 6):
        bad = frozenset(BAD_ANCHORS[b_count])
        groups = defaultdict(list)
        for pair in anchor_pairs:
            if pair not in bad:
                groups[lcm(*pair)].append(pair)
        expected = {2: (60, 45), 3: (60, 45), 4: (57, 42), 5: (55, 40)}[b_count]
        require((sum(map(len, groups.values())), len(groups)) == expected, (b_count, expected))
        worst = (Fraction(0), None)
        for modulus, fixed_pairs in sorted(groups.items()):
            value, support_size, candidate_count, neighbour_count = solve_anchor_group(
                modulus,
                tuple(fixed_pairs),
                b_count,
                adjacency,
                certificate_digest,
            )
            anchor_support_total += support_size
            anchor_candidate_total += candidate_count
            anchor_neighbour_total += neighbour_count
            if value > worst[0]:
                worst = (value, modulus)
        expected_worst = {
            2: (Fraction(1179, 1232), 102),
            3: (Fraction(6147, 6160), 102),
            4: (Fraction(121741, 122760), 1702),
            5: (Fraction(173249, 174240), 1748),
        }[b_count]
        require(worst == expected_worst, (b_count, worst))
        anchor_rows.append((b_count, expected[0], expected[1], worst))

    exception_rows = []
    exception_support_total = 0
    exception_candidate_total = 0
    exception_families_by_b = {}
    for b_count in (4, 5):
        families = exceptional_multisets(b_count, pairs)
        exception_families_by_b[b_count] = families
        expected_count = {4: 5, 5: 20}[b_count]
        require(len(families) == expected_count, (b_count, len(families)))
        groups = defaultdict(list)
        for family in families:
            groups[lcm(*family)].append(family)
        expected_groups = {4: 3, 5: 6}[b_count]
        require(len(groups) == expected_groups, (b_count, len(groups)))
        worst = (Fraction(0), None)
        for modulus, fixed_families in sorted(groups.items()):
            value, support_size, candidate_count = solve_exception_group(
                modulus,
                tuple(fixed_families),
                b_count,
                certificate_digest,
            )
            exception_support_total += support_size
            exception_candidate_total += candidate_count
            if value > worst[0]:
                worst = (value, modulus)
        expected_worst = {
            4: (Fraction(591, 616), 102),
            5: (Fraction(47, 57), 2346),
        }[b_count]
        require(worst == expected_worst, (b_count, worst))
        exception_rows.append((b_count, len(families), len(groups), worst))

    six_count, six_worst, six_divisor_counts = residual_six_families(
        pairs, certificate_digest
    )
    require(six_count == 14940, six_count)
    require(
        six_worst
        == (
            Fraction(1, 2),
            ((21, 33, 39, 77, 91, 143), 3003, 3),
        ),
        six_worst,
    )

    dilation_checks = 0
    for base, speeds in sorted(POSITIVE_WITNESSES.items()):
        require(len(speeds) <= 6, (base, speeds))
        for multiplier in range(1, 11):
            modulus = base * multiplier
            lifted = tuple(multiplier * speed for speed in speeds)
            require(all(speed % modulus != 0 for speed in lifted), (base, multiplier))
            covered = frozenset().union(*(danger_set(modulus, speed) for speed in lifted))
            require(covered == frozenset(range(modulus)), ("positive lift", base, multiplier))
            dilation_checks += 1

    q21_rank = minimum_literal_cover(21)
    q22_rank = minimum_literal_cover(22)
    require(q21_rank[:3] == (8, (1, 2, 3, 4, 5, 7, 8, 10), 8), q21_rank)
    require(q22_rank[:3] == (7, (1, 2, 3, 5, 7, 9, 11), 7), q22_rank)

    q102_orders = divisors(102)[1:]
    q102_owner_loads = tuple(
        (
            owner,
            sum(
                weight * alpha(image_order(102, owner, sheet))
                for sheet, weight in Q102_SHEET_WEIGHTS.items()
            ),
        )
        for owner in q102_orders
    )
    q102_sheet_loads = tuple(
        (
            sheet,
            sum(
                weight * alpha(image_order(102, owner, sheet))
                for owner, weight in Q102_OWNER_WEIGHTS.items()
            ),
        )
        for sheet in q102_orders
    )
    require(sum(Q102_SHEET_WEIGHTS.values()) == 1, "q102 sheet normalization")
    require(sum(Q102_OWNER_WEIGHTS.values()) == 1, "q102 owner normalization")
    require(max(load for _, load in q102_owner_loads) == Fraction(4, 25), q102_owner_loads)
    require(min(load for _, load in q102_sheet_loads) == Fraction(4, 25), q102_sheet_loads)
    require(
        tuple(owner for owner, load in q102_owner_loads if load == Fraction(4, 25))
        == (2, 3, 17, 102),
        q102_owner_loads,
    )
    require(
        tuple(sheet for sheet, load in q102_sheet_loads if load == Fraction(4, 25))
        == (6, 34, 51, 102),
        q102_sheet_loads,
    )

    certificate_hash = certificate_digest.hexdigest()
    semantic.add(("constants", BASES, B_EXCEPTIONS, H_CORRIDOR))
    semantic.add(("graph", admissible_pairs, tuple(sorted(pairs)), expected_maximal))
    semantic.add(("pair", tuple(sorted((key, tuple(value)) for key, value in pair_groups.items())), pair_worst))
    semantic.add(("anchor", tuple(anchor_rows), BAD_ANCHORS))
    semantic.add(("exceptions", tuple(exception_rows), tuple(sorted(exception_families_by_b.items()))))
    semantic.add(("six", six_count, six_worst, six_divisor_counts))
    semantic.add(("certificates", certificate_hash))
    semantic.add(("positive", POSITIVE_WITNESSES, dilation_checks))
    semantic.add(("hostiles", q21_rank, q22_rank, q102_owner_loads, q102_sheet_loads))
    semantic_digest = semantic.hexdigest()
    if EXPECTED_SEMANTIC_DIGEST is not None:
        require(
            semantic_digest == EXPECTED_SEMANTIC_DIGEST,
            (semantic_digest, EXPECTED_SEMANTIC_DIGEST),
        )

    print("THM-3414 FIXED-ZERO SIX-OWNER BASE CLASSIFICATION EXACT COMPANION")
    print(f"source_sha256_lf={lf_hash(source)}")
    print("status=COMPUTER-ASSISTED PROVED global finite-profile atlas; exact rational certificate verification")
    print("scope=fixed source centre zero and transverse cyclic-sheet owners only; no arbitrary common centre, mobile centre, physical-time, or LRC14 conclusion")
    print(f"positive_bases_and_witnesses={tuple(sorted(POSITIVE_WITNESSES.items()))}")
    print(f"positive_dilation_checks={dilation_checks};classification=rank_at_most_six_iff_one_base_divides_q")
    print(f"unit_stratum_b_at_most_one_upper={unit_b_one_bound};therefore_B_owner_count_at_least_2")
    print(
        f"pair_pruning=admissible_{len(admissible_pairs)},excluded_{len(excluded_pairs)},"
        f"residual_{len(pairs)},lcm_groups_{len(pair_groups)};worst={pair_worst[0]}@L{pair_worst[1]}"
    )
    print(
        f"pair_certificate_support_total={pair_support_total};"
        f"finite_candidate_rows_total={pair_candidate_total}"
    )
    print(f"residual_maximal_cliques={expected_maximal};repeatable_vertices={repeatable}")
    print(f"anchor_rows={tuple(anchor_rows)}")
    print(
        f"anchor_certificate_support_total={anchor_support_total};"
        f"ordinary_candidate_rows_total={anchor_candidate_total};"
        f"neighbour_rows_total={anchor_neighbour_total}"
    )
    print(f"exception_rows={tuple(exception_rows)}")
    print(
        f"exception_certificate_support_total={exception_support_total};"
        f"ordinary_candidate_rows_total={exception_candidate_total}"
    )
    print(
        f"six_B_owner_families={six_count};worst_min_profile={six_worst};"
        f"distinct_direct_divisor_witnesses={len(six_divisor_counts)}"
    )
    print(
        f"q21_shared_prime_hostile_rank={q21_rank[0]};witness={q21_rank[1]};"
        f"q22_H_corridor_hostile_rank={q22_rank[0]};witness={q22_rank[1]}"
    )
    print(
        "q102_sharp_fractional_hostile=4/25;"
        "active_owner_orders=(2,3,17,102);tight_sheet_orders=(6,34,51,102)"
    )
    print("global_certificate_margin_at_least=13/6160")
    print(f"certificate_sha256={certificate_hash}")
    print(f"semantic_sha256={semantic_digest}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
