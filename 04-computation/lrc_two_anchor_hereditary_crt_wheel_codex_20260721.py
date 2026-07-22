#!/usr/bin/env python3
"""Exact finite-field and CRT referee for THM-2062."""

from itertools import combinations
from math import gcd


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def determinant(left, right):
    return left[0] * right[1] - left[1] * right[0]


def dot(row, parameter):
    return row[0] * parameter[0] + row[1] * parameter[1]


def gcd_many(values):
    answer = 0
    for value in values:
        answer = gcd(answer, abs(value))
    return answer


def deletion_index(rows, deleted):
    retained = [row for index, row in enumerate(rows) if index != deleted]
    return gcd_many(
        determinant(left, right) for left, right in combinations(retained, 2)
    )


def prime_factors(value):
    answer = set()
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            answer.add(divisor)
            while value % divisor == 0:
                value //= divisor
        divisor += 1
    if value > 1:
        answer.add(value)
    return answer


def radical(value):
    answer = 1
    for prime in prime_factors(value):
        answer *= prime
    return answer


def projective_directions(prime):
    return [(1, slope) for slope in range(prime)] + [(0, 1)]


def bad_deletions_and_directions(rows, prime, indices):
    answer = {}
    for deleted, index in enumerate(indices):
        if index % prime:
            continue
        directions = [
            direction
            for direction in projective_directions(prime)
            if all(
                dot(row, direction) % prime == 0
                for position, row in enumerate(rows)
                if position != deleted
            )
        ]
        require(len(directions) == 1,
                ("unique projective kernel", rows, prime, deleted, directions))
        answer[deleted] = directions[0]
    return answer


def direct_hereditary(rows, parameter):
    return all(
        gcd_many(
            dot(row, parameter)
            for position, row in enumerate(rows)
            if position != deleted
        ) == 1
        for deleted in range(len(rows))
    )


def primitive_direction(row):
    divisor = gcd(abs(row[0]), abs(row[1]))
    direction = (row[0] // divisor, row[1] // divisor)
    if direction[0] < 0 or (direction[0] == 0 and direction[1] < 0):
        direction = (-direction[0], -direction[1])
    return direction


def scalar_multiple(row, direction):
    if direction[0]:
        require(row[0] % direction[0] == 0, ("nonintegral multiple", row, direction))
        value = row[0] // direction[0]
    else:
        require(direction[1] != 0 and row[1] % direction[1] == 0,
                ("nonintegral vertical multiple", row, direction))
        value = row[1] // direction[1]
    require(row == (value * direction[0], value * direction[1]),
            ("wrong line", row, direction))
    return value


def main():
    vectors = [
        (left, right)
        for left in range(-2, 3)
        for right in range(-2, 3)
        if (left, right) != (0, 0)
    ]

    saturated = 0
    all_rank_two = 0
    rank_one_templates = 0
    determinant_checks = 0
    local_prime_checks = 0
    global_direction_checks = 0
    two_direction_primes = 0
    crt_fibres = 0
    crt_residues = 0
    rank_one_parameter_checks = 0

    for rows in combinations(vectors, 4):
        if gcd_many(determinant(left, right) for left, right in combinations(rows, 2)) != 1:
            continue
        saturated += 1
        indices = [deletion_index(rows, deleted) for deleted in range(4)]

        for deleted, index in enumerate(indices):
            if index == 0:
                retained = [row for position, row in enumerate(rows) if position != deleted]
                direction = primitive_direction(retained[0])
                multipliers = [scalar_multiple(row, direction) for row in retained]
                require(gcd_many(multipliers) == 1, ("rank-one multiplier gcd", rows))
                require(abs(determinant(direction, rows[deleted])) == 1,
                        ("rank-one unimodular outlier", rows))
                rank_one_templates += 1
                for first in range(-8, 9):
                    for second in range(-8, 9):
                        if gcd(abs(first), abs(second)) != 1:
                            continue
                        parameter = (first, second)
                        deletion_gcd = gcd_many(dot(row, parameter) for row in retained)
                        require(deletion_gcd == abs(dot(direction, parameter)),
                                ("rank-one exact gcd", rows, parameter))
                        rank_one_parameter_checks += 1
            else:
                for first in range(1, 9):
                    for second in range(-8, 9):
                        if gcd(first, abs(second)) != 1:
                            continue
                        deletion_gcd = gcd_many(
                            dot(row, (first, second))
                            for position, row in enumerate(rows)
                            if position != deleted
                        )
                        require(index % deletion_gcd == 0,
                                ("determinantal divisibility", rows, deleted, first, second))
                        determinant_checks += 1

        primes = {2, 3, 5, 7}
        for index in indices:
            primes.update(prime_factors(index))
        for prime in sorted(primes):
            bad = bad_deletions_and_directions(rows, prime, indices)
            directions = set(bad.values())
            require(len(directions) <= 2,
                    ("too many projective directions", rows, prime, bad))
            if len(directions) == 2:
                require(len(bad) == 2, ("duplicate bad direction", rows, prime, bad))
                deleted_left, deleted_right = sorted(bad)
                for position, row in enumerate(rows):
                    if position not in (deleted_left, deleted_right):
                        require(row[0] % prime == 0 and row[1] % prime == 0,
                                ("two-direction common row", rows, prime, bad))
                require(determinant(rows[deleted_left], rows[deleted_right]) % prime != 0,
                        ("two-direction outliers", rows, prime, bad))
                two_direction_primes += 1
            direct_good = 0
            for first in range(prime):
                for second in range(prime):
                    if first == 0 and second == 0:
                        continue
                    parameter = (first, second)
                    if all(
                        not all(
                            dot(row, parameter) % prime == 0
                            for position, row in enumerate(rows)
                            if position != deleted
                        )
                        for deleted in range(len(rows))
                    ):
                        direct_good += 1
            require(direct_good == (prime - 1) * (prime + 1 - len(directions)),
                    ("global projective count", rows, prime, direct_good, directions))
            global_direction_checks += 1
            local_prime_checks += 1

        if all(indices):
            all_rank_two += 1
            index_product = 1
            for index in indices:
                index_product *= index
            for fixed_n in range(1, 9):
                modulus = radical(fixed_n * index_product)
                direct = sum(
                    gcd(fixed_n, longitudinal) == 1
                    and direct_hereditary(rows, (fixed_n, longitudinal))
                    for longitudinal in range(modulus)
                )
                predicted = 1
                for prime in sorted(prime_factors(modulus)):
                    bad = bad_deletions_and_directions(rows, prime, indices)
                    forbidden = set(bad.values())
                    if fixed_n % prime:
                        nu = sum(direction[0] != 0 for direction in forbidden)
                        local = prime - nu
                    else:
                        local = 0 if (0, 1) in forbidden else prime - 1
                    predicted *= local
                require(direct == predicted,
                        ("CRT product", rows, fixed_n, modulus, direct, predicted))
                crt_fibres += 1
                crt_residues += modulus

        if saturated >= 5000:
            break

    require(saturated == 5000, ("template target", saturated))

    guardrail_rows = ((2, 0), (3, 0), (4, 0), (0, 1))
    require(gcd_many(determinant(x, y) for x, y in combinations(guardrail_rows, 2)) == 1,
            "guardrail saturation")
    require(deletion_index(guardrail_rows, 3) == 0, "guardrail rank one")
    require(abs(dot((1, 0), (1, 2))) == 1, "guardrail affine terminal")
    require(not direct_hereditary(guardrail_rows, (1, 2)),
            "affine terminal alone is not hereditary")

    fibre_death_rows = (
        (0, 1), (1, 1), (2, 0), (2, 2), (2, 4), (4, 2),
    )
    require(gcd_many(determinant(x, y) for x, y in combinations(fibre_death_rows, 2)) == 1,
            "fibre-death saturation")
    fibre_indices = [deletion_index(fibre_death_rows, deleted)
                     for deleted in range(len(fibre_death_rows))]
    fibre_bad = bad_deletions_and_directions(fibre_death_rows, 2, fibre_indices)
    require(set(fibre_bad.values()) == {(1, 0), (1, 1)},
            ("two affine residues at p=2", fibre_bad))
    require(not any(direct_hereditary(fibre_death_rows, (1, m)) for m in range(2)),
            "p=2 fibre death with p not dividing N")

    print("THM-2062 TWO-ANCHOR HEREDITARY CRT-WHEEL AUDIT")
    print("saturated four-row templates checked:", saturated)
    print("all-rank-two templates:", all_rank_two)
    print("rank-one deletion templates/incidences:", rank_one_templates)
    print("determinantal parameter checks:", determinant_checks)
    print("local prime/template checks:", local_prime_checks)
    print("global projective-density checks:", global_direction_checks)
    print("two-direction prime patterns:", two_direction_primes)
    print("fixed-N CRT fibres / residues checked:", crt_fibres, crt_residues)
    print("rank-one parameter checks:", rank_one_parameter_checks)
    print("guardrails: affine line alone insufficient; p=2 affine fibre can die")
    print("carrier: deletion indices + prime-labelled projective annihilators")
    print("tournament fingerprint: at most two vertices, always transitive")
    print("PASS")


if __name__ == "__main__":
    main()
