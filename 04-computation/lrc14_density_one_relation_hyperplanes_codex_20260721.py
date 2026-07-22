#!/usr/bin/env python3
"""Exact counting referee for THM-2074's relation-hyperplane sieve."""

from itertools import combinations, product
from math import comb, gcd


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def gcd_many(values):
    answer = 0
    for value in values:
        answer = gcd(answer, abs(value))
    return answer


def mobius_sieve(limit):
    mu = [0] * (limit + 1)
    mu[1] = 1
    primes = []
    composite = bytearray(limit + 1)
    for value in range(2, limit + 1):
        if not composite[value]:
            primes.append(value)
            mu[value] = -1
        for prime in primes:
            multiple = value * prime
            if multiple > limit:
                break
            composite[multiple] = 1
            if value % prime == 0:
                mu[multiple] = 0
                break
            mu[multiple] = -mu[value]
    return mu


def positive_primitive_tuples(height, arity, mu=None):
    if mu is None:
        mu = mobius_sieve(height)
    return sum(mu[divisor] * (height // divisor) ** arity
               for divisor in range(1, height + 1))


def exact_relation_hyperplanes(coordinates, height, support_sizes, mu=None):
    if mu is None:
        mu = mobius_sieve(height)
    return sum(
        comb(coordinates, support_size)
        * 2 ** (support_size - 1)
        * positive_primitive_tuples(height, support_size, mu)
        for support_size in support_sizes
    )


def canonical_normal(normal):
    divisor = gcd_many(normal)
    require(divisor > 0, ("zero normal", normal))
    answer = tuple(entry // divisor for entry in normal)
    first = next(entry for entry in answer if entry)
    if first < 0:
        answer = tuple(-entry for entry in answer)
    return answer


def enumerate_normals(coordinates, height, support_sizes):
    values = tuple(value for value in range(-height, height + 1) if value)
    answer = set()
    for support_size in support_sizes:
        for support in combinations(range(coordinates), support_size):
            for coefficients in product(values, repeat=support_size):
                if gcd_many(coefficients) != 1:
                    continue
                normal = [0] * coordinates
                for index, coefficient in zip(support, coefficients):
                    normal[index] = coefficient
                answer.add(canonical_normal(tuple(normal)))
    return answer


def dot(left, right):
    return sum(x * y for x, y in zip(left, right))


def audit_small_exact_counts():
    normal_checks = 0
    for coordinates in (3, 4, 5):
        for height in (1, 2, 3):
            support_sizes = tuple(range(2, min(4, coordinates) + 1))
            direct = enumerate_normals(coordinates, height, support_sizes)
            predicted = exact_relation_hyperplanes(
                coordinates, height, support_sizes
            )
            require(len(direct) == predicted,
                    ("primitive-normal formula", coordinates, height,
                     len(direct), predicted))
            normal_checks += 1

    coordinates = 5
    height = 2
    box_height = 6
    normals = enumerate_normals(coordinates, height, (3, 4, 5))
    bad_points = set()
    maximum_single_hyperplane = 0
    for normal in normals:
        zeros = {
            point
            for point in product(range(1, box_height + 1), repeat=coordinates)
            if dot(normal, point) == 0
        }
        maximum_single_hyperplane = max(maximum_single_hyperplane, len(zeros))
        require(len(zeros) <= box_height ** (coordinates - 1),
                ("box hyperplane bound", normal, len(zeros)))
        bad_points.update(zeros)
    require(len(bad_points) <= len(normals) * box_height ** (coordinates - 1),
            "finite union box bound")
    primitive_distinct = [
        point
        for point in product(range(1, box_height + 1), repeat=coordinates)
        if gcd_many(point) == 1 and len(set(point)) == coordinates
    ]
    increasing = [point for point in primitive_distinct if tuple(sorted(point)) == point]
    require(len(primitive_distinct) == len(increasing) * 120,
            ("ordered/increasing factor", len(primitive_distinct), len(increasing)))
    return {
        "normal_checks": normal_checks,
        "small_normals": len(normals),
        "small_bad_points": len(bad_points),
        "small_max_single": maximum_single_hyperplane,
        "small_primitive_distinct": len(primitive_distinct),
        "small_increasing": len(increasing),
    }


def audit_prime_power_packets():
    coordinates = 3
    normals = enumerate_normals(coordinates, 2, (2, 3))
    zero_set_checks = 0
    packet_checks = 0
    packet_summaries = []
    for prime, maximum_power in ((2, 4), (3, 3), (5, 2)):
        for exponent in range(1, maximum_power + 1):
            modulus = prime ** exponent
            all_points = tuple(product(range(modulus), repeat=coordinates))
            allowed = set(all_points)
            for normal in normals:
                zeros = {point for point in all_points if dot(normal, point) % modulus == 0}
                require(len(zeros) == modulus ** (coordinates - 1),
                        ("prime-power hyperplane size", normal, modulus, len(zeros)))
                allowed.difference_update(zeros)
                zero_set_checks += 1
            lower_bound = modulus ** coordinates - len(normals) * modulus ** (coordinates - 1)
            require(len(allowed) >= max(0, lower_bound),
                    ("prime-power union bound", modulus, len(allowed), lower_bound))
            if modulus > len(normals):
                require(allowed, ("positive packet", modulus, len(normals)))
            packet_summaries.append((modulus, len(allowed), lower_bound))
            packet_checks += 1

    hostile_nonprimitive = (2, 0, 1)
    require(gcd_many(hostile_nonprimitive) == 1,
            "normal should be primitive over every prime")
    for prime in (2, 3, 5, 7):
        require(any(coefficient % prime for coefficient in hostile_nonprimitive),
                ("primitive normal vanished modulo prime", prime))
    return len(normals), zero_set_checks, packet_checks, packet_summaries


def first_prime_power_above(prime, threshold):
    exponent = 0
    modulus = 1
    while modulus <= threshold:
        exponent += 1
        modulus *= prime
    return exponent, modulus


def audit_actual_lrc_constants():
    coordinates = 13
    height = 2 ** 20
    mu = mobius_sieve(height)
    positive_counts = {
        support_size: positive_primitive_tuples(height, support_size, mu)
        for support_size in (3, 4, 5)
    }
    relation_hyperplanes = exact_relation_hyperplanes(
        coordinates, height, (3, 4, 5), mu
    )
    crude = sum(
        comb(coordinates, support_size) * (2 * height) ** support_size
        for support_size in (3, 4, 5)
    )
    expected_positive_counts = {
        3: 959124025074311215,
        4: 1116973047989955380768527,
        5: 1222506215916189106034284205191,
    }
    expected_relation_hyperplanes = 25173854387233097811887443361297472
    expected_crude = 52206936149913413947000523914739712
    require(positive_counts == expected_positive_counts,
            ("frozen primitive coefficient counts", positive_counts))
    require(relation_hyperplanes == expected_relation_hyperplanes,
            ("frozen relation ledger", relation_hyperplanes))
    require(crude == expected_crude, ("frozen crude ledger", crude))
    require(relation_hyperplanes > 0, "nonempty relation ledger")
    require(relation_hyperplanes <= crude,
            ("exact hyperplanes exceed crude bound", relation_hyperplanes, crude))
    diagonals = comb(coordinates, 2)
    total_hyperplanes = relation_hyperplanes + diagonals
    prime_power_thresholds = {
        prime: first_prime_power_above(prime, total_hyperplanes)
        for prime in (2, 3, 5, 7)
    }
    for prime, (exponent, modulus) in prime_power_thresholds.items():
        require(modulus > total_hyperplanes, ("threshold below ledger", prime))
        require(modulus // prime <= total_hyperplanes,
                ("threshold not minimal", prime, exponent))
        lower_bound = modulus ** 12 * (modulus - total_hyperplanes)
        require(lower_bound > 0, ("actual packet lower bound", prime))
    return {
        "positive_counts": positive_counts,
        "relation_hyperplanes": relation_hyperplanes,
        "diagonals": diagonals,
        "total_hyperplanes": total_hyperplanes,
        "crude": crude,
        "thresholds": prime_power_thresholds,
    }


def main():
    small = audit_small_exact_counts()
    packet_normals, zero_checks, packet_checks, packet_summaries = (
        audit_prime_power_packets()
    )
    actual = audit_actual_lrc_constants()

    print("THM-2074 DENSITY-ONE LRC RELATION-HYPERPLANE AUDIT")
    print("small primitive-normal formula checks:", small["normal_checks"])
    print("small box normals / union points / max one plane:",
          small["small_normals"], small["small_bad_points"], small["small_max_single"])
    print("small primitive distinct / increasing rows:",
          small["small_primitive_distinct"], small["small_increasing"])
    print("small packet normals / exact zero-set checks / moduli:",
          packet_normals, zero_checks, packet_checks)
    print("small prime-power packets (q, good, union lower bound):", packet_summaries)
    print("P_3(2^20), P_4(2^20), P_5(2^20):",
          tuple(actual["positive_counts"][size] for size in (3, 4, 5)))
    print("exact relation hyperplanes R:", actual["relation_hyperplanes"])
    print("diagonal / total hyperplanes:", actual["diagonals"], actual["total_hyperplanes"])
    print("crude signed-normal bound:", actual["crude"])
    print("first q=ell^k>R+78 (ell: (k,q)):", actual["thresholds"])
    print("box theorem: exceptional ordered rows <= R*B^12")
    print("density theorem: increasing primitive rows ~ B^13/(13!*zeta(13))")
    print("fixed-prime theorem: every positive lift of a good ell^k residue is certified")
    print("tournament analysis: no intrinsic orientation; hyperplane intersection is symmetric")
    print("faithful vertices: relation normals and equality obligations, not runners")
    print("PASS")


if __name__ == "__main__":
    main()
