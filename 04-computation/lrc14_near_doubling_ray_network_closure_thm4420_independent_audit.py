#!/usr/bin/env python3
"""Clean-room exact referee for the THM-4420 near-doubling-ray claim.

This script imports no repository code and does not use the candidate's
carrier constructor when performing its relation-lattice controls.  It uses
the canonical THM-4414 raw projection formula as the definition of the three
degree-zero capacities.
"""

from fractions import Fraction
from hashlib import sha256
from json import dumps
import sys


Q = Fraction
TARGET = Q(6, 77)
CHECKS = 0


class RefereeFailure(RuntimeError):
    pass


def require(predicate, label):
    global CHECKS
    CHECKS += 1
    if not predicate:
        raise RefereeFailure(label)


def strict_integer_cutoff(numerator, denominator=14):
    """Largest nonnegative integer k satisfying denominator*k<numerator."""
    k = numerator // denominator
    if denominator * k == numerator:
        k -= 1
    return k


def eligible_rows(max_third_speed):
    for sigma, first_m in ((1, 5), (-1, 7)):
        m = first_m
        while 2 * m + sigma <= max_third_speed:
            yield m, sigma
            m += 6


def margins(w, carrier):
    result = []
    for i in range(3):
        j, k = (index for index in range(3) if index != i)
        result.append(3 * (w[j] + w[k]) - 14 * abs(carrier[i]))
    return tuple(result)


def projection_terms(w, carrier):
    p = margins(w, carrier)
    q = Q(3, 7 * max(w))
    result = []
    for i in range(3):
        j, k = (index for index in range(3) if index != i)
        result.append(min(q, Q(p[i], 14 * w[j] * w[k])))
    return tuple(result)


def ray_carriers(m, sigma):
    bound_numerator = 3 * (m + 1 if sigma == 1 else m)
    cutoff = strict_integer_cutoff(bound_numerator)
    carriers = {
        (-sigma * k, -2 * k, k)
        for k in range(-cutoff, cutoff + 1)
        if k % 3
    }
    return cutoff, carriers


def ray_cutoff_and_count(m, sigma):
    bound_numerator = 3 * (m + 1 if sigma == 1 else m)
    cutoff = strict_integer_cutoff(bound_numerator)
    return cutoff, 2 * (cutoff - cutoff // 3)


def exhaustive_carriers(m, sigma):
    """Enumerate the support box directly, solving only the original relation."""
    w = (1, m, 2 * m + sigma)
    b1 = strict_integer_cutoff(3 * (w[1] + w[2]))
    b3 = strict_integer_cutoff(3 * (w[0] + w[1]))
    found = set()
    for c1 in range(-b1, b1 + 1):
        for c3 in range(-b3, b3 + 1):
            numerator = -(c1 + w[2] * c3)
            if numerator % m:
                continue
            c2 = numerator // m
            carrier = (c1, c2, c3)
            if any(value % 3 == 0 for value in carrier):
                continue
            if all(value > 0 for value in margins(w, carrier)):
                found.add(carrier)
    return found


def main():
    sys.stdout.reconfigure(newline="\n")
    max_third_speed = 2_000_003
    exact_controls = {
        (5, 1), (11, 1), (17, 1), (23, 1), (41, 1), (83, 1),
        (101, 1), (251, 1),
        (7, -1), (13, -1), (19, -1), (25, -1), (61, -1),
        (103, -1), (145, -1), (253, -1),
    }
    controls_seen = set()
    envelope_equalities = []
    strict_endpoint_rows = []
    semantic = []
    rows = 0

    for m, sigma in eligible_rows(max_third_speed):
        rows += 1
        w = (1, m, 2 * m + sigma)
        require(w[0] < w[1] < w[2], ("sorted", w))
        require(all(value % 2 for value in w), ("odd", w))
        require(all(value % 3 for value in w), ("ternary-unit", w))

        cutoff, carrier_count = ray_cutoff_and_count(m, sigma)

        outside_k = cutoff + 1
        outside = (-sigma * outside_k, -2 * outside_k, outside_k)
        require(min(margins(w, outside)) <= 0,
                ("strict-cutoff", w, cutoff, outside, margins(w, outside)))
        if sigma == 1 and 3 * (m + 1) % 14 == 0:
            require(min(margins(w, outside)) == 0,
                    ("open-endpoint", w, cutoff, margins(w, outside)))
            strict_endpoint_rows.append(w)

        require(11 * carrier_count <= 2 * w[2],
                ("integer-target", w, carrier_count))
        if 11 * carrier_count == 2 * w[2]:
            envelope_equalities.append(w)

        if (m, sigma) in exact_controls:
            controls_seen.add((m, sigma))
            formula_cutoff, carriers = ray_carriers(m, sigma)
            require(formula_cutoff == cutoff and len(carriers) == carrier_count,
                    ("ray-count", w, cutoff, carrier_count, len(carriers)))
            require(all(sum(w[i] * c[i] for i in range(3)) == 0
                        for c in carriers), ("ray-relation", w))
            require(all(all(value > 0 for value in margins(w, c))
                        for c in carriers), ("ray-support", w))
            complete = exhaustive_carriers(m, sigma)
            require(complete == carriers,
                    ("complete-lattice-vs-ray", w, complete ^ carriers))
            require(all(c[1] + 2 * c[2] == 0 for c in complete),
                    ("integral-collapse", w, complete))

            capacities = [Q(0), Q(0), Q(0)]
            physical = Q(0)
            q = Q(3, 7 * w[2])
            for carrier in complete:
                terms = projection_terms(w, carrier)
                for i, term in enumerate(terms):
                    capacities[i] += term
                physical += min(terms)
            count_envelope = carrier_count * q
            require(all(value <= count_envelope for value in capacities),
                    ("count-envelope", w, capacities, count_envelope))
            require(count_envelope <= TARGET,
                    ("target", w, count_envelope))
            if w == (1, 5, 11):
                require(capacities == [TARGET, TARGET, TARGET],
                        ("equality-capacities", w, capacities))
                require(physical == TARGET,
                        ("equality-physical", w, physical))

        semantic.append((w, cutoff, carrier_count))

    require(controls_seen == exact_controls,
            ("missing-controls", exact_controls - controls_seen))
    require(envelope_equalities == [(1, 5, 11)],
            ("envelope-equalities", envelope_equalities))
    require(strict_endpoint_rows[:4] == [(1, 41, 83), (1, 83, 167),
                                         (1, 125, 251), (1, 167, 335)],
            ("endpoint-progression", strict_endpoint_rows[:4]))

    digest = sha256(dumps(semantic, separators=(",", ":")).encode()).hexdigest()
    print("THM-4420 CLEAN-ROOM NEAR-DOUBLING REFEREE")
    print("verdict=PASS")
    print(f"max_third_speed={max_third_speed}; family_rows={rows}")
    print(f"complete_relation_lattice_controls={len(controls_seen)}")
    print("all_three_projection_equalities_only_at=(1,5,11)")
    print("physical_equality_only_at=(1,5,11)")
    print("first_strict_endpoint_rows=" + str(strict_endpoint_rows[:4]))
    print("semantic_sha256=" + digest)
    print(f"checks={CHECKS}")


if __name__ == "__main__":
    main()
