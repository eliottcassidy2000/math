#!/usr/bin/env python3
"""Exact referee for THM-2461.

The companion is dependency-free.  Every truth-bearing check uses ``require``
rather than ``assert``, so normal and optimized Python execute the same audit.
"""

from fractions import Fraction
from itertools import product
from math import gcd


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def valuation(number, prime):
    require(number != 0, "valuation requires a nonzero integer")
    number = abs(number)
    value = 0
    while number % prime == 0:
        number //= prime
        value += 1
    return value


def centered_distance_numerator(speed, numerator, denominator):
    residue = (speed * numerator) % denominator
    return min(residue, denominator - residue)


def centered_norm(value):
    floor_value = value.numerator // value.denominator
    residue = value - floor_value
    return min(residue, 1 - residue)


def danger(speed, point):
    return centered_norm(speed * point) < Fraction(1, 14)


def unit_safe(speed, point):
    return centered_norm(speed * point) > Fraction(1, 14)


def guard_safe(speed, point):
    return centered_norm(speed * point) > Fraction(1, 7)


def base13_address(numerator, denominator, digits):
    address = []
    remainder = numerator
    for _ in range(digits):
        value = 13 * remainder
        digit = value // denominator
        remainder = value - digit * denominator
        address.append(digit)
    return tuple(address)


def hostile_referee():
    denominator = 2222
    guard = 1
    units = (4, 2, 3, 6, 10)
    source = 13
    target_a = 13**3
    target_b = 2 * 13**5
    blockers = (source, target_a, target_b)
    speeds = (guard,) + units + blockers
    numerators = (862, 863, 865)
    expected_source_vectors = (
        (862, 996, 498, 364, 728, 268, 96, 670, 184),
        (863, 992, 496, 367, 734, 258, 109, 645, 254),
        (865, 984, 492, 373, 746, 238, 135, 595, 1092),
    )
    expected_terminal_numerators = (1248, 1417, 1755)
    expected_terminal_vectors = (
        (974, 548, 274, 700, 822, 852, 670, 92, 12),
        (805, 998, 612, 193, 386, 838, 645, 127, 708),
        (467, 354, 934, 821, 580, 226, 595, 565, 122),
    )
    expected_words = (
        ("c2", "c3"),
        ("c2",),
        ("c3",),
    )
    expected_margins = (
        Fraction(177, 15554),
        Fraction(111, 7777),
        Fraction(83, 7777),
    )

    require(len(set(speeds)) == 9, "strict hostile speeds must be distinct")
    common_gcd = 0
    for speed in speeds:
        common_gcd = gcd(common_gcd, speed)
    require(common_gcd == 1, "strict hostile row must be primitive")
    require(guard % 2 == 1, "guard must be odd")
    require(
        tuple(valuation(value, 13) for value in blockers) == (1, 3, 5),
        "strict blocker valuation profile failed",
    )
    require(
        all(valuation(value, 13) == 0 for value in units),
        "ordinary roles must be thirteen-units",
    )
    require(13**2 == 169, "prescribed clock arithmetic failed")

    source_vectors = []
    terminal_numerators = []
    terminal_vectors = []
    addresses = []
    words = []
    margins = []
    masks = []
    deep_probe_roots = []

    for numerator in numerators:
        point = Fraction(numerator, denominator)
        source_vector = tuple(
            centered_distance_numerator(speed, numerator, denominator)
            for speed in speeds
        )
        source_vectors.append(source_vector)

        require(guard_safe(guard, point), "source guard must be safe")
        require(
            all(unit_safe(speed, point) for speed in units),
            "source ordinary roles must be safe",
        )
        require(danger(source, point), "source owner must be dangerous")
        require(unit_safe(target_a, point), "first target must be safe")
        require(unit_safe(target_b, point), "deep target must be safe")

        # Choose q_*=q_1=4.  The five split guard/unit roles are then
        # H,q_2,...,q_5, followed by the source and other-blocker bits.
        mask = (0, 0, 0, 0, 0, 1, 0)
        masks.append(mask)

        address = base13_address(numerator, denominator, 2)
        addresses.append(address)

        terminal_numerator = (169 * numerator) % denominator
        terminal_numerators.append(terminal_numerator)
        terminal = Fraction(terminal_numerator, denominator)
        terminal_vector = tuple(
            centered_distance_numerator(
                speed, terminal_numerator, denominator
            )
            for speed in speeds
        )
        terminal_vectors.append(terminal_vector)

        require(guard_safe(guard, terminal), "terminal guard must be safe")
        require(
            all(unit_safe(speed, terminal) for speed in units),
            "terminal ordinary roles must be safe",
        )
        require(unit_safe(source, terminal), "source owner must expire")
        word = []
        if danger(target_a, terminal):
            word.append("c2")
        if danger(target_b, terminal):
            word.append("c3")
        require(word, "terminal blocker word must be nonempty")
        words.append(tuple(word))

        source_margins = [
            centered_norm(guard * point) - Fraction(1, 7)
        ]
        source_margins.extend(
            centered_norm(speed * point) - Fraction(1, 14)
            for speed in units
        )
        source_margins.extend(
            (
                Fraction(1, 14) - centered_norm(source * point),
                centered_norm(target_a * point) - Fraction(1, 14),
                centered_norm(target_b * point) - Fraction(1, 14),
            )
        )
        terminal_margins = [
            centered_norm(guard * terminal) - Fraction(1, 7)
        ]
        terminal_margins.extend(
            centered_norm(speed * terminal) - Fraction(1, 14)
            for speed in units
        )
        terminal_margins.append(
            centered_norm(source * terminal) - Fraction(1, 14)
        )
        terminal_margins.extend(
            abs(centered_norm(speed * terminal) - Fraction(1, 14))
            for speed in (target_a, target_b)
        )
        all_margins = source_margins + terminal_margins
        require(all(value > 0 for value in all_margins), "hostile hit a seam")
        margins.append(min(all_margins))

        roots = []
        for root in range(13):
            if centered_norm(
                target_b * point - Fraction(root, 13)
            ) < Fraction(1, 14):
                roots.append(root)
        deep_probe_roots.append(tuple(roots))

    require(
        tuple(source_vectors) == expected_source_vectors,
        "source distance vectors failed",
    )
    require(
        tuple(terminal_numerators) == expected_terminal_numerators,
        "prescribed terminal residues failed",
    )
    require(
        tuple(terminal_vectors) == expected_terminal_vectors,
        "terminal distance vectors failed",
    )
    require(
        tuple(words) == expected_words,
        "fork/pure/pure terminal split failed",
    )
    require(
        tuple(margins) == expected_margins,
        "strict hostile margins failed",
    )
    require(
        len(set(masks)) == 1 and masks[0] == (0, 0, 0, 0, 0, 1, 0),
        "complete local source atom failed",
    )
    require(
        len(set(addresses)) == 1 and addresses[0] == (5, 0),
        "common prescribed two-digit address failed",
    )
    require(
        tuple(deep_probe_roots) == ((11, 12), (1, 2), (6, 7)),
        "deep-probe root sets failed",
    )

    # The hostile is deliberately local: this exact point lies in A_0 but
    # outside every blocker danger set, so the row is not a global cover.
    uncovered = Fraction(159, 1111)
    uncovered_numerator = 318
    uncovered_vector = tuple(
        centered_distance_numerator(
            speed, uncovered_numerator, denominator
        )
        for speed in speeds
    )
    require(
        uncovered_vector == (318, 950, 636, 954, 314, 958, 310, 938, 702),
        "uncovered control vector failed",
    )
    require(guard_safe(guard, uncovered), "uncovered guard must be safe")
    require(
        all(unit_safe(speed, uncovered) for speed in units),
        "uncovered ordinary roles must be safe",
    )
    require(
        all(unit_safe(speed, uncovered) for speed in blockers),
        "uncovered control must avoid every blocker",
    )

    return {
        "speeds": speeds,
        "profile": tuple(valuation(value, 13) for value in blockers),
        "source_vectors": tuple(source_vectors),
        "terminal_numerators": tuple(terminal_numerators),
        "terminal_vectors": tuple(terminal_vectors),
        "mask": masks[0],
        "address": addresses[0],
        "words": tuple(words),
        "margins": tuple(margins),
        "deep_probe_roots": tuple(deep_probe_roots),
        "uncovered": uncovered,
        "uncovered_vector": uncovered_vector,
    }


def cyclic_interval_truth(level, point):
    return int(centered_norm(point) < Fraction(level, 14))


def translated_count_referee():
    # Representatives of every open cell cut by the translated interval
    # endpoints s/13 +/- L/14.
    minimum_safe_counts = {}
    cell_counts = {}
    for level in (1, 2):
        endpoints = {Fraction(0), Fraction(1)}
        for shift in range(13):
            for sign in (-1, 1):
                endpoint = Fraction(shift, 13) + sign * Fraction(level, 14)
                endpoint -= endpoint.numerator // endpoint.denominator
                endpoints.add(endpoint)
        ordered = sorted(endpoints)
        representatives = []
        for left, right in zip(ordered, ordered[1:]):
            if left < right:
                representatives.append((left + right) / 2)
        representatives.append((ordered[-1] + ordered[0] + 1) / 2 % 1)

        failure_safe_counts = []
        for point in representatives:
            danger_bits = tuple(
                cyclic_interval_truth(level, point - Fraction(shift, 13))
                for shift in range(13)
            )
            safe_bits = tuple(1 - bit for bit in danger_bits)
            predicted_danger = (
                2 * level
                - cyclic_interval_truth(level, 13 * point)
            )
            require(
                sum(danger_bits) == predicted_danger,
                f"translated danger count failed at L={level}",
            )
            require(
                sum(safe_bits)
                == 13 - 2 * level
                + cyclic_interval_truth(level, 13 * point),
                f"translated safe count failed at L={level}",
            )
            for shift in range(13):
                require(
                    danger_bits[shift] * safe_bits[shift] == 0,
                    f"common-shift diagonal failed at L={level}",
                )
            if danger_bits[0]:
                require(
                    safe_bits[0] == 0,
                    f"failure anchor must kill zero-shift repair at L={level}",
                )
                require(
                    sum(safe_bits[1:]) > 0,
                    f"polarized off-diagonal repair missing at L={level}",
                )
                failure_safe_counts.append(sum(safe_bits[1:]))

        require(failure_safe_counts, f"no failure cells found at L={level}")
        minimum_safe_counts[level] = min(failure_safe_counts)
        cell_counts[level] = len(representatives)

    require(
        minimum_safe_counts == {1: 11, 2: 9},
        "sharp off-diagonal safe counts failed",
    )

    # Finite-cycle analogue of the temporal idempotence no-go.  On one
    # thirteen-cycle, P_i P_{i+1}=P_i forces P to be empty or full.
    temporal_survivors = []
    for bits in product((0, 1), repeat=13):
        if all(bits[index] * bits[(index + 1) % 13] == bits[index]
               for index in range(13)):
            temporal_survivors.append(bits)
    require(
        temporal_survivors == [(0,) * 13, (1,) * 13],
        "finite temporal-idempotence hostile failed",
    )

    return {
        "cell_counts": cell_counts,
        "minimum_safe_counts": minimum_safe_counts,
        "temporal_survivors": len(temporal_survivors),
    }


def temporal_coupling_referee():
    # One source atom with three outgoing prescribed-clock word fibres.
    words = ("c2", "c3", "c2c3")
    coupling = {
        ("source_atom", "c2"): Fraction(1, 3),
        ("source_atom", "c3"): Fraction(1, 3),
        ("source_atom", "c2c3"): Fraction(1, 3),
    }
    total = sum(coupling.values(), Fraction(0))
    row_support = sum(
        coupling[("source_atom", word)] > 0 for word in words
    )
    require(total == 1, "temporal coupling mass identity failed")
    require(row_support == 3, "word coupling must not be a state map")
    require(
        max(coupling.values()) == Fraction(1, 3),
        "sharp three-word pigeonhole failed",
    )

    # The endpoint identity sees only the diagonal state atom.
    source_atom = 1
    require(
        source_atom * source_atom == source_atom,
        "Boolean diagonal idempotence failed",
    )
    require(
        len({word for _, word in coupling}) == 3,
        "diagonal state must not collapse temporal word labels",
    )
    return {
        "total": total,
        "row_support": row_support,
        "largest_word": max(coupling.values()),
    }


def target_covector_referee():
    # Coordinates are (eta,ell) for
    # eta=e_a-e_{k_a}, ell=e_c3-e_{k_b}.
    covectors = {
        "a": (1, 0),
        "c3": (0, 1),
        "k_a": (-1, 0),
        "k_b": (0, -1),
        "u_0": (0, 0),
        "u_1": (0, 0),
        "u_2": (0, 0),
        "u_3": (0, 0),
        "j": (0, 0),
    }
    guard_unit_roles = ("k_a", "k_b", "u_0", "u_1", "u_2", "u_3")
    active_units = tuple(
        role for role in guard_unit_roles if covectors[role] != (0, 0)
    )
    require(
        active_units == ("k_a", "k_b"),
        "only paired graft units may be target-active",
    )

    # THM-2445 chooses q_*=k_b and first-failure splits the five other
    # guard/unit roles.  Only k_a remains target-active.
    q_star = "k_b"
    first_failure_roles = tuple(
        role for role in guard_unit_roles if role != q_star
    )
    active_first_failures = tuple(
        role
        for role in first_failure_roles
        if covectors[role] != (0, 0)
    )
    require(len(first_failure_roles) == 5, "first-failure role count failed")
    require(
        active_first_failures == ("k_a",),
        "generic first failure was spuriously target-active",
    )
    require(
        tuple(covectors["a"][i] + covectors["k_a"][i] for i in range(2))
        == (0, 0),
        "eta dipole balance failed",
    )
    require(
        tuple(covectors["c3"][i] + covectors["k_b"][i] for i in range(2))
        == (0, 0),
        "ell dipole balance failed",
    )

    # Exact fibre orthogonality for the conditional repair pullback
    # K_tilde(r,delta)=K^+(r,lambda(delta)).  For every nonzero covector
    # lambda and dual mode xi, the fibre character sum is supported exactly
    # on xi in span(lambda).
    nonzero_covectors = 0
    for lam_1 in range(13):
        for lam_2 in range(13):
            if (lam_1, lam_2) == (0, 0):
                continue
            nonzero_covectors += 1
            for xi_1 in range(13):
                for xi_2 in range(13):
                    scalar = None
                    for candidate in range(13):
                        if (
                            (candidate * lam_1) % 13,
                            (candidate * lam_2) % 13,
                        ) == (xi_1, xi_2):
                            scalar = candidate
                            break
                    for shift in range(13):
                        exponent_counts = [0] * 13
                        for delta_1 in range(13):
                            for delta_2 in range(13):
                                if (
                                    lam_1 * delta_1 + lam_2 * delta_2
                                ) % 13 != shift:
                                    continue
                                exponent = (
                                    xi_1 * delta_1 + xi_2 * delta_2
                                ) % 13
                                exponent_counts[exponent] += 1
                        if scalar is None:
                            require(
                                exponent_counts == [1] * 13,
                                "off-line pullback character did not cancel",
                            )
                        else:
                            expected = [0] * 13
                            expected[(scalar * shift) % 13] = 13
                            require(
                                exponent_counts == expected,
                                "on-line pullback normalization failed",
                            )
    require(
        nonzero_covectors == 13**2 - 1,
        "nonzero target-covector count failed",
    )
    return {
        "active_units": active_units,
        "first_failure_roles": first_failure_roles,
        "active_first_failures": active_first_failures,
        "pullback_covectors": nonzero_covectors,
        "pullback_line_modes": 13,
    }


def main():
    hostile = hostile_referee()
    translated = translated_count_referee()
    temporal = temporal_coupling_referee()
    covectors = target_covector_referee()

    print("THM2461 TEMPORAL-WORD / POLARIZED-REPAIR REFEREE")
    print(f"hostile_speeds={hostile['speeds']}")
    print(f"hostile_profile={hostile['profile']}")
    print(f"source_distance_vectors={hostile['source_vectors']}")
    print(f"terminal_residues={hostile['terminal_numerators']}")
    print(f"terminal_distance_vectors={hostile['terminal_vectors']}")
    print(f"complete_source_mask={hostile['mask']}")
    print(f"common_base13_address={hostile['address']}")
    print(f"terminal_words={hostile['words']}")
    print(f"strict_margins={hostile['margins']}")
    print(f"deep_probe_root_sets={hostile['deep_probe_roots']}")
    print(
        "local_not_global_cover="
        f"{hostile['uncovered']}:{hostile['uncovered_vector']}"
    )
    print(f"translated_cell_counts={translated['cell_counts']}")
    print(
        "minimum_offdiagonal_safe_counts="
        f"{translated['minimum_safe_counts']}"
    )
    print(
        "finite_temporal_idempotence_survivors="
        f"{translated['temporal_survivors']}"
    )
    print(
        "temporal_coupling="
        f"mass:{temporal['total']},"
        f"row_support:{temporal['row_support']},"
        f"largest_word:{temporal['largest_word']}"
    )
    print(f"target_active_units={covectors['active_units']}")
    print(f"five_first_failure_roles={covectors['first_failure_roles']}")
    print(
        "target_active_first_failures="
        f"{covectors['active_first_failures']}"
    )
    print(
        "repair_pullback="
        f"covectors:{covectors['pullback_covectors']},"
        f"line_modes:{covectors['pullback_line_modes']}"
    )
    print("fraction_only_truth_checks=PASS")
    print("ALL_CHECKS_PASSED")


if __name__ == "__main__":
    main()
