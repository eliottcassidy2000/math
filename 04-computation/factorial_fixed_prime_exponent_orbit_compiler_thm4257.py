#!/usr/bin/env python3
"""Primary exact companion for the fixed-prime exponent-orbit compiler.

This file is intentionally self-contained and uses only Python's standard
library.  Every mathematical gate is optimization-stable: no bare ``assert``
statement carries truth.
"""

from math import isqrt


def require(condition, message):
    if not condition:
        raise RuntimeError("THM-ORBIT PRIMARY FAIL: " + message)


def is_prime(n):
    if n < 2:
        return False
    for d in range(2, isqrt(n) + 1):
        if n % d == 0:
            return False
    return True


def valuation_two(n):
    require(n > 0, "valuation_two received a nonpositive integer")
    value = 0
    while n % 2 == 0:
        value += 1
        n //= 2
    return value


def overlap_masks(b, ell, residue):
    modulus = 1 << ell
    return tuple(
        ((s * residue) % modulus) & (((b - s) * residue) % modulus)
        for s in range(1, (b - 1) // 2 + 1)
    )


def closes(b, ell, residue):
    modulus = 1 << ell
    residue %= modulus
    return residue % 2 == 1 and all(overlap_masks(b, ell, residue))


def closure_residues(b, ell):
    modulus = 1 << ell
    return tuple(r for r in range(1, modulus, 2) if closes(b, ell, r))


def direct_order(unit, ell):
    modulus = 1 << ell
    unit %= modulus
    require(unit % 2 == 1, "direct_order requires an odd unit")
    value = 1
    for order in range(1, modulus + 1):
        value = (value * unit) % modulus
        if value == 1:
            return order
    raise RuntimeError("THM-ORBIT PRIMARY FAIL: order search escaped unit group")


def formula_order(unit, ell):
    modulus = 1 << ell
    unit %= modulus
    require(unit % 2 == 1, "formula_order requires an odd unit")
    if ell == 1 or unit == 1:
        return 1
    if unit % 4 == 1:
        height = valuation_two(unit - 1)
        return 1 if ell <= height else 1 << (ell - height)
    height = valuation_two(unit * unit - 1)
    exponent = max(1, ell - height + 1)
    return 1 << exponent


def coordinate_tables(ell):
    require(ell >= 3, "C2 x cyclic coordinates require ell at least three")
    modulus = 1 << ell
    cyclic_order = 1 << (ell - 2)
    logs = {}
    value = 1
    for exponent in range(cyclic_order):
        require(value not in logs, "power of five repeated prematurely")
        logs[value] = exponent
        value = value * 5 % modulus
    require(value == 1, "power of five failed to close at predicted order")
    require(len(logs) == cyclic_order, "power-of-five subgroup has wrong size")
    return logs


def unit_coordinates(unit, ell, logs=None):
    modulus = 1 << ell
    unit %= modulus
    require(unit % 2 == 1, "unit_coordinates requires an odd unit")
    if logs is None:
        logs = coordinate_tables(ell)
    sign = 0 if unit % 4 == 1 else 1
    positive_part = unit if sign == 0 else (-unit) % modulus
    require(positive_part in logs, "unit missing from signed power-of-five chart")
    return sign, logs[positive_part]


def unit_from_coordinates(sign, exponent, ell):
    modulus = 1 << ell
    value = pow(5, exponent, modulus)
    return value if sign % 2 == 0 else (-value) % modulus


def exponent_classes(b, ell, unit):
    modulus = 1 << ell
    unit %= modulus
    order = direct_order(unit, ell)
    return tuple(k for k in range(order) if closes(b, ell, pow(unit, k, modulus)))


def allowed_child_bits(b, ell, parent):
    """Return exactly the epsilon in {0,1} whose residue child closes."""
    modulus = 1 << ell
    parent %= modulus
    failed = []
    for s in range(1, (b - 1) // 2 + 1):
        left = s * parent % modulus
        right = (b - s) * parent % modulus
        if left & right == 0:
            failed.append(s)
    if not failed:
        return (0, 1)
    demanded = set()
    for s in failed:
        even = s if s % 2 == 0 else b - s
        odd = b - s if s % 2 == 0 else s
        even_high = (even * parent // modulus) % 2
        odd_high = (odd * parent // modulus) % 2
        if even_high != 1:
            return ()
        demanded.add(1 ^ odd_high)
    return tuple(sorted(demanded)) if len(demanded) == 1 else ()


def positive_divisors_of_power_two(order):
    require(order >= 1 and order & (order - 1) == 0, "order is not a power of two")
    values = []
    divisor = 1
    while divisor <= order:
        values.append(divisor)
        divisor *= 2
    return tuple(values)


def minimal_membership_period(classes, order):
    chosen = set(classes)
    for period in positive_divisors_of_power_two(order):
        if all((k in chosen) == ((k + period) % order in chosen) for k in range(order)):
            return period
    raise RuntimeError("THM-ORBIT PRIMARY FAIL: no membership period found")


def maximal_exponent_cosets(classes, order):
    chosen = set(classes)
    maximal = []
    for modulus in positive_divisors_of_power_two(order):
        for residue in range(modulus):
            candidate = {k for k in range(order) if k % modulus == residue}
            if not candidate or not candidate <= chosen:
                continue
            contained = any(
                modulus % old_modulus == 0 and residue % old_modulus == old_residue
                for old_modulus, old_residue in maximal
            )
            if not contained:
                maximal.append((modulus, residue))
    return tuple(maximal)


def negative_collar(b, ell):
    modulus = 1 << ell
    upper = modulus // (2 * (b - 1))
    return tuple((-w) % modulus for w in range(1, upper + 1, 2))


def coordinate_orbit_contains(unit, target, ell, logs=None):
    """Solve the exact two coordinate congruences by a finite exponent cycle."""
    if logs is None:
        logs = coordinate_tables(ell)
    sign_u, log_u = unit_coordinates(unit, ell, logs)
    sign_t, log_t = unit_coordinates(target, ell, logs)
    cyclic_order = 1 << (ell - 2)
    order = direct_order(unit, ell)
    return any(
        sign_u * k % 2 == sign_t and log_u * k % cyclic_order == log_t
        for k in range(order)
    )


def least_nonempty_level(b, maximum=16):
    for ell in range(1, maximum + 1):
        if closure_residues(b, ell):
            return ell
    raise RuntimeError("THM-ORBIT PRIMARY FAIL: least level search exhausted")


def least_universal_level(b):
    ell = 1
    while 1 << (ell - 1) < b:
        ell += 1
    return ell


def format_classes(classes):
    return "{" + ",".join(str(value) for value in classes) + "}"


def run_group_audit():
    coordinate_cells = 0
    order_cells = 0
    minus_one_cells = 0
    for ell in range(3, 12):
        modulus = 1 << ell
        logs = coordinate_tables(ell)
        require(set(logs) == {r for r in range(1, modulus, 4)},
                f"power-of-five chart mismatch at ell={ell}")
        for unit in range(1, modulus, 2):
            sign, exponent = unit_coordinates(unit, ell, logs)
            rebuilt = unit_from_coordinates(sign, exponent, ell)
            require(rebuilt == unit, f"coordinate round trip failed at ell={ell},u={unit}")
            coordinate_cells += 1
            direct = direct_order(unit, ell)
            formula = formula_order(unit, ell)
            require(direct == formula, f"order formula failed at ell={ell},u={unit}")
            order_cells += 1
            orbit = {pow(unit, k, modulus) for k in range(direct)}
            require((modulus - 1 in orbit) == (unit == modulus - 1),
                    f"minus-one orbit criterion failed at ell={ell},u={unit}")
            minus_one_cells += 1
    return coordinate_cells, order_cells, minus_one_cells


def run_orbit_and_lift_audit():
    iff_cells = 0
    lift_parent_cells = 0
    stable_levels = 0
    doubled_levels = 0
    reached_child_histogram = {0: 0, 1: 0, 2: 0}
    for b in range(3, 32, 2):
        for ell in range(2, 10):
            modulus = 1 << ell
            closure = set(closure_residues(b, ell))
            for unit in range(1, modulus, 2):
                order = direct_order(unit, ell)
                classes = set(exponent_classes(b, ell, unit))
                for k in range(order):
                    residue = pow(unit, k, modulus)
                    require((k in classes) == (residue in closure),
                            f"orbit iff failed at b={b},ell={ell},u={unit},k={k}")
                    iff_cells += 1
        for ell in range(2, 9):
            modulus = 1 << ell
            double_modulus = 2 * modulus
            for unit in range(1, double_modulus, 2):
                parent_order = direct_order(unit, ell)
                child_order = direct_order(unit, ell + 1)
                ratio = child_order // parent_order
                require(child_order in (parent_order, 2 * parent_order),
                        f"order lift ratio failed at ell={ell},u={unit}")
                if ratio == 1:
                    stable_levels += 1
                else:
                    doubled_levels += 1
                for k in range(parent_order):
                    parent = pow(unit, k, modulus)
                    allowed = set(allowed_child_bits(b, ell, parent))
                    reached = []
                    accepted = 0
                    for sheet in range(ratio):
                        exponent = k + sheet * parent_order
                        child = pow(unit, exponent, double_modulus)
                        epsilon = (child - parent) // modulus
                        require(epsilon in (0, 1), "child is not a binary lift")
                        reached.append(epsilon)
                        predicted = epsilon in allowed
                        actual = closes(b, ell + 1, child)
                        require(predicted == actual,
                                f"class lift failed at b={b},ell={ell},u={unit},k={k},sheet={sheet}")
                        accepted += int(actual)
                    require(len(set(reached)) == ratio, "reached lift sheets collided")
                    reached_child_histogram[accepted] += 1
                    lift_parent_cells += 1
    return iff_cells, lift_parent_cells, stable_levels, doubled_levels, reached_child_histogram


def run_collar_and_coset_audit():
    collar_residue_cells = 0
    collar_orbit_cells = 0
    collar_hits = 0
    collar_misses = 0
    coset_cells = 0
    for b in range(3, 32, 2):
        for ell in range(3, 10):
            modulus = 1 << ell
            logs = coordinate_tables(ell)
            collar = set(negative_collar(b, ell))
            for target in collar:
                require(closes(b, ell, target),
                        f"negative collar failed at b={b},ell={ell},r={target}")
                collar_residue_cells += 1
            for unit in range(1, modulus, 2):
                order = direct_order(unit, ell)
                orbit = {pow(unit, k, modulus) for k in range(order)}
                direct_hit = bool(orbit & collar)
                coordinate_hit = any(
                    coordinate_orbit_contains(unit, target, ell, logs) for target in collar
                )
                require(direct_hit == coordinate_hit,
                        f"collar congruence criterion failed at b={b},ell={ell},u={unit}")
                collar_hits += int(direct_hit)
                collar_misses += int(not direct_hit)
                collar_orbit_cells += 1

                classes = set(exponent_classes(b, ell, unit))
                closure = set(closure_residues(b, ell))
                for divisor in positive_divisors_of_power_two(order):
                    subgroup = {pow(unit, divisor * j, modulus) for j in range(order // divisor)}
                    for residue in range(divisor):
                        exponent_coset = {k for k in range(order) if k % divisor == residue}
                        image_coset = {pow(unit, residue, modulus) * h % modulus for h in subgroup}
                        require(
                            (exponent_coset <= classes) == (image_coset <= closure),
                            f"coset/AP equivalence failed at b={b},ell={ell},u={unit},d={divisor},a={residue}",
                        )
                        coset_cells += 1
    return (
        collar_residue_cells,
        collar_orbit_cells,
        collar_hits,
        collar_misses,
        coset_cells,
    )


def run_atlas():
    expected_least = {3: 2, 5: 3, 7: 3, 9: 4, 11: 4, 13: 4}
    expected_residues = {
        (3, 2): (3,),
        (5, 3): (5, 7),
        (7, 3): (7,),
        (7, 4): (3, 5, 7, 15),
        (9, 4): (9, 15),
        (9, 5): (9, 15, 19, 25, 29, 31),
        (11, 4): (3, 15),
        (11, 5): (3, 7, 9, 15, 19, 21, 27, 31),
        (13, 4): (5, 15),
        (13, 5): (5, 15, 21, 31),
    }
    for b, expected in expected_least.items():
        require(least_nonempty_level(b) == expected, f"least nonempty level changed for b={b}")
    for key, expected in expected_residues.items():
        require(closure_residues(*key) == expected, f"closure atlas changed at {key}")

    print("LEAST_LEVELS")
    for b in expected_least:
        print(
            f"b={b} least_nonempty={least_nonempty_level(b)} "
            f"least_universal={least_universal_level(b)}"
        )

    print("EXACT_ATLAS")
    atlas_specs = ((3, 2), (5, 3), (7, 3), (7, 4),
                   (9, 4), (9, 5), (11, 4), (11, 5), (13, 4), (13, 5))
    atlas_rows = 0
    empty_rows = 0
    for b, ell in atlas_specs:
        modulus = 1 << ell
        closure = closure_residues(b, ell)
        print(f"b={b} ell={ell} R={format_classes(closure)}")
        for unit in range(1, modulus, 2):
            order = direct_order(unit, ell)
            classes = exponent_classes(b, ell, unit)
            period = minimal_membership_period(classes, order)
            cosets = maximal_exponent_cosets(classes, order)
            coset_text = "{" + ",".join(f"{a}mod{d}" for d, a in cosets) + "}"
            print(
                f"  u={unit:>{len(str(modulus - 1))}} ord={order:<2} "
                f"K={format_classes(classes):<18} period={period:<2} AP={coset_text}"
            )
            atlas_rows += 1
            empty_rows += int(not classes)

    seed_masks = {
        (9, 4, 9): (8, 2, 2, 4),
        (11, 4, 3): (2, 2, 8, 4, 2),
        (13, 4, 5): (4, 2, 2, 4, 8, 2),
    }
    for key, expected in seed_masks.items():
        require(overlap_masks(*key) == expected, f"seed mask changed at {key}")

    spot_rows = {
        (9, 5, 19): (8, (1, 2, 6)),
        (11, 5, 23): (4, (3,)),
        (13, 5, 29): (8, (3, 7)),
    }
    print("FIXED_PRIME_FAMILIES")
    for (b, ell, prime), (expected_order, expected_classes) in spot_rows.items():
        require(is_prime(prime), f"spot base {prime} is not prime")
        require(prime > 2 * b, f"spot prime {prime} misses multiplier 2b")
        order = direct_order(prime, ell)
        classes = exponent_classes(b, ell, prime)
        require((order, classes) == (expected_order, expected_classes),
                f"spot family changed at b={b},p={prime}")
        cosets = maximal_exponent_cosets(classes, order)
        print(
            f"a={2*b} p={prime} ell={ell} ord={order} "
            f"K={format_classes(classes)} AP="
            + "{" + ",".join(f"k={a}mod{d}" for d, a in cosets) + "}"
        )

    collar_b, collar_ell, collar_prime, collar_exponent = 13, 7, 43, 21
    collar_modulus = 1 << collar_ell
    collar_target = pow(collar_prime, collar_exponent, collar_modulus)
    collar_logs = coordinate_tables(collar_ell)
    require(collar_target == 123, "structural collar target changed")
    require(collar_target in negative_collar(collar_b, collar_ell),
            "structural collar target left the negative collar")
    require(unit_coordinates(collar_prime, collar_ell, collar_logs) == (1, 29),
            "structural collar prime coordinates changed")
    require(unit_coordinates(collar_target, collar_ell, collar_logs) == (1, 1),
            "structural collar target coordinates changed")
    require(direct_order(collar_prime, collar_ell) == 32,
            "structural collar prime order changed")
    require(collar_exponent in exponent_classes(collar_b, collar_ell, collar_prime),
            "structural collar exponent is not certified")
    print(
        "COLLAR_ORBIT_FAMILY "
        "a=26 p=43 ell=7 p_coords=(-1,29) collar=-5_coords=(-1,1) "
        "congruences=k_odd,29k=1mod32 solution=k=21mod32"
    )

    hostile_b, hostile_ell, hostile_prime = 9, 4, 23
    hostile_classes = exponent_classes(hostile_b, hostile_ell, hostile_prime)
    require(hostile_classes == (), "empty-orbit hostile acquired a suffix class")
    low_masks = tuple((s * hostile_prime) & ((hostile_b - s) * hostile_prime)
                      for s in range(1, (hostile_b - 1) // 2 + 1))
    high_height = hostile_prime ** 3
    high_masks = tuple((s * high_height) & ((hostile_b - s) * high_height)
                       for s in range(1, (hostile_b - 1) // 2 + 1))
    require(not all(low_masks), "hostile low height unexpectedly closes")
    require(high_masks == (11264, 19456, 3072, 44032), "hostile repair masks changed")
    require(all(high_masks), "hostile high-bit repair failed")
    print(
        "EMPTY_ORBIT_HOSTILE "
        f"b={hostile_b} ell={hostile_ell} p={hostile_prime} K=empty "
        f"k=1_masks={low_masks} k=3_masks={high_masks}"
    )
    return atlas_rows, empty_rows


def main():
    print("FIXED-PRIME EXPONENT-ORBIT COMPILER PRIMARY")
    print("method=direct_bit_products+cyclic_orbits+signed_power_of_five_coordinates")
    coordinate_cells, order_cells, minus_one_cells = run_group_audit()
    iff_cells, lift_cells, stable, doubled, histogram = run_orbit_and_lift_audit()
    collar_cells, collar_orbits, collar_hits, collar_misses, coset_cells = run_collar_and_coset_audit()
    atlas_rows, empty_rows = run_atlas()
    print("AUDIT_COUNTS")
    print(
        f"coordinate_cells={coordinate_cells} order_cells={order_cells} "
        f"minus_one_cells={minus_one_cells}"
    )
    print(
        f"orbit_iff_cells={iff_cells} lift_parent_cells={lift_cells} "
        f"stable_order_rows={stable} doubled_order_rows={doubled} "
        f"accepted_child_histogram={histogram}"
    )
    print(
        f"collar_residue_cells={collar_cells} collar_orbit_cells={collar_orbits} "
        f"collar_hits={collar_hits} collar_misses={collar_misses} "
        f"coset_ap_cells={coset_cells}"
    )
    print(f"atlas_rows={atlas_rows} empty_orbit_rows={empty_rows}")
    print("scope=suffix_certificate_iff_only; suffix_failure_is_not_factorial_failure")
    print("THM-ORBIT PRIMARY PASS")


if __name__ == "__main__":
    main()
