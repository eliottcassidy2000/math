#!/usr/bin/env python3
"""Exact controls for THM-2522.

All arithmetic is integral or ``fractions.Fraction`` arithmetic.  Equal-grid
step vectors represent functions which are constant on the corresponding
half-open circle cells; midpoint evaluations therefore give exact integrals
because every breakpoint used below lies on the grid.
"""

from fractions import Fraction as F


P = 13


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def average(values):
    return sum(values, F(0)) / len(values)


def inner(left, right):
    require(len(left) == len(right), "inner-product dimensions")
    return sum((F(x) * F(y) for x, y in zip(left, right)), F(0)) / len(left)


def perron13(values):
    """Exact P_13 on a step grid whose length is divisible by thirteen."""
    if len(values) == 1:
        return list(map(F, values))
    require(len(values) % P == 0, "Perron grid must have thirteen fibres")
    base = len(values) // P
    return [sum((F(values[k + r * base]) for r in range(P)), F(0)) / P
            for k in range(base)]


def koopman13(values):
    """Exact U f=f(13x), on the once-refined grid."""
    values = list(map(F, values))
    return [values[j % len(values)] for j in range(P * len(values))]


def jumps(values):
    values = list(map(F, values))
    return [values[j] - values[j - 1] for j in range(len(values))]


def replicate_current(current):
    """Unnormalised inverse-image replication of an atomic jump current."""
    current = list(map(F, current))
    base = len(current)
    answer = [F(0) for _ in range(P * base)]
    for k, value in enumerate(current):
        for r in range(P):
            answer[k + r * base] = value
    return answer


def base_tile_value(level, j, grid):
    """Value of 1_[0,1/13)({13^level x}) on cell j of ``grid``."""
    period = grid // (P**level)
    require(period % P == 0, "tile represented below its natural grid")
    return int((j % period) < period // P)


def tower_values(level_weights, highest_level):
    grid = P ** (highest_level + 1)
    return [sum(weight * base_tile_value(level, j, grid)
                for level, weight in level_weights.items())
            for j in range(grid)]


def lift_to_grid(values, levels):
    answer = list(map(F, values))
    for _ in range(levels):
        answer = koopman13(answer)
    return answer


def collision_table(values, level):
    """Direct depth-(level+1) fixed-last-digit table on one fine grid."""
    grid = len(values)
    stalk = P**level
    cover = P * stalk
    require(grid % cover == 0, "collision shift must lie on the fine grid")
    answer = []
    for u in range(P):
        total = F(0)
        for e in range(stalk):
            shift = grid * (u + P * e) // cover
            total += sum((F(values[j]) * F(values[(j + shift) % grid])
                          for j in range(grid)), F(0))
        answer.append(total / (stalk * grid))
    return answer


def reduced_collision_table(level_values):
    grid = len(level_values)
    require(grid % P == 0, "last-digit shift grid")
    return [sum((F(level_values[j])
                 * F(level_values[(j + u * grid // P) % grid])
                 for j in range(grid)), F(0)) / grid
            for u in range(P)]


def tower_controls():
    # q+2(q o T^2) has active digits exactly at levels 0 and 2.  This is
    # the smallest useful hostile to any false monotonicity assertion.
    values = tower_values({0: 1, 2: 2}, highest_level=2)
    levels = [list(map(F, values))]
    while len(levels[-1]) > 1:
        levels.append(perron13(levels[-1]))

    digits = []
    for m in range(len(levels) - 1):
        replicated = koopman13(levels[m + 1])
        digits.append([x - y for x, y in zip(levels[m], replicated)])
    digits.append([F(0)])
    active = [any(value != 0 for value in digit) for digit in digits]
    require(active == [True, False, True, False], "gapped digit spectrum")
    require(all(average(level) == F(3, 13) for level in levels),
            "Perron mean preservation")

    # Endpoint-current descent is exactly current minus replicated child.
    for m in range(3):
        defect = [x - y for x, y in zip(
            jumps(levels[m]), replicate_current(jumps(levels[m + 1])))]
        require(defect == jumps(digits[m]), "replicated endpoint-current recursion")

    # Function recursion, orthogonality, and Parseval telescope.
    reconstructed = lift_to_grid(levels[-1], 3)
    lifted_digits = []
    for m, digit in enumerate(digits[:3]):
        lifted = lift_to_grid(digit, m)
        lifted_digits.append(lifted)
        reconstructed = [x + y for x, y in zip(reconstructed, lifted)]
    require(reconstructed == levels[0], "finite Perron-Wold reconstruction")
    for i in range(len(lifted_digits)):
        for j in range(i + 1, len(lifted_digits)):
            require(inner(lifted_digits[i], lifted_digits[j]) == 0,
                    "distinct valuation digits are orthogonal")
    centered_energy = average([(x - F(3, 13)) ** 2 for x in levels[0]])
    digit_energy = sum((average([x * x for x in digit])
                        for digit in digits[:3]), F(0))
    require(centered_energy == digit_energy, "orthogonal digit energy")

    # A twice-replicated tile has m_*=2, hence first collision depth L=3.
    pure = tower_values({2: 1}, highest_level=2)
    pure_levels = [list(map(F, pure))]
    while len(pure_levels[-1]) > 1:
        pure_levels.append(perron13(pure_levels[-1]))
    pure_digits = [
        [x - y for x, y in zip(pure_levels[m], koopman13(pure_levels[m + 1]))]
        for m in range(3)
    ]
    pure_active = [any(x != 0 for x in digit) for digit in pure_digits]
    require(pure_active == [False, False, True], "sharp first valuation level")

    collision_energies = []
    for m in range(3):
        direct = collision_table(pure, m)
        reduced = reduced_collision_table(pure_levels[m])
        require(direct == reduced, "exact toothpick reduction")
        energy = direct[0] - average(direct)
        require(energy == average([x * x for x in pure_digits[m]]),
                "collision energy equals digit norm")
        require((energy > 0) == pure_active[m], "collision/valuation equivalence")
        collision_energies.append(energy)

    return active, pure_active, collision_energies, centered_energy


def danger(n, numerator, denominator):
    residue = (n * numerator) % denominator
    distance_numerator = min(residue, denominator - residue)
    return 14 * distance_numerator < denominator


def danger_values(n, grid):
    # Midpoints avoid the strict danger boundary.
    return [int(danger(n, 2 * j + 1, 2 * grid)) for j in range(grid)]


def guard(n, numerator, denominator, radius_denominator=7):
    residue = (n * numerator) % denominator
    distance_numerator = min(residue, denominator - residue)
    return radius_denominator * distance_numerator > denominator


def shallow_support_controls():
    unit = 5
    coarse_grid = 14 * unit * P
    fine_grid = P * coarse_grid
    original = danger_values(P * unit, fine_grid)
    shallow = danger_values(unit, coarse_grid)
    require(perron13(original) == list(map(F, shallow)),
            "P_13 D_(13a)=D_a")
    require(original == list(map(int, koopman13(shallow))),
            "standard danger comb is once replicated")

    child = list(map(F, shallow))
    child_digit = [x - y for x, y in zip(child, koopman13(perron13(child)))]
    require(any(x != 0 for x in child_digit), "unit shallow support forces level-one digit")
    table = reduced_collision_table(child)
    energy = table[0] - average(table)
    require(energy == average([x * x for x in child_digit]) and energy > 0,
            "forced depth-two collision energy")
    require(any(value != table[0] for value in table[1:]),
            "rational phase bank is nonconstant")

    # The geometric kernel is empty: no 13-translation orbit lies in D_a.
    orbit_intersection_cells = 0
    shift = coarse_grid // P
    for j in range(coarse_grid):
        if all(shallow[(j + u * shift) % coarse_grid] for u in range(P)):
            orbit_intersection_cells += 1
    require(orbit_intersection_cells == 0, "empty translated-danger intersection")

    # Nonnegativity is load-bearing.  Two opposite branch weights supported
    # in D_(13a) form a nonzero signed function with zero Perron image.
    signed = [F(0) for _ in range(fine_grid)]
    for k, value in enumerate(shallow):
        signed[k] = F(value)
        signed[k + coarse_grid] = -F(value)
    require(any(signed) and all(value == 0 for value in perron13(signed)),
            "signed support hostile")
    for j, value in enumerate(signed):
        if value:
            require(danger(P * unit, 2 * j + 1, 2 * fine_grid),
                    "signed hostile remains in the danger carrier")

    # The unit condition is sharp: a=13 leaves D_a invariant by 1/13.
    nonunit = P
    nonunit_grid = 14 * nonunit * P
    nonunit_child = list(map(F, danger_values(nonunit, nonunit_grid)))
    nonunit_digit = [x - y for x, y in zip(
        nonunit_child, koopman13(perron13(nonunit_child)))]
    require(all(x == 0 for x in nonunit_digit), "nonunit shallow hostile")

    # The actual all-row event has a second carrier: the unit guard core
    # C_H={||Hx||>1/7}.  No complete 1/13 orbit fits in that core, so the
    # level-zero digit is forced.  Intersecting it with D_(13a) exhibits
    # consecutive active levels zero and one on one exact Boolean event.
    guard_unit = 1
    guarded = [int(value and guard(guard_unit, 2 * j + 1, 2 * fine_grid))
               for j, value in enumerate(original)]
    require(any(guarded), "guarded shallow event is nonempty")
    guarded_child = perron13(guarded)
    guarded_digit0 = [x - y for x, y in zip(
        map(F, guarded), koopman13(guarded_child))]
    guarded_grandchild = perron13(guarded_child)
    guarded_digit1 = [x - y for x, y in zip(
        guarded_child, koopman13(guarded_grandchild))]
    require(any(x != 0 for x in guarded_digit0), "unit guard forces level-zero digit")
    require(any(x != 0 for x in guarded_digit1), "shallow danger forces level-one digit")
    guarded_energies = []
    for m, digit in enumerate((guarded_digit0, guarded_digit1)):
        direct = collision_table(guarded, m)
        reduced = reduced_collision_table((list(map(F, guarded)), guarded_child)[m])
        require(direct == reduced, "guarded direct/reduced collision")
        energy_m = direct[0] - average(direct)
        require(energy_m == average([x * x for x in digit]) and energy_m > 0,
                "guarded consecutive digit energy")
        guarded_energies.append(energy_m)

    guard_shift = fine_grid // P
    guard_orbit_cells = 0
    for j in range(fine_grid):
        if all(guard(guard_unit, 2 * ((j + u * guard_shift) % fine_grid) + 1,
                     2 * fine_grid)
               for u in range(P)):
            guard_orbit_cells += 1
    require(guard_orbit_cells == 0, "unit guard has no complete root orbit")

    # Both guard hypotheses are sharp.  If 13 divides H, C_H is itself
    # 1/13-invariant.  If the excluded central radius is below 1/26, the
    # half-grid orbit (2u+1)/26 fits wholly in the carrier.
    nonunit_guard = [F(int(guard(P, 2 * j + 1, 2 * fine_grid)))
                     for j in range(fine_grid)]
    nonunit_guard_digit = [x - y for x, y in zip(
        nonunit_guard, koopman13(perron13(nonunit_guard)))]
    require(any(nonunit_guard) and all(x == 0 for x in nonunit_guard_digit),
            "nonunit guard hostile")
    require(all(52 * min((2 * u + 1) % 26, 26 - ((2 * u + 1) % 26)) > 26
                for u in range(P)), "small-radius half-grid hostile")

    return (energy, table, orbit_intersection_cells, guarded_energies,
            guard_orbit_cells)


def interval_zero(numerator, denominator):
    return 14 * (numerator % denominator) < denominator


def interval_top(numerator, denominator):
    return 14 * (numerator % denominator) >= 13 * denominator


def weighted_owner_table(delay):
    # F=1_[0,1/14), G=1_[13/14,1); collision depth L=1 is fixed.
    multiplier = P ** (delay + 1)
    grid = 14 * multiplier
    denominator = 2 * grid
    counts = []
    for u in range(P):
        count = 0
        for j in range(grid):
            numerator = 2 * j + 1
            if not interval_zero(numerator, denominator):
                continue
            shifted_num = P * numerator + u * denominator
            shifted_den = P * denominator
            if not interval_zero(shifted_num, shifted_den):
                continue
            if interval_top(multiplier * numerator, denominator):
                count += 1
        counts.append(F(count, grid))
    return counts


def unweighted_interval_table():
    grid = 14 * P
    denominator = 2 * grid
    answer = []
    for u in range(P):
        count = 0
        for j in range(grid):
            numerator = 2 * j + 1
            shifted_num = P * numerator + u * denominator
            shifted_den = P * denominator
            if interval_zero(numerator, denominator) and interval_zero(shifted_num, shifted_den):
                count += 1
        answer.append(F(count, grid))
    return answer


def owner_delay_controls():
    unweighted = unweighted_interval_table()
    drift = unweighted[0] - average(unweighted)
    require(drift > 0, "unweighted collision digit")
    rho = F(1, 14)
    first_positive = None
    rows = []
    for delay in range(4):
        table = weighted_owner_table(delay)
        weighted_drift = table[0] - average(table)
        multiplier = P ** (delay + 1)
        # B=1, Var(F)=Var(G)=2 gives the aggregate covariance bound 2/(3K).
        require(abs(weighted_drift - rho * drift) <= F(2, 3 * multiplier),
                "late-owner BV covariance")
        nonconstant = any(value != table[0] for value in table[1:])
        if nonconstant and first_positive is None:
            first_positive = delay
        rows.append((delay, weighted_drift, nonconstant))
    require(all(value == 0 for value in weighted_owner_table(0)),
            "undelayed owner hostile")
    require(first_positive is not None and rows[-1][2], "eventual owner landing")
    return drift, first_positive, rows


def conductor_horizon_controls():
    # A sparse coprime current at endpoint conductor 13^2*5.
    nu = 2
    modulus = 5
    current = [F(1), F(-1), F(0), F(0), F(0)]
    order = 4  # ord_5(13)
    labelled = []
    recursion_checks = 0
    crt_checks = 0
    previous = None
    for m in range(nu, nu + 2 * order + 1):
        unit = pow(P, m - nu, modulus)
        inv_unit = pow(unit, -1, modulus)
        mu = [current[(inv_unit * v) % modulus] / (P**m)
              for v in range(modulus)]
        labelled.append([value * (P**m) for value in mu])
        if previous is not None:
            inv_p = pow(P, -1, modulus)
            expected = [previous[(inv_p * v) % modulus] / P
                        for v in range(modulus)]
            require(mu == expected, "late current recursion")
            recursion_checks += modulus
        previous = mu

        # Every nonzero last digit meets every coprime-current residue by CRT.
        cover = P * modulus
        for a in range(1, P):
            residues = []
            for b in range(modulus):
                candidates = [r for r in range(1, cover)
                              if r % P == a and (unit * r) % modulus == b]
                require(len(candidates) == 1, "CRT ladder address")
                residues.append(candidates[0])
                crt_checks += 1
            require(len(set(residues)) == modulus, "complete phase ladder")

    require(labelled[:order] == labelled[order:2 * order],
            "phase-bank labels have ord_d(13) horizon period")
    require(sum(value * value for value in current) == 2,
            "coprime current energy control")
    return order, recursion_checks, crt_checks, labelled[:order]


def main():
    active, pure_active, collision_energies, centered_energy = tower_controls()
    (shallow_energy, shallow_table, empty_cells, guarded_energies,
     guard_empty_cells) = shallow_support_controls()
    owner_drift, first_owner_delay, owner_rows = owner_delay_controls()
    period, recursion_checks, crt_checks, labels = conductor_horizon_controls()

    print("THM-2522 intrinsic collision-depth companion: PASS")
    print(f"gapped_digit_pattern={active}; sharp_first_pattern={pure_active}")
    print("collision_energies=" + ",".join(map(str, collision_energies)))
    print(f"orthogonal_centered_energy={centered_energy}")
    print(f"depth_two_shallow_energy={shallow_energy}; empty_orbit_cells={empty_cells}")
    print("depth_two_collision_table=" + ",".join(map(str, shallow_table)))
    print("guarded_depth_one_two_energies=" + ",".join(map(str, guarded_energies))
          + f"; guard_orbit_cells={guard_empty_cells}")
    print(f"unweighted_owner_control_drift={owner_drift}; first_positive_owner_delay={first_owner_delay}")
    print("owner_delay_rows=" + ";".join(
        f"R{delay}:D={drift}:bank={int(bank)}" for delay, drift, bank in owner_rows))
    print(f"horizon_period={period}; current_recursion_checks={recursion_checks}; crt_checks={crt_checks}")
    print("horizon_labels=" + ";".join(
        ",".join(map(str, row)) for row in labels))
    print("collision depth is m+1; late owner delay changes no ancestry address")


if __name__ == "__main__":
    main()
