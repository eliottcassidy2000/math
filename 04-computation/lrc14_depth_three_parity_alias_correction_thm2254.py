#!/usr/bin/env python3
"""Exact audit of the unit-13 parity alias in the all-equal (3,4,5) branch.

The 169-step kernel has two stochastic square roots related by swapping the
two equal-mass annuli A and B.  Exact 13-root counting selects the root with
centered eigenvalue -1/13 in both centered directions.  The other root gives
the two formerly reported (but invalid) hostile-control closures.
"""

from fractions import Fraction
from itertools import product


Q = Fraction
TARGET = Q(961, 6930)
PI = (Q(1, 7), Q(1, 7), Q(5, 7))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def state(x):
    """Return A=D_1, B=E_1\\D_1, or C=E_1^c off null endpoints."""
    x %= 1
    distance = min(x, 1 - x)
    if distance < Q(1, 14):
        return 0
    if distance < Q(1, 7):
        return 1
    return 2


def guard(x):
    """Indicator of C_1={x: ||x||>1/7}, off null endpoints."""
    x %= 1
    return int(min(x, 1 - x) > Q(1, 7))


def derive_transition():
    """Derive the 13-step state law by exact midpoint-cell enumeration."""
    cell_count = 14 * 13
    counts = [[0] * 3 for _ in range(3)]
    for index in range(cell_count):
        z = Q(2 * index + 1, 2 * cell_count)
        counts[state(z)][state(13 * z)] += 1
    transition = tuple(
        tuple(Q(value, sum(row)) for value in row) for row in counts
    )
    return tuple(tuple(row) for row in counts), transition


def derive_reverse_root_counts():
    """Count the states of all 13 roots, pointwise on each future state."""
    values = [set() for _ in range(3)]
    for index in range(14):
        future = Q(2 * index + 1, 28)
        counts = [0, 0, 0]
        for root in range(13):
            counts[state((future + root) / 13)] += 1
        values[state(future)].add(tuple(counts))
    require(
        all(len(bucket) == 1 for bucket in values),
        "reverse root law is not constant on a coarse future state",
    )
    return tuple(next(iter(bucket)) for bucket in values)


def derive_guard_weights(multiplier):
    """Compute E[1_{C_H} | ux=z] for u=multiplier*H, multiplier 1 or 6."""
    values = [set() for _ in range(3)]
    for index in range(14):
        z = Q(2 * index + 1, 28)
        weight = Q(
            sum(guard((z + root) / multiplier) for root in range(multiplier)),
            multiplier,
        )
        values[state(z)].add(weight)
    require(
        all(len(bucket) == 1 for bucket in values),
        f"guard weight is not state-measurable for multiplier {multiplier}",
    )
    return tuple(next(iter(bucket)) for bucket in values)


def multiply(left, right):
    return tuple(
        tuple(
            sum(left[i][k] * right[k][j] for k in range(len(right)))
            for j in range(len(right[0]))
        )
        for i in range(len(left))
    )


def matvec(matrix, vector):
    return tuple(
        sum(matrix[i][j] * vector[j] for j in range(len(vector)))
        for i in range(len(matrix))
    )


def rowvec(vector, matrix):
    return tuple(
        sum(vector[i] * matrix[i][j] for i in range(len(vector)))
        for j in range(len(matrix[0]))
    )


def cylinder_value(transition, guard_weights):
    """Evaluate E[G W_1 W_3 (W_2+(1-W_0)/169)] exactly."""
    total = Q()
    leading = Q()
    repair = Q()
    law_mass = Q()
    for word in product(range(3), repeat=6):
        probability = PI[word[0]]
        for left, right in zip(word, word[1:]):
            probability *= transition[left][right]
        law_mass += probability
        danger = tuple(int(value == 0) for value in word)
        windows = tuple(int(any(danger[t : t + 3])) for t in range(4))
        positive = windows[1] * windows[3]
        leading_atom = positive * windows[2]
        repair_atom = positive * (1 - windows[0])
        weight = guard_weights[word[0]]
        leading += probability * weight * leading_atom
        repair += probability * weight * repair_atom
        total += probability * weight * (
            leading_atom + Q(repair_atom, 169)
        )
    require(law_mass == 1, "six-state word law is not normalized")
    return total, leading, repair


def direct_grid_cylinders(weight_specs):
    """Independently integrate the length-six cylinder on its exact grid."""
    cell_count = 14 * 13**5
    modulus = 2 * cell_count
    accumulators = [[0, 0] for _ in weight_specs]
    for index in range(cell_count):
        residue = 2 * index + 1
        dangers = []
        initial_state = None
        for step in range(6):
            distance_numerator = min(residue, modulus - residue)
            if step == 0:
                if 14 * distance_numerator < modulus:
                    initial_state = 0
                elif 7 * distance_numerator < modulus:
                    initial_state = 1
                else:
                    initial_state = 2
            dangers.append(int(14 * distance_numerator < modulus))
            residue = (13 * residue) % modulus
        windows = tuple(int(any(dangers[t : t + 3])) for t in range(4))
        positive = windows[1] * windows[3]
        leading_atom = positive * windows[2]
        repair_atom = positive * (1 - windows[0])
        for accumulator, (weight_numerators, _) in zip(
            accumulators, weight_specs
        ):
            guard_numerator = weight_numerators[initial_state]
            accumulator[0] += guard_numerator * leading_atom
            accumulator[1] += guard_numerator * repair_atom
    values = []
    for (leading, repair), (_, weight_denominator) in zip(
        accumulators, weight_specs
    ):
        values.append(
            Q(
                169 * leading + repair,
                weight_denominator * 169 * cell_count,
            )
        )
    return cell_count, tuple(tuple(pair) for pair in accumulators), tuple(values)


def format_matrix(matrix):
    return "; ".join(
        "(" + ",".join(str(value) for value in row) + ")" for row in matrix
    )


def main():
    counts, correct = derive_transition()
    reverse_counts = derive_reverse_root_counts()
    require(
        counts == ((2, 4, 20), (4, 2, 20), (20, 20, 90)),
        "wrong exact 182-cell transition counts",
    )
    require(
        reverse_counts == ((1, 2, 10), (2, 1, 10), (2, 2, 9)),
        "wrong pointwise reverse 13-root counts",
    )
    require(
        correct
        == (
            (Q(1, 13), Q(2, 13), Q(10, 13)),
            (Q(2, 13), Q(1, 13), Q(10, 13)),
            (Q(2, 13), Q(2, 13), Q(9, 13)),
        ),
        "wrong normalized 13-step transition",
    )
    require(rowvec(PI, correct) == PI, "stationary law failed")

    identity = (
        (Q(1), Q(0), Q(0)),
        (Q(0), Q(1), Q(0)),
        (Q(0), Q(0), Q(1)),
    )
    projector = (PI, PI, PI)
    rho = Q(-1, 13)
    spectral_form = tuple(
        tuple(
            rho * identity[i][j] + (1 - rho) * projector[i][j]
            for j in range(3)
        )
        for i in range(3)
    )
    require(correct == spectral_form, "rank-one spectral form failed")

    involution = (
        (Q(0), Q(1), Q(0)),
        (Q(1), Q(0), Q(0)),
        (Q(0), Q(0), Q(1)),
    )
    twisted = multiply(involution, correct)
    require(twisted == multiply(correct, involution), "J and P do not commute")
    require(
        twisted
        == (
            (Q(2, 13), Q(1, 13), Q(10, 13)),
            (Q(1, 13), Q(2, 13), Q(10, 13)),
            (Q(2, 13), Q(2, 13), Q(9, 13)),
        ),
        "wrong parity-twisted root",
    )
    expected_square = (
        (Q(25, 169), Q(24, 169), Q(120, 169)),
        (Q(24, 169), Q(25, 169), Q(120, 169)),
        (Q(24, 169), Q(24, 169), Q(121, 169)),
    )
    require(multiply(correct, correct) == expected_square, "P^2 failed")
    require(multiply(twisted, twisted) == expected_square, "(JP)^2 failed")

    antisymmetric = (Q(1), Q(-1), Q(0))
    symmetric_centered = (Q(5), Q(5), Q(-2))
    for vector in (antisymmetric, symmetric_centered):
        require(sum(PI[i] * vector[i] for i in range(3)) == 0, "not centered")
        require(
            matvec(correct, vector) == tuple(rho * x for x in vector),
            "P centered eigenlaw failed",
        )
    require(
        matvec(twisted, antisymmetric)
        == tuple(-rho * x for x in antisymmetric),
        "twisted antisymmetric sign failed",
    )
    require(
        matvec(twisted, symmetric_centered)
        == tuple(rho * x for x in symmetric_centered),
        "twisted symmetric-centered sign failed",
    )

    bridge_paths = tuple(
        tuple(
            tuple(correct[i][k] * correct[k][j] * 169 for k in range(3))
            for j in range(3)
        )
        for i in range(3)
    )
    require(
        bridge_paths
        == (
            ((Q(1), Q(4), Q(20)), (Q(2), Q(2), Q(20)), (Q(10), Q(20), Q(90))),
            ((Q(2), Q(2), Q(20)), (Q(4), Q(1), Q(20)), (Q(20), Q(10), Q(90))),
            ((Q(2), Q(4), Q(18)), (Q(4), Q(2), Q(18)), (Q(20), Q(20), Q(81))),
        ),
        "wrong exact midpoint bridge path weights",
    )
    bridge_phi = tuple(
        tuple(
            (bridge_paths[i][j][0] - bridge_paths[i][j][1])
            / sum(bridge_paths[i][j])
            for j in range(3)
        )
        for i in range(3)
    )
    bridge_danger = tuple(
        tuple(
            bridge_paths[i][j][0] / sum(bridge_paths[i][j])
            for j in range(3)
        )
        for i in range(3)
    )
    require(
        bridge_phi
        == (
            (Q(-3, 25), Q(0), Q(-1, 12)),
            (Q(0), Q(3, 25), Q(1, 12)),
            (Q(-1, 12), Q(1, 12), Q(0)),
        ),
        "wrong signed midpoint bridge",
    )
    require(
        bridge_danger
        == (
            (Q(1, 25), Q(1, 12), Q(1, 12)),
            (Q(1, 12), Q(4, 25), Q(1, 6)),
            (Q(1, 12), Q(1, 6), Q(20, 121)),
        ),
        "wrong midpoint-danger bridge",
    )
    twisted_bridge_paths = tuple(
        tuple(
            tuple(twisted[i][k] * twisted[k][j] * 169 for k in range(3))
            for j in range(3)
        )
        for i in range(3)
    )
    twisted_bridge_phi = tuple(
        tuple(
            (
                twisted_bridge_paths[i][j][0]
                - twisted_bridge_paths[i][j][1]
            )
            / sum(twisted_bridge_paths[i][j])
            for j in range(3)
        )
        for i in range(3)
    )
    require(
        twisted_bridge_phi
        == tuple(tuple(-value for value in row) for row in bridge_phi),
        "twisted midpoint bridge does not negate the Z/2 coordinate",
    )

    # Lump to danger A versus safe B union C, using stationary conditional
    # weights 1/6 and 5/6 inside the safe state.
    correct_safe_to_danger = (
        PI[1] * correct[1][0] + PI[2] * correct[2][0]
    ) / (PI[1] + PI[2])
    twisted_safe_to_danger = (
        PI[1] * twisted[1][0] + PI[2] * twisted[2][0]
    ) / (PI[1] + PI[2])
    require(correct[0][0] == Q(1, 13), "correct A-to-A lump failed")
    require(correct_safe_to_danger == Q(2, 13), "correct safe lump failed")
    require(twisted[0][0] == Q(2, 13), "twisted A-to-A control failed")
    require(
        twisted_safe_to_danger == Q(11, 78),
        "twisted safe-to-A control failed",
    )

    weights_h = derive_guard_weights(1)
    weights_6h = derive_guard_weights(6)
    require(weights_h == (Q(0), Q(0), Q(1)), "u=H guard weights failed")
    require(
        weights_6h == (Q(5, 6), Q(5, 6), Q(4, 6)),
        "u=6H guard weights failed",
    )

    controls = (("u=H", weights_h), ("u=6H", weights_6h))
    expected = {
        ("correct", "u=H"): Q(28460, 199927),
        ("correct", "u=6H"): Q(734515, 5198102),
        ("twisted", "u=H"): Q(354710, 2599051),
        ("twisted", "u=6H"): Q(1057220, 7797153),
    }
    values = {}
    print("counts", counts)
    print("reverse-root-counts", reverse_counts)
    print("P", format_matrix(correct))
    print("JP", format_matrix(twisted))
    print("P2=(JP)2", format_matrix(expected_square))
    print("stationary", PI, "centered-eigenvalue-P", rho)
    print("midpoint-bridge-phi", format_matrix(bridge_phi))
    print("midpoint-bridge-danger", format_matrix(bridge_danger))
    print(
        "danger-lump P",
        correct[0][0],
        correct_safe_to_danger,
        "danger-lump JP",
        twisted[0][0],
        twisted_safe_to_danger,
    )
    print("guard-weights u=H", weights_h, "u=6H", weights_6h)
    for matrix_name, matrix in (("correct", correct), ("twisted", twisted)):
        for control_name, weights in controls:
            value, leading, repair = cylinder_value(matrix, weights)
            require(
                value == expected[(matrix_name, control_name)],
                f"wrong cylinder for {matrix_name}, {control_name}",
            )
            values[(matrix_name, control_name)] = value
            print(
                matrix_name,
                control_name,
                "value",
                value,
                "value-target",
                value - TARGET,
                "leading",
                leading,
                "repair-raw",
                repair,
            )

    require(values[("correct", "u=H")] > TARGET, "correct u=H should fail")
    require(values[("correct", "u=6H")] > TARGET, "correct u=6H should fail")
    require(values[("twisted", "u=H")] < TARGET, "twisted u=H should pass")
    require(values[("twisted", "u=6H")] < TARGET, "twisted u=6H should pass")
    correction_h = (
        values[("correct", "u=H")] - values[("twisted", "u=H")]
    )
    correction_6h = (
        values[("correct", "u=6H")] - values[("twisted", "u=6H")]
    )
    require(correction_h == Q(15270, 2599051), "u=H correction failed")
    require(correction_6h == Q(89105, 15594306), "u=6H correction failed")
    print("corrections", correction_h, correction_6h)

    unguarded, _, _ = cylinder_value(correct, (Q(1), Q(1), Q(1)))
    require(unguarded == Q(514705, 2599051), "unguarded value failed")
    decorrelated = Q(5, 7) * unguarded
    require(
        decorrelated == Q(2573525, 18193357),
        "decorrelated value failed",
    )
    require(decorrelated > TARGET, "decorrelated hostile control should fail")
    print(
        "unguarded",
        unguarded,
        "decorrelated",
        decorrelated,
        "decorrelated-target",
        decorrelated - TARGET,
    )
    direct_cells, direct_counts, direct_values = direct_grid_cylinders(
        (
            ((0, 0, 1), 1),
            ((5, 5, 4), 6),
            ((1, 1, 1), 1),
        )
    )
    require(
        direct_values
        == (
            values[("correct", "u=H")],
            values[("correct", "u=6H")],
            unguarded,
        ),
        "independent exact-grid census disagrees with Markov evaluation",
    )
    print(
        "direct-grid-cells",
        direct_cells,
        "integer-leading-repair",
        direct_counts,
        "values",
        direct_values,
    )
    print("AUDIT PASS")


if __name__ == "__main__":
    main()
