#!/usr/bin/env python3
"""Dependency-free exact referee for the product-torus part of THM-2534.

This script checks the finite algebra behind the transverse target--root
boundary construction.  It deliberately does *not* reprove the live LRC
input imported from proved THM-2533.  Its finite scope is conditional on
those distinct mixed gains and their nonzero linear Fourier witnesses: one
transverse product-torus translation produces a signed atomwise exclusive
Boolean boundary, and one nonnegative integer sum retains all the witnesses.

All checks use integers, fractions, or exact cyclotomic normal forms.  No
``assert`` statement is used, so normal and ``python -O`` executions have
identical force and output.
"""

from fractions import Fraction
from itertools import combinations, product


ROOT_P = 13
SOURCE_P = 7
MODULUS = ROOT_P * SOURCE_P
checks = 0


def require(condition: bool, message: str) -> None:
    """Raise on a failed audit in both normal and optimized Python."""

    global checks
    checks += 1
    if not condition:
        raise RuntimeError("FAILED: " + message)


def inverse(value: int, modulus: int) -> int:
    return pow(value % modulus, -1, modulus)


# ---------------------------------------------------------------------------
# Exact Q(zeta_13) normal forms.
#
# A group-ring coefficient vector c_0,...,c_12 evaluates at zeta_13.  The
# single relation 1+z+...+z^12=0 gives the unique degree-at-most-11 form
# (c_0-c_12,...,c_11-c_12).  Integer coefficients suffice for this referee.


def cyclo13_reduce(coefficients: tuple[int, ...]) -> tuple[int, ...]:
    require(len(coefficients) == ROOT_P, "cyclotomic input length")
    tail = coefficients[-1]
    return tuple(value - tail for value in coefficients[:-1])


def cyclo13_root(exponent: int) -> tuple[int, ...]:
    coefficients = [0] * ROOT_P
    coefficients[exponent % ROOT_P] = 1
    return cyclo13_reduce(tuple(coefficients))


def cyclo13_add(
    left: tuple[int, ...], right: tuple[int, ...]
) -> tuple[int, ...]:
    return tuple(a + b for a, b in zip(left, right, strict=True))


def cyclo13_sub(
    left: tuple[int, ...], right: tuple[int, ...]
) -> tuple[int, ...]:
    return tuple(a - b for a, b in zip(left, right, strict=True))


def cyclo13_scale(value: tuple[int, ...], scalar: int) -> tuple[int, ...]:
    return tuple(scalar * coefficient for coefficient in value)


def cyclo13_mul(
    left: tuple[int, ...], right: tuple[int, ...]
) -> tuple[int, ...]:
    coefficients = [0] * ROOT_P
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            coefficients[(i + j) % ROOT_P] += a * b
    return cyclo13_reduce(tuple(coefficients))


CYCLO13_ZERO = (0,) * (ROOT_P - 1)
CYCLO13_ONE = cyclo13_root(0)

cyclotomic_product_checks = 0
for first_exponent, second_exponent in product(range(ROOT_P), repeat=2):
    require(
        cyclo13_mul(
            cyclo13_root(first_exponent), cyclo13_root(second_exponent)
        )
        == cyclo13_root(first_exponent + second_exponent),
        "exact cyclotomic monomial product",
    )
    cyclotomic_product_checks += 1


# The characteristic multiplier of
#
#   (P_(tau,h) f)(s,r) = f(s+h,r+tau)
#
# on the target--root character zeta^(b s-a r) is
# zeta^(a tau-b h).  At prime order it vanishes after subtracting one
# exactly on the characteristic/tangent line a tau=b h.
characteristic_checks = 0
mixed_characteristic_checks = 0
tangent_count = 0
transverse_count = 0
for a, b, tau, h in product(range(ROOT_P), repeat=4):
    exponent = (a * tau - b * h) % ROOT_P
    multiplier = cyclo13_sub(cyclo13_root(exponent), CYCLO13_ONE)
    require(
        (multiplier == CYCLO13_ZERO) == (exponent == 0),
        "prime cyclotomic characteristic kernel",
    )
    characteristic_checks += 1

for a, b, tau, h in product(range(1, ROOT_P), repeat=4):
    exponent = (a * tau - b * h) % ROOT_P
    tangent = h == (a * inverse(b, ROOT_P) * tau) % ROOT_P
    require((exponent == 0) == tangent, "mixed tangent iff characteristic zero")
    if tangent:
        tangent_count += 1
    else:
        transverse_count += 1
    mixed_characteristic_checks += 1

require(tangent_count == 12**3, "mixed tangent census")
require(transverse_count == 12**4 - 12**3, "mixed transverse census")


# Audit the translation law on every delta basis vector of the 13 x 13
# product torus, every target--root character, and every nonzero direction.
# Linearity then proves the multiplier identity for every table.
basis_translation_checks = 0
for source_target, root in product(range(ROOT_P), repeat=2):
    for a, b in product(range(ROOT_P), repeat=2):
        old_exponent = (b * source_target - a * root) % ROOT_P
        for tau, h in product(range(1, ROOT_P), repeat=2):
            # P sends a delta at (s,r) to the delta at (s-h,r-tau).
            new_exponent = (
                b * (source_target - h) - a * (root - tau)
            ) % ROOT_P
            expected = (old_exponent + a * tau - b * h) % ROOT_P
            require(new_exponent == expected, "delta-basis product-torus multiplier")
            basis_translation_checks += 1


# Every triple of distinct nonzero gains leaves exactly nine nonzero target
# shifts transverse to all three gains, for each fixed nonzero root step.
gain_triples = tuple(combinations(range(1, ROOT_P), 3))
transverse_triple_checks = 0
transverse_gain_multiplier_checks = 0
for gains in gain_triples:
    for tau in range(1, ROOT_P):
        forbidden = {(gain * tau) % ROOT_P for gain in gains}
        allowed = tuple(h for h in range(1, ROOT_P) if h not in forbidden)
        require(len(forbidden) == 3, "distinct gain tangent directions")
        require(len(allowed) == 9, "three gains leave nine transverse shifts")
        for b in range(1, ROOT_P):
            for h in allowed:
                for gain in gains:
                    a = (gain * b) % ROOT_P
                    exponent = (a * tau - b * h) % ROOT_P
                    require(exponent != 0, "common transverse shift detects every gain")
                    transverse_gain_multiplier_checks += 1
        transverse_triple_checks += 1


# ---------------------------------------------------------------------------
# Pointwise Boolean boundary split and invariant-factor retention.

exclusive_boolean_checks = 0
for old, new in product((0, 1), repeat=2):
    forward = new * (1 - old)
    backward = old * (1 - new)
    require(new - old == forward - backward, "exclusive Boolean difference")
    require(forward in (0, 1) and backward in (0, 1), "Boolean boundary values")
    require(forward * backward == 0, "forward/backward boundaries are disjoint")
    exclusive_boolean_checks += 1

invariant_factor_checks = 0
for common, old_tail, new_tail in product((0, 1), repeat=3):
    old = common * old_tail
    new = common * new_tail
    forward = new * (1 - old)
    backward = old * (1 - new)
    require(
        forward == common * new_tail * (1 - old_tail),
        "forward boundary retains invariant Boolean factor",
    )
    require(
        backward == common * old_tail * (1 - new_tail),
        "backward boundary retains invariant Boolean factor",
    )
    invariant_factor_checks += 1


# ---------------------------------------------------------------------------
# Common nonnegative packing scalar.
#
# If u_i-v_i is nonzero, the equation u_i+t v_i=0 excludes at most one t.
# Hence r rows cannot cover r+1 prescribed scalars.  We exhaust every exact
# rational pair (u,v) in [-3,3]^2 with u-v != 0 for r=1,2,3.  The tight
# controls u_i=-i, v_i=1 force the first r candidates to fail.


coefficient_pairs = tuple(
    (Fraction(u), Fraction(v))
    for u, v in product(range(-3, 4), repeat=2)
    if u != v
)
require(len(coefficient_pairs) == 42, "packing coefficient-pair census")

packing_case_counts: dict[int, int] = {}
packing_choice_histograms: dict[int, dict[int, int]] = {}
for row_count in range(1, 4):
    candidates = tuple(Fraction(t) for t in range(1, row_count + 2))
    case_count = 0
    histogram = {t: 0 for t in range(1, row_count + 2)}
    for rows in product(coefficient_pairs, repeat=row_count):
        good = tuple(
            t
            for t in candidates
            if all(u + t * v != 0 for u, v in rows)
        )
        require(bool(good), "r+1 scalar packing leaves a common survivor")
        choice = int(good[0])
        histogram[choice] += 1
        case_count += 1
    packing_case_counts[row_count] = case_count
    packing_choice_histograms[row_count] = histogram

    tight_rows = tuple((Fraction(-i), Fraction(1)) for i in range(1, row_count + 1))
    tight_good = tuple(
        t
        for t in candidates
        if all(u + t * v != 0 for u, v in tight_rows)
    )
    require(
        tight_good == (Fraction(row_count + 1),),
        "r+1 candidate bound is sharp on rational controls",
    )

require(
    packing_case_counts == {1: 42, 2: 42**2, 3: 42**3},
    "packing exhaustive case counts",
)

# Repeat the one-row exclusion audit on genuinely non-real cyclotomic
# controls.  The subsequent forbidden-value census isolates the universal
# pigeonhole step from the chosen coefficient sample.
cyclotomic_values = {
    CYCLO13_ZERO,
    CYCLO13_ONE,
    cyclo13_scale(CYCLO13_ONE, -1),
}
for exponent in range(ROOT_P):
    root = cyclo13_root(exponent)
    cyclotomic_values.add(root)
    cyclotomic_values.add(cyclo13_scale(root, -1))
    cyclotomic_values.add(cyclo13_add(CYCLO13_ONE, root))

cyclotomic_pair_packing_checks = 0
for u, v in product(tuple(cyclotomic_values), repeat=2):
    if cyclo13_sub(u, v) == CYCLO13_ZERO:
        continue
    killed = tuple(
        t
        for t in range(1, 5)
        if cyclo13_add(u, cyclo13_scale(v, t)) == CYCLO13_ZERO
    )
    require(len(killed) <= 1, "one cyclotomic row excludes at most one scalar")
    cyclotomic_pair_packing_checks += 1

for row_count in range(1, 4):
    candidates = tuple(range(1, row_count + 2))
    possible_bad_values = (None,) + candidates
    for bad_rows in product(possible_bad_values, repeat=row_count):
        survivors = tuple(t for t in candidates if t not in bad_rows)
        require(bool(survivors), "r abstract bad values cannot cover r+1 candidates")


# ---------------------------------------------------------------------------
# Exact Galois index covariance and the 216-incidence count.


def crt_unit(unit7: int, unit13: int) -> int:
    candidates = tuple(
        value
        for value in range(1, MODULUS)
        if value % SOURCE_P == unit7 % SOURCE_P
        and value % ROOT_P == unit13 % ROOT_P
    )
    require(len(candidates) == 1, "CRT unit uniqueness")
    return candidates[0]


galois_units = tuple(
    (unit7, unit13, crt_unit(unit7, unit13))
    for unit7 in range(1, SOURCE_P)
    for unit13 in range(1, ROOT_P)
)
require(len(galois_units) == 72, "Gal(Q(zeta_91)/Q) size")
require(len({unit for _, _, unit in galois_units}) == 72, "CRT Galois units distinct")
require(
    all(inverse(unit, MODULUS) for _, _, unit in galois_units),
    "every CRT representative is a unit",
)


def phase91(
    kappa: int,
    target_colour: int,
    root_colour: int,
    source_difference: int,
    target_difference: int,
    root_difference: int,
) -> int:
    return (
        13 * kappa * source_difference
        + 7 * (target_colour * target_difference + root_colour * root_difference)
    ) % MODULUS


# A rational overlap coefficient is fixed by Galois.  On every phase basis
# monomial relevant to the three control gains, the CRT automorphism scales
# (kappa,b,a) exactly as claimed.
control_gains = (1, 2, 4)
galois_basis_checks = 0
for source_difference in range(SOURCE_P):
    for target_difference, root_difference in product(range(ROOT_P), repeat=2):
        for gain in control_gains:
            base_phase = phase91(
                1,
                1,
                gain,
                source_difference,
                target_difference,
                root_difference,
            )
            for unit7, unit13, unit91 in galois_units:
                transformed_phase = phase91(
                    unit7,
                    unit13,
                    (unit13 * gain) % ROOT_P,
                    source_difference,
                    target_difference,
                    root_difference,
                )
                require(
                    transformed_phase == (unit91 * base_phase) % MODULUS,
                    "Galois covariance on a rational phase basis",
                )
                galois_basis_checks += 1


galois_orbit_checks = 0
for gains in gain_triples:
    total_orbit: set[tuple[int, int, int]] = set()
    for gain in gains:
        gain_orbit = {
            (
                unit7,
                unit13,
                (unit13 * gain) % ROOT_P,
            )
            for unit7, unit13, _ in galois_units
        }
        require(len(gain_orbit) == 72, "one gain has 72 Galois incidences")
        require(
            all(
                root_colour * inverse(target_colour, ROOT_P) % ROOT_P == gain
                for _, target_colour, root_colour in gain_orbit
            ),
            "Galois action preserves the gain ratio",
        )
        require(total_orbit.isdisjoint(gain_orbit), "distinct gains have disjoint orbits")
        total_orbit.update(gain_orbit)
    require(len(total_orbit) == 216, "three gains give 216 incidences")
    galois_orbit_checks += 1


# ---------------------------------------------------------------------------
# Sharp Boolean simultaneous-selection hostile.
#
# On the two-bit cube use L_1=x_1, L_2=x_2, L_3=x_1-x_2.  Each functional
# is live, but their zero hyperplanes cover every Boolean point.  Thus the
# abstract hypotheses cannot select one Boolean atom retaining all three.
# A nonnegative integer packing can: (1,0)+(1,1)=(2,1).


boolean_cube = tuple(product((0, 1), repeat=2))


def hostile_functionals(point: tuple[int, int]) -> tuple[int, int, int]:
    first, second = point
    return first, second, first - second


hostile_values = {point: hostile_functionals(point) for point in boolean_cube}
require(
    all(any(value == 0 for value in values) for values in hostile_values.values()),
    "three zero hyperplanes cover the Boolean square",
)
require(
    all(any(values[index] != 0 for values in hostile_values.values()) for index in range(3)),
    "each hostile functional is individually live",
)
require(
    not any(all(value != 0 for value in values) for values in hostile_values.values()),
    "no one Boolean selector retains all three incidences",
)
packed_point = (2, 1)
require(
    hostile_functionals(packed_point) == (2, 1, 1),
    "nonnegative multiplicity packing retains all three incidences",
)
require(
    tuple(a + b for a, b in zip((1, 0), (1, 1), strict=True)) == packed_point,
    "hostile packing is a sum of Boolean selectors",
)


# The familiar matched positive-flux equality needs a weight invariant under
# the same shift.  This finite orbit audit records that scope boundary rather
# than silently treating a target-dependent sidecar as common.
flux_invariant_checks = 0
flux_mismatch_witnesses = 0
for mask in range(1 << ROOT_P):
    values = tuple((mask >> index) & 1 for index in range(ROOT_P))
    shifted = tuple(values[(index + 1) % ROOT_P] for index in range(ROOT_P))
    for constant_weight in (Fraction(0), Fraction(1), Fraction(5, 3)):
        boundary = sum(
            constant_weight * values[index] * (1 - shifted[index])
            for index in range(ROOT_P)
        )
        square = Fraction(1, 2) * sum(
            constant_weight * (shifted[index] - values[index]) ** 2
            for index in range(ROOT_P)
        )
        require(boundary == square, "matched flux for invariant weight")
        flux_invariant_checks += 1

    noninvariant_weight = tuple(Fraction(int(index == 0)) for index in range(ROOT_P))
    boundary = sum(
        noninvariant_weight[index] * values[index] * (1 - shifted[index])
        for index in range(ROOT_P)
    )
    square = Fraction(1, 2) * sum(
        noninvariant_weight[index] * (shifted[index] - values[index]) ** 2
        for index in range(ROOT_P)
    )
    if boundary != square:
        flux_mismatch_witnesses += 1

require(flux_mismatch_witnesses > 0, "noninvariant weight produces a flux mismatch")


packing_summary = ",".join(
    f"r{row_count}:{packing_case_counts[row_count]}"
    for row_count in sorted(packing_case_counts)
)
packing_histogram = ";".join(
    "r"
    + str(row_count)
    + "="
    + ",".join(
        f"{choice}:{packing_choice_histograms[row_count][choice]}"
        for choice in sorted(packing_choice_histograms[row_count])
    )
    for row_count in sorted(packing_choice_histograms)
)

print("THM-2534 transverse product-torus exact referee")
print("status=CONDITIONAL-ALGEBRA-VERIFIED")
print(f"cyclotomic_product_checks={cyclotomic_product_checks}")
print(f"characteristic_checks={characteristic_checks}")
print(f"mixed_characteristic_checks={mixed_characteristic_checks}")
print(f"tangent_count={tangent_count} transverse_count={transverse_count}")
print(f"basis_translation_checks={basis_translation_checks}")
print(f"gain_triples={len(gain_triples)} transverse_triple_checks={transverse_triple_checks}")
print(f"transverse_gain_multiplier_checks={transverse_gain_multiplier_checks}")
print(f"exclusive_boolean_checks={exclusive_boolean_checks}")
print(f"invariant_factor_checks={invariant_factor_checks}")
print(f"packing_cases={packing_summary}")
print(f"packing_first_choice_histogram={packing_histogram}")
print(f"cyclotomic_pair_packing_checks={cyclotomic_pair_packing_checks}")
print("packing_candidates=t_in_{1,...,r+1} sharp_for_r=1,2,3")
print(f"galois_group_size={len(galois_units)}")
print(f"galois_basis_checks={galois_basis_checks}")
print(f"galois_gain_triples={galois_orbit_checks} incidences_per_triple=216")
print("boolean_simultaneous_selection_hostile=L1=x1,L2=x2,L3=x1-x2")
print("boolean_common_atom=IMPOSSIBLE nonnegative_pack=(2,1) values=(2,1,1)")
print(f"flux_invariant_checks={flux_invariant_checks}")
print(f"flux_noninvariant_mismatch_witnesses={flux_mismatch_witnesses}")
print("live_three_gain_input=IMPORTED_FROM_THM2533_PROVED_DEPENDENCY")
print("semantic_arrival_or_owner_loop=NOT_VERIFIED")
print("anchor_and_zero_row_preservation=NOT_VERIFIED")
print(f"checks={checks}")
print("all_checks=PASS")
