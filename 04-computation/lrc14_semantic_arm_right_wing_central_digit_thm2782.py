#!/usr/bin/env python3
"""Finite-exact semantic-arm attachment of the THM-2771 right wing.

Candidate theorem companion.  It rebuilds the physical interval objects from
the canonical carrier constructors and proves three explicit target-root-one
semantic witnesses inside the inherited root-zero gauge, their raw contents
and positive weighted masses, the clock-one Q3 arm contrast, its intrinsic
mod-13 Bockstein decoding, and the exact augmentation-quotient inverse.

It also keeps three coordinates separate:

* n in Z/13^6: the physical semantic address;
* (sigma,tau) in F_13^2: lawful present-packet target-action labels;
* R_0 in V=F_13^2: an allocated endpoint origin.

After the standard mod-169 digit section, the physical +13 address step has
the same coordinate formula as the Heisenberg central direction.  This is
not an allocated-endpoint identification.  Exact descent, cospan, and target
label controls record the three missing naturality statements.

No floating point, randomized calculation, Python assertion, or canon edit is
used.
"""

from fractions import Fraction
from math import gcd
from pathlib import Path
import hashlib
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))

import lrc14_fully_marked_root_zero_target_profile_thm2749 as marked
import lrc14_root_zero_clutch_mayer_vietoris_wing_shear_thm2751 as wing
import lrc14_relative_present_semantic_lift_probe_20260728 as relative


P = 13
Q = 7
ADDRESS_MODULUS = P**6
T = marked.T
TAU = Fraction(7, ADDRESS_MODULUS)
ARM_STEP = 13
ARM_OFFSETS = (0, ARM_STEP, 2 * ARM_STEP)
EXPECTED_T = 297836897838480
EXPECTED_SHIFT = 431933040
EXPECTED_Q_RADIUS = Fraction(1, 100360982066072)
EXPECTED_BASES = {1: 3454614, 2: 1393828, 3: 708364}
EXPECTED_CONTENTS = {
    1: 790161473087466480,
    2: 790135314376327920,
    3: 790135314376327920,
}
EXPECTED_PRIMITIVE = {
    1: (
        (0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
        (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1),
        (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1),
    ),
    2: (
        (0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1),
        (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1),
        (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1),
    ),
    3: (
        (0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1),
        (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1),
        (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1),
    ),
}
K_BETA = (3, 5, 5, 5, 7, 12, 2, 8, 2, 9, 8, 11, 0)
EXPECTED_E1_Q3_BETA_DECODED = (
    11, 9, 2, 12, 10, 9, 10, 4, 3, 5, 6, 1, 12
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def frac(value):
    return value - value.numerator // value.denominator


def valuation(value, prime):
    require(value != 0, "valuation of zero requested")
    result = 0
    value = abs(value)
    while value % prime == 0:
        result += 1
        value //= prime
    return result


def cyclic_convolution(left, right, modulus=P):
    return tuple(
        sum(
            left[index] * right[(degree - index) % modulus]
            for index in range(modulus)
        )
        for degree in range(modulus)
    )


def cyclic_convolution_mod(left, right, modulus=P):
    return tuple(
        value % modulus for value in cyclic_convolution(left, right, modulus)
    )


def q3_projection(rows):
    """Integral arm augmentation quotient: 3*v_i-sum_j(v_j)."""
    return tuple(
        tuple(
            3 * rows[arm][target]
            - sum(rows[other][target] for other in range(3))
            for target in range(P)
        )
        for arm in range(3)
    )


def weighted_mass(pieces):
    return sum(
        (right - left) * weight for left, right, weight in pieces
    )


def build_context():
    module, _prefixes, _a, _b, rails, present, _starts = (
        marked.lift.m.core.build_carrier_data()
    )
    delayed = marked.marked_prefixes(
        module,
        marked.private.build_pair_prefixes(module),
        marked.two.deepest_fork(module),
    )
    source = marked.two.exclusive_source(module, 3)
    clock_comb = tuple(
        module.make_comb(module.C1, 182, 26 * e - 13, 26 * e + 13)
        for e in range(Q)
    )
    source_weight, _target_weight, rail_common = marked.rail_data(
        rails, marked.RAIL
    )
    return (
        module, delayed, present, source, clock_comb,
        source_weight, rail_common,
    )


def carrier_objects(module, source, clock_comb, rail_common, e, sigma, tau):
    """Full lawful (sigma,tau) present section and corrected cofiber objects."""
    raw = tuple(marked.two.intersect_sorted(source, clock_comb[e]))
    raw = tuple(module.subtract_comb(
        raw, module.W[1], 182,
        -14 * sigma - 13, -14 * sigma + 13,
    ))
    raw = tuple(module.subtract_comb(
        raw, module.C2, 182,
        14 * sigma - 13, 14 * sigma + 13,
    ))
    raw = tuple(module.subtract_comb(
        raw, module.W[2], 182,
        -14 * tau - 13, -14 * tau + 13,
    ))
    raw = tuple(module.subtract_comb(
        raw, module.C3, 182,
        14 * tau - 13, 14 * tau + 13,
    ))
    a = marked.intersect(rail_common, raw)
    b = marked.intersect(
        rail_common, marked.shift_union(raw, marked.SHIFT)
    )
    common = marked.intersect(a, b)
    left = wing.difference(a, b)
    right = wing.difference(b, a)
    require(
        marked.merge(common + left) == marked.merge(a)
        and marked.merge(common + right) == marked.merge(b),
        f"cofiber decomposition failed at {(e, sigma, tau)}",
    )
    return {
        "A": a,
        "B": b,
        "M": common,
        "L": left,
        "R": right,
        # Push pulled target objects into the target coordinate.
        "R+": marked.shift_union(right, -marked.SHIFT),
    }


def arm_source_center(base, arm):
    return frac(relative.Z + (base + ARM_STEP * arm) * TAU)


def arm_target_center(base, arm):
    return frac(arm_source_center(base, arm) + TAU)


def target_cylinder(base, arm):
    center = arm_target_center(base, arm)
    half_width = relative.Q_RADIUS * T
    coordinate = center * T
    left = coordinate - half_width
    right = coordinate + half_width
    require(0 < left < right < T, "selected target cylinder crossed seam")
    return ((left, right),)


def predecessor_carry(center):
    coordinate = center * T
    root_coordinate = coordinate * marked.private.S % T
    return (P * root_coordinate) // T


def local_cell(module, delayed, present, source_weight, objects,
               base, arm, tau, terminal_clock=1, target_side=True):
    """Physical mass and actual delayed carry-six coefficient on one arm."""
    source_safe = marked.complement(present[terminal_clock, 7])
    target_safe = marked.complement(marked.shift_union(
        present[terminal_clock, 7], marked.SHIFT
    ))
    carrier = marked.intersect(objects["R"], source_safe)
    carrier = marked.intersect(carrier, target_safe)
    if target_side:
        carrier = marked.shift_union(carrier, -marked.SHIFT)
        root_bounds = (1, 13)
        carry = 6
    else:
        root_bounds = (169, 181)
        carry = 12
    weighted = marked.restrict_weighted(source_weight, carrier)
    weighted = marked.private.old.intersect_weighted_comb(
        weighted, module.C3, 182, *root_bounds
    )
    weighted = marked.restrict_weighted(
        weighted, target_cylinder(base, arm)
    )
    mass = weighted_mass(weighted)
    coefficient = marked.private.delayed_carry_pair(
        weighted, delayed[terminal_clock], {}
    )[carry][1]
    if coefficient:
        center = arm_target_center(base, arm)
        sigma_labels, tau_labels = relative.full_target_labels(center)
        require(
            relative.semantic.semantic_record(center) == (3, (1, 2)),
            "positive arm left the E3/Q_(3,{1,2}) semantic word",
        )
        require(
            0 in sigma_labels and tau in tau_labels,
            "positive arm lost its lawful (sigma=0,tau) target label",
        )
        require(
            predecessor_carry(center) == 6,
            "positive arm lost actual delayed carry six",
        )
    return mass, coefficient


def witness_table(context, e, base):
    (
        module, delayed, present, source, clocks,
        source_weight, rail_common,
    ) = context
    masses = [[Fraction(0) for _tau in range(P)] for _arm in range(3)]
    coefficients = [[0 for _tau in range(P)] for _arm in range(3)]
    for tau in range(P):
        objects = carrier_objects(
            module, source, clocks, rail_common, e, 0, tau
        )
        for arm in range(3):
            mass, coefficient = local_cell(
                module, delayed, present, source_weight,
                objects, base, arm, tau
            )
            masses[arm][tau] = mass
            coefficients[arm][tau] = coefficient
    return tuple(map(tuple, masses)), tuple(map(tuple, coefficients))


def cospan_side_cell(context, channel, base, arm, tau, target_side):
    """Evaluate one canonical source- or target-endpoint coefficient."""
    (
        module, delayed, present, source, clocks,
        source_weight, rail_common,
    ) = context
    objects = carrier_objects(
        module, source, clocks, rail_common, 1, 0, tau
    )
    source_safe = marked.complement(present[1, 7])
    target_safe = marked.complement(marked.shift_union(
        present[1, 7], marked.SHIFT
    ))
    carrier = marked.intersect(objects[channel], source_safe)
    carrier = marked.intersect(carrier, target_safe)
    if target_side:
        carrier = marked.shift_union(carrier, -marked.SHIFT)
        root_bounds = (1, 13)
        carry = 6
        center = arm_target_center(base, arm)
    else:
        root_bounds = (169, 181)
        carry = 12
        center = arm_source_center(base, arm)
    half_width = relative.Q_RADIUS * T
    cylinder = ((
        center * T - half_width,
        center * T + half_width,
    ),)
    weighted = marked.restrict_weighted(source_weight, carrier)
    weighted = marked.private.old.intersect_weighted_comb(
        weighted, module.C3, 182, *root_bounds
    )
    weighted = marked.restrict_weighted(weighted, cylinder)
    return (
        weighted_mass(weighted),
        marked.private.delayed_carry_pair(
            weighted, delayed[1], {}
        )[carry][1],
    )


def omega_digits(address):
    """n=v+13w on Omega=Z/169, with v,w in F_13."""
    residue = address % (P**2)
    v = residue % P
    w = residue // P
    return v, w


def omega_from_digits(v, w):
    return (v % P) + P * (w % P)


def omega_x(address):
    """Endpoint transvection coordinate: (v,w)->(v,w+v)."""
    return 14 * address % (P**2)


def omega_y(address):
    """Carry-suppressed low-digit shift: (v,w)->(v+1,w)."""
    v, w = omega_digits(address)
    return omega_from_digits(v + 1, w)


def omega_z(address):
    """Central digit shift: (v,w)->(v,w+1)."""
    return (address + P) % (P**2)


def arm_metadata(base):
    return tuple(
        (
            base + ARM_OFFSETS[arm],
            (base + ARM_OFFSETS[arm]) % P,
            arm_source_center(base, arm),
            arm_target_center(base, arm),
            relative.semantic.semantic_record(arm_source_center(base, arm)),
            relative.semantic.semantic_record(arm_target_center(base, arm)),
            predecessor_carry(arm_target_center(base, arm)),
        )
        for arm in range(3)
    )


def main():
    require(
        (T, marked.SHIFT, relative.Q_RADIUS)
        == (EXPECTED_T, EXPECTED_SHIFT, EXPECTED_Q_RADIUS),
        "canonical scale, target shift, or semantic radius changed",
    )
    require(
        ADDRESS_MODULUS == 4826809
        and ADDRESS_MODULUS // gcd(ADDRESS_MODULUS, ARM_STEP) == P**5,
        "semantic address-step order changed",
    )
    require(
        ARM_STEP * TAU == Fraction(91, ADDRESS_MODULUS),
        "semantic arm translation changed",
    )

    context = build_context()
    witness_rows = []
    witness_data = {}
    for e in (1, 2, 3):
        base = EXPECTED_BASES[e]
        metadata = arm_metadata(base)
        require(
            all(
                target_record == (3, (1, 2))
                and carry == 6
                and residue == 7
                for (
                    _address, residue, _source, _target,
                    source_record, target_record, carry,
                ) in metadata
            ),
            f"clock-{e} target witness lost semantic/carry/root-gauge typing",
        )
        masses, raw = witness_table(context, e, base)
        nonzero = tuple(value for row in raw for value in row if value)
        content = 0
        for value in nonzero:
            content = gcd(content, abs(value))
        require(
            content == EXPECTED_CONTENTS[e]
            and valuation(content, 13) == 1,
            f"clock-{e} raw content or v13 changed",
        )
        primitive = tuple(
            tuple(value // content for value in row) for row in raw
        )
        require(
            primitive == EXPECTED_PRIMITIVE[e],
            f"clock-{e} primitive arm table changed",
        )
        require(
            all(
                mass == Fraction(content, ADDRESS_MODULUS)
                for arm in range(3)
                for tau, mass in enumerate(masses[arm])
                if primitive[arm][tau]
            )
            and all(
                mass == 0
                for arm in range(3)
                for tau, mass in enumerate(masses[arm])
                if not primitive[arm][tau]
            ),
            f"clock-{e} positive mass/content law changed",
        )

        bockstein_unit = (content // 13) % 13
        beta = tuple(
            tuple(bockstein_unit * value % 13 for value in row)
            for row in primitive
        )
        decoded = tuple(
            cyclic_convolution_mod(row, K_BETA) for row in beta
        )
        q3_decoded = tuple(
            tuple(value % 13 for value in row)
            for row in q3_projection(decoded)
        )
        if e == 1:
            require(
                q3_decoded[0] == EXPECTED_E1_Q3_BETA_DECODED
                and all(q3_decoded[0])
                and q3_decoded[1] == q3_decoded[2]
                and all(
                    (q3_decoded[0][target]
                     + q3_decoded[1][target]
                     + q3_decoded[2][target]) % 13 == 0
                    for target in range(P)
                ),
                "clock-one decoded Q3 Bockstein charge changed",
            )
        witness_rows.append((
            e, base, metadata, content, bockstein_unit,
            Fraction(content, ADDRESS_MODULUS),
            primitive, q3_decoded,
        ))
        witness_data[e] = (masses, raw, primitive)

    # The clock-one exceptional tau=12 column is uniform across the arms.
    # Q3 kills it and leaves W=z^3+...+z^11 on arm zero.
    w = tuple(1 if 3 <= exponent <= 11 else 0 for exponent in range(P))
    inverse = tuple(1 if exponent in (2, 6, 10) else 0 for exponent in range(P))
    product = cyclic_convolution(w, inverse)
    require(
        product == (3,) + (2,) * 12,
        "W inverse identity W*(z^2+z^6+z^10)=delta_0+2N changed",
    )
    require(
        tuple(
            EXPECTED_PRIMITIVE[1][0][target]
            - EXPECTED_PRIMITIVE[1][1][target]
            for target in range(P)
        ) == w,
        "clock-one arm contrast stopped being W",
    )

    # The three printed arms are a genuine central-Z segment after the exact
    # endpoint-origin quotient Omega=Z/169.  All three clock witnesses reduce
    # to 85=7+13*6, and n->n+13 fixes v while incrementing w.
    require(
        all(base % (P**2) == 85 for base in EXPECTED_BASES.values())
        and omega_digits(85) == (7, 6),
        "three clock witnesses lost their common endpoint-origin residue",
    )
    for address in range(P**2):
        v, w_digit = omega_digits(address)
        require(
            omega_digits(omega_x(address)) == (v, (w_digit + v) % P)
            and omega_digits(omega_y(address))
            == ((v + 1) % P, w_digit)
            and omega_digits(omega_z(address))
            == (v, (w_digit + 1) % P),
            "Omega action formulas changed",
        )
        # X Y = Z Y X, equivalent to [X,Y]=Z in the theorem convention.
        require(
            omega_x(omega_y(address))
            == omega_z(omega_y(omega_x(address))),
            "Omega Heisenberg commutator changed",
        )

    central_metadata = tuple(
        (
            j,
            (EXPECTED_BASES[1] + ARM_STEP * j) % (P**2),
            omega_digits(EXPECTED_BASES[1] + ARM_STEP * j),
            relative.semantic.semantic_record(
                arm_target_center(EXPECTED_BASES[1], j)
            ),
            predecessor_carry(arm_target_center(EXPECTED_BASES[1], j)),
        )
        for j in range(P)
    )
    require(
        tuple(row[2] for row in central_metadata)
        == tuple((7, (6 + j) % P) for j in range(P))
        and all(
            semantic_record == (3, (1, 2)) and carry == 6
            for _j, _address, _digits, semantic_record, carry
            in central_metadata
        ),
        "clock-one central lift lost endpoint/semantic/carry typing",
    )

    # Extend the exact coefficient table from the first three arms to one
    # thirteen-point lift of the central cycle.  Reuse the first three rows.
    (
        module, delayed, present, source, clocks,
        source_weight, rail_common,
    ) = context
    central_raw = [
        list(witness_data[1][1][arm]) for arm in range(3)
    ] + [[0 for _tau in range(P)] for _j in range(3, P)]
    central_mass = [
        list(witness_data[1][0][arm]) for arm in range(3)
    ] + [[Fraction(0) for _tau in range(P)] for _j in range(3, P)]
    for tau in range(P):
        objects = carrier_objects(
            module, source, clocks, rail_common, 1, 0, tau
        )
        for j in range(3, P):
            mass, coefficient = local_cell(
                module, delayed, present, source_weight,
                objects, EXPECTED_BASES[1], j, tau
            )
            central_mass[j][tau] = mass
            central_raw[j][tau] = coefficient
    central_raw = tuple(tuple(row) for row in central_raw)
    central_mass = tuple(tuple(row) for row in central_mass)
    central_primitive = tuple(
        tuple(value // EXPECTED_CONTENTS[1] for value in row)
        for row in central_raw
    )
    expected_central = (
        EXPECTED_PRIMITIVE[1][0],
        *((0,) * 12 + (1,) for _j in range(1, 7)),
        *((0,) * 13 for _j in range(7, 13)),
    )
    require(
        central_primitive == expected_central
        and all(
            central_mass[j][tau]
            == (
                Fraction(EXPECTED_CONTENTS[1], ADDRESS_MODULUS)
                if central_primitive[j][tau] else 0
            )
            for j in range(P) for tau in range(P)
        ),
        "full central-lift coefficient table changed",
    )

    # The target decoder makes every one of the thirteen target columns
    # nonconstant along this central lift.  Therefore each column has all
    # twelve primitive complex central characters.  The stronger descent to
    # Omega nevertheless fails: the next physical lift point j=13 has the
    # same Omega address as j=0 but zero coefficient at tau=3 and tau=12.
    central_beta_integer = tuple(
        tuple((EXPECTED_CONTENTS[1] // 13) * value for value in row)
        for row in central_primitive
    )
    central_decoded_integer = tuple(
        cyclic_convolution(row, K_BETA) for row in central_beta_integer
    )
    require(
        all(
            len({
                central_decoded_integer[j][target] for j in range(P)
            }) > 1
            for target in range(P)
        ),
        "a chosen integer-lift decoded target column became central-flat",
    )
    next_lift = {}
    for tau in (3, 12):
        objects = carrier_objects(
            module, source, clocks, rail_common, 1, 0, tau
        )
        next_lift[tau] = local_cell(
            module, delayed, present, source_weight,
            objects, EXPECTED_BASES[1], 13, tau
        )
    require(
        next_lift == {3: (0, 0), 12: (0, 0)}
        and omega_digits(EXPECTED_BASES[1] + 13 * ARM_STEP)
        == omega_digits(EXPECTED_BASES[1]),
        "central quotient-descent hostile changed",
    )

    # Consecutive base addresses are the low-digit Y direction, not Z.
    # Their target/source cylinders coincide literally, but the two canonical
    # endpoint coefficient conventions alternate: carry-12 source survives on
    # 6714->6715, while carry-6 target survives on 6715->6716.
    cospan_rows = []
    for arm in range(3):
        require(
            arm_target_center(6714, arm)
            == arm_source_center(6715, arm)
            and arm_target_center(6715, arm)
            == arm_source_center(6716, arm),
            "successive-base source/target cylinder identity changed",
        )
        edge_one = (
            cospan_side_cell(context, "R", 6714, arm, 3, True),
            cospan_side_cell(context, "M", 6715, arm, 3, False),
        )
        edge_two = (
            cospan_side_cell(context, "R", 6715, arm, 12, True),
            cospan_side_cell(context, "L", 6716, arm, 12, False),
        )
        require(
            edge_one == (
                (0, 0),
                (
                    Fraction(EXPECTED_CONTENTS[1], ADDRESS_MODULUS),
                    EXPECTED_CONTENTS[1],
                ),
            )
            and edge_two == (
                (
                    Fraction(EXPECTED_CONTENTS[1], ADDRESS_MODULUS),
                    EXPECTED_CONTENTS[1],
                ),
                (0, 0),
            ),
            "successive-base cospan carry boundary changed",
        )
        cospan_rows.append((arm, edge_one, edge_two))

    # The nearest naive identification of an arm shift with a lawful target
    # action already fails on exact physical support.  At (e,tau)=(1,3),
    # arm zero is positive at sigma=0, while arm one is empty after either
    # sigma->sigma+1 or tau->tau+1.
    (
        module, delayed, present, source, clocks,
        source_weight, rail_common,
    ) = context
    base = EXPECTED_BASES[1]
    sigma0_tau3 = carrier_objects(
        module, source, clocks, rail_common, 1, 0, 3
    )
    sigma1_tau3 = carrier_objects(
        module, source, clocks, rail_common, 1, 1, 3
    )
    sigma0_tau4 = carrier_objects(
        module, source, clocks, rail_common, 1, 0, 4
    )
    m00, c00 = local_cell(
        module, delayed, present, source_weight,
        sigma0_tau3, base, 0, 3
    )
    m_sigma, c_sigma = local_cell(
        module, delayed, present, source_weight,
        sigma1_tau3, base, 1, 3
    )
    m_tau, c_tau = local_cell(
        module, delayed, present, source_weight,
        sigma0_tau4, base, 1, 4
    )
    require(
        m00 > 0 and c00 > 0
        and (m_sigma, c_sigma) == (0, 0)
        and (m_tau, c_tau) == (0, 0),
        "naive target-label/arm-shift hostile changed",
    )

    # Coordinate-only Heisenberg check.  In the normalized frame
    # q=(1,0), t=(0,1), an endpoint origin R_0=(w,v) has Delta=v.
    # The centre maps R_0 to R_0+q=(w+1,v), order 13.  On Delta=1 this
    # endpoint is L=R_0+q.  No semantic/address datum enters this action.
    heisenberg_rows = []
    for v in range(P):
        cycle = tuple(((w + 1) % P, v) for w in range(P))
        require(len(set(cycle)) == P, "central endpoint cycle lost order 13")
        if v == 1:
            require(
                all(
                    ((w + 1) % P, v)
                    == (((w, v)[0] + 1) % P, (w, v)[1])
                    for w in range(P)
                ),
                "Delta-one endpoint edge changed",
            )
        heisenberg_rows.append((v, cycle[0], cycle[-1]))

    print("THM-2782 SEMANTIC RIGHT-ARM / BOCKSTEIN ATTACHMENT")
    print("status=THM-2782 FINITE-EXACT candidate; no LRC conclusion")
    print(
        "universe="
        f"T={T}; address=Z/{ADDRESS_MODULUS}; "
        f"source_center={{Z+7n/{ADDRESS_MODULUS}}}; "
        f"target_center=source+7/{ADDRESS_MODULUS}; "
        f"radius={relative.Q_RADIUS}"
    )
    print(
        "carrier_formula="
        "E3 intersect clock_e; subtract "
        "D_q1(-sigma/13),D_c2(+sigma/13),"
        "D_q2(-tau/13),D_c3(+tau/13); "
        "R=(T_tau^-1 A)\\A; push R to target; "
        "retain c3-root1,D6-Q_(3,{1,2}),carry6"
    )
    for (
        e, base, metadata, content, bockstein_unit, mass,
        primitive, q3_decoded,
    ) in witness_rows:
        print(
            f"clock={e} base_n={base} base_mod13={base % 13} "
            f"arm_metadata={metadata}"
        )
        print(
            f"clock={e} raw_content={content} v13=1 "
            f"(content/13 mod13)={bockstein_unit} "
            f"positive_mass_each_cell={mass} "
            f"positive_cells={sum(sum(row) for row in primitive)}"
        )
        print(f"clock={e} primitive_arm_by_tau={primitive}")
        print(f"clock={e} decoded_Q3_mod13={q3_decoded}")
    print(
        "clock1_arm_contrast_W="
        f"{w}; inverse={inverse}; W*inverse={product}=delta0+2*N13"
    )
    print(
        "clock1_actual_Bockstein_Q3_arm0_all13="
        f"{EXPECTED_E1_Q3_BETA_DECODED}; nonzero=13/13"
    )
    print(
        "arm_address_step="
        f"n->n+13; physical_shift=91/{ADDRESS_MODULUS}; "
        f"full_address_order={P**5}; Omega=Z/169 order=13"
    )
    print(
        "endpoint_central_Z="
        "fixed target q and Delta; R0->R0+q; "
        "in frame q=(1,0),t=(0,1): (v,w)->(v,w+1); order=13; "
        "at Delta=1, Z(R0)=L"
    )
    print(
        "Omega_Heisenberg="
        "n=v+13w; X:n->14n; Y:(v,w)->(v+1,w) "
        "(carry-suppressed); Z:n->n+13; X*Y=Z*Y*X"
    )
    print(
        f"clock1_central_lift_metadata={central_metadata}"
    )
    print(
        f"clock1_central_lift_primitive_by_(w,tau)={central_primitive}"
    )
    print(
        "clock1_central_lift_decoded_columns="
        "13/13 nonconstant for the printed integer lift of K_beta, hence "
        "every such column has all 12 primitive complex Z-characters"
    )
    print(
        "central_descent_hostile="
        f"j=0 Omega={omega_digits(EXPECTED_BASES[1])} "
        f"raw_tau3={central_raw[0][3]} raw_tau12={central_raw[0][12]}; "
        f"j=13 same_Omega="
        f"{omega_digits(EXPECTED_BASES[1] + 13 * ARM_STEP)} "
        f"cells={next_lift}"
    )
    print(
        "successive_base_Y_cospan="
        "target_cylinder(n)=source_cylinder(n+1) exactly; "
        f"rows={tuple(cospan_rows)}"
    )
    print(
        "naive_shift_hostile="
        f"(sigma,tau,arm)=(0,3,0) mass={m00},coeff={c00}; "
        f"(1,3,1)=({m_sigma},{c_sigma}); "
        f"(0,4,1)=({m_tau},{c_tau})"
    )
    print(
        "typing_verdict=positive physical semantic-arm attachment to the "
        "right-wing coefficient Bockstein.  Modulo 169 the arm step is "
        "exactly the endpoint central Z direction, and the three printed "
        "arms are one same-semantic/carry central segment.  The coefficient "
        "does not descend around the full quotient cycle; Y cospan endpoint "
        "coefficients alternate zero/nonzero; both naive lawful-target-label "
        "identifications fail."
    )
    print(
        "scope=central-Z segment and decoded central charge, not a full "
        "H13-equivariant endpoint allocation: carry-suppressed Y is not a "
        "physical translation, the high address digit prevents Omega "
        "descent, and full factor/current covariance is absent.  No "
        "THM2542 C7-chart identification, root-deck intertwiner, or row "
        "exclusion."
    )
    print(
        "label_hostile=Omega low digit v=7, predecessor carry=6, and the "
        "THM2542 marker/root-deck transition are three distinct typed "
        "coordinates; the numerical congruence 7=-6 mod13 identifies none "
        "of them."
    )
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
