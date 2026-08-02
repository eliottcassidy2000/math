#!/usr/bin/env python3
"""Exact hostile companion for candidate THM-3065.

Only integer and rational arithmetic is used.  The analytic Selberg--Schur
identity is proved in the theorem; this companion independently checks its
sign consequence on arbitrary offset sets, the elementary one-sided
alternant formula, the Gregory--Newton wall, and the THM-3062 carrier ledger.
"""

from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from math import comb, factorial


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def rising(x, length):
    out = Fraction(1)
    for step in range(length):
        out *= x + step
    return out


def falling(x, length):
    out = Fraction(1)
    for step in range(length):
        out *= x - step
    return out


def determinant(matrix):
    """Fraction-free Bareiss determinant over Q."""
    a = [[Fraction(entry) for entry in row] for row in matrix]
    size = len(a)
    if size == 0:
        return Fraction(1)
    sign = 1
    previous = Fraction(1)
    for column in range(size - 1):
        pivot = next((row for row in range(column, size) if a[row][column]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            a[column], a[pivot] = a[pivot], a[column]
            sign = -sign
        value = a[column][column]
        for row in range(column + 1, size):
            for j in range(column + 1, size):
                a[row][j] = (
                    a[row][j] * value - a[row][column] * a[column][j]
                ) / previous
        previous = value
    return sign * a[-1][-1]


def sign(value):
    return int(value > 0) - int(value < 0)


def moment(a, b, index):
    return rising(b, index) / rising(a, index)


def universal_factor(a, b, order):
    out = Fraction(1)
    for j in range(order):
        out *= rising(a - b, j)
    return out


def vandermonde(values):
    out = Fraction(1)
    for i, left in enumerate(values):
        for right in values[i + 1 :]:
            out *= right - left
    return out


def generalized_minor(a, b, rows, columns):
    return determinant(
        [[moment(a, b, row + column) for column in columns] for row in rows]
    )


def one_sided_formula(a, b, rows, start):
    order = len(rows)
    out = vandermonde(rows) * universal_factor(a, b, order)
    for row in rows:
        out *= moment(a, b, row + start)
        out /= rising(a + row + start, order - 1)
    return out


def contiguous_formula(a, b, start, order):
    out = Fraction(1)
    for j in range(order):
        out *= factorial(j) * rising(b, start + j) * rising(a - b, j)
        out /= rising(a, start + order + j - 1)
    return out


def inventory_moment(base, inventory, index):
    out = Fraction(1)
    for shift, exponent in inventory.items():
        factor = rising(base + shift, index)
        if exponent >= 0:
            out *= factor**exponent
        else:
            out /= factor ** (-exponent)
    return out


def consolidate(entries):
    out = {}
    for shift, exponent in entries:
        out[shift] = out.get(shift, 0) + exponent
    return {shift: exponent for shift, exponent in out.items() if exponent}


def ordered_prefixes(inventory):
    total = 0
    out = []
    for shift in sorted(inventory):
        total += inventory[shift]
        out.append((shift, total))
    return out


def canonical_mesh_moment(base, inventory, index):
    nodes = sorted(inventory)
    prefixes = ordered_prefixes(inventory)
    require(all(value >= 0 for _, value in prefixes), "negative mesh prefix")
    out = rising(base + nodes[-1], index) ** prefixes[-1][1]
    for position in range(len(nodes) - 1):
        left = nodes[position]
        right = nodes[position + 1]
        flow = prefixes[position][1]
        out *= (
            rising(base + left, index) / rising(base + right, index)
        ) ** flow
    return out


def actual_carrier_inventory(gap):
    base_entries = [
        (Fraction(0), 26),
        (Fraction(1, 4), 6),
        (Fraction(1, 3), 8),
        (Fraction(1, 2), 6),
        (Fraction(2, 3), 8),
        (Fraction(3, 4), 6),
        (Fraction(1), -14),
    ]
    # W_(C+h)/W_C, after removal of its positive C-independent scale, is
    # 26 transfers 0->h and 20 transfers 1->h+1 on the shape line.
    return consolidate(
        base_entries
        + [
            (Fraction(0), -26),
            (Fraction(gap), 26),
            (Fraction(1), -20),
            (Fraction(gap + 1), 20),
        ]
    )


def feed(digest, *items):
    digest.update("|".join(str(item) for item in items).encode("ascii"))
    digest.update(b"\n")


print("THM-3065 RECIPROCAL BETA GENERALIZED HANKEL WALL")
digest = sha256()

# Gregory--Newton is the elementary source of the integer walls:
# Delta^k m_n = d(d-1)...(d-k+1) m_n/(a+n)_k, d=b-a.
parameters = [
    ("forward", Fraction(5, 4), Fraction(1, 3)),
    ("forward", Fraction(7, 3), Fraction(5, 4)),
    ("forward", Fraction(9, 2), Fraction(3, 2)),
    ("equal", Fraction(1, 2), Fraction(1, 2)),
    ("equal", Fraction(5, 3), Fraction(5, 3)),
]
for position, gap in enumerate(
    (
        Fraction(1, 4),
        Fraction(3, 4),
        Fraction(5, 4),
        Fraction(7, 4),
        Fraction(9, 4),
        Fraction(11, 4),
        Fraction(13, 4),
    )
):
    a = (Fraction(1, 2), Fraction(4, 3), Fraction(7, 3))[position % 3]
    parameters.append(("nonintegral", a, a + gap))
for gap in range(1, 5):
    a = (Fraction(2, 3), Fraction(5, 3))[gap % 2]
    parameters.append(("integer", a, a + gap))

newton_cells = 0
integer_terminations = 0
for category, a, b in parameters:
    d = b - a
    for index in range(9):
        for order in range(8):
            direct = sum(
                (-1) ** (order - j)
                * comb(order, j)
                * moment(a, b, index + j)
                for j in range(order + 1)
            )
            closed = falling(d, order) * moment(a, b, index)
            closed /= rising(a + index, order)
            require(direct == closed, "Gregory--Newton identity failed")
            if category == "integer" and order > int(d):
                require(direct == 0, "integer Newton termination failed")
                integer_terminations += 1
            feed(digest, "N", category, a, b, index, order, direct)
            newton_cells += 1
print(
    f"gregory_newton_cells={newton_cells} "
    f"integer_termination_cells={integer_terminations}"
)

# Full arbitrary-offset hostile bank.  The determinant itself is computed;
# only its sign/rank is compared with the universal product.
minor_counts = {name: 0 for name in ("forward", "equal", "nonintegral", "integer")}
for category, a, b in parameters:
    for order in range(2, 6):
        offset_sets = tuple(combinations(range(7), order))
        expected = sign(universal_factor(a, b, order))
        for rows in offset_sets:
            for columns in offset_sets:
                value = generalized_minor(a, b, rows, columns)
                require(sign(value) == expected, "arbitrary-offset signature failed")
                feed(digest, "G", category, a, b, rows, columns, value)
                minor_counts[category] += 1
print(
    "generalized_minor_cells="
    f"{sum(minor_counts.values())} "
    + " ".join(f"{key}={minor_counts[key]}" for key in minor_counts)
)

# The elementary alternant and its contiguous specialization are equality
# checks, not merely sign checks.
alternant_parameters = parameters[0:3] + parameters[5:10] + parameters[12:16]
alternant_cells = 0
contiguous_cells = 0
for category, a, b in alternant_parameters:
    for order in range(1, 6):
        for rows in combinations(range(7), order):
            for start in range(4):
                columns = tuple(start + j for j in range(order))
                direct = generalized_minor(a, b, rows, columns)
                closed = one_sided_formula(a, b, rows, start)
                require(direct == closed, "one-sided alternant formula failed")
                feed(digest, "A", category, a, b, rows, start, direct)
                alternant_cells += 1
    for order in range(1, 8):
        for start in range(6):
            rows = tuple(start + j for j in range(order))
            columns = tuple(range(order))
            direct = generalized_minor(a, b, rows, columns)
            closed = contiguous_formula(a, b, start, order)
            require(direct == closed, "contiguous formula failed")
            contiguous_cells += 1
print(
    f"one_sided_alternant_cells={alternant_cells} "
    f"contiguous_formula_cells={contiguous_cells}"
)

# Explicitly cross every positive integer gap wall from both sides, on a
# nonconsecutive row/column pair.  The factor at d=m has multiplicity
# max(0,r-m-1) in the universal product.
wall_cells = 0
wall_table = []
for integer_gap in range(1, 5):
    for order in range(2, 8):
        rows = tuple(2 * j for j in range(order))
        columns = tuple(3 * j + (j % 2) for j in range(order))
        row = []
        for d in (
            Fraction(7 * integer_gap - 1, 7),
            Fraction(integer_gap),
            Fraction(7 * integer_gap + 1, 7),
        ):
            a = Fraction(5, 4)
            b = a + d
            value = generalized_minor(a, b, rows, columns)
            expected = sign(universal_factor(a, b, order))
            require(sign(value) == expected, "integer-wall hostile failed")
            row.append(sign(value))
            feed(digest, "W", integer_gap, order, d, value)
            wall_cells += 1
        multiplicity = max(0, order - integer_gap - 1)
        require((row[1] == 0) == (multiplicity > 0), "wall multiplicity/rank mismatch")
        if multiplicity:
            require(row[0] * row[2] == (-1) ** multiplicity,
                    "wall crossing parity failed")
        wall_table.append((integer_gap, order, tuple(row), multiplicity))
print(
    f"integer_wall_cells={wall_cells} "
    f"wall_table_digest={sha256(repr(wall_table).encode('ascii')).hexdigest()}"
)

# The THM-3062 carrier, after removal of positive constant/geometric scales,
# first has this rational mesh inventory.  Its prefixes are the exact cut
# flows.  The actual fixed terminal gap contributes two forward transfers.
base_inventory = consolidate(
    [
        (Fraction(0), 26),
        (Fraction(1, 4), 6),
        (Fraction(1, 3), 8),
        (Fraction(1, 2), 6),
        (Fraction(2, 3), 8),
        (Fraction(3, 4), 6),
        (Fraction(1), -14),
    ]
)
base_prefixes = tuple(value for _, value in ordered_prefixes(base_inventory))
require(base_prefixes == (26, 32, 40, 46, 54, 60, 46),
        "four-slot base carrier prefix changed")
mesh_cells = 0
actual_gap_cells = 0
for base in (Fraction(1, 2), Fraction(1), Fraction(5, 2)):
    for index in range(9):
        direct = inventory_moment(base, base_inventory, index)
        canonical = canonical_mesh_moment(base, base_inventory, index)
        require(direct == canonical, "rational-mesh carrier flow failed")
        mesh_cells += 1
    for gap in range(1, 6):
        inventory = actual_carrier_inventory(gap)
        prefixes = ordered_prefixes(inventory)
        require(all(value >= 0 for _, value in prefixes), "actual carrier prefix failed")
        require(prefixes[-1][1] == 46, "actual carrier terminal prefix changed")
        for index in range(8):
            direct = inventory_moment(base, inventory, index)
            canonical = canonical_mesh_moment(base, inventory, index)
            require(direct == canonical, "actual carrier canonical flow failed")
            actual_gap_cells += 1
print(
    f"mesh_prefix_cells={mesh_cells} base_prefixes={','.join(map(str, base_prefixes))} "
    f"actual_gap_cells={actual_gap_cells} terminal_prefix=46"
)

# Finite-exact scout only: reciprocal actual carriers exhibit the checkerboard
# Hankel signature in this hostile bank.  The theorem promotes only order two;
# closure of that higher signature under these Hadamard products remains open.
reciprocal_cells = 0
for base in (Fraction(1), Fraction(3, 2)):
    for gap in range(1, 5):
        inventory = actual_carrier_inventory(gap)
        values = [Fraction(1, inventory_moment(base, inventory, index)) for index in range(13)]
        for order in range(2, 6):
            expected = (-1) ** (order * (order - 1) // 2)
            offset_sets = tuple(combinations(range(7), order))
            for rows in offset_sets:
                for columns in offset_sets:
                    value = determinant(
                        [[values[row + column] for column in columns] for row in rows]
                    )
                    require(sign(value) == expected,
                            "finite reciprocal-carrier signature changed")
                    reciprocal_cells += 1
print(
    f"reciprocal_carrier_minor_cells={reciprocal_cells} "
    "finite_signature=(-1)^binom(order,2)"
)

# A nonpositive dual prefix is enough for every order-two sign, but not for
# checkerboard sign regularity.  This exact zero-cut hostile is the first
# failed lift.  Two nearest strict-prefix repairs are positive controls only.
dual_shapes = (Fraction(1, 7), Fraction(2, 3), Fraction(20))
dual_hostile = (-1, 1, -1)
require(tuple(value for _, value in ordered_prefixes(
    dict(zip(dual_shapes, dual_hostile)))) == (-1, 0, -1),
        "dual-prefix hostile ledger changed")
dual_values = []
for index in range(14):
    value = Fraction(1)
    for shape, exponent in zip(dual_shapes, dual_hostile):
        factor = rising(shape, index)
        value = value * factor**exponent if exponent >= 0 else value / factor ** (-exponent)
    dual_values.append(value)
minimal_rows = (0, 1, 2)
minimal_columns = (0, 1, 2)
minimal_determinant = determinant(
    [[dual_values[row + column] for column in minimal_columns] for row in minimal_rows]
)
require(minimal_determinant == Fraction(4914161, 84138683904000),
        "minimal dual-prefix order-three hostile changed")

order_four_inventory = (-3, 3, -1)
order_four_values = []
for index in range(14):
    value = Fraction(1)
    for shape, exponent in zip(dual_shapes, order_four_inventory):
        factor = rising(shape, index)
        value = value * factor**exponent if exponent >= 0 else value / factor ** (-exponent)
    order_four_values.append(value)
hostile_rows = (0, 1, 2, 7)
hostile_columns = (0, 1, 2, 6)
dual_determinant = determinant(
    [[order_four_values[row + column] for column in hostile_columns] for row in hostile_rows]
)
expected_dual_determinant = Fraction(
    -89813384398251034317167637460105582061544597308528939184258027657887067231528473035487,
    2712795677868954759419954836150485676184464443881644193874894868030111422768506341215436800000000000000,
)
require(dual_determinant == expected_dual_determinant,
        "dual-prefix order-four hostile changed")
dual_h2_cells = 0
for rows in combinations(range(7), 2):
    for columns in combinations(range(7), 2):
        value = determinant(
            [[dual_values[row + column] for column in columns] for row in rows]
        )
        require(value < 0, "dual-prefix order-two survivor failed")
        dual_h2_cells += 1

strict_repair_cells = 0
for strict_inventory in ((-2, 1, -1), (-3, 2, -1)):
    prefixes = tuple(value for _, value in ordered_prefixes(
        dict(zip(dual_shapes, strict_inventory))))
    require(all(value < 0 for value in prefixes), "strict-prefix control is not strict")
    values = []
    for index in range(13):
        value = Fraction(1)
        for shape, exponent in zip(dual_shapes, strict_inventory):
            factor = rising(shape, index)
            value = value * factor**exponent if exponent >= 0 else value / factor ** (-exponent)
        values.append(value)
    for order in range(2, 6):
        expected = (-1) ** (order * (order - 1) // 2)
        offset_sets = tuple(combinations(range(7), order))
        for rows in offset_sets:
            for columns in offset_sets:
                value = determinant(
                    [[values[row + column] for column in columns] for row in rows]
                )
                require(sign(value) == expected,
                        "strict dual-prefix finite control changed")
                strict_repair_cells += 1
feed(digest, "D", dual_shapes, dual_hostile, minimal_rows, minimal_columns,
     minimal_determinant, order_four_inventory, hostile_rows, hostile_columns,
     dual_determinant)
print(
    f"dual_prefix_hostile_signs=H3:{sign(minimal_determinant)},H4:{sign(dual_determinant)} "
    f"h2_survivor_cells={dual_h2_cells} strict_repair_cells={strict_repair_cells}"
)

print(f"exact_value_digest={digest.hexdigest()}")
print("scope=single_ratio_all_orders;carrier_reciprocal_higher_orders_finite_only")
print("all_exact_checks=PASS")
