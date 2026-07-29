#!/usr/bin/env python3
"""Exact source-endpoint-mask/coboundary audit above THM-2825.

Canonical THM-2859 companion.  The THM-2825 physical companion builds all
169 endpoint present sets as
large interval unions.  This scratch audit instead evaluates the defining
periodic inequalities directly on each selected interval.  It first matches
the complete stored R/M1/M2 source-mask and translation censuses, then
decorates every vertex of all 587 cofiber-rooted half-step paths.

Three deliberately different notions of "coboundary" are separated:

1. additive path defect: D(x->y)=1_E(y)-1_E(x), which is tautologically the
   coboundary of the full endpoint-mask 0-cochain;
2. address-translation derivative: D=(T_delta-I)phi on F_13^2, which has
   zero sum on every delta-line and in particular zero total augmentation;
3. translation gauge: E(y)=T_delta E(x), a prerequisite for a V-valued edge
   gauge and hence for a V-valued vertex potential along the forest.

All arithmetic is exact.  No Python ``assert`` statement is used.
"""

from collections import Counter, defaultdict
from hashlib import sha256
from pathlib import Path
import ast
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


PRIMARY_NAME = "lrc14_nearest_half_step_common_right_collar_thm2825.py"
PRIMARY_SHA256 = (
    "bd9ffe7f6815b5c563bd483c300118fbdd683f3d9303babbab7912e031747c9a"
)
primary_payload = (COMP / PRIMARY_NAME).read_bytes().replace(b"\r\n", b"\n")
require(
    sha256(primary_payload).hexdigest() == PRIMARY_SHA256,
    "pinned THM-2825 primary companion changed",
)


import lrc14_nearest_half_step_common_right_collar_thm2825 as collar


P = 13
H = collar.H
T = collar.copies.T
LENGTH = collar.copies.LENGTH
SHIFT = collar.copies.SHIFT
ENDPOINT = collar.copies.endpoint_base
KEYS = tuple(ENDPOINT.KEYS)
KEY_INDEX = {key: index for index, key in enumerate(KEYS)}
W = tuple(ENDPOINT.endpoint.W)
PATTERN = dict(ENDPOINT.PAT_E3)

require(
    KEYS == tuple((a, b) for a in range(P) for b in range(P)),
    "endpoint-address enumeration changed",
)
require(
    PATTERN
    == {
        0: "gout",
        1: "out",
        2: "out",
        3: "out",
        4: "out",
        5: "out",
        6: "out",
        7: "out",
        8: "in",
    },
    "endpoint E3 pattern changed",
)
require(
    W == (1, 14, 27, 40, 53, 66, 13, 13**3, 2 * 13**5),
    "endpoint speed row changed",
)


def periodic_parameters(index, ell_value):
    """Return period, base, width for the defining danger/present comb."""
    if PATTERN[index] == "gout":
        denominator = 91
        low = -13 - 7 * ell_value
        high = 13 - 7 * ell_value
    else:
        denominator = 182
        low = -13 - 14 * ell_value
        high = 13 - 14 * ell_value
    unit = T // (denominator * W[index])
    period = denominator * unit
    base = (low % denominator) * unit
    width = (high - low) * unit
    return period, base, width


def nearby_boundaries(left, right, period, base, width):
    """Interior start/stop boundaries of one periodic half-open comb."""
    first = (left - base) // period - 2
    last = (right - base) // period + 2
    rows = []
    for number in range(first, last + 1):
        start = base + number * period
        stop = start + width
        if left < start < right:
            rows.append(start)
        if left < stop < right:
            rows.append(stop)
    return tuple(rows)


def midpoint_in_comb(left, right, period, base, width):
    """Membership of the exact midpoint, represented with denominator two."""
    twice_x = left + right
    twice_period = 2 * period
    twice_base = 2 * base
    quotient = (twice_x - twice_base) // twice_period
    residue = twice_x - twice_base - quotient * twice_period
    return 0 <= residue < 2 * width


def address_indicator(interval, address):
    """Exact 0/1/partial indicator without materializing a present-set union."""
    left, right = interval
    ell = ENDPOINT.REPS[address]
    factor_data = tuple(
        periodic_parameters(index, ell[index]) for index in range(len(W))
    )
    cuts = {left, right}
    for period, base, width in factor_data:
        cuts.update(nearby_boundaries(left, right, period, base, width))
    cuts = tuple(sorted(cuts))
    statuses = []
    for start, stop in zip(cuts, cuts[1:]):
        deep_period, deep_base, deep_width = factor_data[8]
        live = midpoint_in_comb(
            start, stop, deep_period, deep_base, deep_width
        )
        if live:
            for index in range(8):
                period, base, width = factor_data[index]
                if midpoint_in_comb(start, stop, period, base, width):
                    live = False
                    break
        statuses.append(live)
    if all(statuses):
        return 1
    if not any(statuses):
        return 0
    return 2


mask_cache = {}


def endpoint_mask(interval):
    if interval not in mask_cache:
        row = tuple(address_indicator(interval, address) for address in KEYS)
        require(2 not in row, f"analytic endpoint mask cuts interval {interval}")
        mask_cache[interval] = row
    return mask_cache[interval]


def translate_mask(mask, delta):
    da, db = delta
    return tuple(
        mask[KEY_INDEX[((a + da) % P, (b + db) % P)]]
        for a, b in KEYS
    )


translation_orbit_cache = {}


def translation_orbit(mask):
    if mask not in translation_orbit_cache:
        translation_orbit_cache[mask] = {
            translate_mask(mask, delta): []
            for delta in KEYS
        }
        for delta in KEYS:
            translation_orbit_cache[mask][translate_mask(mask, delta)].append(
                delta
            )
    return translation_orbit_cache[mask]


def translation_deltas(source, target):
    return tuple(translation_orbit(source).get(target, ()))


CYCLIC_ORDERS = ("13a+q", "a+13q")
cyclic_translation_cache = {}


def cyclic_encode(address, order):
    a, q = address
    if order == "13a+q":
        return 13 * a + q
    if order == "a+13q":
        return a + 13 * q
    raise RuntimeError(f"unknown cyclic digit order {order}")


def cyclic_decode(number, order):
    low = number % 13
    high = number // 13
    if order == "13a+q":
        return high, low
    if order == "a+13q":
        return low, high
    raise RuntimeError(f"unknown cyclic digit order {order}")


def cyclic_translate_mask(mask, shift, order):
    return tuple(
        mask[KEY_INDEX[cyclic_decode(
            (cyclic_encode(address, order) + shift) % (P**2),
            order,
        )]]
        for address in KEYS
    )


def cyclic_translation_deltas(source, target, order):
    key = (source, target, order)
    if key not in cyclic_translation_cache:
        if sum(source) != sum(target):
            answer = ()
        else:
            answer = tuple(
                shift
                for shift in range(P**2)
                if cyclic_translate_mask(source, shift, order) == target
            )
        cyclic_translation_cache[key] = answer
    return cyclic_translation_cache[key]


def support(mask):
    return frozenset(
        address for address, value in zip(KEYS, mask) if value
    )


factor_cache = {}


def cartesian_factor(mask):
    if mask in factor_cache:
        return factor_cache[mask]
    points = support(mask)
    if not points:
        answer = (frozenset(), frozenset(), True)
    else:
        first = frozenset(a for a, _b in points)
        second = frozenset(b for _a, b in points)
        answer = (
            first,
            second,
            points
            == frozenset((a, b) for a in first for b in second),
        )
    factor_cache[mask] = answer
    return answer


PROJECTIVE_DIRECTIONS = tuple(
    [(1, slope) for slope in range(P)] + [(0, 1)]
)


derivative_direction_cache = {}
cyclic_derivative_cache = {}
derivative_direction_mod2_cache = {}
cyclic_derivative_mod2_cache = {}


def translation_derivative_directions(defect):
    """Directions d for which defect=(T_d-I)phi over the integers.

    On each order-13 d-line this is equivalent to zero line sum.
    """
    if defect in derivative_direction_cache:
        return derivative_direction_cache[defect]
    values = {key: defect[index] for index, key in enumerate(KEYS)}
    valid = []
    for da, db in PROJECTIVE_DIRECTIONS:
        seen = set()
        good = True
        for origin in KEYS:
            if origin in seen:
                continue
            line = tuple(
                ((origin[0] + k * da) % P, (origin[1] + k * db) % P)
                for k in range(P)
            )
            seen.update(line)
            if sum(values[point] for point in line) != 0:
                good = False
                break
        if good:
            valid.append((da, db))
    answer = tuple(valid)
    derivative_direction_cache[defect] = answer
    return answer


def translation_derivative_directions_mod2(defect):
    """Directions d for which defect=(T_d-I)phi over F_2."""
    parity_defect = tuple(value % 2 for value in defect)
    if parity_defect in derivative_direction_mod2_cache:
        return derivative_direction_mod2_cache[parity_defect]
    values = {
        key: parity_defect[index] for index, key in enumerate(KEYS)
    }
    valid = []
    for da, db in PROJECTIVE_DIRECTIONS:
        seen = set()
        good = True
        for origin in KEYS:
            if origin in seen:
                continue
            line = tuple(
                ((origin[0] + k * da) % P, (origin[1] + k * db) % P)
                for k in range(P)
            )
            seen.update(line)
            if sum(values[point] for point in line) % 2:
                good = False
                break
        if good:
            valid.append((da, db))
    answer = tuple(valid)
    derivative_direction_mod2_cache[parity_defect] = answer
    return answer


def cyclic_derivative_shifts(defect, order):
    """Nonzero C_169 shifts for which defect=(T_shift-I)phi over Z."""
    key = (defect, order)
    if key in cyclic_derivative_cache:
        return cyclic_derivative_cache[key]
    if sum(defect) != 0:
        answer = ()
    elif not any(defect):
        answer = tuple(range(1, P**2))
    else:
        values = {
            cyclic_encode(address, order): defect[index]
            for index, address in enumerate(KEYS)
        }
        valid = []
        for shift in range(1, P**2):
            seen = set()
            good = True
            for origin in range(P**2):
                if origin in seen:
                    continue
                orbit = []
                cursor = origin
                while cursor not in seen:
                    seen.add(cursor)
                    orbit.append(cursor)
                    cursor = (cursor + shift) % (P**2)
                if sum(values[point] for point in orbit) != 0:
                    good = False
                    break
            if good:
                valid.append(shift)
        answer = tuple(valid)
    cyclic_derivative_cache[key] = answer
    return answer


def cyclic_derivative_shifts_mod2(defect, order):
    """Nonzero C_169 shifts admitting a derivative potential over F_2."""
    parity_defect = tuple(value % 2 for value in defect)
    key = (parity_defect, order)
    if key in cyclic_derivative_mod2_cache:
        return cyclic_derivative_mod2_cache[key]
    values = {
        cyclic_encode(address, order): parity_defect[index]
        for index, address in enumerate(KEYS)
    }
    valid = []
    for shift in range(1, P**2):
        seen = set()
        good = True
        for origin in range(P**2):
            if origin in seen:
                continue
            orbit = []
            cursor = origin
            while cursor not in seen:
                seen.add(cursor)
                orbit.append(cursor)
                cursor = (cursor + shift) % (P**2)
            if sum(values[point] for point in orbit) % 2:
                good = False
                break
        if good:
            valid.append(shift)
    answer = tuple(valid)
    cyclic_derivative_mod2_cache[key] = answer
    return answer


def mod2_direction_potential(defect, direction):
    """One normalized phi with phi(x+d)+phi(x)=defect(x) over F_2."""
    parity_defect = {
        key: defect[index] % 2 for index, key in enumerate(KEYS)
    }
    da, db = direction
    potential = {}
    seen = set()
    for origin in KEYS:
        if origin in seen:
            continue
        orbit = tuple(
            ((origin[0] + k * da) % P, (origin[1] + k * db) % P)
            for k in range(P)
        )
        seen.update(orbit)
        potential[orbit[0]] = 0
        for point, successor in zip(orbit, orbit[1:]):
            potential[successor] = (
                potential[point] + parity_defect[point]
            ) % 2
        require(
            (
                potential[orbit[-1]]
                + potential[orbit[0]]
            ) % 2
            == parity_defect[orbit[-1]],
            "mod-two derivative potential failed to close",
        )
    require(
        all(
            (
                potential[((a + da) % P, (b + db) % P)]
                + potential[(a, b)]
            ) % 2
            == parity_defect[(a, b)]
            for a, b in KEYS
        ),
        "mod-two derivative potential failed pointwise",
    )
    return tuple(potential[address] for address in KEYS)


def defect_row(source, target):
    return tuple(b - a for a, b in zip(source, target))


def relation_type(source, target):
    left = support(source)
    right = support(target)
    if left == right:
        return "equal"
    if left < right:
        return "birth"
    if right < left:
        return "death"
    return "crossing"


def axis_type(source, target):
    a0, b0, product0 = cartesian_factor(source)
    a1, b1, product1 = cartesian_factor(target)
    require(product0 and product1, "non-Cartesian endpoint mask")
    if source == target:
        return "equal"
    if not a0 or not b0 or not a1 or not b1:
        return "empty_boundary"
    changed = (a0 != a1, b0 != b1)
    if changed == (True, False):
        return "first_axis"
    if changed == (False, True):
        return "second_axis"
    return "both_axes"


def mask_description(mask):
    first, second, product = cartesian_factor(mask)
    return (
        tuple(sorted(first)),
        tuple(sorted(second)),
        len(support(mask)),
        product,
    )


def path_geometry():
    (
        _module,
        _rails,
        _present,
        details,
        full_module,
        e3,
        clocks,
        _q_pairs,
        _delayed,
        _source_weight,
        _target_weight,
        _rail_common,
    ) = collar.copies.physical_setup()
    paths = []
    for clock in range(7):
        for sigma in collar.COMMON_S:
            for target in collar.COMMON_T:
                _source, _target, common, right = collar.cell_objects(
                    details,
                    full_module,
                    e3,
                    clocks,
                    clock,
                    sigma,
                    target,
                )
                if not right:
                    continue
                common_by_left = {piece[0]: piece for piece in common}
                for root_index, root in enumerate(right):
                    vertices = [root]
                    cursor = root[0] + H
                    while cursor in common_by_left:
                        vertices.append(common_by_left[cursor])
                        cursor += H
                    require(
                        len(vertices) >= 3,
                        "rooted endpoint path lost its two collars",
                    )
                    paths.append(
                        ((clock, sigma, target), root_index, tuple(vertices))
                    )
    require(len(paths) == 587, "rooted path census changed")
    return tuple(paths)


def edge_audit(paths, jump):
    hamming = Counter()
    augmentation = Counter()
    relations = Counter()
    axes = Counter()
    translation_counts = Counter()
    derivative_direction_counts = Counter()
    derivative_direction_mod2_counts = Counter()
    cyclic_translation_counts = {
        order: Counter() for order in CYCLIC_ORDERS
    }
    cyclic_derivative_counts = {
        order: Counter() for order in CYCLIC_ORDERS
    }
    cyclic_derivative_mod2_counts = {
        order: Counter() for order in CYCLIC_ORDERS
    }
    pair_types = Counter()
    total = 0
    first_weight_obstruction = None
    first_small_obstruction = None
    first_equal_weight_nontranslation = None
    first_nonzero_translation = None
    first_mod2_potential = None
    solvable = 0

    for cell, root_index, vertices in paths:
        masks = tuple(endpoint_mask(piece[:2]) for piece in vertices)
        for source_index in range(len(vertices) - jump):
            source_piece = vertices[source_index]
            target_piece = vertices[source_index + jump]
            source_mask = masks[source_index]
            target_mask = masks[source_index + jump]
            defect = defect_row(source_mask, target_mask)
            distance = sum(a != b for a, b in zip(source_mask, target_mask))
            aug = sum(defect)
            deltas = translation_deltas(source_mask, target_mask)
            directions = translation_derivative_directions(defect)
            directions_mod2 = translation_derivative_directions_mod2(defect)
            cyclic_deltas = {
                order: cyclic_translation_deltas(
                    source_mask, target_mask, order
                )
                for order in CYCLIC_ORDERS
            }
            cyclic_derivatives = {
                order: cyclic_derivative_shifts(defect, order)
                for order in CYCLIC_ORDERS
            }
            cyclic_derivatives_mod2 = {
                order: cyclic_derivative_shifts_mod2(defect, order)
                for order in CYCLIC_ORDERS
            }
            relation = relation_type(source_mask, target_mask)
            axis = axis_type(source_mask, target_mask)
            total += 1
            hamming[distance] += 1
            augmentation[aug] += 1
            relations[relation] += 1
            axes[axis] += 1
            translation_counts[len(deltas)] += 1
            derivative_direction_counts[len(directions)] += 1
            derivative_direction_mod2_counts[len(directions_mod2)] += 1
            for order in CYCLIC_ORDERS:
                cyclic_translation_counts[order][
                    len(cyclic_deltas[order])
                ] += 1
                cyclic_derivative_counts[order][
                    len(cyclic_derivatives[order])
                ] += 1
                cyclic_derivative_mod2_counts[order][
                    len(cyclic_derivatives_mod2[order])
                ] += 1
            pair_types[
                (
                    mask_description(source_mask),
                    mask_description(target_mask),
                    distance,
                    aug,
                    relation,
                    axis,
                    len(deltas),
                    directions,
                    tuple(
                        (order, len(cyclic_deltas[order]))
                        for order in CYCLIC_ORDERS
                    ),
                    tuple(
                        (order, len(cyclic_derivatives[order]))
                        for order in CYCLIC_ORDERS
                    ),
                    len(directions_mod2),
                    tuple(
                        (order, len(cyclic_derivatives_mod2[order]))
                        for order in CYCLIC_ORDERS
                    ),
                )
            ] += 1
            witness = (
                cell,
                root_index,
                source_index,
                source_piece[:2],
                target_piece[:2],
                mask_description(source_mask),
                mask_description(target_mask),
                tuple(sorted(support(target_mask) ^ support(source_mask))),
                deltas,
                directions,
            )
            if aug and first_weight_obstruction is None:
                first_weight_obstruction = witness
            if (
                aug
                and (
                    first_small_obstruction is None
                    or distance
                    < len(first_small_obstruction[7])
                )
            ):
                first_small_obstruction = witness
            if not aug and source_mask != target_mask and not deltas:
                if first_equal_weight_nontranslation is None:
                    first_equal_weight_nontranslation = witness
            if deltas:
                solvable += 1
                nonzero = tuple(delta for delta in deltas if delta != (0, 0))
                if nonzero and source_mask != target_mask:
                    if first_nonzero_translation is None:
                        first_nonzero_translation = witness
            if (
                source_mask != target_mask
                and directions_mod2
                and first_mod2_potential is None
            ):
                potential = mod2_direction_potential(
                    defect, directions_mod2[0]
                )
                first_mod2_potential = (
                    witness,
                    directions_mod2[0],
                    tuple(
                        address
                        for address, value in zip(KEYS, potential)
                        if value
                    ),
                )

    return {
        "jump": jump,
        "total": total,
        "hamming": hamming,
        "augmentation": augmentation,
        "relations": relations,
        "axes": axes,
        "translation_counts": translation_counts,
        "derivative_direction_counts": derivative_direction_counts,
        "derivative_direction_mod2_counts":
            derivative_direction_mod2_counts,
        "cyclic_translation_counts": cyclic_translation_counts,
        "cyclic_derivative_counts": cyclic_derivative_counts,
        "cyclic_derivative_mod2_counts":
            cyclic_derivative_mod2_counts,
        "pair_types": pair_types,
        "solvable": solvable,
        "first_weight_obstruction": first_weight_obstruction,
        "first_small_obstruction": first_small_obstruction,
        "first_equal_weight_nontranslation":
            first_equal_weight_nontranslation,
        "first_nonzero_translation": first_nonzero_translation,
        "first_mod2_potential": first_mod2_potential,
    }


def matrix_multiply(left, right):
    a, b, c, d = left
    e, f, g, h = right
    return (
        (a * e + b * g) % P,
        (a * f + b * h) % P,
        (c * e + d * g) % P,
        (c * f + d * h) % P,
    )


def matrix_add(left, right):
    return tuple((a + b) % P for a, b in zip(left, right))


def matrix_vector(matrix, vector):
    a, b, c, d = matrix
    x, y = vector
    return ((a * x + b * y) % P, (c * x + d * y) % P)


IDENTITY2 = (1, 0, 0, 1)
ZERO2 = (0, 0, 0, 0)


def matrix_power_sum(matrix, exponent):
    """Return A^n and I+A+...+A^(n-1), exactly."""
    power = IDENTITY2
    total = ZERO2
    for _index in range(exponent):
        total = matrix_add(total, power)
        power = matrix_multiply(power, matrix)
    return power, total


def agl_fixed_power_count(exponent):
    """Count affine maps x->Ax+b whose exponent-th power is identity."""
    total = 0
    linear_parts = 0
    for a in range(P):
        for b in range(P):
            for c in range(P):
                for d in range(P):
                    if (a * d - b * c) % P == 0:
                        continue
                    matrix = (a, b, c, d)
                    power, geometric_sum = matrix_power_sum(
                        matrix, exponent
                    )
                    if power != IDENTITY2:
                        continue
                    linear_parts += 1
                    total += sum(
                        matrix_vector(geometric_sum, vector) == (0, 0)
                        for vector in KEYS
                    )
    return linear_parts, total


def validate_collar(paths):
    endpoint_source_census = {
        "R": Counter(),
        "M1": Counter(),
        "M2": Counter(),
    }
    pair_census = {"h": Counter(), "2h": Counter()}
    for _cell, _root_index, vertices in paths:
        masks = {
            "R": endpoint_mask(vertices[0][:2]),
            "M1": endpoint_mask(vertices[1][:2]),
            "M2": endpoint_mask(vertices[2][:2]),
        }
        for name, mask in masks.items():
            endpoint_source_census[name][(mask.count(0), mask.count(1))] += 1
        for label, name in (("h", "M1"), ("2h", "M2")):
            left = masks["R"]
            right = masks[name]
            pair_census[label][
                (
                    sum(a != b for a, b in zip(left, right)),
                    translation_deltas(left, right),
                )
            ] += 1

    expected_source = {
        "R": Counter({
            (59, 110): 21,
            (69, 100): 28,
            (70, 99): 16,
            (79, 90): 333,
            (88, 81): 108,
            (169, 0): 81,
        }),
        "M1": Counter({
            (59, 110): 28,
            (69, 100): 259,
            (70, 99): 95,
            (79, 90): 124,
            (88, 81): 81,
        }),
        "M2": Counter({
            (59, 110): 28,
            (69, 100): 259,
            (70, 99): 95,
            (79, 90): 124,
            (88, 81): 81,
        }),
    }
    expected_pair = Counter({
        (0, ((0, 0),)): 74,
        (9, ()): 187,
        (10, ()): 245,
        (81, ()): 81,
    })
    require(
        endpoint_source_census == expected_source,
        f"analytic source-mask census missed stored physical audit: "
        f"{endpoint_source_census}",
    )
    require(
        pair_census["h"] == expected_pair
        and pair_census["2h"] == expected_pair,
        f"analytic pair census missed stored physical audit: {pair_census}",
    )
    return endpoint_source_census, pair_census


def path_orbit_audit(paths):
    all_one_orbit = 0
    root_orbit_prefix = Counter()
    orbit_count = Counter()
    translation_solvable_h_edges = 0
    vertices = 0
    for _cell, _root_index, pieces in paths:
        masks = tuple(endpoint_mask(piece[:2]) for piece in pieces)
        canonical = tuple(
            min(translation_orbit(mask)) for mask in masks
        )
        distinct = len(set(canonical))
        orbit_count[distinct] += 1
        if distinct == 1:
            all_one_orbit += 1
        prefix = 1
        while (
            prefix < len(masks)
            and canonical[prefix] == canonical[0]
        ):
            prefix += 1
        root_orbit_prefix[prefix - 1] += 1
        vertices += len(masks)
        translation_solvable_h_edges += sum(
            bool(translation_deltas(left, right))
            for left, right in zip(masks, masks[1:])
        )
    # Removing every unsolved edge from a forest gives this many components.
    translation_gauge_components = (
        vertices - translation_solvable_h_edges
    )
    return {
        "vertices": vertices,
        "all_one_orbit": all_one_orbit,
        "orbit_count": orbit_count,
        "root_orbit_prefix": root_orbit_prefix,
        "translation_solvable_h_edges": translation_solvable_h_edges,
        "translation_gauge_components": translation_gauge_components,
    }


TARGETED_CARRY_DELTAS = (
    (0, 8),
    (0, 9),
    (0, 4),
    (8, 0),
    (9, 0),
    (4, 0),
)
DISTINGUISHED_I = (142004992589460, 142005019034340)
DISTINGUISHED_HORN_CELLS = frozenset(
    (1, sigma, target)
    for sigma in (0, 3, 8, 9, 12)
    for target in (5, 6, 9, 10)
)


def nonlocal_translation_audit(paths):
    """Enumerate every translated mask pair on one rooted physical path.

    The adjacent audit only sees identity translations.  This pass retains
    every ordered pair i<j, records the gap j-i, and then asks whether the
    formal F_13 identities 8+9=4 in either coordinate are realized by an
    actual composable triple of masks on a single path.
    """
    total_pairs = 0
    solvable_pairs = 0
    nonzero_pairs = 0
    delta_census = Counter()
    delta_gap_census = Counter()
    targeted_census = Counter()
    targeted_gap_census = Counter()
    targeted_path_sets = defaultdict(set)
    first_targeted = {}
    carry_triangles = Counter()
    first_carry_triangle = {}
    composition_failures = []

    carry_patterns = (
        ((0, 8), (0, 9), (0, 4)),
        ((0, 9), (0, 8), (0, 4)),
        ((8, 0), (9, 0), (4, 0)),
        ((9, 0), (8, 0), (4, 0)),
    )

    for path_number, (cell, root_index, pieces) in enumerate(paths):
        masks = tuple(endpoint_mask(piece[:2]) for piece in pieces)
        labelled_pairs = {}
        outgoing = defaultdict(lambda: defaultdict(list))
        incoming = defaultdict(lambda: defaultdict(list))
        for source_index in range(len(masks)):
            for target_index in range(source_index + 1, len(masks)):
                total_pairs += 1
                deltas = translation_deltas(
                    masks[source_index], masks[target_index]
                )
                if not deltas:
                    continue
                require(
                    len(deltas) == 1,
                    "a within-path translated pair acquired a nontrivial "
                    "mask stabilizer",
                )
                delta = deltas[0]
                gap = target_index - source_index
                solvable_pairs += 1
                delta_census[delta] += 1
                delta_gap_census[(delta, gap)] += 1
                labelled_pairs[(source_index, target_index)] = delta
                outgoing[source_index][delta].append(target_index)
                incoming[target_index][delta].append(source_index)
                if delta == (0, 0):
                    continue
                nonzero_pairs += 1
                if delta in TARGETED_CARRY_DELTAS:
                    targeted_census[delta] += 1
                    targeted_gap_census[(delta, gap)] += 1
                    targeted_path_sets[delta].add(path_number)
                    if delta not in first_targeted:
                        first_targeted[delta] = (
                            path_number,
                            cell,
                            root_index,
                            source_index,
                            target_index,
                            gap,
                            pieces[source_index][:2],
                            pieces[target_index][:2],
                            mask_description(masks[source_index]),
                            mask_description(masks[target_index]),
                        )

        # Translation labels must compose whenever all three physical pair
        # maps exist.  Search the two orientations of 8,9 in each coordinate.
        for first, second, direct in carry_patterns:
            for middle_index in range(len(masks)):
                for source_index in incoming[middle_index].get(first, ()):
                    for target_index in outgoing[middle_index].get(second, ()):
                        observed = labelled_pairs.get(
                            (source_index, target_index)
                        )
                        if observed != direct:
                            composition_failures.append(
                                (
                                    path_number,
                                    source_index,
                                    middle_index,
                                    target_index,
                                    first,
                                    second,
                                    direct,
                                    observed,
                                )
                            )
                            continue
                        pattern = (first, second, direct)
                        carry_triangles[pattern] += 1
                        if pattern not in first_carry_triangle:
                            first_carry_triangle[pattern] = (
                                path_number,
                                cell,
                                root_index,
                                source_index,
                                middle_index,
                                target_index,
                                pieces[source_index][:2],
                                pieces[middle_index][:2],
                                pieces[target_index][:2],
                            )

    require(
        not composition_failures,
        f"within-path translation labels failed composition: "
        f"{composition_failures[:3]}",
    )
    nonzero_delta_census = Counter({
        delta: count
        for delta, count in delta_census.items()
        if delta != (0, 0)
    })
    target_rows = tuple(
        (
            delta,
            targeted_census[delta],
            len(targeted_path_sets[delta]),
            tuple(
                sorted(
                    (gap, count)
                    for (row_delta, gap), count
                    in targeted_gap_census.items()
                    if row_delta == delta
                )
            ),
            first_targeted.get(delta),
        )
        for delta in TARGETED_CARRY_DELTAS
    )
    return {
        "total_pairs": total_pairs,
        "solvable_pairs": solvable_pairs,
        "nonzero_pairs": nonzero_pairs,
        "delta_census": delta_census,
        "nonzero_delta_census": nonzero_delta_census,
        "delta_gap_census": delta_gap_census,
        "target_rows": target_rows,
        "carry_triangles": carry_triangles,
        "first_carry_triangle": first_carry_triangle,
    }


def distinguished_i_translation_audit(paths):
    """Check the selected I section and its long-range Z^8 reference."""
    selected = []
    for cell, root_index, pieces in paths:
        if (
            cell in DISTINGUISHED_HORN_CELLS
            and len(pieces) > 2
            and pieces[2][:2] == DISTINGUISHED_I
        ):
            selected.append((cell, root_index, pieces))
    require(
        len(selected) == 20
        and {cell for cell, _root, _pieces in selected}
        == DISTINGUISHED_HORN_CELLS,
        "distinguished I horn-section census changed",
    )

    rows = []
    long_cells = []
    normalized_hits = []
    for cell, root_index, pieces in sorted(selected):
        base = endpoint_mask(pieces[2][:2])
        translated = tuple(
            (target_index, translation_deltas(
                base, endpoint_mask(pieces[target_index][:2])
            ))
            for target_index in range(3, len(pieces))
            if translation_deltas(
                base, endpoint_mask(pieces[target_index][:2])
            )
        )
        hits = tuple(
            (target_index, deltas)
            for target_index, deltas in translated
            if any(delta[1] == 1 for delta in deltas)
        )
        normalized_hits.extend(
            (cell, root_index, target_index, deltas)
            for target_index, deltas in hits
        )
        if len(pieces) > 68:
            deltas_68 = translation_deltas(
                base, endpoint_mask(pieces[68][:2])
            )
            require(
                deltas_68 == ((0, 8),),
                f"selected I lost its step-68 Z^8 reference at {cell}",
            )
            long_cells.append(cell)
        rows.append(
            (
                cell,
                root_index,
                len(pieces),
                tuple(
                    target_index
                    for target_index, deltas in translated
                    if deltas == ((0, 0),)
                ),
                tuple(
                    target_index
                    for target_index, deltas in translated
                    if deltas == ((0, 8),)
                ),
                hits,
            )
        )
    require(
        tuple(sorted(long_cells))
        == tuple(
            sorted(
                (1, sigma, target)
                for sigma in (0, 3, 12)
                for target in (5, 6, 9, 10)
            )
        )
        and not normalized_hits,
        "distinguished I long-carry/normalized-carry profile changed",
    )
    return {
        "rows": tuple(rows),
        "long_cells": tuple(sorted(long_cells)),
        "normalized_hits": tuple(normalized_hits),
    }


def allocation_shifted_i_audit():
    """Compare full endpoint masks on I+q(T/13), source and target."""
    rows = []
    source_masks = {}
    target_masks = {}
    for allocation in (0, 3, 7, 11):
        interval = (
            DISTINGUISHED_I[0] + allocation * T // P,
            DISTINGUISHED_I[1] + allocation * T // P,
        )
        target_interval = (
            interval[0] + SHIFT,
            interval[1] + SHIFT,
        )
        source_mask = endpoint_mask(interval)
        target_mask = endpoint_mask(target_interval)
        require(
            source_mask == target_mask,
            f"allocation q{allocation} source/target endpoint masks differ",
        )
        source_masks[allocation] = source_mask
        target_masks[allocation] = target_mask
        rows.append(
            (
                allocation,
                interval,
                target_interval,
                sum(source_mask),
                sum(target_mask),
                mask_description(source_mask),
            )
        )
    require(
        tuple(sum(source_masks[q]) for q in (0, 3, 7, 11))
        == (81, 90, 0, 81)
        and not translation_deltas(source_masks[0], source_masks[11])
        and not translation_deltas(target_masks[0], target_masks[11]),
        "allocation-shifted I endpoint boundary changed",
    )
    return {
        "rows": tuple(rows),
        "q0_q11_source_deltas":
            translation_deltas(source_masks[0], source_masks[11]),
        "q0_q11_target_deltas":
            translation_deltas(target_masks[0], target_masks[11]),
    }


def hinge_witness(paths):
    matches = []
    for cell, root_index, pieces in paths:
        if cell != (1, 0, 5):
            continue
        root_mask = endpoint_mask(pieces[0][:2])
        second_mask = endpoint_mask(pieces[2][:2])
        if sum(root_mask) == 0 and sum(second_mask) == 81:
            matches.append(
                (
                    cell,
                    root_index,
                    pieces[0][:2],
                    pieces[2][:2],
                    mask_description(root_mask),
                    mask_description(second_mask),
                )
            )
    require(len(matches) == 1, f"hinge witness multiplicity changed: {matches}")
    return matches[0]


def main():
    paths = path_geometry()
    endpoint_source_census, pair_census = validate_collar(paths)
    one_step = edge_audit(paths, 1)
    two_step = edge_audit(paths, 2)
    path_orbits = path_orbit_audit(paths)
    nonlocal_translations = nonlocal_translation_audit(paths)
    distinguished_i = distinguished_i_translation_audit(paths)
    allocation_shifted_i = allocation_shifted_i_audit()
    hinge = hinge_witness(paths)
    agl_13 = agl_fixed_power_count(13)
    agl_169 = agl_fixed_power_count(169)
    require(
        agl_13 == agl_169
        and agl_13 == (169, 28_561),
        f"AGL(2,13) p-primary exponent census changed: "
        f"{agl_13}, {agl_169}",
    )

    unique_masks = set()
    rectangle_types = Counter()
    for _cell, _root_index, pieces in paths:
        for piece in pieces:
            mask = endpoint_mask(piece[:2])
            unique_masks.add(mask)
            rectangle_types[mask_description(mask)] += 1
    require(
        all(cartesian_factor(mask)[2] for mask in unique_masks),
        "a full source endpoint mask is not Cartesian",
    )
    rectangle_size_types = Counter()
    vertex_weight_census = Counter()
    for description, count in rectangle_types.items():
        first, second, weight, _product = description
        rectangle_size_types[len(first), len(second), weight] += count
        vertex_weight_census[weight] += count

    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    require(
        not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
        "scratch audit contains an executable Python assert",
    )

    print("THM-2859 FULL SOURCE-ENDPOINT COBOUNDARY AUDIT")
    print(
        f"primary_sha256={PRIMARY_SHA256};"
        f"paths={len(paths)};path_vertices={path_orbits['vertices']};"
        f"physical_interval_masks_cached={len(mask_cache)};"
        f"unique_full_masks={len(unique_masks)}"
    )
    print(
        "independent_evaluator="
        "direct periodic inequalities, no materialized 169 present-set unions;"
        f"stored_source_census={tuple((name, tuple(sorted(rows.items()))) for name, rows in endpoint_source_census.items())};"
        f"stored_pair_census={tuple((name, tuple(sorted(rows.items(), key=repr))) for name, rows in pair_census.items())};"
        "verdict=exact match"
    )
    rectangle_digest = sha256(
        repr(tuple(sorted(rectangle_types.items(), key=repr))).encode("ascii")
    ).hexdigest()
    print(
        "full_mask_structure="
        "every mask is an exact Cartesian rectangle A_x*B_x;"
        f"rectangle_size_types="
        f"{tuple(sorted(rectangle_size_types.items()))};"
        f"vertex_weight_census={tuple(sorted(vertex_weight_census.items()))};"
        f"full_rectangle_census_sha256={rectangle_digest}"
    )
    for result in (one_step, two_step):
        label = f"{result['jump']}h"
        print(
            f"edges_{label}={result['total']};"
            f"hamming={tuple(sorted(result['hamming'].items()))};"
            f"augmentation={tuple(sorted(result['augmentation'].items()))};"
            f"relations={tuple(sorted(result['relations'].items()))};"
            f"axis_types={tuple(sorted(result['axes'].items()))}"
        )
        print(
            f"translation_gauge_{label}="
            f"solvable:{result['solvable']}/{result['total']};"
            f"delta_count_histogram="
            f"{tuple(sorted(result['translation_counts'].items()))};"
            f"first_nonzero_translation="
            f"{result['first_nonzero_translation']}"
        )
        print(
            f"cyclic_C169_translation_{label}="
            f"{tuple((order, tuple(sorted(result['cyclic_translation_counts'][order].items()))) for order in CYCLIC_ORDERS)};"
            "orders encode n=13a+q and n=a+13q"
        )
        print(
            f"address_derivative_{label}="
            "D=(T_delta-I)phi possible-direction-count histogram:"
            f"{tuple(sorted(result['derivative_direction_counts'].items()))};"
            f"first_equal_weight_nontranslation="
            f"{result['first_equal_weight_nontranslation']}"
        )
        mod2_valid = sum(
            count
            for direction_count, count
            in result["derivative_direction_mod2_counts"].items()
            if direction_count
        )
        print(
            f"address_derivative_mod2_{label}="
            f"possible_direction_count_histogram:"
            f"{tuple(sorted(result['derivative_direction_mod2_counts'].items()))};"
            f"valid_edges={mod2_valid}/{result['total']};"
            f"maximal_mod2_derivative_subforest_components="
            f"{path_orbits['vertices'] - mod2_valid};"
            f"first_nonzero_potential={result['first_mod2_potential']}"
        )
        print(
            f"cyclic_C169_derivative_{label}="
            f"{tuple((order, tuple(sorted(result['cyclic_derivative_counts'][order].items()))) for order in CYCLIC_ORDERS)}"
        )
        print(
            f"cyclic_C169_derivative_mod2_{label}="
            f"{tuple((order, tuple(sorted(result['cyclic_derivative_mod2_counts'][order].items()))) for order in CYCLIC_ORDERS)}"
        )
        print(
            f"minimal_weight_obstruction_{label}="
            f"{result['first_small_obstruction']}"
        )
        compact_types = Counter()
        for row, count in result["pair_types"].items():
            (
                source_description,
                target_description,
                distance,
                augmentation,
                relation,
                axis,
                delta_count,
                directions,
                cyclic_delta_counts,
                cyclic_derivative_counts,
                mod2_direction_count,
                cyclic_mod2_derivative_counts,
            ) = row
            compact_types[
                source_description[2],
                target_description[2],
                distance,
                augmentation,
                relation,
                axis,
                delta_count,
                len(directions),
                cyclic_delta_counts,
                cyclic_derivative_counts,
                mod2_direction_count,
                cyclic_mod2_derivative_counts,
            ] += count
        compact_digest = sha256(
            repr(tuple(sorted(compact_types.items(), key=repr))).encode(
                "ascii"
            )
        ).hexdigest()
        full_pair_digest = sha256(
            repr(
                tuple(sorted(result["pair_types"].items(), key=repr))
            ).encode("ascii")
        ).hexdigest()
        print(
            f"pair_type_digests_{label}="
            f"compact:{compact_digest},full:{full_pair_digest};"
            f"compact_type_count={len(compact_types)};"
            f"full_type_count={len(result['pair_types'])}"
        )
    print(
        "path_translation_orbits="
        f"all_one_orbit:{path_orbits['all_one_orbit']}/{len(paths)};"
        f"distinct_orbits_per_path="
        f"{tuple(sorted(path_orbits['orbit_count'].items()))};"
        f"root_solvable_prefix_edges="
        f"{tuple(sorted(path_orbits['root_orbit_prefix'].items()))};"
        f"solvable_h_edges="
        f"{path_orbits['translation_solvable_h_edges']};"
        f"maximal_gauge_subforest_components="
        f"{path_orbits['translation_gauge_components']}"
    )
    nonlocal_gap_digest = sha256(
        repr(
            tuple(
                sorted(
                    nonlocal_translations["delta_gap_census"].items(),
                    key=repr,
                )
            )
        ).encode("ascii")
    ).hexdigest()
    print(
        "within_path_translation_pairs="
        f"all_ordered_i_lt_j:{nonlocal_translations['total_pairs']};"
        f"solvable:{nonlocal_translations['solvable_pairs']};"
        f"nonzero:{nonlocal_translations['nonzero_pairs']};"
        f"nonzero_delta_census="
        f"{tuple(sorted(nonlocal_translations['nonzero_delta_census'].items()))};"
        f"delta_gap_census_sha256={nonlocal_gap_digest}"
    )
    compact_target_rows = tuple(
        (
            delta,
            count,
            path_count,
            (gaps[0][0], gaps[-1][0]) if gaps else (),
            sha256(repr(gaps).encode("ascii")).hexdigest(),
            witness,
        )
        for delta, count, path_count, gaps, witness
        in nonlocal_translations["target_rows"]
    )
    print(
        "targeted_8_9_4_translation_rows="
        f"{compact_target_rows}"
    )
    print(
        "composable_8_9_vs_4_triangles="
        f"{tuple(sorted(nonlocal_translations['carry_triangles'].items()))};"
        f"first_witnesses="
        f"{tuple(sorted(nonlocal_translations['first_carry_triangle'].items()))};"
        "scope=physical endpoint-mask translations on one rooted collar path"
    )
    distinguished_profile_census = Counter(
        (length, identity_steps, z8_steps, normalized_hits)
        for (
            _cell,
            _root_index,
            length,
            identity_steps,
            z8_steps,
            normalized_hits,
        ) in distinguished_i["rows"]
    )
    distinguished_profile_digest = sha256(
        repr(distinguished_i["rows"]).encode("ascii")
    ).hexdigest()
    print(
        "distinguished_I_nonlocal_carry="
        f"I={DISTINGUISHED_I};"
        f"long_step68_Z8_cells={distinguished_i['long_cells']};"
        f"normalized_Z_hits={distinguished_i['normalized_hits']};"
        f"profile_census="
        f"{tuple(sorted(distinguished_profile_census.items()))};"
        f"complete_profile_sha256={distinguished_profile_digest};"
        "verdict=generator-valued split-plane Z^8 source-mask reference on "
        "12 long horn paths, but no selected-I Z=(0,1) reference and no "
        "physical C13 action"
    )
    print(
        "allocation_shifted_I_endpoint_boundary="
        f"rows={allocation_shifted_i['rows']};"
        f"q0_to_q11_source_deltas="
        f"{allocation_shifted_i['q0_q11_source_deltas']};"
        f"q0_to_q11_target_deltas="
        f"{allocation_shifted_i['q0_q11_target_deltas']};"
        "verdict=masses q0,q3,q7,q11 are 81,90,0,81 on both source "
        "and target; equal mass does not make q11 a q0 translate"
    )
    print(
        "hinge_clock1_s0_t5="
        f"{hinge};"
        "verdict=empty source mask -> 81-point M2 mask, impossible under "
        "every permutation of the 169 addresses"
    )
    print(
        "affine_group_boundary="
        f"AGL2_maps_with_f13_identity={agl_13};"
        f"AGL2_maps_with_f169_identity={agl_169};"
        "there is no affine-plane element of order 169, while cyclic "
        "translation by one on either digit encoding has order 169; "
        "all affine relabelings preserve mask mass, so the nonconstant "
        "edge obstructions survive AGL(2,13), Aff(C169), coordinate "
        "reversal, and both digit orders"
    )
    print(
        "COCHAIN_LEDGER="
        "ambient additive path cochain is exact with potential "
        "Phi(vertex)=its full 169-bit source endpoint mask; "
        "translation-gauge exactness requires edgewise orbit equivalence; "
        "over Z or F13, address-translation derivative exactness requires "
        "zero line sums and every nonzero augmentation obstructs; over F2 "
        "the ten-point slices acquire an exact one-axis potential, but "
        "odd nine/81-point defects remain and no new mask conjugacy follows"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
