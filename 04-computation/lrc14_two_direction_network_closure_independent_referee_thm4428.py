#!/usr/bin/env python3
"""Clean-room exact referee for the THM-4428 two-direction closure.

The all-height part is elementary once THM-4386 Lemma 2.1 is inherited:
in a multi-direction support each primitive live direction has l1 norm at
least 16 and hence maximum coefficient at least 7.  This verifier audits all
remaining algebraic implications, exhausts the complete c<55 head, and uses
three mutually different finite representations:

* a simultaneous owner-interval sweep;
* direct relation-lattice boxes solved in the xy and yz coordinates; and
* literal shifted-sheet contact intervals for the three THM-4414 networks.

No repository module or incoming verifier is imported.  All proof gates use
``require`` and therefore remain active under optimized Python.
"""

from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations, product
from math import gcd
import sys


TARGET = Q(6, 77)
TAIL_START = 55
CHECKS = 0
LITERAL_CHECKS = 0


class RefereeFailure(RuntimeError):
    pass


def require(condition, payload):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RefereeFailure(payload)


def literal_require(condition, payload):
    global LITERAL_CHECKS
    LITERAL_CHECKS += 1
    if not condition:
        raise RefereeFailure(payload)


def strict_floor(numerator, denominator):
    """Largest integer strictly less than numerator/denominator."""
    return (numerator - 1) // denominator


def cross(left, right):
    return (
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    )


def primitive_direction(carrier):
    divisor = gcd(gcd(abs(carrier[0]), abs(carrier[1])), abs(carrier[2]))
    require(divisor > 0, ("nonzero carrier", carrier))
    ray = tuple(value // divisor for value in carrier)
    first = next(value for value in ray if value)
    return ray if first > 0 else tuple(-value for value in ray)


def carrier_bounds(w):
    return tuple(strict_floor(3 * (sum(w) - w[index]), 14) for index in range(3))


def live_gate(w, carrier, bounds):
    return (
        all(abs(carrier[index]) <= bounds[index] for index in range(3))
        and all(value % 3 for value in carrier)
        and sum(w[index] * carrier[index] for index in range(3)) == 0
    )


def lattice_box_xy(w):
    """Solve the relation for z after independently ranging x and y."""
    bx, by, bz = carrier_bounds(w)
    result = set()
    for x in range(-bx, bx + 1):
        for y in range(-by, by + 1):
            numerator = -(w[0] * x + w[1] * y)
            if numerator % w[2]:
                continue
            carrier = (x, y, numerator // w[2])
            if live_gate(w, carrier, (bx, by, bz)):
                result.add(carrier)
    return result


def lattice_box_yz(w):
    """Independent coordinate order: solve the relation for x from y,z."""
    bx, by, bz = carrier_bounds(w)
    result = set()
    for y in range(-by, by + 1):
        for z in range(-bz, bz + 1):
            numerator = -(w[1] * y + w[2] * z)
            if numerator % w[0]:
                continue
            carrier = (numerator // w[0], y, z)
            if live_gate(w, carrier, (bx, by, bz)):
                result.add(carrier)
    return result


@lru_cache(maxsize=None)
def danger_intervals(speed):
    """Open owner intervals, with endpoints numerator/(14*speed)."""
    return tuple(
        (max(0, 14 * owner - 3), min(14 * speed, 14 * owner + 3), owner)
        for owner in range(speed + 1)
    )


def owner_interval_sweep(w):
    """Enumerate physical components before translating them to carriers."""
    rows = tuple(danger_intervals(speed) for speed in w)
    cursors = [0, 0, 0]
    result = set()
    while all(cursors[index] < len(rows[index]) for index in range(3)):
        current = tuple(rows[index][cursors[index]] for index in range(3))

        lower_num, lower_den = current[0][0], w[0]
        upper_num, upper_den = current[0][1], w[0]
        for index in (1, 2):
            if current[index][0] * lower_den > lower_num * w[index]:
                lower_num, lower_den = current[index][0], w[index]
            if current[index][1] * upper_den < upper_num * w[index]:
                upper_num, upper_den = current[index][1], w[index]

        if lower_num * upper_den < upper_num * lower_den:
            owners = tuple(current[index][2] for index in range(3))
            labels = tuple((-w[index] * owners[index]) % 3 for index in range(3))
            if set(labels) == {0, 1, 2}:
                carrier = cross(w, owners)
                require(carrier not in result, ("duplicate sweep carrier", w, owners, carrier))
                result.add(carrier)

        earliest_num, earliest_den = current[0][1], w[0]
        for index in (1, 2):
            if current[index][1] * earliest_den < earliest_num * w[index]:
                earliest_num, earliest_den = current[index][1], w[index]
        for index in range(3):
            if current[index][1] * earliest_den == earliest_num * w[index]:
                cursors[index] += 1
    return result


def raw_projections(w, carriers):
    """THM-4414 E_i, equivalently THM-4422 S_i, in coordinate order."""
    q = Q(3, 7 * w[2])
    totals = [Q(0), Q(0), Q(0)]
    for carrier in carriers:
        for omitted in range(3):
            pair = [index for index in range(3) if index != omitted]
            numerator = 3 * (w[pair[0]] + w[pair[1]]) - 14 * abs(carrier[omitted])
            require(numerator > 0, ("positive raw roof", w, carrier, omitted, numerator))
            totals[omitted] += min(
                q, Q(numerator, 14 * w[pair[0]] * w[pair[1]])
            )
    return tuple(totals)


def circular_piece(left, right, period):
    if left < 0:
        return ((0, right), (period + left, period))
    if right > period:
        return ((left, period), (0, right - period))
    return ((left, right),)


def shifted_sheet_intervals(w):
    """Literal six-sheet intervals over one common integer denominator."""
    denominator = 42 * w[0] * w[1] * w[2]
    sheets = {}
    for index, speed in enumerate(w):
        scale = denominator // (42 * speed)
        for label in range(3):
            pieces = []
            for owner in range(speed):
                residue = (3 * owner - speed * label) % (3 * speed)
                center = 14 * residue * scale
                pieces.extend(circular_piece(center - 3 * scale, center + 3 * scale, denominator))
            sheets[index, label] = tuple(sorted(pieces))
    return denominator, sheets


def interval_intersections(left, right):
    first = second = 0
    result = []
    while first < len(left) and second < len(right):
        start = max(left[first][0], right[second][0])
        stop = min(left[first][1], right[second][1])
        if start < stop:
            result.append((start, stop))
        left_stop, right_stop = left[first][1], right[second][1]
        first += left_stop <= right_stop
        second += right_stop <= left_stop
    return tuple(result)


def literal_sheet_networks(w):
    """Reconstruct all three edgewise-minimum networks from sheet contacts."""
    denominator, sheets = shifted_sheet_intervals(w)
    network_totals = []
    physical_totals = []
    for omitted in range(3):
        pair_indices = [index for index in range(3) if index != omitted]
        network = physical = 0
        for labels in permutations(range(3)):
            pair_components = interval_intersections(
                sheets[pair_indices[0], labels[pair_indices[0]]],
                sheets[pair_indices[1], labels[pair_indices[1]]],
            )
            opposite = sheets[omitted, labels[omitted]]
            first = second = 0
            while first < len(pair_components) and second < len(opposite):
                start = max(pair_components[first][0], opposite[second][0])
                stop = min(pair_components[first][1], opposite[second][1])
                if start < stop:
                    network += min(
                        pair_components[first][1] - pair_components[first][0],
                        opposite[second][1] - opposite[second][0],
                    )
                    physical += stop - start
                pair_stop, opposite_stop = pair_components[first][1], opposite[second][1]
                first += pair_stop <= opposite_stop
                second += opposite_stop <= pair_stop
        network_totals.append(Q(network, denominator))
        physical_totals.append(Q(physical, denominator))
    literal_require(len(set(physical_totals)) == 1, ("literal physical grouping", w, physical_totals))
    return tuple(network_totals), physical_totals[0]


def ray_partition(carriers):
    result = {}
    for carrier in carriers:
        direction = primitive_direction(carrier)
        result.setdefault(direction, set()).add(carrier)
    return result


def direction_cutoff(w, direction):
    return min(
        Q(3 * (sum(w) - w[index]), 14 * abs(direction[index]))
        for index in range(3)
    )


def exact_ray(w, direction):
    cutoff = direction_cutoff(w, direction)
    last = strict_floor(cutoff.numerator, cutoff.denominator)
    points = {
        tuple(sign * multiplier * coordinate for coordinate in direction)
        for multiplier in range(1, last + 1)
        if multiplier % 3
        for sign in (-1, 1)
    }
    return cutoff, last, points


def analytic_gate_audit():
    """Finite residue classifications and exact constant arithmetic."""
    # Full-support F_3 kernel vectors form one unoriented ray for every unit w.
    residue_cases = 0
    for wbar in product((1, 2), repeat=3):
        kernel = tuple(
            vector for vector in product((1, 2), repeat=3)
            if sum(wbar[index] * vector[index] for index in range(3)) % 3 == 0
        )
        require(len(kernel) == 2, ("F3 full-support kernel", wbar, kernel))
        require(kernel[1] == tuple((-value) % 3 for value in kernel[0]),
                ("F3 opposite rays", wbar, kernel))
        residue_cases += 1

    # If a ternary-unit coefficient vector has M<7 and even l1>14, it cannot
    # exist: magnitude 6 is forbidden, and three magnitudes at most 5 have
    # even sum at most 14.  Exhaust the complete signed small universe.
    small_vectors = 0
    for vector in product(tuple(value for value in range(-6, 7) if value and value % 3), repeat=3):
        norm = sum(abs(value) for value in vector)
        if norm % 2 == 0:
            small_vectors += 1
            require(norm <= 14, ("M<7 even-norm classification", vector, norm))

    # For K=3t+j the exact signed nonmultiples-of-three count obeys
    # 3N <= 4K+4.  Checking one representative for each residue is complete.
    for last in range(3):
        signed_count = 2 * (last - last // 3)
        require(3 * signed_count <= 4 * last + 4,
                ("three-block residue bound", last, signed_count))

    tail_at_55 = Q(12, 343) + Q(16, 7 * TAIL_START)
    require(tail_at_55 == Q(1444, 18865), ("tail arithmetic", tail_at_55))
    require(TARGET == Q(1470, 18865), ("target common denominator", TARGET))
    require(tail_at_55 < TARGET, ("strict tail threshold", tail_at_55, TARGET))
    return residue_cases, small_vectors, tail_at_55


def audit_two_direction_row(w, carriers, finite_head):
    partition = ray_partition(carriers)
    require(len(partition) == 2, ("exactly two primitive unoriented directions", w, partition))
    directions = tuple(sorted(partition))
    maxima = []
    counts = []
    for direction in directions:
        require(all(value % 3 for value in direction), ("direction ternary units", w, direction))
        norm = sum(abs(value) for value in direction)
        require(norm % 2 == 0, ("odd-speed parity", w, direction, norm))
        # THM-4386 Lemma 2.1 supplies norm>14 in a multi-direction support.
        require(norm >= 16, ("inherited short-relation obstruction", w, direction, norm))
        maximum = max(abs(value) for value in direction)
        require(maximum >= 7, ("maximum coefficient at least seven", w, direction, maximum))
        maxima.append(maximum)

        cutoff, last, predicted = exact_ray(w, direction)
        require(partition[direction] == predicted,
                ("complete residue-deleted multiplier list", w, direction,
                 partition[direction] ^ predicted))
        exact_count = 2 * (last - last // 3)
        require(exact_count == len(predicted), ("exact ray count", w, direction, exact_count))
        require(cutoff < Q(3 * w[2], 7 * maximum),
                ("narrow-coordinate cutoff", w, direction, cutoff, maximum))
        require(Q(exact_count) < Q(4, 3) * cutoff + Q(4, 3),
                ("strict residue-deleted count bound", w, direction, cutoff, exact_count))
        counts.append(exact_count)

    u, v = directions
    determinant = cross(u, v)
    require(determinant != (0, 0, 0), ("distinct directions independent", w, directions))
    require(all(value % 3 == 0 for value in determinant),
            ("owner residue determinant", w, directions, determinant))
    require(determinant[2] % w[2] == 0, ("integral determinant multiplier", w, determinant))
    multiplier = determinant[2] // w[2]
    require(multiplier and multiplier % 3 == 0,
            ("nonzero multiplier divisible by three", w, directions, multiplier))
    require(determinant == tuple(multiplier * value for value in w),
            ("oriented cross product", w, directions, determinant, multiplier))
    require(2 * maxima[0] * maxima[1] >= 3 * w[2],
            ("height product bound", w, maxima))
    require((maxima[0] - 7) * (maxima[1] - 7) >= 0,
            ("reciprocal factorization hypotheses", w, maxima))
    reciprocal = Q(1, maxima[0]) + Q(1, maxima[1])
    reciprocal_bound = Q(1, 7) + Q(14, 3 * w[2])
    require(reciprocal <= reciprocal_bound,
            ("reciprocal determinant bound", w, maxima, reciprocal, reciprocal_bound))

    require(sum(counts) == len(carriers), ("two disjoint ray counts", w, counts, len(carriers)))
    values = raw_projections(w, carriers)
    q = Q(3, 7 * w[2])
    intermediate = Q(12, 49) * reciprocal + Q(8, 7 * w[2])
    envelope = Q(12, 343) + Q(16, 7 * w[2])
    for index, value in enumerate(values):
        require(value <= q * len(carriers), ("common cap pays projection", w, index, value))
        require(value < intermediate, ("two-ray count envelope", w, index, value, intermediate))
        require(value < envelope, ("determinant reciprocal envelope", w, index, value, envelope))
        if finite_head or w[2] >= TAIL_START:
            require(value < TARGET, ("strict network target", w, index, value))
    return values, directions, multiplier, tuple(maxima), tuple(counts)


def eligible_head():
    values = tuple(value for value in range(1, TAIL_START, 2) if value % 3)
    return tuple(sorted((
        w for w in combinations(values, 3)
        if gcd(gcd(w[0], w[1]), w[2]) == 1
    ), key=lambda w: (w[2], w[0], w[1])))


def main():
    sys.stdout.reconfigure(newline="\n")
    residue_cases, small_vectors, tail_value = analytic_gate_audit()
    rows = eligible_head()
    require(len(rows) == 814, ("complete finite-head universe", len(rows)))

    selected = 0
    literal_rows = 0
    first_two_direction = None
    leader = (Q(0), None, None, None)
    all_rows_digest = sha256()
    selected_digest = sha256()
    three_direction_control = None

    for w in rows:
        sweep = owner_interval_sweep(w)
        xy = lattice_box_xy(w)
        yz = lattice_box_yz(w)
        require(sweep == xy == yz, ("three-way carrier equality", w, sweep ^ xy, xy ^ yz))
        rays = ray_partition(sweep)
        values = raw_projections(w, sweep)
        all_rows_digest.update(repr((w, tuple(sorted(sweep)), tuple(sorted(rays)), values)).encode("ascii"))

        if len(rays) == 3 and w == (19, 23, 29):
            control_directions = tuple(sorted(rays))
            require(cross(control_directions[0], control_directions[1]) != (0, 0, 0),
                    ("three-direction control has rank at least two", control_directions))
            require(all(sum(w[index] * direction[index] for index in range(3)) == 0
                        for direction in control_directions),
                    ("three-direction control lies in rank-two kernel", w, control_directions))
            three_direction_control = (w, 2, 3, control_directions, values)
        if len(rays) != 2:
            continue

        selected += 1
        if first_two_direction is None:
            first_two_direction = w
        values, directions, multiplier, maxima, counts = audit_two_direction_row(w, sweep, True)
        literal_values, physical = literal_sheet_networks(w)
        literal_rows += 1
        require(values == literal_values,
                ("raw projection equals literal sheet network", w, values, literal_values))
        selected_digest.update(
            repr((w, tuple(sorted(sweep)), values, directions, multiplier, maxima, counts, physical)).encode("ascii")
        )
        candidate = (max(values), w, values, directions)
        if candidate > leader:
            leader = candidate

    require(selected == 192, ("exactly-two-direction finite-head count", selected))
    require(first_two_direction == (17, 23, 25),
            ("first two-direction row", first_two_direction))
    require(leader == (
        Q(114, 2233),
        (11, 23, 29),
        (Q(2, 203), Q(114, 2233), Q(1634, 51359)),
        ((8, 5, -7), (11, -4, -1)),
    ), ("finite-head global projection leader", leader))
    require(three_direction_control is not None, "three-direction control present")
    require(three_direction_control[1:3] == (2, 3),
            ("direction count is not rank", three_direction_control))

    tail_controls = []
    for w in ((1, 11, 55), (1, 5, 101), (5, 49, 251)):
        sweep = owner_interval_sweep(w)
        xy = lattice_box_xy(w)
        yz = lattice_box_yz(w)
        require(sweep == xy == yz, ("tail three-way carrier equality", w))
        values, directions, multiplier, maxima, counts = audit_two_direction_row(w, sweep, False)
        literal_values, physical = literal_sheet_networks(w)
        literal_rows += 1
        require(values == literal_values,
                ("tail raw/literal projection equality", w, values, literal_values))
        tail_controls.append((w, len(sweep), directions, multiplier, values, physical))

    print("THM-4428 TWO-DIRECTION NETWORK CLOSURE -- CLEAN-ROOM INDEPENDENT REFEREE")
    print("status=PASS; theorem_scope=complete support with exactly two primitive unoriented directions")
    print("notation=E_i(THM-4414)=S_i(THM-4422)=raw coordinate-i network projection")
    print("analytic_dependency=THM-4386 Lemma 2.1 short full-support ternary-unit relation rigidity")
    print("analytic_chain=l1>=16; M_u,M_v>=7; 3c<=2M_uM_v; E_i<12/343+16/(7c)")
    print("tail_start=55; tail_envelope_at_55=%s; target=%s" % (tail_value, TARGET))
    print("residue_kernel_cases=%d; signed_small_vector_cases=%d" % (residue_cases, small_vectors))
    print("finite_head=max_speed<55; universe_rows=%d; three_representation_rows=%d" % (len(rows), len(rows)))
    print("finite_head_exactly_two_direction_rows=%d; literal_network_rows=%d" % (selected, literal_rows - 3))
    print("finite_head_first_two_direction=%s" % (first_two_direction,))
    print("finite_head_leader=value:%s,triple:%s,E:%s,directions:%s" % leader)
    print("direction_not_rank_control=triple:%s,rank:%d,direction_count:%d,directions:%s,E:%s" %
          three_direction_control)
    for control in tail_controls:
        print("tail_control=triple:%s,N:%d,directions:%s,k:%d,E:%s,physical:%s" % control)
    print("all_head_rows_sha256=" + all_rows_digest.hexdigest())
    print("selected_rows_sha256=" + selected_digest.hexdigest())
    print("checks=%d; literal_checks=%d; literal_network_rows_total=%d" %
          (CHECKS, LITERAL_CHECKS, literal_rows))
    print("verdict=PASS; three_or_more_directions=OPEN; entry=OPEN; synchronization=OPEN; LRC14=OPEN")


if __name__ == "__main__":
    main()
