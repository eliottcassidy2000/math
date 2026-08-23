#!/usr/bin/env python3
"""Independent hostile audit of the scratch F12/F13 reset-observer bridge.

The reset graph is rebuilt directly from the two multiplicity boxes and
L1-distance descent.  Response rows are rebuilt from the older THM-3238
product-Gamma coefficient source with a fresh partition/upset constructor;
the bridge probe itself and the THM-3286 executable are not imported.
"""

from collections import Counter
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE_PATH = ROOT / "04-computation/gmc_complete_physical_bank_unique_reset_thm3238.py"


def load_base():
    spec = spec_from_file_location("thm3238_direct_audit_source", BASE_PATH)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


base = load_base()


@lru_cache(maxsize=None)
def partitions(total, maximum=None):
    if total == 0:
        return ((),)
    bound = total if maximum is None else min(total, maximum)
    return tuple(
        (head,) + tail
        for head in range(bound, 0, -1)
        for tail in partitions(total - head, head)
    )


@lru_cache(maxsize=None)
def coarsens(fine, coarse):
    """Return whether coarse can be obtained by merging parts of fine."""
    if sum(fine) != sum(coarse) or len(fine) < len(coarse):
        return False
    pieces = tuple(sorted(fine, reverse=True))
    bins = tuple(sorted(coarse, reverse=True))

    @lru_cache(maxsize=None)
    def place(index, capacities):
        if index == len(pieces):
            return not any(capacities)
        piece = pieces[index]
        tried = set()
        for slot, capacity in enumerate(capacities):
            if capacity < piece or capacity in tried:
                continue
            tried.add(capacity)
            changed = list(capacities)
            changed[slot] -= piece
            changed.sort(reverse=True)
            if place(index + 1, tuple(changed)):
                return True
        return False

    return place(0, bins)


UPSETS = tuple(
    frozenset(
        shape for shape in partitions(degree)
        if any(coarsens(generator, shape) for generator in generators)
    )
    for degree, _expected, generators, _factor in base.CERTIFICATE
)
if tuple(map(len, UPSETS)) != tuple(row[1] for row in base.CERTIFICATE):
    raise RuntimeError("independent upset reconstruction failed")


def row_values(a_value, b_value, state):
    vectors = base.coefficient_vectors(1, base.BANK, a_value, b_value, state)
    return tuple(
        sum(vectors[entry[0]][shape] for shape in upset)
        for entry, upset in zip(base.CERTIFICATE, UPSETS)
    )


FACES = {
    "F12": {
        "parameters": (1, 2),
        "capacities": (4, 3, 2, 1, 1, 0, 0, 0),
        "reset": (1, 2, 1, 1, 1, 0, 0, 0),
    },
    "F13": {
        "parameters": (1, 3),
        "capacities": (4, 3, 2, 2, 2, 1, 1, 1),
        "reset": (1, 0, 2, 1, 1, 1, 1, 1),
    },
}
N0 = (3, 4, 5)
N1 = (1, 3, 4, 5)


def as_counts(state):
    counts = Counter(state)
    return tuple(counts[root] for root in range(1, 9))


def as_state(counts):
    return tuple(
        root for root, count in enumerate(counts, 1) for _ in range(count)
    )


def l1(left, right):
    return sum(abs(a - b) for a, b in zip(left, right))


def neighbours(face, state):
    data = FACES[face]
    counts = as_counts(state)
    old_distance = l1(counts, data["reset"])
    answer = []
    for root, capacity in enumerate(data["capacities"]):
        for step in (-1, 1):
            changed = list(counts)
            changed[root] += step
            if not 0 <= changed[root] <= capacity or not any(changed):
                continue
            if l1(changed, data["reset"]) == old_distance - 1:
                answer.append(as_state(changed))
    # Preserve the independently chosen root order.  The candidate digest
    # serializes this order, even though the underlying graph is unordered.
    return tuple(answer)


def face_profile(face):
    a_value, b_value = FACES[face]["parameters"]
    poles = base.reduced_poles(1, base.BANK, a_value, b_value)[0]
    reset = tuple(sorted(base.residual_roots(
        1, base.dominant_row(1), a_value, b_value
    )))
    return tuple(sorted(Counter(poles).items())), reset


def local_record(face, state):
    a_value, b_value = FACES[face]["parameters"]
    targets = neighbours(face, state)
    source = row_values(a_value, b_value, state)
    witness_fibres = []
    for index in range(22):
        witnesses = []
        for target in targets:
            delta = row_values(a_value, b_value, target)[index] - source[index]
            if delta > 0:
                witnesses.append((target, delta))
        if witnesses:
            witness_fibres.append((index + 1, tuple(witnesses)))
    return {
        "face": face,
        "state": state,
        "neighbours": targets,
        "availability": tuple(row for row, _ in witness_fibres),
        "witness_fibres": tuple(witness_fibres),
    }


def digest(value):
    return sha256(repr(value).encode("ascii")).hexdigest()


def graph_audit():
    capacities = FACES["F12"]["capacities"]
    states = tuple(sorted((
        as_state(counts)
        for counts in product(*(range(capacity + 1) for capacity in capacities))
        if any(counts)
    ), key=lambda state: (len(state), state)))
    graph = []
    sinks = []
    histogram = Counter()
    for state in states:
        right = set(neighbours("F13", state))
        common = tuple(
            target for target in neighbours("F12", state) if target in right
        )
        histogram[len(common)] += 1
        if common:
            graph.extend((state, target) for target in common)
        else:
            sinks.append(state)
    return states, tuple(graph), tuple(sinks), tuple(sorted(histogram.items()))


def complex_add(left, right):
    return left[0] + right[0], left[1] + right[1]


def complex_mul(left, right):
    return (
        left[0] * right[0] - left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


def edge_integral_by_coefficients(start, end):
    """Integrate (a+bt)^2 conjugate(b) dt coefficient by coefficient."""
    direction = (end[0] - start[0], end[1] - start[1])
    conjugate_direction = (direction[0], -direction[1])
    coefficients = (
        complex_mul(start, start),
        (2 * complex_mul(start, direction)[0],
         2 * complex_mul(start, direction)[1]),
        complex_mul(direction, direction),
    )
    integral_square = (Fraction(0), Fraction(0))
    for degree, coefficient in enumerate(coefficients):
        integral_square = complex_add(
            integral_square,
            (coefficient[0] / (degree + 1), coefficient[1] / (degree + 1)),
        )
    return complex_mul(conjugate_direction, integral_square)


def main():
    expected_profiles = {
        "F12": (((1, 4), (2, 3), (3, 2), (4, 1), (5, 1)),
                (1, 2, 2, 3, 4, 5)),
        "F13": (((1, 4), (2, 3), (3, 2), (4, 2), (5, 2),
                  (6, 1), (7, 1), (8, 1)),
                (1, 3, 3, 4, 5, 6, 7, 8)),
    }
    if {face: face_profile(face) for face in FACES} != expected_profiles:
        raise RuntimeError("face profile mismatch")

    records = {
        (face, state): local_record(face, state)
        for face in ("F12", "F13") for state in (N0, N1)
    }
    deltas = {}
    for face in ("F12", "F13"):
        a_value, b_value = FACES[face]["parameters"]
        before = row_values(a_value, b_value, N0)
        after = row_values(a_value, b_value, N1)
        deltas[face] = tuple(b - a for a, b in zip(before, after))
    positive = {
        face: tuple(i + 1 for i, value in enumerate(delta) if value > 0)
        for face, delta in deltas.items()
    }
    equal_positive = tuple(
        (i + 1, j + 1, left)
        for i, left in enumerate(deltas["F12"]) if left > 0
        for j, right in enumerate(deltas["F13"]) if right > 0 and left == right
    )

    states, graph, sinks, histogram = graph_audit()

    vertices = (
        (Fraction(-1), Fraction(-1)),
        (Fraction(2), Fraction(-1)),
        (Fraction(-1), Fraction(2)),
    )
    increments = tuple(edge_integral_by_coefficients(
        vertices[index], vertices[(index + 1) % 3]
    ) for index in range(3))
    offsets = (Fraction(0), Fraction(0)), increments[0], complex_add(increments[0], increments[1])
    holonomy = complex_add(offsets[-1], increments[-1])

    print("profiles=" + repr(expected_profiles))
    print("neighbours=" + repr(tuple(
        (face, state, neighbours(face, state))
        for face in ("F12", "F13") for state in (N0, N1)
    )))
    print("positive_rows=" + repr(positive))
    print("availability_sizes=" + repr(tuple(
        (key, len(records[key]["availability"])) for key in sorted(records)
    )))
    print("equal_positive_delta_relation=" + repr(equal_positive))
    print("shared_row_intersections=" + repr(tuple(
        (state, tuple(sorted(
            set(records[("F12", state)]["availability"])
            & set(records[("F13", state)]["availability"])
        ))) for state in (N0, N1)
    )))
    print("common_graph=" + repr((
        len(states), len(graph), histogram, sinks, digest(graph)
    )))
    print("shared_record_digest=" + digest(tuple(
        (key, records[key]) for key in sorted(records)
    )))
    print("boundary=" + repr((increments, offsets, holonomy)))
    print("boundary_digest=" + digest((vertices, increments, offsets, holonomy)))
    print("all_independent_checks=PASS")


if __name__ == "__main__":
    main()
