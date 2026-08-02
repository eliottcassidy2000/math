#!/usr/bin/env python3
"""Exact cover census for THM-3119.

The proof of factorial conjugacy and stochastic same-label deletion is
symbolic.  This companion exhausts every integer-partition cover through
source degree nine, verifies the conjugacy coefficient by coefficient, and
checks the induced refinement transport by exact integer max flow.  It also
records the first raw block-deletion hostiles on the sign representation.
"""

from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from math import factorial, prod


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


@lru_cache(maxsize=None)
def partitions(total, largest=None):
    if total == 0:
        return ((),)
    if largest is None or largest > total:
        largest = total
    return tuple(
        (first,) + tail
        for first in range(largest, 0, -1)
        for tail in partitions(total - first, first)
    )


def weight(shape):
    return prod(factorial(part) for part in shape)


def lower(shape, part):
    answer = list(shape)
    answer.remove(part)
    if part > 1:
        answer.append(part - 1)
    return tuple(sorted(answer, reverse=True))


def block_deletion(shape):
    counts = Counter(shape)
    return {
        lower(shape, part): multiplicity
        for part, multiplicity in counts.items()
    }


def labelled_deletion(shape):
    counts = Counter(shape)
    return {
        lower(shape, part): part * multiplicity
        for part, multiplicity in counts.items()
    }


def covers(shape):
    """All types obtained by merging exactly two blocks."""

    answer = set()
    for first in range(len(shape)):
        for second in range(first + 1, len(shape)):
            merged = [shape[index] for index in range(len(shape))
                      if index not in (first, second)]
            merged.append(shape[first] + shape[second])
            answer.add(tuple(sorted(merged, reverse=True)))
    return tuple(sorted(answer, reverse=True))


@lru_cache(maxsize=None)
def coarsens(fine, coarse):
    if sum(fine) != sum(coarse) or len(coarse) > len(fine):
        return False
    if not coarse:
        return not fine
    target = coarse[0]
    remaining_targets = coarse[1:]
    remainders = set()

    def choose(start, subtotal, selected):
        if subtotal == target:
            remainders.add(tuple(sorted(
                (fine[index] for index in range(len(fine))
                 if index not in selected),
                reverse=True,
            )))
            return
        if subtotal > target:
            return
        previous = None
        for index in range(start, len(fine)):
            if fine[index] == previous:
                continue
            previous = fine[index]
            choose(index + 1, subtotal + fine[index], selected + (index,))

    choose(0, 0, ())
    return any(coarsens(remainder, remaining_targets)
               for remainder in remainders)


def refinement_flow(delta):
    """Check that a zero-mass vector is a positive fine-to-coarse boundary."""

    require(sum(delta.values()) == 0, "deletion difference lost total mass")
    positive = tuple((shape, value) for shape, value in delta.items()
                     if value > 0)
    negative = tuple((shape, -value) for shape, value in delta.items()
                     if value < 0)
    demand = sum(value for _, value in negative)
    source = ("source",)
    sink = ("sink",)
    adjacency = defaultdict(list)
    capacity = {}

    def add_edge(left, right, value):
        adjacency[left].append(right)
        adjacency[right].append(left)
        capacity[left, right] = Fraction(value)
        capacity[right, left] = Fraction(0)

    for shape, value in positive:
        add_edge(source, ("positive", shape), value)
    for shape, value in negative:
        add_edge(("negative", shape), sink, value)
    for high, _ in positive:
        for low, _ in negative:
            if coarsens(low, high):
                add_edge(("positive", high), ("negative", low), demand)

    flow = Fraction(0)
    while True:
        parent = {source: None}
        queue = [source]
        for vertex in queue:
            for neighbour in adjacency[vertex]:
                if neighbour not in parent and capacity[vertex, neighbour] > 0:
                    parent[neighbour] = vertex
                    queue.append(neighbour)
            if sink in parent:
                break
        if sink not in parent:
            break
        increment = Fraction(demand)
        vertex = sink
        while parent[vertex] is not None:
            increment = min(increment, capacity[parent[vertex], vertex])
            vertex = parent[vertex]
        vertex = sink
        while parent[vertex] is not None:
            previous = parent[vertex]
            capacity[previous, vertex] -= increment
            capacity[vertex, previous] += increment
            vertex = previous
        flow += increment
    require(flow == demand, "same-label deletion violated refinement order")
    return len(positive), len(negative)


def difference(left, right):
    answer = defaultdict(int)
    for shape, value in left.items():
        answer[shape] += value
    for shape, value in right.items():
        answer[shape] -= value
    return {shape: value for shape, value in answer.items() if value}


def sign_gap_scalar(shape):
    """Scalar of Lbar_shape on the sign representation."""

    return 0 if all(part == 1 for part in shape) else 1


records = []
total_shapes = total_lowerings = total_covers = total_equal = 0
total_flow_positive = total_flow_negative = 0
raw_sign_negative = raw_sign_zero = raw_sign_positive = 0
minimum_weight_ratio = None
EXPECTED = {
    2: (2, 2, 1, 1, 0, 1, 0),
    3: (3, 4, 2, 0, 0, 1, 1),
    4: (5, 7, 5, 0, 2, 2, 1),
    5: (7, 12, 9, 0, 6, 2, 1),
    6: (11, 19, 17, 0, 14, 2, 1),
    7: (15, 30, 28, 0, 25, 2, 1),
    8: (22, 45, 47, 0, 44, 2, 1),
    9: (30, 67, 73, 0, 70, 2, 1),
}

for degree in range(2, 10):
    degree_shapes = partitions(degree)
    degree_lowerings = degree_covers = degree_equal = 0
    degree_raw_negative = degree_raw_zero = degree_raw_positive = 0
    for shape in degree_shapes:
        block = block_deletion(shape)
        labelled = labelled_deletion(shape)
        require(sum(labelled.values()) == degree,
                "labelled deletion is not mass preserving")
        for part, multiplicity in Counter(shape).items():
            target = lower(shape, part)
            require(
                Fraction(weight(shape), weight(target)) == part,
                "factorial diagonal intertwining failed",
            )
            require(
                Fraction(weight(shape) * multiplicity, weight(target))
                == labelled[target],
                "conjugated block deletion is not labelled deletion",
            )
            degree_lowerings += 1

        for coarse in covers(shape):
            degree_covers += 1
            require(coarsens(shape, coarse), "cover is not a coarsening")
            ratio = Fraction(weight(coarse), weight(shape))
            require(ratio >= 2, "proper merge did not increase factorial weight")
            minimum_weight_ratio = (ratio if minimum_weight_ratio is None
                                    else min(minimum_weight_ratio, ratio))

            labelled_delta = difference(
                labelled_deletion(coarse), labelled_deletion(shape))
            if not labelled_delta:
                degree_equal += 1
            positive, negative = refinement_flow(labelled_delta)
            total_flow_positive += positive
            total_flow_negative += negative

            raw_delta = difference(
                block_deletion(coarse), block_deletion(shape))
            # Raw block deletion has total mass -1 on every proper cover.
            require(sum(raw_delta.values()) == -1,
                    "raw block-deletion cover mass drift")
            sign_scalar = sum(
                coefficient * sign_gap_scalar(target)
                for target, coefficient in raw_delta.items()
            )
            if sign_scalar < 0:
                degree_raw_negative += 1
            elif sign_scalar == 0:
                degree_raw_zero += 1
            else:
                degree_raw_positive += 1

    total_shapes += len(degree_shapes)
    total_lowerings += degree_lowerings
    total_covers += degree_covers
    total_equal += degree_equal
    raw_sign_negative += degree_raw_negative
    raw_sign_zero += degree_raw_zero
    raw_sign_positive += degree_raw_positive
    require(
        (len(degree_shapes), degree_lowerings, degree_covers, degree_equal,
         degree_raw_negative, degree_raw_zero, degree_raw_positive)
        == EXPECTED[degree],
        "degree cover census drift",
    )
    records.append(
        f"N{degree}:partitions={len(degree_shapes)}:lowerings={degree_lowerings}:"
        f"covers={degree_covers}:labelled_equal={degree_equal}:"
        f"raw_sign=neg{degree_raw_negative}/zero{degree_raw_zero}/"
        f"pos{degree_raw_positive}"
    )

# The first nonzero raw failures live in degree four.
raw_31 = difference(block_deletion((4,)), block_deletion((3, 1)))
raw_22 = difference(block_deletion((4,)), block_deletion((2, 2)))
require(raw_31 == {(2, 1): -1}, "(3,1)->(4) hostile drift")
require(raw_22 == {(3,): 1, (2, 1): -2}, "(2,2)->(4) hostile drift")
require(sum(value * sign_gap_scalar(shape)
            for shape, value in raw_31.items()) == -1,
        "(3,1)->(4) sign hostile disappeared")
require(sum(value * sign_gap_scalar(shape)
            for shape, value in raw_22.items()) == -1,
        "(2,2)->(4) sign hostile disappeared")
require(
    (total_shapes, total_lowerings, total_covers, total_equal,
     total_flow_positive, total_flow_negative,
     raw_sign_negative, raw_sign_zero, raw_sign_positive,
     minimum_weight_ratio)
    == (95, 186, 182, 1, 340, 325, 161, 14, 7, Fraction(2)),
    "global cover census drift",
)

print("weights=w_lambda=product_i(lambda_i!):unique_up_to_global_scale")
print("intertwining=W_N^-1*block_delete*W_(N+1)=labelled_delete")
print("labelled_delete=delete_same_uniform_label:refinement_stochastic")
print("carrier=Kbar_lambda=w_lambda*Lbar_lambda:coarsening_monotone")
for record in records:
    print(record)
print(
    f"totals=partitions:{total_shapes}:lowerings:{total_lowerings}:"
    f"covers:{total_covers}:labelled_equal:{total_equal}:"
    f"flow_support:{total_flow_positive}/{total_flow_negative}:"
    f"raw_sign:neg{raw_sign_negative}/zero{raw_sign_zero}/"
    f"pos{raw_sign_positive}:min_weight_ratio:{minimum_weight_ratio}"
)
print("first_raw_hostiles=(3,1)->(4):-Lbar_(2,1);(2,2)->(4):Lbar_3-2Lbar_(2,1)")
print("all_exact_checks=PASS")
