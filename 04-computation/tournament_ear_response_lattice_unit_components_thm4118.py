#!/usr/bin/env python3
"""Exact referee for THM-4118 using the Start/End/Q ear boundary."""

from collections import deque
from functools import reduce
from hashlib import sha256
from itertools import combinations, permutations
import json
from math import gcd


EXPECTED_SEMANTIC = "9b1eb26c20a92b7220f1534169d466c551ceec03188d1d438b2c76e9ac03589b"


def decode(code, n):
    adjacency = [0] * n
    bit = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (code >> bit) & 1:
                adjacency[i] |= 1 << j
            else:
                adjacency[j] |= 1 << i
            bit += 1
    return tuple(adjacency)


def encode(adjacency):
    code = 0
    bit = 0
    for i in range(len(adjacency)):
        for j in range(i + 1, len(adjacency)):
            code |= ((adjacency[i] >> j) & 1) << bit
            bit += 1
    return code


def relabel(adjacency, new_to_old):
    n = len(adjacency)
    out = [0] * n
    for i in range(n):
        for j in range(n):
            if adjacency[new_to_old[i]] & (1 << new_to_old[j]):
                out[i] |= 1 << j
    return tuple(out)


def is_strong(adjacency):
    n = len(adjacency)
    full = (1 << n) - 1
    for root in range(n):
        seen = 1 << root
        frontier = seen
        while frontier:
            nxt = 0
            todo = frontier
            while todo:
                bit = todo & -todo
                todo -= bit
                nxt |= adjacency[bit.bit_length() - 1]
            nxt &= ~seen
            seen |= nxt
            frontier = nxt
        if seen != full:
            return False
    return True


def boundary(adjacency):
    n = len(adjacency)
    size = 1 << n
    full = size - 1
    starts = [[0] * n for _ in range(size)]
    ends = [[0] * n for _ in range(size)]
    for v in range(n):
        starts[1 << v][v] = 1
        ends[1 << v][v] = 1
    for mask in range(1, size):
        if mask & (mask - 1) == 0:
            continue
        for v in range(n):
            bit = 1 << v
            if not (mask & bit):
                continue
            rest = mask ^ bit
            todo = rest
            while todo:
                other_bit = todo & -todo
                todo -= other_bit
                other = other_bit.bit_length() - 1
                if adjacency[v] & other_bit:
                    starts[mask][v] += starts[rest][other]
                if adjacency[other] & bit:
                    ends[mask][v] += ends[rest][other]
    q = [[0] * n for _ in range(n)]
    for left_mask in range(1, full):
        right_mask = full ^ left_mask
        left_todo = left_mask
        while left_todo:
            left_bit = left_todo & -left_todo
            left_todo -= left_bit
            left = left_bit.bit_length() - 1
            left_count = ends[left_mask][left]
            if not left_count:
                continue
            right_todo = right_mask
            while right_todo:
                right_bit = right_todo & -right_todo
                right_todo -= right_bit
                right = right_bit.bit_length() - 1
                q[left][right] += left_count * starts[right_mask][right]
    return tuple(starts[full]), tuple(ends[full]), tuple(tuple(row) for row in q)


def responses(data):
    starts, ends, q = data
    n = len(starts)
    values = []
    for mask in range(1 << n):
        value = 0
        for i in range(n):
            value += starts[i] if (mask >> i) & 1 else ends[i]
        for i in range(n):
            if (mask >> i) & 1:
                continue
            for j in range(n):
                if (mask >> j) & 1:
                    value += q[i][j]
        values.append(value)
    return tuple(values)


def gcd_all(values):
    return reduce(gcd, (abs(x) for x in values), 0)


def coefficients(values, q):
    n = len(q)
    h = values[0]
    linear = tuple(values[1 << i] - h for i in range(n))
    curvature = tuple(tuple(
        0 if i == j else q[i][j] + q[j][i] for j in range(n)
    ) for i in range(n))
    flat = linear + tuple(
        curvature[i][j] for i, j in combinations(range(n), 2)
    )
    return h, linear, curvature, gcd_all(flat)


def verify_quadratic(values, h, linear, curvature):
    n = len(linear)
    for mask, actual in enumerate(values):
        predicted = h
        chosen = [i for i in range(n) if (mask >> i) & 1]
        predicted += sum(linear[i] for i in chosen)
        predicted -= sum(curvature[i][j] for i, j in combinations(chosen, 2))
        assert predicted == actual


def neighbors(mask, n):
    for i in range(n):
        yield mask ^ (1 << i)
    inside = [i for i in range(n) if (mask >> i) & 1]
    outside = [i for i in range(n) if not ((mask >> i) & 1)]
    for i in inside:
        for j in outside:
            yield mask ^ (1 << i) ^ (1 << j)


def unit_components(values):
    n = len(values).bit_length() - 1
    full = (1 << n) - 1
    unseen = set(range(1, full))
    components = []
    while unseen:
        root = min(unseen)
        unseen.remove(root)
        component = {root}
        queue = deque([root])
        while queue:
            mask = queue.popleft()
            for nxt in neighbors(mask, n):
                if nxt in unseen and abs(values[nxt] - values[mask]) in (0, 2):
                    unseen.remove(nxt)
                    component.add(nxt)
                    queue.append(nxt)
        image = {values[mask] for mask in component}
        interval = set(range(min(image), max(image) + 1, 2))
        assert image == interval
        components.append((tuple(sorted(component)), min(image), max(image)))
    return tuple(components)


def row(code, n):
    adjacency = decode(code, n)
    data = boundary(adjacency)
    values = responses(data)
    h, linear, curvature, lattice = coefficients(values, data[2])
    verify_quadratic(values, h, linear, curvature)
    response_gcd = gcd_all(tuple(values[mask] - h for mask in range(1, (1 << n) - 1)))
    assert response_gcd == lattice
    assert all(value % 2 == 1 for value in values)
    assert lattice % 2 == 0
    image = tuple(sorted(set(values[1:-1])))
    components = unit_components(values)
    solid = image == tuple(range(image[0], image[-1] + 1, 2))
    spans = any(lo == image[0] and hi == image[-1] for _, lo, hi in components)
    largest = max(((lo, hi) for _, lo, hi in components),
                  key=lambda x: x[1] - x[0])
    return {
        "H": h, "values": values, "image": image, "lattice": lattice,
        "solid": solid, "components": components, "spans": spans,
        "largest": largest,
    }


def unlabeled_representatives(n):
    pair_count = n * (n - 1) // 2
    seen = set()
    representatives = []
    relabelings = tuple(permutations(range(n)))
    for code in range(1 << pair_count):
        if code in seen:
            continue
        adjacency = decode(code, n)
        orbit = {encode(relabel(adjacency, p)) for p in relabelings}
        seen.update(orbit)
        representatives.append(min(orbit))
    return tuple(representatives)


def main():
    census = {}
    for n, expected in ((3, (2, 2, 2, 2)),
                        (4, (24, 24, 24, 24)),
                        (5, (544, 0, 544, 0))):
        strong = solid = lattice_two = spanning = 0
        for code in range(1 << (n * (n - 1) // 2)):
            if not is_strong(decode(code, n)):
                continue
            record = row(code, n)
            strong += 1
            solid += record["solid"]
            lattice_two += record["lattice"] == 2
            spanning += record["spans"]
        census[n] = (strong, solid, lattice_two, spanning)
        assert census[n] == expected

    named = {code: row(code, 5) for code in (8, 759, 1015)}
    assert named[8]["image"] == named[1015]["image"] == (
        15, 17, 19, 23, 25, 27, 29, 33, 37, 41
    )
    assert named[759]["image"] == (15, 17, 19, 23, 29, 31, 33, 37, 43)
    assert all(named[code]["lattice"] == 2 for code in named)
    assert len(named[8]["components"]) == 10 and named[8]["largest"] == (15, 19)
    assert len(named[759]["components"]) == 7 and named[759]["largest"] == (29, 33)
    assert relabel(decode(8, 5), (4, 3, 2, 1, 0)) == decode(1015, 5)

    chain = ((8, (16, 18)), (8, (5, 21)), (759, (5, 9, 17)))
    chain_intervals = []
    for code, masks in chain:
        record = named[code]
        containing = [component for component in record["components"]
                      if set(masks) <= set(component[0])]
        assert len(containing) == 1
        chain_intervals.append((containing[0][1], containing[0][2]))
    assert tuple(chain_intervals) == ((23, 25), (27, 29), (29, 33))
    assert all(b <= c + 2 for (_, b), (c, _) in zip(chain_intervals, chain_intervals[1:]))

    reps6 = unlabeled_representatives(6)
    assert len(reps6) == 56
    strong6 = [code for code in reps6 if is_strong(decode(code, 6))]
    records6 = [(code, row(code, 6)) for code in strong6]
    histogram = {}
    for _, record in records6:
        histogram[record["lattice"]] = histogram.get(record["lattice"], 0) + 1
    assert len(records6) == 35 and histogram == {2: 34, 6: 1}
    hostile = row(30571, 6)
    assert hostile["H"] == 45 and hostile["lattice"] == 6
    assert hostile["image"] == (93, 111, 123, 129, 135, 141, 153, 159, 189)
    hostile_canonical = min(
        encode(relabel(decode(30571, 6), p)) for p in permutations(range(6))
    )
    exceptional = next(code for code, record in records6 if record["lattice"] == 6)
    assert hostile_canonical == exceptional
    assert sum(record["solid"] for _, record in records6) == 0

    ledger = {
        "labelled_census_n3_n5": census,
        "code8_image": named[8]["image"],
        "code8_components_largest": (len(named[8]["components"]), named[8]["largest"]),
        "code759_image": named[759]["image"],
        "code759_components_largest": (len(named[759]["components"]), named[759]["largest"]),
        "selected_component_intervals": tuple(chain_intervals),
        "order6_unlabeled_strong": len(records6),
        "order6_lattice_histogram": histogram,
        "code30571": (hostile["H"], hostile["lattice"], hostile["image"]),
    }
    semantic = sha256(json.dumps(
        ledger, sort_keys=True, separators=(",", ":")
    ).encode()).hexdigest()
    if EXPECTED_SEMANTIC is not None:
        assert semantic == EXPECTED_SEMANTIC

    print("status=PASS")
    print("implementation=Start/End/Q boundary plus exhaustive labelled and orbit census")
    print(f"labelled_census_n3_n5={census}")
    print(f"code8=d2;image={named[8]['image']};components={len(named[8]['components'])};largest={named[8]['largest']}")
    print(f"code759=d2;image={named[759]['image']};components={len(named[759]['components'])};largest={named[759]['largest']}")
    print(f"selected_component_intervals={tuple(chain_intervals)};certified_union=(23,33)")
    print(f"order6_unlabeled=56;strong={len(records6)};lattice_histogram={histogram};solid=0")
    print(f"code30571=H{hostile['H']};d{hostile['lattice']};image={hostile['image']}")
    print(f"semantic_sha256={semantic}")


if __name__ == "__main__":
    main()
