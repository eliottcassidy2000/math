#!/usr/bin/env python3
"""Clean-room THM-4118 audit by literal child construction.

No Start/End/Q data and no primary code are imported.  Every ear child is
constructed, its Hamiltonian paths are counted by subset DP, and unit-state
components are recovered by disjoint-set union.
"""

from functools import reduce
from hashlib import sha256
from itertools import combinations, permutations
import json
from math import gcd


EXPECTED_SEMANTIC = "9b1eb26c20a92b7220f1534169d466c551ceec03188d1d438b2c76e9ac03589b"


def tournament(code, n):
    out = [0] * n
    position = 0
    for i in range(n):
        for j in range(i + 1, n):
            if code & (1 << position):
                out[i] |= 1 << j
            else:
                out[j] |= 1 << i
            position += 1
    return tuple(out)


def tournament_code(out):
    code = 0
    position = 0
    for i in range(len(out)):
        for j in range(i + 1, len(out)):
            if out[i] & (1 << j):
                code |= 1 << position
            position += 1
    return code


def permute(out, order):
    changed = [0] * len(out)
    for i, old_i in enumerate(order):
        for j, old_j in enumerate(order):
            if out[old_i] & (1 << old_j):
                changed[i] |= 1 << j
    return tuple(changed)


def strong(out):
    n = len(out)
    target = (1 << n) - 1
    for initial in range(n):
        reached = 1 << initial
        while True:
            expanded = reached
            todo = reached
            while todo:
                bit = todo & -todo
                todo -= bit
                expanded |= out[bit.bit_length() - 1]
            if expanded == reached:
                break
            reached = expanded
        if reached != target:
            return False
    return True


def hamiltonian_paths(out):
    n = len(out)
    ways = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        ways[1 << v][v] = 1
    for mask in range(1, 1 << n):
        if mask & (mask - 1) == 0:
            continue
        for last in range(n):
            bit = 1 << last
            if not (mask & bit):
                continue
            previous_mask = mask ^ bit
            total = 0
            todo = previous_mask
            while todo:
                previous_bit = todo & -todo
                todo -= previous_bit
                previous = previous_bit.bit_length() - 1
                if out[previous] & bit:
                    total += ways[previous_mask][previous]
            ways[mask][last] = total
    return sum(ways[-1])


def ear_child(base, signature):
    n = len(base)
    child = list(base) + [0]
    new_bit = 1 << n
    for old in range(n):
        if signature & (1 << old):
            child[n] |= 1 << old
        else:
            child[old] |= new_bit
    return tuple(child)


def ear_values(base):
    return tuple(
        hamiltonian_paths(ear_child(base, signature))
        for signature in range(1 << len(base))
    )


def gcd_values(values):
    return reduce(gcd, (abs(x) for x in values), 0)


def extract_quadratic(base, values):
    n = len(base)
    h = hamiltonian_paths(base)
    assert h == values[0] == values[-1]
    linear = tuple(values[1 << i] - h for i in range(n))
    curvature = {}
    for i, j in combinations(range(n), 2):
        curvature[i, j] = (
            values[1 << i] + values[1 << j]
            - h - values[(1 << i) | (1 << j)]
        )
    for signature, actual in enumerate(values):
        chosen = [i for i in range(n) if signature & (1 << i)]
        predicted = h + sum(linear[i] for i in chosen)
        predicted -= sum(curvature[min(i, j), max(i, j)]
                         for i, j in combinations(chosen, 2))
        assert predicted == actual
    lattice = gcd_values(linear + tuple(curvature.values()))
    response_lattice = gcd_values(
        tuple(values[s] - h for s in range(1, (1 << n) - 1))
    )
    assert lattice == response_lattice
    return h, lattice


class DSU:
    def __init__(self, items):
        self.parent = {item: item for item in items}

    def find(self, item):
        while self.parent[item] != item:
            self.parent[item] = self.parent[self.parent[item]]
            item = self.parent[item]
        return item

    def union(self, first, second):
        a = self.find(first)
        b = self.find(second)
        if a != b:
            self.parent[b] = a


def state_components(values):
    n = len(values).bit_length() - 1
    full = (1 << n) - 1
    domain = tuple(range(1, full))
    domain_set = set(domain)
    dsu = DSU(domain)
    for mask in domain:
        candidates = [mask ^ (1 << i) for i in range(n)]
        inside = [i for i in range(n) if mask & (1 << i)]
        outside = [j for j in range(n) if not (mask & (1 << j))]
        candidates += [mask ^ (1 << i) ^ (1 << j)
                       for i in inside for j in outside]
        for other in candidates:
            if other in domain_set and abs(values[other] - values[mask]) in (0, 2):
                dsu.union(mask, other)
    groups = {}
    for mask in domain:
        groups.setdefault(dsu.find(mask), set()).add(mask)
    result = []
    for masks in groups.values():
        image = {values[mask] for mask in masks}
        lo, hi = min(image), max(image)
        assert image == set(range(lo, hi + 1, 2))
        result.append((tuple(sorted(masks)), lo, hi))
    return tuple(sorted(result))


def analyze(code, n):
    base = tournament(code, n)
    values = ear_values(base)
    h, lattice = extract_quadratic(base, values)
    assert all(value % 2 for value in values)
    components = state_components(values)
    image = tuple(sorted(set(values[1:-1])))
    solid = image == tuple(range(image[0], image[-1] + 1, 2))
    spans = any(lo == image[0] and hi == image[-1] for _, lo, hi in components)
    largest = max(((lo, hi) for _, lo, hi in components),
                  key=lambda interval: interval[1] - interval[0])
    return {
        "H": h, "lattice": lattice, "values": values, "image": image,
        "solid": solid, "spans": spans, "components": components,
        "largest": largest,
    }


def orbit_representatives(n):
    unseen = set(range(1 << (n * (n - 1) // 2)))
    representatives = []
    orders = tuple(permutations(range(n)))
    while unseen:
        code = min(unseen)
        base = tournament(code, n)
        orbit = {tournament_code(permute(base, order)) for order in orders}
        unseen.difference_update(orbit)
        representatives.append(min(orbit))
    return tuple(representatives)


def main():
    census = {}
    for n, expected in ((3, (2, 2, 2, 2)),
                        (4, (24, 24, 24, 24)),
                        (5, (544, 0, 544, 0))):
        totals = [0, 0, 0, 0]
        for code in range(1 << (n * (n - 1) // 2)):
            if not strong(tournament(code, n)):
                continue
            record = analyze(code, n)
            totals[0] += 1
            totals[1] += record["solid"]
            totals[2] += record["lattice"] == 2
            totals[3] += record["spans"]
        census[n] = tuple(totals)
        assert census[n] == expected

    named = {code: analyze(code, 5) for code in (8, 759, 1015)}
    assert named[8]["image"] == named[1015]["image"] == (
        15, 17, 19, 23, 25, 27, 29, 33, 37, 41
    )
    assert named[759]["image"] == (15, 17, 19, 23, 29, 31, 33, 37, 43)
    assert tuple(named[code]["lattice"] for code in (8, 759, 1015)) == (2, 2, 2)
    assert (len(named[8]["components"]), named[8]["largest"]) == (10, (15, 19))
    assert (len(named[759]["components"]), named[759]["largest"]) == (7, (29, 33))

    selected = ((8, (16, 18)), (8, (5, 21)), (759, (5, 9, 17)))
    intervals = []
    for code, masks in selected:
        found = [component for component in named[code]["components"]
                 if set(masks) <= set(component[0])]
        assert len(found) == 1
        intervals.append(found[0][1:])
    assert tuple(intervals) == ((23, 25), (27, 29), (29, 33))

    representatives = orbit_representatives(6)
    assert len(representatives) == 56
    strong_codes = [code for code in representatives if strong(tournament(code, 6))]
    records = [(code, analyze(code, 6)) for code in strong_codes]
    histogram = {}
    for _, record in records:
        histogram[record["lattice"]] = histogram.get(record["lattice"], 0) + 1
    assert len(records) == 35 and histogram == {2: 34, 6: 1}
    assert sum(record["solid"] for _, record in records) == 0
    hostile = analyze(30571, 6)
    assert (hostile["H"], hostile["lattice"], hostile["image"]) == (
        45, 6, (93, 111, 123, 129, 135, 141, 153, 159, 189)
    )
    hostile_orbit = {
        tournament_code(permute(tournament(30571, 6), order))
        for order in permutations(range(6))
    }
    exceptional = next(code for code, record in records if record["lattice"] == 6)
    assert exceptional in hostile_orbit

    ledger = {
        "labelled_census_n3_n5": census,
        "code8_image": named[8]["image"],
        "code8_components_largest": (len(named[8]["components"]), named[8]["largest"]),
        "code759_image": named[759]["image"],
        "code759_components_largest": (len(named[759]["components"]), named[759]["largest"]),
        "selected_component_intervals": tuple(intervals),
        "order6_unlabeled_strong": len(records),
        "order6_lattice_histogram": histogram,
        "code30571": (hostile["H"], hostile["lattice"], hostile["image"]),
    }
    semantic = sha256(json.dumps(
        ledger, sort_keys=True, separators=(",", ":")
    ).encode()).hexdigest()
    if EXPECTED_SEMANTIC is not None:
        assert semantic == EXPECTED_SEMANTIC

    print("status=PASS")
    print("implementation=literal ear-child construction plus Hamiltonian subset DP and DSU")
    print(f"labelled_census_n3_n5={census}")
    print(f"code8=d2;image={named[8]['image']};components={len(named[8]['components'])};largest={named[8]['largest']}")
    print(f"code759=d2;image={named[759]['image']};components={len(named[759]['components'])};largest={named[759]['largest']}")
    print(f"selected_component_intervals={tuple(intervals)};certified_union=(23,33)")
    print(f"order6_unlabeled=56;strong={len(records)};lattice_histogram={histogram};solid=0")
    print(f"code30571=H{hostile['H']};d{hostile['lattice']};image={hostile['image']}")
    print(f"semantic_sha256={semantic}")


if __name__ == "__main__":
    main()
