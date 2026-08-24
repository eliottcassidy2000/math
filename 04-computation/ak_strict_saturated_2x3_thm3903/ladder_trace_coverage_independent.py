#!/usr/bin/env python3
"""Independent combinatorial coverage audit for saturated [2,3] traces.

This does not repeat the Groebner elimination.  It independently proves that
the three multiround families fed to that elimination exhaust every trace
not already covered by the one- and two-round arguments, and recomputes all
family/headroom counts without importing any other scratch module.
"""

from itertools import combinations, combinations_with_replacement, product
import hashlib
import json
import sys

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


VERTICES = tuple(range(6))
EDGES = ((0, 3), (1, 4), (2, 5), (0, 1), (1, 2), (3, 4), (4, 5))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def has_cycle(live):
    live = set(live)
    parent = {v: v for v in live}

    def root(v):
        while parent[v] != v:
            parent[v] = parent[parent[v]]
            v = parent[v]
        return v

    for u, v in EDGES:
        if u not in live or v not in live:
            continue
        ru, rv = root(u), root(v)
        if ru == rv:
            return True
        parent[ru] = rv
    return False


def ordered_partitions(items, blocks):
    items = tuple(items)
    for assignment in product(range(blocks), repeat=len(items)):
        if set(assignment) != set(range(blocks)):
            continue
        yield tuple(frozenset(items[i] for i, value in enumerate(assignment)
                              if value == block) for block in range(blocks))


def forest_compatible_trace(waves):
    live = set(VERTICES)
    # Before each of the final two waves, a forest is allowed.  Earlier, a
    # forest would leave more than the single additional round permitted by
    # THM-3875 and hence cannot realize this trace.
    for j, wave in enumerate(waves):
        live.difference_update(wave)
        if j <= len(waves) - 3 and not has_cycle(live):
            return False
    return True


def headroom_profiles(sizes):
    result = []
    for profile in product(range(5), repeat=len(sizes)):
        if profile[-1] != sizes[-1]:
            continue
        if any(profile[j] < sizes[j] for j in range(len(sizes))):
            continue
        if any(not 0 <= profile[j] - profile[j + 1] <= sizes[j]
               for j in range(len(sizes) - 1)):
            continue
        result.append(profile)
    return tuple(result)


def canonical_counts():
    corner3 = rung3 = corner4 = 0
    # Fix A as first singleton.
    remaining = (1, 2, 3, 4, 5)
    for size in range(1, 5):
        for second in combinations(remaining, size):
            sizes = (1, size, 5 - size)
            corner3 += len(headroom_profiles(sizes))
    # Fix AD as first outer rung.
    remaining = (1, 2, 4, 5)
    for size in range(1, 4):
        for second in combinations(remaining, size):
            sizes = (2, size, 4 - size)
            rung3 += len(headroom_profiles(sizes))
    # Fix A, then D alone; split the remaining C4 in two nonempty waves.
    cell = (1, 2, 4, 5)
    for size in range(1, 4):
        for third in combinations(cell, size):
            sizes = (1, 1, size, 4 - size)
            corner4 += len(headroom_profiles(sizes))
    return corner3, rung3, corner4


def main():
    gates = 0
    require(has_cycle(VERTICES), "full ladder lost its cycles")
    require(not has_cycle({1, 2, 4}), "forest detector failed")
    gates += 2

    compatible = {}
    first_sets = {}
    for q in range(3, 7):
        traces = tuple(w for w in ordered_partitions(VERTICES, q)
                       if forest_compatible_trace(w))
        compatible[q] = len(traces)
        first_sets[q] = tuple(sorted(set(tuple(sorted(w[0])) for w in traces)))
    require(compatible == {3: 148, 4: 56, 5: 0, 6: 0},
            "multiround trace coverage changed")
    require(first_sets[3] == ((0,), (0, 3), (2,), (2, 5), (3,), (5,)),
            "three-round first-wave classification changed")
    require(first_sets[4] == ((0,), (2,), (3,), (5,)),
            "four-round first-wave classification changed")
    gates += 3

    # Four-round traces must delete the other vertex of the outer rung alone
    # in wave two.
    mates = {0: 3, 3: 0, 2: 5, 5: 2}
    for waves in ordered_partitions(VERTICES, 4):
        if forest_compatible_trace(waves):
            first = next(iter(waves[0]))
            require(waves[1] == frozenset((mates[first],)),
                    "unexpected four-round second wave")
            gates += 1

    corner3, rung3, corner4 = canonical_counts()
    seed_multisets = len(tuple(combinations_with_replacement(VERTICES, 3)))
    require((seed_multisets, corner3, rung3, corner4) == (56, 70, 60, 80),
            "canonical family counts changed")
    require((seed_multisets * rung3, seed_multisets * corner3,
             seed_multisets * corner4) == (3360, 3920, 4480),
            "Groebner job universe changed")
    gates += 2

    semantic = {
        "compatible_ordered_traces": compatible,
        "first_sets": first_sets,
        "canonical_profiles_per_seed": (rung3, corner3, corner4),
        "seed_multisets": seed_multisets,
        "exact_jobs": (3360, 3920, 4480),
    }
    digest = hashlib.sha256(json.dumps(
        semantic, sort_keys=True, separators=(",", ":"), default=str
    ).encode()).hexdigest()
    print("AK_STRICT_2x3_MULTIROUND_COVERAGE_INDEPENDENT_20260823")
    print("status=INDEPENDENT_COMBINATORIAL_COVERAGE_PASS;NOT_GROEBNER_REPLICATION")
    print(f"forest_compatible_ordered_traces={compatible}")
    print("canonical_first_families=(outer_rung_3,outer_corner_3,outer_corner_4)")
    print(f"profiles_per_seed=(rung3:{rung3},corner3:{corner3},corner4:{corner4})")
    print(f"seed_multisets={seed_multisets};exact_jobs=(3360,3920,4480)")
    print(f"semantic_sha256={digest}")
    print(f"gates={gates}")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
