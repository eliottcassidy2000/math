#!/usr/bin/env python3
"""Independent hostile audit of origin/main THM-3756.

This script does not import the candidate implementation.  It starts from
primitive Euclid pairs and the standard three Berggren matrices on ordered
triples, then compares those objects with the proposed odd-root chart.
"""

from collections import Counter
from hashlib import sha256
from math import gcd, isqrt


MAX_N = 1200
MAX_SHELL_R = 1999
MAX_CONE_R = 1200
MAX_BOUNDARY_S = 9999
TREE_DEPTH = 8
HEAP_DEPTH = 10

BERGGREN = {
    "L": ((1, -2, 2), (2, -1, 2), (2, -2, 3)),
    "M": ((1, 2, 2), (2, 1, 2), (2, 2, 3)),
    "R": ((-1, 2, 2), (-2, 1, 2), (-2, 2, 3)),
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def odd_root(rank):
    return 2 * rank - 1


def triple(pair):
    r, s = pair
    require(r > s >= 1, ("cone domain", pair))
    q, d = odd_root(r), odd_root(s)
    return q * d, (q * q - d * d) // 2, (q * q + d * d) // 2


def inverse_triple(value):
    a, b, c = value
    require(a > 0 and b > 0 and c > 0 and a * a + b * b == c * c,
            ("positive triple domain", value))
    require(a % 2 == 1 and b % 2 == 0, ("ordered legs", value))
    q, d = isqrt(b + c), isqrt(c - b)
    require((q * q, d * d) == (b + c, c - b), ("square inverse", value))
    require(q % 2 == d % 2 == 1, ("odd inverse roots", value))
    return (q + 1) // 2, (d + 1) // 2


def mat_vec(matrix, vector):
    return tuple(sum(coefficient * entry for coefficient, entry in zip(row, vector))
                 for row in matrix)


def affine_children(pair):
    r, s = pair
    return {
        "L": (r + 2 * s - 1, s),
        "M": (2 * r + s - 1, r),
        "R": (2 * r - s, r),
    }


def inverse_cone(pair):
    r, s = pair
    require(r > s >= 1, ("inverse cone domain", pair))
    if r == 3 * s - 1:
        return None
    if r >= 3 * s:
        answer = "L", (r - 2 * s + 1, s)
    elif 2 * s <= r <= 3 * s - 2:
        answer = "M", (s, r - 2 * s + 1)
    elif s < r <= 2 * s - 1:
        answer = "R", (s, 2 * s - r)
    else:
        raise RuntimeError(("uncovered cone point", pair))
    require(answer[1][0] < r and answer[1][0] > answer[1][1] >= 1,
            ("invalid strict parent", pair, answer))
    require(affine_children(answer[1])[answer[0]] == pair,
            ("parent round trip", pair, answer))
    return answer


def content(pair):
    return gcd(odd_root(pair[0]), odd_root(pair[1]))


def triangular(n):
    return n * (n + 1) // 2


def ambient_address(pair):
    r, s = pair
    return triangular(r - 2) + s


def heap_children(address):
    return tuple(3 * address + digit for digit in (1, 2, 3))


def main():
    # Ordered-leg hostile and inverse endpoints.
    require(inverse_triple((3, 4, 5)) == (2, 1), "root inverse")
    require(isqrt(3 + 5) ** 2 != 3 + 5, "leg-swap hostile disappeared")
    require(triple((2, 1)) == (3, 4, 5), "root chart")

    # A disjoint source census: every primitive opposite-parity Euclid pair.
    euclid_count = 0
    euclid_pairs = set()
    euclid_triples = set()
    matrix_checks = 0
    for n in range(2, MAX_N):
        for m in range(1, n):
            if gcd(m, n) != 1 or (m + n) % 2 != 1:
                continue
            value = (n * n - m * m, 2 * m * n, n * n + m * m)
            pair = (m + n + 1) // 2, (n - m + 1) // 2
            require(odd_root(pair[0]) == m + n, ("outer root", m, n, pair))
            require(odd_root(pair[1]) == n - m, ("inner root", m, n, pair))
            require(gcd(odd_root(pair[0]), odd_root(pair[1])) == 1,
                    ("chart coprimality", m, n, pair))
            require(triple(pair) == value, ("Euclid/chart mismatch", m, n))
            require(inverse_triple(value) == pair, ("lossless inverse", value))
            for label, child in affine_children(pair).items():
                require(
                    triple(child) == mat_vec(BERGGREN[label], value),
                    ("matrix/affine branch mismatch", label, pair),
                )
                require(inverse_cone(child) == (label, pair),
                        ("child inverse mismatch", label, pair))
                matrix_checks += 1
            euclid_count += 1
            euclid_pairs.add(pair)
            euclid_triples.add(value)
    require(euclid_count == len(euclid_pairs) == len(euclid_triples),
            "Euclid chart collision")

    # Exact shell counts, including q=3 and both prime/composite controls.
    shell_total = 0
    shell_digest_rows = []
    for r in range(2, MAX_SHELL_R + 1):
        q = odd_root(r)
        odd_choices = tuple(d for d in range(1, q, 2) if gcd(d, q) == 1)
        all_reduced = sum(gcd(d, q) == 1 for d in range(1, q))
        require(2 * len(odd_choices) == all_reduced,
                ("phi/2 shell failure", r, q, len(odd_choices), all_reduced))
        require(odd_choices[0] == 1 and odd_choices[-1] <= q - 2,
                ("shell endpoint failure", r, odd_choices[:1], odd_choices[-1:]))
        shell_total += len(odd_choices)
        if r <= 30:
            shell_digest_rows.append((r, q, len(odd_choices)))
    require(shell_digest_rows[0] == (2, 3, 1), "r=2 endpoint shell")
    require(next(row for row in shell_digest_rows if row[0] == 5) == (5, 9, 3),
            "first composite shell")

    # Every inverse interval boundary and every cone point through MAX_CONE_R.
    cone_counts = Counter()
    cone_round_trips = 0
    for r in range(2, MAX_CONE_R + 1):
        for s in range(1, r):
            pair = (r, s)
            step = inverse_cone(pair)
            g = content(pair)
            value = triple(pair)
            require(gcd(gcd(value[0], value[1]), value[2]) == g * g,
                    ("exact square content", pair, value, g))
            if step is None:
                cone_counts["root"] += 1
                require(r == 3 * s - 1 and g == 2 * s - 1,
                        ("root line/content", pair, g))
                require(value == (3 * g * g, 4 * g * g, 5 * g * g),
                        ("scaled root", pair, value, g))
            else:
                label, previous = step
                cone_counts[label] += 1
                require(content(previous) == g, ("inverse changed content", pair))
                cone_round_trips += 1
    for s in range(1, MAX_BOUNDARY_S + 1):
        cases = []
        if s + 1 <= 2 * s - 1:
            cases.append((s + 1, "R"))
            cases.append((2 * s - 1, "R"))
        if 2 * s <= 3 * s - 2:
            cases.append((2 * s, "M"))
            cases.append((3 * s - 2, "M"))
        cases.extend(((3 * s - 1, None), (3 * s, "L")))
        for r, expected in cases:
            result = inverse_cone((r, s))
            require((None if result is None else result[0]) == expected,
                    ("inverse interval boundary", r, s, expected, result))

    # Fixed-content copies: distinct 3^h levels and unique inverse words.
    forest_controls = {}
    for g in (1, 3, 5, 9, 15, 101, 9999):
        root = ((3 * g + 1) // 2, (g + 1) // 2)
        require(content(root) == g, ("root content", g, root))
        level = {root}
        seen = set()
        for depth in range(TREE_DEPTH + 1):
            require(len(level) == 3**depth and not (level & seen),
                    ("forest level collision", g, depth, len(level)))
            for pair in level:
                require(content(pair) == g, ("forest content drift", g, pair))
                q0 = odd_root(pair[0]) // g
                d0 = odd_root(pair[1]) // g
                normalized = ((q0 + 1) // 2, (d0 + 1) // 2)
                require(triple(pair) == tuple(g * g * entry for entry in triple(normalized)),
                        ("component scaling", g, pair, normalized))
                cursor = pair
                while cursor != root:
                    step = inverse_cone(cursor)
                    require(step is not None, ("premature component root", g, cursor))
                    cursor = step[1]
                require(cursor == root, ("wrong component root", g, pair, cursor))
            seen.update(level)
            if depth < TREE_DEPTH:
                level = {
                    child
                    for pair in level
                    for child in affine_children(pair).values()
                }
        forest_controls[g] = len(seen)

    # Ambient and selected-shell addresses are independently packed.
    ambient = {
        ambient_address((r, s))
        for r in range(2, MAX_CONE_R + 1)
        for s in range(1, r)
    }
    require(ambient == set(range(1, triangular(MAX_CONE_R - 1) + 1)),
            "ambient address packing")
    require(ambient_address((5, 2)) == 8 and content((5, 2)) == 3,
            "first ambient primitive hole")

    selected = []
    selected_pairs = []
    offset = 0
    for r in range(2, MAX_SHELL_R + 1):
        allowed_s = [s for s in range(1, r) if gcd(odd_root(r), odd_root(s)) == 1]
        for local_index, s in enumerate(allowed_s, 1):
            selected.append(offset + local_index)
            selected_pairs.append((r, s))
        offset += len(allowed_s)
    require(selected == list(range(1, offset + 1)), "selected address packing")
    require(len(selected_pairs) == len(set(selected_pairs)), "selected pair collision")

    # THM-3382 heap recursion: finite words through depth ten map to one interval.
    heap_levels = [{0}]
    heap_seen = set()
    for depth in range(HEAP_DEPTH + 1):
        level = heap_levels[depth]
        require(len(level) == 3**depth and not (level & heap_seen),
                ("heap level collision", depth, len(level)))
        heap_seen.update(level)
        if depth < HEAP_DEPTH:
            heap_levels.append({child for address in level for child in heap_children(address)})
    require(heap_seen == set(range((3 ** (HEAP_DEPTH + 1) - 1) // 2)),
            "heap address interval")

    # Scope hostiles: degree/order address loss and non-odd-square scale.
    require(triple((3, 1)) == (5, 12, 13), "rank-three left hostile")
    require(triple((3, 2)) == (15, 8, 17), "rank-three right hostile")
    require(
        tuple(child[0] for child in affine_children((3, 1)).values())
        != tuple(child[0] for child in affine_children((3, 2)).values()),
        "rank unexpectedly determines child signatures",
    )
    scaled_two = (6, 8, 10)
    require(isqrt(scaled_two[1] + scaled_two[2]) ** 2 != sum(scaled_two[1:]),
            "arbitrary scale entered odd-root chart")

    semantic = (
        euclid_count,
        matrix_checks,
        shell_total,
        tuple(shell_digest_rows),
        tuple(sorted(cone_counts.items())),
        cone_round_trips,
        tuple(sorted(forest_controls.items())),
        len(ambient),
        len(selected),
        tuple(len(level) for level in heap_levels),
        max(heap_seen),
    )
    print("source=independent_Euclid_pairs+standard_3x3_Berggren_matrices")
    print("ordered_root_and_swap_hostile=PASS")
    print("euclid_pairs_n_lt_1200=" + repr(euclid_count))
    print("matrix_affine_child_checks=" + repr(matrix_checks))
    print("shells_r_le_1999=" + repr((shell_total, tuple(shell_digest_rows))))
    print("cone_r_le_1200=" + repr((tuple(sorted(cone_counts.items())), cone_round_trips)))
    print("inverse_boundaries_s_le_9999=PASS")
    print("fixed_content_tree_controls=" + repr(tuple(sorted(forest_controls.items()))))
    print("ambient_address_interval=" + repr((1, max(ambient), len(ambient))))
    print("selected_address_interval=" + repr((1, len(selected), len(selected))))
    print("heap_through_depth_10=" + repr((len(heap_seen), min(heap_seen), max(heap_seen))))
    print("scope_hostiles=rank_collision;child_signature_loss;scale_2_outside")
    print("semantic_sha256=" + sha256(repr(semantic).encode("ascii")).hexdigest())
    print("RESULT=PROMOTE;ARTIFACT=THM-3756_INDEPENDENT_AUDIT")


if __name__ == "__main__":
    main()
