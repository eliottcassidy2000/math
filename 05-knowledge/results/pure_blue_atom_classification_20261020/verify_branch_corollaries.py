#!/usr/bin/env python3
"""Independent finite controls for the proved counting/geometry corollaries."""

from functools import lru_cache
from itertools import permutations
import json
from pathlib import Path


@lru_cache(None)
def words(weight):
    if weight < 0:
        return ()
    if weight == 0:
        return ((),)
    return tuple((letter,) + rest for letter in (1, 3) for rest in words(weight - letter))


def counts(limit):
    a = [0] * (limit + 3)
    a[0] = 1
    for j in range(1, len(a)):
        a[j] = a[j - 1] + (a[j - 3] if j >= 3 else 0)
    def get(j):
        return a[j] if j >= 0 else 0
    p = []
    for n in range(limit + 1):
        m = n // 2
        p.append(get(m) if n % 2 == 0 else get(m) + get(m - 1) + get(m - 2))
    for n, value in enumerate(p):
        assert value == (p[n - 2] if n >= 2 else 0) + (p[n - 6] if n >= 6 else 0) + int(n in (0, 1, 3, 5))
        if n >= 4:
            m = n // 2
            L = get(m) if n % 2 == 0 else get(m) + get(m - 1)
            J = 0 if n % 2 == 0 else get(m - 2)
            beta = 1 << (((n - 1) ** 2) // 4 - 1)
            assert beta >= L + 3 * J
    for m in range((limit - 4) // 2 + 1):
        assert p[2 * m + 1] == p[2 * m + 4]
    return p


def count_paths_and_aut(A):
    n = len(A)
    perms = list(permutations(range(n)))
    hp = sum(all(A[u][v] for u, v in zip(p, p[1:])) for p in perms)
    aut = sum(all(A[i][j] == A[p[i]][p[j]] for i in range(n) for j in range(n)) for p in perms)
    return hp, aut


def stack_A_C5_Aop(word):
    blocks = list(word) + [5] + list(reversed(word))
    n = sum(blocks)
    A = [[False] * n for _ in range(n)]
    start = 0
    for size in blocks:
        end = start + size
        for i in range(start, end):
            for j in range(start, end):
                A[i][j] = 0 < (j - i) % size <= size // 2
            for j in range(end, n):
                A[i][j] = True
        start = end
    return A


def complement(A, p):
    n = len(A)
    assert len(p) == n and sorted(p) == list(range(n))
    assert all(A[u][v] for u, v in zip(p, p[1:]))
    protected = {frozenset((u, v)) for u, v in zip(p, p[1:])}
    return [
        [i != j and (A[i][j] if frozenset((i, j)) in protected else not A[i][j]) for j in range(n)]
        for i in range(n)
    ]


def run():
    p = counts(100)
    bijection_rows = []
    for m in range(13):
        encoded = []
        for center, suffix in ((1, (1, 1)), (3, (3,)), (5, (3, 1))):
            for left in words(m - (center - 1) // 2):
                encoded.append(left + suffix)
        assert len(encoded) == len(set(encoded))
        assert set(encoded) == set(words(m + 2))
        bijection_rows.append({"source_order": 2 * m + 1, "target_order": 2 * m + 4, "classes": len(encoded)})
    center_paths = [(0, 1, 2, 3, 4), (0, 2, 3, 4, 1), (0, 2, 4, 1, 3)]
    geometry_rows = []
    for r in range(13):
        for w in words(r):
            A = stack_A_C5_Aop(w)
            n = len(A)
            center = set(range(r, r + 5))
            invariants = []
            for cp in center_paths:
                path = list(range(r)) + [r + i for i in cp] + list(range(r + 5, n))
                B = complement(A, path)
                if r:
                    assert {i for i, row in enumerate(B) if sum(row) == r + 2} == center
                    assert max(map(sum, B[:r])) <= r
                    assert min(map(sum, B[r + 5:])) >= r + 4
                C = [[B[r + i][r + j] for j in range(5)] for i in range(5)]
                invariants.append(count_paths_and_aut(C))
                # Literal score update, including both endpoint corrections.
                for i in range(n):
                    expected = n - 1 - sum(A[i]) + int(i == path[0]) - int(i == path[-1])
                    assert sum(B[i]) == expected
            assert invariants == [(11, 1), (15, 3), (13, 1)]
            assert len(set(hp for hp, _ in invariants)) == 3
            geometry_rows.append({"half_word": list(w), "order": n, "partner_center_invariants": invariants})
    # Independently count the staircase-reflection orbits.
    blue_dimensions = []
    for n in range(3, 51):
        tiles = {(i, j) for i in range(n) for j in range(i + 2, n)}
        seen = set()
        dim = 0
        for i, j in tiles:
            if (i, j) not in seen:
                seen.update(((i, j), (n - 1 - j, n - 1 - i)))
                dim += 1
        assert seen == tiles
        assert dim == (n - 1) ** 2 // 4
        blue_dimensions.append({"n": n, "blue_dimension": dim})
    return {
        "all_checks_passed": True,
        "counts_through_order": 100,
        "pure_blue_counts": p,
        "center_to_tail_bijections": bijection_rows,
        "score_shell_checks": geometry_rows,
        "score_shell_pure_blue_classes": len(geometry_rows),
        "score_shell_complement_partners": 3 * len(geometry_rows),
        "reflection_orbit_checks": blue_dimensions,
    }


if __name__ == "__main__":
    result = run()
    out = Path(__file__).with_name("branch_corollaries_verification.json")
    out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if not isinstance(v, list)}))
    print(f"Report written to {out}")
