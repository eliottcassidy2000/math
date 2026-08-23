#!/usr/bin/env python3
"""Independent wreath-recursive audit of Rule 30 level transitivity.

This implementation does not import the canonical companion.  It constructs
the finite tree permutations directly from the three wreath recursions and
audits the self-replicating stabilizer plus forward-monoid orbits.
"""

from __future__ import annotations

import collections
import hashlib
import json
import sys


sys.stdout.reconfigure(newline="\n")

MAX_DEPTH = 14
EXPECTED_SEMANTIC_SHA256 = "2654e5dd189e0bc70eedc82fdc9a57b810af0122c47314606042e2d96296a613"


def require(condition: bool, label: object) -> None:
    if not condition:
        raise RuntimeError(f"requirement failed: {label}")


def build_actions(max_depth: int) -> dict[int, tuple[tuple[int, ...], ...]]:
    actions = {0: ((0,), (0,), (0,))}
    for depth in range(1, max_depth + 1):
        old_a, old_b, old_c = actions[depth - 1]
        size = 1 << depth
        a = [0] * size
        b = [0] * size
        c = [0] * size
        for value in range(size):
            bit = value & 1
            suffix = value >> 1
            if bit == 0:
                a[value] = old_a[suffix] << 1
                b[value] = 1 | (old_c[suffix] << 1)
                c[value] = 1 | (old_a[suffix] << 1)
            else:
                a[value] = 1 | (old_b[suffix] << 1)
                b[value] = old_b[suffix] << 1
                c[value] = old_b[suffix] << 1
        rows = (tuple(a), tuple(b), tuple(c))
        for label, row in zip("ABC", rows):
            require(sorted(row) == list(range(size)), ("permutation", depth, label))
        actions[depth] = rows
    return actions


def inverse(row: tuple[int, ...]) -> tuple[int, ...]:
    answer = [0] * len(row)
    for source, target in enumerate(row):
        answer[target] = source
    return tuple(answer)


def forward_orbit(a: tuple[int, ...], b: tuple[int, ...], start: int):
    queue = collections.deque((start,))
    distance = {start: 0}
    while queue:
        value = queue.popleft()
        for row in (a, b):
            target = row[value]
            if target not in distance:
                distance[target] = distance[value] + 1
                queue.append(target)
    return distance


def permutation_digest(row: tuple[int, ...], depth: int) -> str:
    width = max(1, (depth + 7) // 8)
    payload = b"".join(value.to_bytes(width, "little") for value in row)
    return hashlib.sha256(payload).hexdigest()


def main() -> None:
    actions = build_actions(MAX_DEPTH)
    identity_rows = []
    orbit_rows = []

    for depth in range(1, MAX_DEPTH + 1):
        a, b, _ = actions[depth]
        inv_b = inverse(b)
        conjugate = tuple(b[a[inv_b[value]]] for value in range(1 << depth))
        if depth >= 2:
            small_a, small_b, _ = actions[depth - 1]
            for suffix in range(1 << (depth - 1)):
                value = suffix << 1
                require(a[value] == (small_a[suffix] << 1),
                        ("A section", depth, suffix))
                require(conjugate[value] == (small_b[suffix] << 1),
                        ("conjugate section", depth, suffix))
        require(all((conjugate[value] ^ value) % 2 == 0
                    for value in range(1 << depth)),
                ("root stabilizer", depth))
        identity_rows.append((depth, permutation_digest(conjugate, depth)))

        starts = (0,) + tuple(1 << exponent for exponent in range(1, depth))
        summaries = []
        for start in starts:
            distance = forward_orbit(a, b, start)
            require(len(distance) == 1 << depth,
                    ("forward monoid orbit", depth, start, len(distance)))
            active_endpoints = {value for value in distance if (value ^ start) & 1}
            require(active_endpoints == set(range(1, 1 << depth, 2)),
                    ("active endpoint half", depth, start))
            summaries.append((start, max(distance.values()), sum(distance.values())))
        orbit_rows.append((depth, tuple(summaries)))

    semantic = (MAX_DEPTH, tuple(identity_rows), tuple(orbit_rows))
    semantic_sha256 = hashlib.sha256(
        json.dumps(semantic, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    require(semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
            "frozen semantic transcript")

    print("RULE30_UNIVERSAL_STATIC_CLOSURE_INDEPENDENT_20260823")
    print("status=INDEPENDENT_FINITE_EXACT_AUDIT;direct_wreath_recursion")
    print(f"depths=1..{MAX_DEPTH};starts={sum(len(rows) for _, rows in orbit_rows)}")
    print("zero_eccentricities=" + repr(tuple(
        (depth, rows[0][1]) for depth, rows in orbit_rows
    )).replace(" ", ""))
    print("identity=(B^-1*A*B)|0=B;positive_AB_orbit=full_level")
    print(f"semantic_sha256={semantic_sha256}")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
