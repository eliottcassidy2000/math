#!/usr/bin/env python3
"""Construct and independently validate the prefix proof's certificates.

No external dependencies, Hamiltonian-cycle search, or supplied census is used.
All arithmetic and certificate checks are exact integers.
"""

import argparse
from collections import Counter
from itertools import permutations
import json
from math import gcd, lcm
from pathlib import Path
import random


def permutation_order(p):
    seen = set()
    ans = 1
    for i in range(len(p)):
        if i in seen:
            continue
        j, length = i, 0
        while j not in seen:
            seen.add(j)
            j = p[j]
            length += 1
        ans = lcm(ans, length)
    assert len(seen) == len(p)
    return ans


def compose(p, q):
    return [p[q[i]] for i in range(len(p))]


def inverse(p):
    q = [None] * len(p)
    for i, j in enumerate(p):
        q[j] = i
    return q


def signs_to_set(n, mask):
    return {
        d if (mask >> (d - 1)) & 1 else n - d
        for d in range(1, (n + 1) // 2)
    }


def find_prefix(n, S):
    """Return proof-family witness or exact three-pattern residual."""
    m = n // 2
    assert n >= 9 and n % 2 and 1 in S
    if 2 not in S:
        if 3 in S:
            return [2, 0, 3, 1] + list(range(4, n)), "F1"
        for k in range(4, m + 1):
            if k in S and k - 2 not in S:
                q = list(range(n))
                q[0], q[k - 1] = q[k - 1], q[0]
                return q, f"F2(k={k})"
        assert S & set(range(1, m + 1)) == {1}
        return None, "A"
    for k in range(3, m):
        if k not in S:
            if k % 2:
                q = [1, k, 0] + list(range(2, k)) + list(range(k + 1, n))
            else:
                q = [1] + list(range(3, k + 1)) + [0, 2] + list(range(k + 1, n))
            return q, f"F3(k={k},{'odd' if k % 2 else 'even'})"
    assert set(range(1, m)) <= S
    return None, "C" if m in S else "B"


def construct(n, S):
    """Return a global HP word and the coordinate-space prefix permutation."""
    a = 1
    route = []
    for _ in range(3):
        assert gcd(a, n) == 1 and a in S
        shifted = {d for d in range(1, n) if (a * d) % n in S}
        q, family = find_prefix(n, shifted)
        route.append({"scale": a, "family": family})
        if q is not None:
            return [(a * x) % n for x in q], q, route
        a = (a * (-2 if family == "A" else 2)) % n
    raise AssertionError(("proof did not terminate", n, sorted(S), route))


def check_path(n, S, p):
    assert len(p) == n and sorted(p) == list(range(n))
    assert all((b - a) % n in S for a, b in zip(p, p[1:]))
    reverse = {p[i]: p[n - 1 - i] for i in range(n)}
    for u in range(n):
        for v in range(n):
            if u != v and (v - u) % n in S:
                if (reverse[u] - reverse[v]) % n not in S:
                    return [u, v], [reverse[v], reverse[u]]
    return None


def verify_certificate(n, S, p, q=None):
    failure = check_path(n, S, p)
    assert failure is not None, ("unexpected anti-reversal", n, sorted(S), p)
    if q is not None:
        m = n // 2
        assert sorted(q) == list(range(n))
        assert all(q[i] == i for i in range(m, n))
        assert permutation_order(q) % 2 == 0
        r = list(reversed(range(n)))
        # Independent full-array validation of the algebra in the proof.
        rqinv_r = compose(compose(r, inverse(q)), r)
        product = compose(q, rqinv_r)
        lhs = compose(compose(compose(q, r), inverse(q)), r)
        assert lhs == product
        supp1 = {i for i in range(n) if q[i] != i}
        supp2 = {i for i in range(n) if rqinv_r[i] != i}
        assert supp1.isdisjoint(supp2)
        assert permutation_order(product) == permutation_order(q)
        # Independently test the proof's prescribed, search-free violated arc.
        scale = p[q.index(1)]
        if len(supp1) == 2:
            # F2: use the actual oriented arc on {0,k-1}.
            u, v = 0, (scale * q[0]) % n
            if (v - u) % n not in S:
                u, v = v, u
        else:
            # F1 and either F3 parity: the base arc 0 -> 1.
            u, v = 0, scale
        rev = dict(zip(p, reversed(p)))
        assert (v - u) % n in S
        assert (rev[u] - rev[v]) % n not in S
        failure = ([u, v], [rev[v], rev[u]])
    return {
        "n": n,
        "S": sorted(S),
        "path": p,
        "failed_arc": failure[0],
        "absent_reversed_image_arc": failure[1],
    }


SEVEN = [
    ({1, 4, 5}, [0, 1, 2, 6, 4, 5, 3], [0, 5]),
    ({1, 2, 4}, [0, 1, 2, 4, 5, 6, 3], [0, 2]),
    ({1, 3, 5}, [0, 1, 2, 5, 3, 6, 4], [2, 0]),
    ({1, 2, 3}, [0, 1, 3, 4, 6, 2, 5], [0, 2]),
]


def check_seven():
    rows = []
    for S, p, pair in SEVEN:
        row = verify_certificate(7, S, p)
        rev = dict(zip(p, reversed(p)))
        u, v = pair
        assert (v - u) % 7 in S
        assert (rev[u] - rev[v]) % 7 not in S
        row["text_table_failed_arc"] = pair
        rows.append(row)
    assert {tuple(sorted(S)) for S, _, _ in SEVEN} == {
        tuple(sorted(signs_to_set(7, mask))) for mask in range(1, 8, 2)
    }
    return rows


def check_small_positive():
    rows = []
    for n in (3, 5):
        for mask in range(1, 1 << (n // 2), 2):
            S = signs_to_set(n, mask)
            paths = [
                p for p in permutations(range(n))
                if all((b - a) % n in S for a, b in zip(p, p[1:]))
            ]
            assert all(check_path(n, S, p) is None for p in paths)
            assert len(paths) == (3 if n == 3 else 15)
            rows.append({"n": n, "S": sorted(S), "HP": len(paths), "pure": True})
    return rows


def run(max_exhaustive, large_max):
    report = {
        "scope": "normalized odd circulants: 1 in S; proof-backed constructive certificate",
        "positive_controls": check_small_positive(),
        "order_seven": check_seven(),
        "exhaustive": [],
        "larger_checks": [],
    }
    for n in range(9, max_exhaustive + 1, 2):
        counts = Counter()
        maximum_steps = 0
        total = 0
        for mask in range(1, 1 << (n // 2), 2):
            S = signs_to_set(n, mask)
            p, q, route = construct(n, S)
            verify_certificate(n, S, p, q)
            counts[route[-1]["family"]] += 1
            maximum_steps = max(maximum_steps, len(route) - 1)
            total += 1
        row = {
            "n": n,
            "normalized_presentations": total,
            "excluded": total,
            "max_unit_reindexings": maximum_steps,
            "certificate_families": dict(counts),
        }
        report["exhaustive"].append(row)
        print(json.dumps(row), flush=True)

    rng = random.Random(20261020)
    for n in range(9, large_max + 1, 2):
        m = n // 2
        # Three worst-case residuals, both checkerboard signs, and random signs.
        sets = [
            {1} | set(range(m + 1, 2 * m)),
            set(range(1, m)) | {m + 1},
            set(range(1, m + 1)),
            signs_to_set(n, sum(1 << (d - 1) for d in range(1, m + 1, 2))),
            signs_to_set(n, rng.getrandbits(m) | 1),
        ]
        for S in sets:
            p, q, route = construct(n, S)
            verify_certificate(n, S, p, q)
        report["larger_checks"].append({"n": n, "structured_and_random_presentations": len(sets)})
    report["large_check_total"] = sum(
        row["structured_and_random_presentations"] for row in report["larger_checks"]
    )
    report["total_exhaustive_n_ge_9"] = sum(
        row["normalized_presentations"] for row in report["exhaustive"]
    )
    report["all_checks_passed"] = True
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-exhaustive", type=int, default=31)
    parser.add_argument("--large-max", type=int, default=501)
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("verification.json"))
    args = parser.parse_args()
    report = run(args.max_exhaustive, args.large_max)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(f"All checks passed; report written to {args.output}")
