#!/usr/bin/env python3
"""Exact finite checks for the norm-one forest and Berggren absorber dictionary.

Standard library only. The infinite theorem needs the descent proof; these
checks independently audit its algebra and bounded enumeration consequences.
"""
import argparse
import hashlib
import json
from math import gcd, isqrt
from pathlib import Path

M = {
    "A": ((1, -2, 2), (2, -1, 2), (2, -2, 3)),
    "B": ((1, 2, 2), (2, 1, 2), (2, 2, 3)),
    "C": ((-1, 2, 2), (-2, 1, 2), (-2, 2, 3)),
}
INV = {
    "A": ((1, 2, -2), (-2, -1, 2), (-2, -2, 3)),
    "B": ((1, 2, -2), (2, 1, -2), (-2, -2, 3)),
    "C": ((-1, -2, 2), (2, 1, -2), (-2, -2, 3)),
}
I = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
J = ((1, 0, 0), (0, 1, 0), (0, 0, -1))
ROOTS = ((1, 1, 1), (2, 1, 2), (1, 2, 2))
PPT_ROOT = (3, 4, 5)


def apply(m, v):
    return tuple(sum(x * y for x, y in zip(row, v)) for row in m)


def multiply(a, b):
    return tuple(tuple(sum(a[i][k] * b[k][j] for k in range(3))
                       for j in range(3)) for i in range(3))


def norm(v):
    a, b, c = v
    return a * a + b * b - c * c


def parent(v):
    a, b, c = v
    r, s, h = a + 2 * b - 2 * c, 2 * a + b - 2 * c, 3 * c - 2 * a - 2 * b
    assert r and s and (r > 0 or s > 0) and 0 < h < c
    label = "A" if s < 0 else "C" if r < 0 else "B"
    p = apply(INV[label], v)
    assert min(p) > 0 and norm(p) == norm(v)
    assert sum(min(apply(inv, v)) > 0 for inv in INV.values()) == 1
    assert apply(M[label], p) == v
    return label, p


def address(v, root):
    reverse = []
    while v != root:
        label, v = parent(v)
        reverse.append(label)
    return "".join(reversed(reverse))


def from_address(word, root):
    for label in word:  # Words are chronological: AAAB means B*A*A*A.
        root = apply(M[label], root)
    return root


def bounded_tree(root, cap):
    seen, stack = {root}, [root]
    while stack:
        p = stack.pop()
        children = [apply(m, p) for m in M.values()]
        assert len(set(children)) == 3
        for label, v in zip(M, children):
            assert min(v) > 0 and norm(v) == norm(root) and v[2] > p[2]
            assert parent(v) == (label, p)
            if v[2] <= cap:
                assert v not in seen, ("collision", v)
                seen.add(v)
                stack.append(v)
    return seen


def brute(cap, q):
    result = set()
    for c in range(1, cap + 1, 2):
        for a in range(1, c + 1, 2):
            b2 = c * c + q - a * a
            b = isqrt(b2)
            if b > 0 and b * b == b2:
                if q == 1 and b % 2 == 1:
                    result.add((a, b, c))
                if q == 0 and b % 2 == 0 and gcd(a, b) == 1:
                    result.add((a, b, c))
    return result


def brute_normone_forest(cap):
    result = set()
    for c in range(1, cap + 1):
        for a in range(1, c + 1):
            b2 = c * c + 1 - a * a
            b = isqrt(b2)
            if b > 0 and b * b == b2:
                result.add((a, b, c))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--height", type=int, default=2001)
    parser.add_argument("--depth", type=int, default=8)
    args = parser.parse_args()
    assert args.height >= 7 and args.depth >= 0
    for label, m in M.items():
        transpose = tuple(zip(*m))
        assert multiply(multiply(transpose, J), m) == J
        assert multiply(INV[label], m) == multiply(m, INV[label]) == I
    sectors = [bounded_tree(root, args.height) for root in ROOTS]
    odd = sectors[0]
    forest = set().union(*sectors)
    assert sum(map(len, sectors)) == len(forest)
    assert forest == brute_normone_forest(args.height)
    for root, sector in zip(ROOTS, sectors):
        assert all(tuple(z % 2 for z in v) == tuple(z % 2 for z in root)
                   for v in sector)
    ppt = bounded_tree(PPT_ROOT, args.height)
    assert odd == brute(args.height, 1)
    assert ppt == brute(args.height, 0)
    rays = {(1, c, c) for c in range(1, args.height + 1, 2)}
    rays |= {(c, 1, c) for c in range(1, args.height + 1, 2)}
    assert {v for v in odd if min(v[:2]) == 1} == rays
    all_rays = {(1, c, c) for c in range(1, args.height + 1)}
    all_rays |= {(c, 1, c) for c in range(1, args.height + 1)}
    assert {v for v in forest if min(v[:2]) == 1} == all_rays
    assert {v for v in forest if v[0] + 2 * v[1] - 2 * v[2] == 0} == {(2, 1, 2)}
    assert {v for v in forest if 2 * v[0] + v[1] - 2 * v[2] == 0} == {(1, 2, 2)}
    assert {v for v in forest if 3 * v[2] - 2 * v[0] - 2 * v[1] <= 0} == set(ROOTS)
    # Full word-level correspondence to a fixed depth, independent of height.
    roots = ROOTS + (PPT_ROOT,)
    layer = [("", roots)]
    seen = [set() for _ in roots]
    for depth in range(args.depth + 1):
        following = []
        assert len(layer) == 3 ** depth
        for word, nodes in layer:
            for i, (root, v) in enumerate(zip(roots, nodes)):
                assert v not in seen[i]
                seen[i].add(v)
                assert address(v, root) == word
                assert tuple(z % 2 for z in v) == tuple(z % 2 for z in root)
                assert gcd(gcd(v[0], v[1]), v[2]) == 1
            if depth < args.depth:
                following.extend((word + label, tuple(apply(m, v) for v in nodes))
                                 for label, m in M.items())
        layer = following
    # Exhaust all distinct ordered positive edges with t <= 10, for each
    # nondegenerate norm-one coefficient triple in the height-bounded set.
    pairs = [(x, t * t - x, t) for t in range(2, 11)
             for x in range(1, t * t) if 2 * x != t * t]
    tests = odd_tests = 0
    for a, b, c in forest:
        if a == 1 or b == 1:
            continue
        for x, y, t in pairs:
            u, v = a * a * t * t - x, b * b * t * t - y
            assert (u + x, x + y, y + v, u + v) == ((a * t) ** 2, t * t, (b * t) ** 2, (c * t) ** 2)
            assert min(u, v) > max(x, y) and len({u, x, y, v}) == 4
            tests += 1
            odd_tests += a % 2 == b % 2 == c % 2 == 1
    expected = {"": ((3, 4, 5), (1, 1, 1)),
                "AAA": ((9, 40, 41), (1, 7, 7)),
                "AAAB": ((171, 140, 221), (29, 23, 37)),
                "B": ((21, 20, 29), (5, 5, 7)),
                "BB": ((119, 120, 169), (29, 29, 41))}
    for word, (zero, one) in expected.items():
        assert from_address(word, (3, 4, 5)) == zero
        assert from_address(word, (1, 1, 1)) == one
    even_seed = from_address("AAA", (2, 1, 2))
    assert even_seed == (8, 31, 32)
    for k in range(32):
        assert from_address("A" * k, (2, 1, 2)) == (2 * k + 2, 2 * (k + 1) ** 2 - 1, 2 * (k + 1) ** 2)
        assert from_address("C" * k, (1, 2, 2)) == (2 * (k + 1) ** 2 - 1, 2 * k + 2, 2 * (k + 1) ** 2)
    assert apply(M["A"], (2, 1, 2)) == (4, 7, 8)
    assert apply(M["C"], (1, 2, 2)) == (7, 4, 8)
    diagonal = {v for v in forest if v[0] == v[1]}
    central, v = set(), (1, 1, 1)
    while v[2] <= args.height:
        central.add(v)
        v = apply(M["B"], v)
    assert diagonal == central
    result = {"status": "PASS", "arithmetic": "exact Python integers; standard library only",
              "verifier_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "height_cap": args.height, "odd_norm_one_nodes": len(odd),
              "parity_oriented_primitive_pythagorean_nodes": len(ppt),
              "degenerate_boundary_nodes": len(rays),
              "nondegenerate_absorber_coefficients": len(odd - rays),
              "complete_address_depth": args.depth, "addresses_checked_per_tree": len(seen[0]),
              "norm_one_forest_nodes": len(forest),
              "norm_one_forest_tree_counts": {str(root): len(sector) for root, sector in zip(ROOTS, sectors)},
              "all_degenerate_boundary_nodes": len(all_rays),
              "all_nondegenerate_absorber_coefficients": len(forest - all_rays),
              "ordered_input_edges_per_coefficient": len(pairs), "absorber_checks": tests,
              "odd_sector_absorber_checks": odd_tests,
              "diagonal_nodes": sorted(diagonal, key=lambda v: v[2]),
              "seed_addresses": expected, "AAA_even_odd_even_seed": even_seed,
              "scope": "Finite algebra/enumeration audit only; infinite classification follows from the descent proof."}
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
