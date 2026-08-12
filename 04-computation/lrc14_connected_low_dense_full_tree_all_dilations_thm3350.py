#!/usr/bin/env python3
"""Exact dense-low full-tree compiler for the reflected LRC(14) atlas.

This THM-3350 referee does three logically separate things.

1. It freezes the universal singleton-debt maximum needed by the structural
   connected-low argument.  Every six-tuple of distinct positive levels has
   a rank permutation in {1,...,6}, and replacing each level by its rank can
   only increase the debt.  Hence the 649*6! finite rank scan is universal.
2. It checks the two intrinsic connected-low shapes having only one high
   complete-graph edge.  These are the only shapes not covered merely by
   putting two high edges in a spanning tree.
3. For every body, every labelling of either exceptional shape, and every
   common dilation, it certifies a deterministic maximum-homogeneous-limit
   spanning tree.  A general (not only low-channel) midpoint/shear estimate
   supplies the tail and direct rational interval geometry checks the head.

It is the tracked exact companion for the dense all-dilation branch of
THM-3350.  Set LRC_REPO only to override the repository inferred from this
file's canonical location.
"""

from __future__ import annotations

import hashlib
import os
from collections import Counter
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from itertools import permutations
from math import gcd
from pathlib import Path


ROOT = Path(
    os.environ.get("LRC_REPO", Path(__file__).resolve().parents[1])
).resolve()
SELECTOR = ROOT / "04-computation/lrc14_j7_reflected_low_two_star_limit_selector_scout_20260812.py"
EXPECTED_SELECTOR_SHA256 = "cd8b08087f0f7e1e0c0c7d0be673629c7c2702c170c5c1e771e1d76df1d3cd1c"

DENSE_SHAPES = (
    (1, 2, 3, 4, 6, 12),
    (2, 3, 4, 6, 8, 12),
)
EXPECTED = {
    "body_count": 649,
    "debt_rows": 467280,
    "dense_assignments": 934560,
    "threshold_hist": ((1, 934335), (2, 223), (3, 2)),
    "head_checks": 227,
    "head_failures": 0,
    "selected_tree_types": 145763,
    "debt_semantic_sha256": "e919f2bc07fae564d518cad9547f639c2693fbf723d3bc96dbc1b4c962f2cba7",
    "assignment_semantic_sha256": "9d185440b228327d438fa70cabbc5474fe40613abed0d728bf915f6ff1a78fca",
    "head_semantic_sha256": "9a955420daa11e583e92384288e16506b346335d580032cf95f1a37e71a7c053",
    "global_debt_max": F(186636088362, 11773143757375),
    "global_debt_body": (1, 2, 3, 4, 6, 12),
    "global_debt_levels": (6, 5, 4, 3, 2, 1),
    "weakest_head_margin": F(1359521303047, 13562237063235),
    "weakest_head_body": (1, 2, 3, 4, 6, 12),
    "weakest_head_levels": (4, 3, 2, 1, 12, 6),
    "weakest_head_tree": ((1, 4), (0, 3), (3, 4), (1, 2), (0, 5)),
    "weakest_limit": (
        F(148335185300354078659, 1924191320819535630720),
        (2, 3, 4, 8, 10, 12),
        1680,
        900,
        (1, 2, 3, 4, 6, 12),
        (6, 4, 2, 1, 12, 3),
        ((1, 5), (0, 3), (3, 4), (0, 1), (2, 4)),
        F(3163, 40320),
    ),
    "max_threshold_rows": (
        (
            3,
            (1, 2, 3, 4, 6, 12),
            168,
            90,
            (1, 2, 3, 4, 6, 12),
            (4, 12, 2, 6, 3, 1),
            ((4, 5), (1, 4), (1, 2), (3, 5), (0, 4)),
            F(14171, 114240),
        ),
        (
            3,
            (1, 2, 3, 4, 6, 12),
            168,
            90,
            (1, 2, 3, 4, 6, 12),
            (12, 6, 4, 2, 3, 1),
            ((3, 4), (4, 5), (2, 4), (0, 5), (1, 5)),
            F(355261, 2242240),
        ),
    ),
}


def require(condition: bool, detail: object) -> None:
    if not condition:
        raise RuntimeError(detail)


def sha256_lf(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def load(path: Path):
    spec = spec_from_file_location("connected_low_dense_base", path)
    require(spec is not None and spec.loader is not None, path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


require(SELECTOR.is_file(), ("missing selector", SELECTOR))
require(
    sha256_lf(SELECTOR) == EXPECTED_SELECTOR_SHA256,
    ("selector hash", sha256_lf(SELECTOR), EXPECTED_SELECTOR_SHA256),
)
M = load(SELECTOR)

EDGES = tuple((i, j) for i in range(6) for j in range(i + 1, 6))


def primitive_pair(p: int, q: int) -> tuple[int, int, int]:
    d = gcd(p, q)
    return d, p // d, q // d


def is_low(p: int, q: int) -> bool:
    _, P, Q = primitive_pair(p, q)
    return P + Q <= 7


def gamma(level: int) -> F:
    return F(1, 2) if level % 2 == 0 else F(level * level + 1, 2 * level * level)


def all_channel_error(
    scale: int,
    ruler: int,
    e: int,
    p: int,
    f: int,
    q: int,
) -> F:
    """General midpoint/shear bound for one common-dilation channel.

    If p=dP,q=dQ and h=sd, the first two terms are the exact jump-grid
    contraction bounds.  The last term uses TV(Phi')<=4/(PQ) per circle and
    allows the determinant path to wrap floor(|C|/L)+1 times.
    """
    d, P, Q = primitive_pair(p, q)
    h = scale * d
    C = Q * e - P * f
    wraps = abs(C) // ruler + 1
    return (
        gamma(P) * F(e * P, h * ruler * P - e)
        + gamma(Q) * F(f * Q, h * ruler * Q - f)
        + F(abs(C) * wraps, 2 * h * h * ruler * P * Q)
    )


def debt(body, ruler: int, levels, scale: int) -> F:
    return sum(
        (F(e, 7 * (scale * q * ruler - e)) for e, q in zip(body, levels)),
        F(0),
    )


class DSU:
    def __init__(self, n: int):
        self.parent = list(range(n))

    def find(self, x: int) -> int:
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def join(self, x: int, y: int) -> bool:
        x, y = self.find(x), self.find(y)
        if x == y:
            return False
        self.parent[y] = x
        return True


def maximum_tree(weights) -> tuple[tuple[tuple[int, int], ...], F]:
    """Kruskal maximum tree; exact ties take lexicographically largest edge."""
    dsu = DSU(6)
    chosen = []
    total = F(0)
    for weight, edge in sorted(
        ((weights[edge], edge) for edge in EDGES), reverse=True
    ):
        if dsu.join(*edge):
            chosen.append(edge)
            total += weight
            if len(chosen) == 5:
                break
    require(len(chosen) == 5, ("disconnected complete graph", weights))
    return tuple(chosen), total


def direct_pair_mass(ruler: int, cell: int, e: int, p: int, f: int, q: int) -> F:
    """Literal rational interval intersection, valid at every ratio."""
    first = M.R.reflected_level_arcs(ruler, e, p, cell)
    second = M.R.reflected_level_arcs(ruler, f, q, cell)
    answer = M.interval_intersection_mass(first, second)
    # On the domain where the optimized floor-moment engine accepts the row,
    # demand an independent exact equality.
    try:
        fast = M.physical_mass(ruler, cell, e, p, f, q)
    except RuntimeError:
        pass
    else:
        require(answer == fast, ("direct/fast mismatch", ruler, cell, e, p, f, q, answer, fast))
    return answer


def tree_error(scale, ruler, body, levels, tree) -> F:
    return sum(
        (
            all_channel_error(
                scale,
                ruler,
                body[i],
                levels[i],
                body[j],
                levels[j],
            )
            for i, j in tree
        ),
        F(0),
    )


def physical_tree_credit(scale, ruler, cell, body, levels, tree) -> F:
    return sum(
        (
            direct_pair_mass(
                ruler,
                cell,
                body[i],
                scale * levels[i],
                body[j],
                scale * levels[j],
            )
            for i, j in tree
        ),
        F(0),
    )


def tail_threshold(limit, ruler, body, levels, tree) -> int:
    for scale in range(1, 1000):
        if limit - tree_error(scale, ruler, body, levels, tree) > debt(
            body, ruler, levels, scale
        ):
            return scale
    raise RuntimeError(("threshold >=1000", body, ruler, levels, tree, limit))


def body_cell(body, ruler):
    cell, *_ = M.body_geometry(body, ruler)
    return cell


def main() -> None:
    bodies = M.MS.body_universe()
    require(len(bodies) == EXPECTED["body_count"], len(bodies))

    dense_profile = tuple(
        (
            shape,
            tuple(edge for edge in EDGES if is_low(shape[edge[0]], shape[edge[1]])),
            tuple(edge for edge in EDGES if not is_low(shape[edge[0]], shape[edge[1]])),
        )
        for shape in DENSE_SHAPES
    )
    require(tuple(len(row[1]) for row in dense_profile) == (14, 14), dense_profile)
    require(tuple(len(row[2]) for row in dense_profile) == (1, 1), dense_profile)

    # Universal debt rank reduction: q_i >= rank(q_i) coordinatewise.
    max_debt = None
    max_debt_row = None
    debt_rows = 0
    debt_semantic = hashlib.sha256()
    rank_levels = tuple(permutations(range(1, 7)))
    for body, ruler in bodies:
        for levels in rank_levels:
            value = debt(body, ruler, levels, 1)
            row = (value, body, ruler, levels)
            if max_debt_row is None or row > max_debt_row:
                max_debt_row = row
            max_debt = value if max_debt is None or value > max_debt else max_debt
            debt_rows += 1
            debt_semantic.update((repr((body, ruler, levels, value)) + "\n").encode())

    require(debt_rows == EXPECTED["debt_rows"], debt_rows)
    require(max_debt == EXPECTED["global_debt_max"], max_debt)
    require(max_debt_row[1] == EXPECTED["global_debt_body"], max_debt_row)
    require(max_debt_row[3] == EXPECTED["global_debt_levels"], max_debt_row)
    require(F(2, 105) > max_debt, (F(2, 105), max_debt))

    threshold_hist = Counter()
    assignment_count = head_checks = failures = 0
    weakest_head = None
    weakest_limit = None
    max_threshold_rows = []
    assignment_semantic = hashlib.sha256()
    head_semantic = hashlib.sha256()
    selected_tree_types = set()

    for body, ruler in bodies:
        cell = body_cell(body, ruler)
        for shape in DENSE_SHAPES:
            perms = tuple(permutations(shape))
            pair_cache = {}
            for edge in EDGES:
                i, j = edge
                for p in shape:
                    for q in shape:
                        if p != q:
                            pair_cache[edge, p, q] = M.homogeneous_limit(
                                ruler, cell, body[i], p, body[j], q
                            )

            for levels in perms:
                weights = {
                    edge: pair_cache[edge, levels[edge[0]], levels[edge[1]]]
                    for edge in EDGES
                }
                tree, limit = maximum_tree(weights)
                limit_margin = limit - debt(body, ruler, levels, 1)
                require(limit_margin > 0, ("homogeneous tree/debt", body, levels, tree, limit_margin))
                limit_row = (limit_margin, body, ruler, cell, shape, levels, tree, limit)
                if weakest_limit is None or limit_row < weakest_limit:
                    weakest_limit = limit_row
                threshold = tail_threshold(limit, ruler, body, levels, tree)
                threshold_hist[threshold] += 1
                selected_tree_types.add(tree)
                assignment_count += 1
                assignment_semantic.update(
                    (
                        repr((body, ruler, cell, shape, levels, tree, limit, threshold))
                        + "\n"
                    ).encode()
                )
                if not max_threshold_rows or threshold > max_threshold_rows[0][0]:
                    max_threshold_rows = [(threshold, body, ruler, cell, shape, levels, tree, limit)]
                elif threshold == max_threshold_rows[0][0]:
                    max_threshold_rows.append((threshold, body, ruler, cell, shape, levels, tree, limit))

                for scale in range(1, threshold):
                    credit = physical_tree_credit(
                        scale, ruler, cell, body, levels, tree
                    )
                    invoice = debt(body, ruler, levels, scale)
                    margin = credit - invoice
                    head_checks += 1
                    if margin <= 0:
                        failures += 1
                    row = (
                        margin,
                        body,
                        ruler,
                        cell,
                        shape,
                        levels,
                        tree,
                        scale,
                        credit,
                        invoice,
                        limit,
                        threshold,
                    )
                    if weakest_head is None or row < weakest_head:
                        weakest_head = row
                    head_semantic.update((repr(row[1:]) + "\n").encode())

    actual_hist = tuple(sorted(threshold_hist.items()))
    require(assignment_count == EXPECTED["dense_assignments"], assignment_count)
    require(actual_hist == EXPECTED["threshold_hist"], actual_hist)
    require(head_checks == EXPECTED["head_checks"], head_checks)
    require(failures == EXPECTED["head_failures"], failures)
    require(weakest_head is not None, "no finite heads")
    require(weakest_head[0] == EXPECTED["weakest_head_margin"], weakest_head)
    require(weakest_head[1] == EXPECTED["weakest_head_body"], weakest_head)
    require(weakest_head[5] == EXPECTED["weakest_head_levels"], weakest_head)
    require(weakest_head[6] == EXPECTED["weakest_head_tree"], weakest_head)
    require(weakest_limit == EXPECTED["weakest_limit"], weakest_limit)
    require(len(selected_tree_types) == EXPECTED["selected_tree_types"], len(selected_tree_types))
    require(tuple(max_threshold_rows) == EXPECTED["max_threshold_rows"], max_threshold_rows)
    require(
        debt_semantic.hexdigest() == EXPECTED["debt_semantic_sha256"],
        debt_semantic.hexdigest(),
    )
    require(
        assignment_semantic.hexdigest() == EXPECTED["assignment_semantic_sha256"],
        assignment_semantic.hexdigest(),
    )
    require(
        head_semantic.hexdigest() == EXPECTED["head_semantic_sha256"],
        head_semantic.hexdigest(),
    )

    print("LRC14 CONNECTED-LOW DENSE FULL-TREE COMPILER")
    print("status=VERIFIED-EXACT finite universe + analytic all-channel tail")
    print(f"selector_sha256={sha256_lf(SELECTOR)}")
    print(f"bodies={len(bodies)}")
    print(f"dense_profiles={dense_profile}")
    print(f"universal_debt_rank_rows={debt_rows}")
    print(f"universal_debt_max_row={max_debt_row}")
    print(f"two_high_floor={F(2,105)};strict_gap={F(2,105)-max_debt}")
    print(f"debt_semantic_sha256={debt_semantic.hexdigest()}")
    print(f"dense_assignments={assignment_count}")
    print(f"selected_tree_types={len(selected_tree_types)}")
    print(f"tail_threshold_hist={actual_hist}")
    print(f"weakest_homogeneous_limit_margin_row={weakest_limit}")
    print(f"head_checks={head_checks};failures={failures}")
    print(f"max_threshold_rows={tuple(max_threshold_rows)}")
    print(f"weakest_head_row={weakest_head}")
    print(f"assignment_semantic_sha256={assignment_semantic.hexdigest()}")
    print(f"head_semantic_sha256={head_semantic.hexdigest()}")
    print("conclusion=dense-low exceptional shapes close at every common dilation")


if __name__ == "__main__":
    main()
