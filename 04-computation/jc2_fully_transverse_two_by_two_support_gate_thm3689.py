#!/usr/bin/env python3
"""Exact companion for THM-3689.

The computation is degree-unbounded.  It enumerates the finite set of bucket
partitions of the eight nonzero contributions to Jac(P,Q)-1 in the fully
transverse two-by-two support chart, solves every exponent-equality system over
QQ, and prints the precise reason that every system leaves the chart.
"""

from __future__ import annotations

import ast
import hashlib
from collections import Counter
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def set_partitions_without_singletons(items: tuple[int, ...]):
    """Restricted-growth enumeration; every set partition occurs once."""

    blocks: list[list[int]] = []

    def visit(index: int):
        if index == len(items):
            if all(len(block) >= 2 for block in blocks):
                yield tuple(tuple(block) for block in blocks)
            return
        value = items[index]
        for block_index in range(len(blocks)):
            blocks[block_index].append(value)
            yield from visit(index + 1)
            blocks[block_index].pop()
        blocks.append([value])
        yield from visit(index + 1)
        blocks.pop()

    yield from visit(0)


def main() -> None:
    ax, ay, cx, cy, bx, by, dx, dy = variables = sp.symbols(
        "ax ay cx cy bx by dx dy"
    )

    # P=x+p_A X^A+p_C X^C and Q=y+q_B X^B+q_D X^D.
    # The labels are P_A,P_C,Q_B,Q_D,J_AB,J_AD,J_CB,J_CD.
    buckets = (
        (ax - 1, ay),
        (cx - 1, cy),
        (bx, by - 1),
        (dx, dy - 1),
        (ax + bx - 1, ay + by - 1),
        (ax + dx - 1, ay + dy - 1),
        (cx + bx - 1, cy + by - 1),
        (cx + dx - 1, cy + dy - 1),
    )

    partitions = list(set_partitions_without_singletons(tuple(range(8))))
    profile_counts = Counter(tuple(sorted(map(len, part))) for part in partitions)
    gate(len(partitions) == 715, "wrong no-singleton partition count")
    gate(
        profile_counts
        == Counter({(2, 3, 3): 280, (2, 2, 4): 210, (4, 4): 35,
                    (2, 6): 28, (3, 5): 56, (2, 2, 2, 2): 105, (8,): 1}),
        "wrong partition profile ledger",
    )

    exits: Counter[str] = Counter()
    semantic_rows: list[str] = []

    for part in partitions:
        equations = []
        for block in part:
            base = buckets[block[0]]
            for label in block[1:]:
                equations.extend(
                    (buckets[label][0] - base[0], buckets[label][1] - base[1])
                )

        matrix, rhs = sp.linear_eq_to_matrix(equations, variables)
        solution_set = sp.linsolve((matrix, rhs), variables)
        gate(solution_set is not sp.EmptySet, f"unexpected inconsistent partition {part}")
        solution = tuple(next(iter(solution_set)))

        A, C = solution[0:2], solution[2:4]
        B, D = solution[4:6], solution[6:8]

        # This priority is part of the frozen consequence ledger.  Each exit
        # contradicts a defining hypothesis of the chart.
        if all(sp.expand(A[k] - C[k]) == 0 for k in range(2)):
            reason = "equal_P_supports"
        elif all(sp.expand(B[k] - D[k]) == 0 for k in range(2)):
            reason = "equal_Q_supports"
        elif sp.expand(A[0] + A[1]) == 1:
            reason = "first_P_exponent_is_linear"
        elif sp.expand(C[0] + C[1]) == 1:
            reason = "second_P_exponent_is_linear"
        else:
            raise RuntimeError(f"unclassified partition {part}: {solution}")

        exits[reason] += 1
        semantic_rows.append(f"{part}:{reason}:{solution}")

    expected_exits = Counter(
        {
            "equal_P_supports": 676,
            "equal_Q_supports": 23,
            "first_P_exponent_is_linear": 8,
            "second_P_exponent_is_linear": 8,
        }
    )
    gate(exits == expected_exits, "wrong exit-reason ledger")

    # A positive Keller control immediately outside the two-by-two chart.
    x, y, alpha, delta = sp.symbols("x y alpha delta")
    P = x + alpha * y**2
    Q = y + delta * P**2
    jacobian = sp.expand(sp.diff(P, x) * sp.diff(Q, y) - sp.diff(P, y) * sp.diff(Q, x))
    gate(jacobian == 1, "triangular 1-by-3 positive control failed")
    gate(
        sp.expand(Q - (y + delta * (x + alpha * y**2) ** 2)) == 0,
        "positive control normal form failed",
    )

    # A fully transverse hostile has eight nonzero formal contributions and
    # visibly exposes a singleton bucket.
    hostile = {
        ax: 2,
        ay: 1,
        cx: 1,
        cy: 3,
        bx: 1,
        by: 2,
        dx: 3,
        dy: 1,
    }
    hostile_buckets = [tuple(sp.expand(z).subs(hostile) for z in bucket) for bucket in buckets]
    multiplicities = Counter(hostile_buckets)
    gate(any(value == 1 for value in multiplicities.values()), "hostile singleton missing")
    determinants = (
        ax * by - ay * bx,
        ax * dy - ay * dx,
        cx * by - cy * bx,
        cx * dy - cy * dx,
    )
    gate(
        all(sp.expand(det).subs(hostile) != 0 for det in determinants),
        "hostile is not fully transverse",
    )

    source = Path(__file__).read_text(encoding="utf-8")
    tree = ast.parse(source)
    gate(
        not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
        "inactive Python assert found",
    )

    semantic = "\n".join(semantic_rows).encode("utf-8")
    digest = hashlib.sha256(semantic).hexdigest()

    print("THM-3689 exact companion -- fully transverse 2x2 sparse-support gate")
    print("universe=no-singleton partitions of 8 labelled Jacobian contributions")
    print(f"partitions={len(partitions)}")
    print("profiles=" + ",".join(f"{profile}:{profile_counts[profile]}" for profile in sorted(profile_counts)))
    print("exits=" + ",".join(f"{key}:{exits[key]}" for key in sorted(exits)))
    print(f"semantic_sha256={digest}")
    print(f"CHECKS={CHECKS}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
