#!/usr/bin/env python3
"""Exact finite controls for THM-2330.

This is not a knot enumerator.  It checks the theorem's two conical
word-metric hostile controls, the partition-lattice potential identities,
and one finite min-plus convolution bank.
"""

from __future__ import annotations

from functools import lru_cache
from itertools import combinations, product


Vector = tuple[int, ...]
Block = frozenset[int]
Partition = tuple[Block, ...]


def add(a: Vector, b: Vector) -> Vector:
    return tuple(x + y for x, y in zip(a, b))


def sub(a: Vector, b: Vector) -> Vector:
    return tuple(x - y for x, y in zip(a, b))


def scale(k: int, a: Vector) -> Vector:
    return tuple(k * x for x in a)


def basis(r: int, i: int) -> Vector:
    return tuple(int(j == i) for j in range(r))


def l1(z: Vector) -> int:
    return sum(abs(x) for x in z)


def shortcut_norm(z: Vector) -> int:
    """Generators +/-e_i cost 1 and +/-(1,...,1) cost r-1."""
    r = len(z)
    # A minimizer lies between the coordinates and zero; this wider exact
    # bank is deliberately redundant.
    radius = max([abs(x) for x in z] + [1]) + 2
    return min(
        abs(k) * (r - 1) + sum(abs(x - k) for x in z)
        for k in range(-radius, radius + 1)
    )


def target_shortcut_norm(z: Vector) -> int:
    """Generators +/-e_i and +/-(1,1,-1), all of cost 1."""
    assert len(z) == 3
    v = (1, 1, -1)
    radius = max([abs(x) for x in z] + [1]) + 2
    return min(
        abs(k) + l1(sub(z, scale(k, v)))
        for k in range(-radius, radius + 1)
    )


def canonical_partition(blocks: list[set[int] | Block]) -> Partition:
    frozen = [frozenset(block) for block in blocks]
    return tuple(sorted(frozen, key=lambda b: (min(b), len(b), tuple(b))))


@lru_cache(maxsize=None)
def set_partitions(items: tuple[int, ...]) -> tuple[Partition, ...]:
    if not items:
        return ((),)
    first, *rest = items
    out: set[Partition] = set()
    for pi in set_partitions(tuple(rest)):
        out.add(canonical_partition([{first}, *pi]))
        for j in range(len(pi)):
            blocks = list(pi)
            blocks[j] = blocks[j] | {first}
            out.add(canonical_partition(blocks))
    return tuple(sorted(out, key=lambda p: (-len(p), tuple(tuple(b) for b in p))))


def refines(pi: Partition, rho: Partition) -> bool:
    return all(any(block <= coarse for coarse in rho) for block in pi)


def covers(pi: Partition, rho: Partition) -> bool:
    return refines(pi, rho) and len(pi) == len(rho) + 1


def block_vector(block: Block, r: int) -> Vector:
    ans = (0,) * r
    for i in block:
        ans = add(ans, basis(r, i))
    return ans


def weak_compositions(n: int, parts: int) -> tuple[tuple[int, ...], ...]:
    if parts == 1:
        return ((n,),)
    out = []
    for first in range(n + 1):
        for tail in weak_compositions(n - first, parts - 1):
            out.append((first, *tail))
    return tuple(out)


def vector_factorizations(target: Vector, parts: int) -> tuple[tuple[Vector, ...], ...]:
    per_coordinate = [weak_compositions(n, parts) for n in target]
    out = []
    for coordinate_splits in product(*per_coordinate):
        vectors = tuple(
            tuple(coordinate_splits[j][p] for j in range(len(target)))
            for p in range(parts)
        )
        out.append(vectors)
    return tuple(out)


def lift_cost(
    starts: tuple[Vector, ...],
    target: Vector,
    norm,
) -> int:
    return min(
        sum(norm(sub(start, end)) for start, end in zip(starts, ends))
        for ends in vector_factorizations(target, len(starts))
    )


def partition_lift_cost(r: int, pi: Partition, target: Vector, norm) -> int:
    starts = tuple(block_vector(block, r) for block in pi)
    return lift_cost(starts, target, norm)


def test_target_hostile() -> dict[str, int]:
    e1, e2, e3 = (basis(3, i) for i in range(3))
    starts = (e1, e2)
    total = add(e1, e2)
    root_l1 = lift_cost(starts, (0, 0, 0), l1)
    root_short = lift_cost(starts, (0, 0, 0), target_shortcut_norm)
    target_l1 = lift_cost(starts, e3, l1)
    target_short = lift_cost(starts, e3, target_shortcut_norm)
    base_l1 = l1(sub(total, e3))
    base_short = target_shortcut_norm(sub(total, e3))

    assert [l1(e1), l1(e2), l1(total)] == [1, 1, 2]
    assert [
        target_shortcut_norm(e1),
        target_shortcut_norm(e2),
        target_shortcut_norm(total),
    ] == [1, 1, 2]
    assert (root_l1, root_short) == (2, 2)
    assert (target_l1, target_short) == (3, 3)
    assert (base_l1, base_short) == (3, 1)
    assert (target_l1 - base_l1, target_short - base_short) == (0, 2)
    return {
        "root_lambda_l1": root_l1,
        "root_lambda_shortcut": root_short,
        "target_lambda_l1": target_l1,
        "target_lambda_shortcut": target_short,
        "target_omega_l1": target_l1 - base_l1,
        "target_omega_shortcut": target_short - base_short,
    }


def test_arity_and_partition_lattice(r: int) -> dict[str, object]:
    items = tuple(range(r))
    partitions = set_partitions(items)
    zero = (0,) * r
    g = (1,) * r

    for s in range(r + 1):
        for subset in combinations(items, s):
            z = tuple(int(i in subset) for i in items)
            expected = s if s < r else r - 1
            assert shortcut_norm(z) == expected
            assert l1(z) == s

    lam_l1 = {
        pi: partition_lift_cost(r, pi, zero, l1)
        for pi in partitions
    }
    lam_short = {
        pi: partition_lift_cost(r, pi, zero, shortcut_norm)
        for pi in partitions
    }
    onehat = next(pi for pi in partitions if len(pi) == 1)
    assert all(value == r for value in lam_l1.values())
    assert lam_short[onehat] == r - 1
    assert all(
        value == r
        for pi, value in lam_short.items()
        if pi != onehat
    )

    cover_drops = []
    for pi in partitions:
        for rho in partitions:
            if covers(pi, rho):
                assert lam_l1[rho] <= lam_l1[pi]
                assert lam_short[rho] <= lam_short[pi]
                cover_drops.append(lam_short[pi] - lam_short[rho])

    # Every comparable triple satisfies the exact potential/coboundary law.
    triples = 0
    for pi in partitions:
        for rho in partitions:
            if not refines(pi, rho):
                continue
            for tau in partitions:
                if refines(rho, tau):
                    lhs = lam_short[pi] - lam_short[tau]
                    rhs = (lam_short[pi] - lam_short[rho]) + (
                        lam_short[rho] - lam_short[tau]
                    )
                    assert lhs == rhs
                    triples += 1

    # The prescribed start is already a factorization of g.
    assert all(
        partition_lift_cost(r, pi, g, shortcut_norm) == 0
        for pi in partitions
    )

    return {
        "r": r,
        "partition_count": len(partitions),
        "cover_count": len(cover_drops),
        "positive_cover_drops": sum(drop > 0 for drop in cover_drops),
        "max_cover_drop": max(cover_drops),
        "comparable_triples": triples,
        "proper_arity_agreement": r - 1,
        "terminal_root_omega_l1": lam_l1[next(pi for pi in partitions if len(pi) == r)]
        - l1(g),
        "terminal_root_omega_shortcut": lam_short[
            next(pi for pi in partitions if len(pi) == r)
        ]
        - shortcut_norm(g),
    }


def test_convolution() -> dict[str, int]:
    r = 3
    x_starts = (basis(r, 0), basis(r, 1))
    z_starts = (basis(r, 2),)
    all_starts = x_starts + z_starts
    target_bank = tuple(product(range(3), repeat=r))
    checks = 0

    for target in target_bank:
        lhs = lift_cost(all_starts, target, shortcut_norm)
        rhs = min(
            lift_cost(x_starts, a, shortcut_norm)
            + lift_cost(z_starts, sub(target, a), shortcut_norm)
            for a in product(*(range(n + 1) for n in target))
        )
        assert lhs == rhs
        checks += 1

    return {
        "target_count": len(target_bank),
        "convolution_checks": checks,
    }


def main() -> None:
    print("THM-2330 exact conical word-metric controls")
    print()
    print("target-conditioned hostile")
    for key, value in test_target_hostile().items():
        print(f"{key}={value}")

    print()
    print("arbitrary-arity partition controls")
    for r in range(3, 7):
        row = test_arity_and_partition_lattice(r)
        print(" ".join(f"{key}={value}" for key, value in row.items()))

    print()
    print("disjoint-packet min-plus convolution")
    for key, value in test_convolution().items():
        print(f"{key}={value}")

    print()
    print("all exact controls passed")


if __name__ == "__main__":
    main()
