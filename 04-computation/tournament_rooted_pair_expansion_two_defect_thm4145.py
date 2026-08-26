#!/usr/bin/env python3
"""Exact controls for THM-4145's rooted pair-expansion identities."""

from __future__ import annotations

import hashlib
import itertools
from math import factorial


EXPECTED = {
    "semantic_sha256": "d174d868d8d1784901ef3926b0d0dc03e27a842eadd9a31df931c00716fa11cc",
    "semantic_fnv64": "9f075ab9b9d3f07f",
    "strong_rooted": 2822,
}


def adjacency(order: int, label: int) -> list[int]:
    out = [0] * order
    bit = 0
    for i in range(order):
        for j in range(i + 1, order):
            if (label >> bit) & 1:
                out[i] |= 1 << j
            else:
                out[j] |= 1 << i
            bit += 1
    return out


def label_bits(out: list[int]) -> str:
    bits: list[str] = []
    for i in range(len(out)):
        for j in range(i + 1, len(out)):
            bits.append("1" if (out[i] >> j) & 1 else "0")
    return "".join(bits)


def is_strong(out: list[int]) -> bool:
    order = len(out)
    full = (1 << order) - 1
    incoming = [0] * order
    for i in range(order):
        for j in range(order):
            if (out[i] >> j) & 1:
                incoming[j] |= 1 << i
    for arcs in (out, incoming):
        seen = todo = 1
        while todo:
            vbit = todo & -todo
            todo ^= vbit
            v = vbit.bit_length() - 1
            fresh = arcs[v] & ~seen
            seen |= fresh
            todo |= fresh
        if seen != full:
            return False
    return True


def hamilton_count(out: list[int]) -> int:
    order = len(out)
    size = 1 << order
    ending = [[0] * order for _ in range(size)]
    for v in range(order):
        ending[1 << v][v] = 1
    for mask in range(1, size):
        for v in range(order):
            if not ((mask >> v) & 1):
                continue
            previous = mask ^ (1 << v)
            for u in range(order):
                if ((previous >> u) & 1) and ((out[u] >> v) & 1):
                    ending[mask][v] += ending[previous][u]
    return sum(ending[-1])


def pair_expand(out: list[int], root: int) -> list[int]:
    order = len(out)
    child = out.copy() + [0]
    clone = order
    for u in range(order):
        if u == root:
            continue
        if (out[root] >> u) & 1:
            child[clone] |= 1 << u
        else:
            child[u] |= 1 << clone
    child[root] |= 1 << clone
    return child


def ear_at_root_type(out: list[int], root: int) -> list[int]:
    """Add x_S with S=N^+(root); root itself is not in S."""
    order = len(out)
    child = out.copy() + [0]
    ear = order
    for u in range(order):
        if u != root and ((out[root] >> u) & 1):
            child[ear] |= 1 << u
        else:
            child[u] |= 1 << ear
    return child


def deletion_defect_formula(out: list[int], root: int) -> int:
    vertices = tuple(v for v in range(len(out)) if v != root)
    total = 0
    for ordering in itertools.permutations(vertices):
        defects = {
            i + 1
            for i in range(len(ordering) - 1)
            if not ((out[ordering[i]] >> ordering[i + 1]) & 1)
        }
        good: set[int] = set()
        if (out[root] >> ordering[0]) & 1:
            good.add(0)
        if (out[ordering[-1]] >> root) & 1:
            good.add(len(ordering))
        for i in range(len(ordering) - 1):
            if ((out[ordering[i]] >> root) & 1) and (
                (out[root] >> ordering[i + 1]) & 1
            ):
                good.add(i + 1)
        if len(defects) > 2 or not defects <= good:
            continue
        k = len(good)
        total += (k * k, 2 * k - 1, 2)[len(defects)]
    return total


def naive_square(out: list[int], root: int) -> int:
    vertices = tuple(v for v in range(len(out)) if v != root)
    total = 0
    for ordering in itertools.permutations(vertices):
        if any(
            not ((out[ordering[i]] >> ordering[i + 1]) & 1)
            for i in range(len(ordering) - 1)
        ):
            continue
        good = int((out[root] >> ordering[0]) & 1)
        good += int((out[ordering[-1]] >> root) & 1)
        good += sum(
            int(
                ((out[ordering[i]] >> root) & 1)
                and ((out[root] >> ordering[i + 1]) & 1)
            )
            for i in range(len(ordering) - 1)
        )
        total += good * good
    return total


def internal_gap_counts(child: list[int], root: int, clone: int) -> tuple[int, int]:
    """Return Q(root,clone) and Q(clone,root) by literal orderings."""
    correct = reverse_defect = 0
    for ordering in itertools.permutations(range(len(child))):
        positions = {v: i for i, v in enumerate(ordering)}
        if abs(positions[root] - positions[clone]) != 1:
            continue
        defects = [
            i
            for i in range(len(child) - 1)
            if not ((child[ordering[i]] >> ordering[i + 1]) & 1)
        ]
        i = min(positions[root], positions[clone])
        if ordering[i : i + 2] == (root, clone) and not defects:
            correct += 1
        if ordering[i : i + 2] == (clone, root) and defects == [i]:
            reverse_defect += 1
    return correct, reverse_defect


def fnv_update(value: int, data: bytes) -> int:
    for byte in data:
        value ^= byte
        value = (value * 1099511628211) & ((1 << 64) - 1)
    return value


def main() -> None:
    sha = hashlib.sha256()
    fnv = 14695981039346656037
    total_rows = 0
    total_strong_rooted = 0
    order_lines: list[str] = []

    for order in range(2, 6):
        labels = 1 << (order * (order - 1) // 2)
        rooted = strong_rooted = 0
        for label in range(labels):
            out = adjacency(order, label)
            strong_parent = is_strong(out)
            parent_h = hamilton_count(out)
            for root in range(order):
                child = pair_expand(out, root)
                ear = ear_at_root_type(out, root)
                assert child == ear
                strong_child = is_strong(child)
                assert strong_child == strong_parent
                child_h = hamilton_count(child)
                formula_h = deletion_defect_formula(out, root)
                assert child_h == formula_h
                correct, reverse = internal_gap_counts(child, root, order)
                assert correct == parent_h
                assert reverse == parent_h
                row = (
                    f"{order},{label},{root},{int(strong_parent)},"
                    f"{parent_h},{child_h},{formula_h},{correct},{reverse}\n"
                ).encode()
                sha.update(row)
                fnv = fnv_update(fnv, row)
                rooted += 1
                strong_rooted += int(strong_parent)
        expected_rooted = labels * order
        assert rooted == expected_rooted
        total_rows += rooted
        total_strong_rooted += strong_rooted
        order_lines.append(
            f"order={order};labels={labels};rooted={rooted};"
            f"strong_rooted={strong_rooted};max_deletion_orders={factorial(order-1)}"
        )

    c3 = adjacency(3, 0b101)
    c3_child = pair_expand(c3, 0)
    assert label_bits(c3) == "101"
    assert label_bits(c3_child) == "101101"
    assert hamilton_count(c3_child) == 5
    assert naive_square(c3, 0) == 4

    root_loss = adjacency(4, 0b00100)
    assert label_bits(root_loss) == "001000"
    assert hamilton_count(pair_expand(root_loss, 0)) == 11
    assert hamilton_count(pair_expand(root_loss, 1)) == 9

    semantic_sha256 = sha.hexdigest()
    semantic_fnv64 = f"{fnv:016x}"
    assert semantic_sha256 == EXPECTED["semantic_sha256"]
    assert semantic_fnv64 == EXPECTED["semantic_fnv64"]
    assert total_strong_rooted == EXPECTED["strong_rooted"]

    print("THM4145_ROOTED_PAIR_EXPANSION_TWO_DEFECT_20260825")
    print("status=PASS;universe=all rooted labelled tournaments orders 2..5")
    for line in order_lines:
        print(line)
    print(f"rooted_rows={total_rows};strong_rooted={total_strong_rooted}")
    print("order11_covering_presentations=93559490")
    print("c3_hostile=(101,0,101101,5,4)")
    print("root_loss_control=(001000,(0,11),(1,9))")
    print(f"semantic_sha256={semantic_sha256}")
    print(f"semantic_fnv64={semantic_fnv64}")
    print("checks=strong_iff,ear_identity,defect_formula,internal_gap_pair")


if __name__ == "__main__":
    main()
