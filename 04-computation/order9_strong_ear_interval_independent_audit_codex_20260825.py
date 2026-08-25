#!/usr/bin/env python3
"""Independent audit of the order-nine strong-ear interval certificate.

This script does not import the generator.  It decodes every parent/child mask,
reconstructs the child from the parent and cut, tests strong connectivity by
two reachability searches, and recomputes H with a memoized last-vertex
recurrence.  Small brute-permutation controls calibrate the recurrence.
"""

from functools import lru_cache
from hashlib import sha256
from itertools import permutations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = (
    ROOT
    / "05-knowledge/results/order9_strong_ear_interval_certificates_codex_20260825.tsv"
)


def require(condition, label):
    if not condition:
        raise RuntimeError("CHECK FAILED: " + label)


def decode(code, n):
    adj = [0] * n
    bit = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (code >> bit) & 1:
                adj[i] |= 1 << j
            else:
                adj[j] |= 1 << i
            bit += 1
    require(code < (1 << bit), "code has no out-of-range bits")
    return tuple(adj)


def encode(adj):
    code = 0
    bit = 0
    for i in range(len(adj)):
        for j in range(i + 1, len(adj)):
            if (adj[i] >> j) & 1:
                code |= 1 << bit
            bit += 1
    return code


def extend(parent, cut):
    n = len(parent)
    child = list(parent) + [0]
    for i in range(n):
        if (cut >> i) & 1:
            child[n] |= 1 << i
        else:
            child[i] |= 1 << n
    return tuple(child)


def reverse(adj):
    n = len(adj)
    out = [0] * n
    for i in range(n):
        for j in range(n):
            if i != j and (adj[i] >> j) & 1:
                out[j] |= 1 << i
    return tuple(out)


def reaches_all(adj):
    full = (1 << len(adj)) - 1
    seen = 1
    stack = [0]
    while stack:
        v = stack.pop()
        unseen = adj[v] & (full ^ seen)
        while unseen:
            bit = unseen & -unseen
            unseen ^= bit
            seen |= bit
            stack.append(bit.bit_length() - 1)
    return seen == full


def is_strong(adj):
    return reaches_all(adj) and reaches_all(reverse(adj))


def h_count_recursive(adj):
    n = len(adj)
    full = (1 << n) - 1

    @lru_cache(None)
    def ending(mask, last):
        last_bit = 1 << last
        if mask == last_bit:
            return 1
        previous_mask = mask ^ last_bit
        return sum(
            ending(previous_mask, previous)
            for previous in range(n)
            if (previous_mask >> previous) & 1
            and (adj[previous] >> last) & 1
        )

    return sum(ending(full, last) for last in range(n))


def h_count_brute(adj):
    return sum(
        all((adj[path[i]] >> path[i + 1]) & 1 for i in range(len(adj) - 1))
        for path in permutations(range(len(adj)))
    )


def controls():
    transitive6 = tuple(
        sum(1 << j for j in range(i + 1, 6)) for i in range(6)
    )
    c3 = (1 << 1, 1 << 2, 1 << 0)
    require(h_count_recursive(transitive6) == h_count_brute(transitive6) == 1,
            "transitive brute control")
    require(not is_strong(transitive6), "transitive connectivity hostile")
    require(h_count_recursive(c3) == h_count_brute(c3) == 3,
            "C3 brute control")
    require(is_strong(c3), "C3 connectivity control")
    return "transitive6:H=1/non-strong;C3:H=3/strong"


def main():
    print("ORDER-NINE STRONG-EAR INTERVAL INDEPENDENT AUDIT")
    print("controls=" + controls())
    payload = CERTIFICATE.read_bytes()
    lines = payload.decode().splitlines()
    require(
        lines[0]
        == "H\tparent_code\tcut\tchild_code\tH_parent\tterminal\tinitial\tinternal",
        "certificate header",
    )
    expected = tuple(range(85, 2882, 2))
    require(len(lines) - 1 == len(expected) == 1399, "certificate row count")
    targets = []
    parent_values = set()
    for row_index, line in enumerate(lines[1:]):
        (
            target,
            parent_code,
            cut,
            child_code,
            claimed_parent_h,
            terminal,
            initial,
            internal,
        ) = map(int, line.split("\t"))
        require(target == expected[row_index], "ordered target interval")
        require(0 < cut < 255, "nonconstant ear cut")
        parent = decode(parent_code, 8)
        child = decode(child_code, 9)
        require(encode(parent) == parent_code, "parent decode/encode round trip")
        require(encode(child) == child_code, "child decode/encode round trip")
        require(extend(parent, cut) == child, "child reconstructs from parent and cut")
        require(is_strong(parent), "strong parent")
        require(is_strong(child), "strong child")
        parent_h = h_count_recursive(parent)
        child_h = h_count_recursive(child)
        require(parent_h == claimed_parent_h, "parent H certificate")
        require(child_h == target, "independent child H certificate")
        require(terminal + initial + internal == target, "ear components total")
        targets.append(target)
        parent_values.add(parent_h)

    require(tuple(targets) == expected, "all odd values in the interval")
    print("certificate_sha256=" + sha256(payload).hexdigest())
    print("verified_rows=1399")
    print("verified_interval=all odd H in [85,2881]")
    print("distinct_parent_H_values=" + str(len(parent_values)))
    print("parent_H_range=" + repr((min(parent_values), max(parent_values))))
    print("strong_parent_failures=0;strong_child_failures=0;H_failures=0")
    print("conclusion=INDEPENDENTLY VERIFIED FINITE-EXACT")


if __name__ == "__main__":
    main()
