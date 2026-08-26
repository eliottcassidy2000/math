#!/usr/bin/env python3
"""Independent small-order referee for THM-4177.

The input is a concatenated gentourng -q stream for parent orders 3 through 6.
This implementation deliberately uses Python integers, direct disjoint-edge
enumeration for D, and separate tuple/list endpoint tables.
"""

from __future__ import annotations

import sys
from collections import defaultdict
from itertools import combinations


def need(condition: bool, message: str) -> None:
    """Enforce an audit condition even when Python runs with -O."""
    if not condition:
        raise RuntimeError(message)


def parse_bits(bits: str) -> tuple[int, tuple[int, ...]]:
    edge_count = len(bits)
    n = 0
    while n * (n - 1) // 2 < edge_count:
        n += 1
    need(
        n * (n - 1) // 2 == edge_count and 3 <= n <= 6,
        f"invalid tournament bitstring length {edge_count}",
    )
    out = [0] * n
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if bits[k] == "1":
                out[i] |= 1 << j
            else:
                need(bits[k] == "0", f"invalid tournament bit {bits[k]!r}")
                out[j] |= 1 << i
            k += 1
    return n, tuple(out)


def extend(n: int, out: tuple[int, ...], z: int) -> tuple[int, tuple[int, ...]]:
    child = list(out) + [0]
    x = n
    for i in range(n):
        if (z >> i) & 1:
            child[x] |= 1 << i
        else:
            child[i] |= 1 << x
    return n + 1, tuple(child)


def has_sink(out: tuple[int, ...]) -> bool:
    return any(mask == 0 for mask in out)


def has_source(n: int, out: tuple[int, ...]) -> bool:
    return any(mask.bit_count() == n - 1 for mask in out)


def is_strong(n: int, out: tuple[int, ...]) -> bool:
    def visit(reverse: bool) -> int:
        seen = 1
        todo = [0]
        while todo:
            u = todo.pop()
            for v in range(n):
                arc = ((out[v] >> u) & 1) if reverse else ((out[u] >> v) & 1)
                if arc and not ((seen >> v) & 1):
                    seen |= 1 << v
                    todo.append(v)
        return seen

    return visit(False) == (1 << n) - 1 and visit(True) == (1 << n) - 1


def is_prime(n: int, out: tuple[int, ...]) -> bool:
    full = (1 << n) - 1
    for module in range(1, full):
        if module.bit_count() < 2:
            continue
        ok = True
        for v in range(n):
            if (module >> v) & 1:
                continue
            relation = out[v] & module
            if relation not in (0, module):
                ok = False
                break
        if ok:
            return False
    return True


def endpoint_data(n: int, out: tuple[int, ...]):
    size = 1 << n
    incoming = [0] * n
    for u in range(n):
        for v in range(n):
            if (out[u] >> v) & 1:
                incoming[v] |= 1 << u
    ending = [[0] * n for _ in range(size)]
    starting = [[0] * n for _ in range(size)]
    for v in range(n):
        ending[1 << v][v] = 1
        starting[1 << v][v] = 1
    for mask in range(1, size):
        for v in range(n):
            if not ((mask >> v) & 1):
                continue
            rest = mask ^ (1 << v)
            ending[mask][v] += sum(
                ending[rest][u]
                for u in range(n)
                if ((rest >> u) & 1) and ((incoming[v] >> u) & 1)
            )
            starting[mask][v] += sum(
                starting[rest][u]
                for u in range(n)
                if ((rest >> u) & 1) and ((out[v] >> u) & 1)
            )

    def before(mask: int, v: int) -> int:
        if mask == 0:
            return 1
        return sum(
            ending[mask][u]
            for u in range(n)
            if ((mask >> u) & 1) and ((out[u] >> v) & 1)
        )

    def after(mask: int, v: int) -> int:
        if mask == 0:
            return 1
        return sum(
            starting[mask][u]
            for u in range(n)
            if ((mask >> u) & 1) and ((out[v] >> u) & 1)
        )

    full = size - 1
    exposed = [[0] * n for _ in range(n)]
    capacity = [[0] * n for _ in range(n)]
    for x in range(n):
        for y in range(n):
            if x == y:
                continue
            rem = full ^ (1 << x) ^ (1 << y)
            left = rem
            while True:
                exposed[x][y] += before(left, x) * after(rem ^ left, y)
                if left == 0:
                    break
                left = (left - 1) & rem
    for x in range(n):
        for y in range(x + 1, n):
            capacity[x][y] = capacity[y][x] = exposed[x][y] + exposed[y][x]
    subset_h = [sum(row) for row in ending]
    subset_h[0] = 1
    return {
        "H": sum(ending[full]),
        "start": tuple(starting[full]),
        "end": tuple(ending[full]),
        "Q": tuple(tuple(row) for row in exposed),
        "c": tuple(tuple(row) for row in capacity),
        "subset_H": tuple(subset_h),
    }


def packet(n: int, out: tuple[int, ...], capacity) -> tuple[int, int]:
    degree = [0] * n
    signed = [0] * n
    edges = []
    for i in range(n):
        for j in range(i + 1, n):
            value = capacity[i][j]
            edges.append((i, j, value))
            degree[i] += value
            degree[j] += value
            if (out[i] >> j) & 1:
                signed[i] += value
                signed[j] -= value
            else:
                signed[i] -= value
                signed[j] += value
    C = sum(d * h for d, h in zip(degree, signed))
    D = 0
    for index, (i, j, value) in enumerate(edges):
        for u, v, other in edges[index + 1 :]:
            if len({i, j, u, v}) == 4:
                D += value * other
    return C, D


def root_split(n: int, out: tuple[int, ...], child_out: tuple[int, ...], cap):
    old = tuple(tuple(cap[i][j] if i < n and j < n else 0 for j in range(n)) for i in range(n))
    C, D = packet(n, out, old)
    degree = [0] * n
    signed = [0] * n
    W = 0
    for i in range(n):
        for j in range(i + 1, n):
            value = old[i][j]
            W += value
            degree[i] += value
            degree[j] += value
            if (out[i] >> j) & 1:
                signed[i] += value
                signed[j] -= value
            else:
                signed[i] -= value
                signed[j] += value
    A = 0
    S = 0
    for i in range(n):
        a = cap[i][n]
        sigma = 1 if ((child_out[i] >> n) & 1) else -1
        D += a * (W - degree[i])
        C += a * signed[i] + sigma * a * degree[i] + sigma * a * a
        A += a
        S += sigma * a
    C -= A * S
    return C, D


def odd_path_capacity(n: int, out: tuple[int, ...], subset_h) -> tuple[tuple[int, ...], ...]:
    full = (1 << n) - 1
    answer = [[0] * n for _ in range(n)]

    def walk(start: int, v: int, mask: int, length: int) -> None:
        if length and length % 2:
            i, j = sorted((start, v))
            answer[i][j] += 2 * subset_h[full ^ mask]
            answer[j][i] = answer[i][j]
        for u in range(n):
            if not ((mask >> u) & 1) and ((out[v] >> u) & 1):
                walk(start, u, mask | (1 << u), length + 1)

    for start in range(n):
        walk(start, start, 1 << start, 0)
    return tuple(tuple(row) for row in answer)


def audit_arbitrary_tensor_root_split() -> int:
    failures = 0
    checks = 0
    for n in (3, 4, 5, 6):
        out = [0] * (n + 1)
        for i in range(n + 1):
            for j in range(i + 1, n + 1):
                if (17 * i + 31 * j + n) % 3:
                    out[i] |= 1 << j
                else:
                    out[j] |= 1 << i
        cap = [[0] * (n + 1) for _ in range(n + 1)]
        for i in range(n + 1):
            for j in range(i + 1, n + 1):
                cap[i][j] = cap[j][i] = 1 + ((13 * i + 29 * j + 7 * n) % 23)
        direct = packet(n + 1, tuple(out), tuple(tuple(row) for row in cap))
        split = root_split(n, tuple(out[:n]), tuple(out), tuple(tuple(row) for row in cap))
        checks += 1
        failures += direct != split
    need(failures == 0, f"arbitrary-tensor root-split failures: {failures}")
    return checks


def audit_generic_tensor_hostile() -> tuple[int, int, int]:
    n, out = parse_bits("100111")
    need(
        n == 4 and not has_sink(out) and not has_source(n, out) and is_strong(n, out),
        "generic-tensor hostile carrier failed its structural controls",
    )
    values = (1, 1, 1, 1, 2, 3)  # 01,02,03,12,13,23
    capacity = [[0] * n for _ in range(n)]
    for value, (i, j) in zip(values, combinations(range(n), 2)):
        capacity[i][j] = capacity[j][i] = value
    C, D = packet(n, out, tuple(tuple(row) for row in capacity))
    need(
        (C, D, D + 2 * C) == (-4, 6, -2),
        f"generic-tensor hostile changed: {(C, D, D + 2 * C)}",
    )
    return C, D, D + 2 * C


def main() -> None:
    stats = {n: defaultdict(int) for n in range(3, 7)}
    arbitrary_checks = audit_arbitrary_tensor_root_split()
    hostile = audit_generic_tensor_hostile()
    for bits in sys.stdin.read().split():
        n, out = parse_bits(bits)
        parent = endpoint_data(n, out)
        s = stats[n]
        s["parents"] += 1
        s["strong_parents"] += is_strong(n, out)
        s["prime_parents"] += is_prime(n, out)
        s["sink_parents"] += has_sink(out)
        s["source_parents"] += has_source(n, out)
        path_cap = odd_path_capacity(n, out, parent["subset_H"])
        for i in range(n):
            for j in range(i + 1, n):
                s["path_checks"] += 1
                s["path_fail"] += path_cap[i][j] != parent["c"][i][j]
        plus = []
        minus = []
        for z in range(1 << n):
            child_n, child_out = extend(n, out, z)
            child = endpoint_data(child_n, child_out)
            C, D = packet(child_n, child_out, child["c"])
            gp, gm = D + 2 * C, D - 2 * C
            plus.append(gp)
            minus.append(gm)
            s["presentations"] += 1
            s["plus_neg"] += gp < 0
            s["plus_zero"] += gp == 0
            s["minus_neg"] += gm < 0
            s["minus_zero"] += gm == 0
            sink = has_sink(child_out)
            source = has_source(child_n, child_out)
            s["sink_children"] += sink
            s["source_children"] += source
            s["nosink_plus_nonpositive"] += (not sink) and gp <= 0
            s["nosource_minus_nonpositive"] += (not source) and gm <= 0
            s["strong_children"] += is_strong(child_n, child_out)
            s["prime_children"] += is_prime(child_n, child_out)
            predicted_split = root_split(n, out, child_out, child["c"])
            s["split_checks"] += 1
            s["split_fail"] += predicted_split != (C, D)
            for i in range(n):
                predicted = parent["start"][i] + parent["end"][i]
                predicted += sum(
                    parent["Q"][i][j] if ((z >> j) & 1) else parent["Q"][j][i]
                    for j in range(n)
                    if j != i
                )
                s["root_checks"] += 1
                s["root_fail"] += predicted != child["c"][i][n]
        full = (1 << n) - 1
        plus_best = min(plus[1:])
        minus_best = min(minus[:full])
        sink_parent = has_sink(out)
        source_parent = has_source(n, out)
        s["plus_exist_disagree"] += (plus_best <= 0) != sink_parent
        s["minus_exist_disagree"] += (minus_best <= 0) != source_parent
        s["plus_corner_disagree"] += (plus_best < plus[0]) != sink_parent
        s["minus_corner_disagree"] += (minus_best < minus[full]) != source_parent

    print("THM4177_INDEPENDENT_V1")
    print(f"arbitrary_tensor_root_split_checks={arbitrary_checks} failures=0")
    print(f"generic_tensor_hostile label=100111 weights=1,1,1,1,2,3 C={hostile[0]} D={hostile[1]} Gplus={hostile[2]}")
    for n in range(3, 7):
        s = stats[n]
        print(
            f"q={n} parents={s['parents']} presentations={s['presentations']} "
            f"strong_parents={s['strong_parents']} prime_parents={s['prime_parents']} "
            f"sink_parents={s['sink_parents']} source_parents={s['source_parents']}"
        )
        print(
            f"  algebra root_formula={s['root_checks']}/{s['root_fail']} "
            f"root_split={s['split_checks']}/{s['split_fail']} "
            f"path_edges={s['path_checks']}/{s['path_fail']}"
        )
        print(
            f"  intrinsic plus_neg={s['plus_neg']} plus_zero={s['plus_zero']} "
            f"sink_children={s['sink_children']} nosink_plus_nonpositive={s['nosink_plus_nonpositive']} "
            f"minus_neg={s['minus_neg']} minus_zero={s['minus_zero']} "
            f"source_children={s['source_children']} nosource_minus_nonpositive={s['nosource_minus_nonpositive']}"
        )
        print(
            f"  child_controls strong={s['strong_children']} prime={s['prime_children']}"
        )
        print(
            "  parent_equivalences "
            f"plus_nonzero_nonpositive_iff_sink_disagree={s['plus_exist_disagree']} "
            f"minus_nonfull_nonpositive_iff_source_disagree={s['minus_exist_disagree']} "
            f"plus_corner_min_iff_sink_disagree={s['plus_corner_disagree']} "
            f"minus_corner_min_iff_source_disagree={s['minus_corner_disagree']}"
        )


if __name__ == "__main__":
    main()
