#!/usr/bin/env python3
"""Exact, deterministic controls for the binary endpoint/cycle construction."""
import hashlib
import itertools
import json


def require(value, message):
    if not value:
        raise RuntimeError(message)


def build(h, m=None, forward=False):
    if m is None:
        m = (h - 1).bit_length()
    require(m >= 0 and 1 <= h <= 1 << m, ("invalid target", h, m))
    n, z, x = m + 4, m + 2, m + 3
    rows = [0] * n
    def arc(i, j):
        rows[i] |= 1 << j
    for i in range(z):
        for j in range(i + 1, z):
            arc(i, j)
    arc(z, 0)
    arc(1, z)
    for j in range(2, z):
        if ((h - 1) >> (j - 2)) & 1:
            arc(j, z)
        else:
            arc(z, j)
    arc(x, 1)
    for i in range(n - 1):
        if i != 1:
            arc(i, x)
    if forward:
        rows[z] ^= 1 << x
        rows[x] ^= 1 << z
    return tuple(rows)


def check_tournament(rows):
    n = len(rows)
    for i in range(n):
        require(rows[i] >= 0 and rows[i] < (1 << n) and not rows[i] >> i & 1,
                ("invalid row", i, rows))
        for j in range(i):
            require(((rows[i] >> j) ^ (rows[j] >> i)) & 1,
                    ("not a tournament", i, j, rows))


def endpoint_dp(rows, root):
    n = len(rows)
    dp = [[0] * n for _ in range(1 << n)]
    dp[1 << root][root] = 1
    for mask in range(1, 1 << n):
        for v, ways in enumerate(dp[mask]):
            if not ways:
                continue
            rem = rows[v] & ~mask
            while rem:
                bit = rem & -rem
                w = bit.bit_length() - 1
                dp[mask | bit][w] += ways
                rem -= bit
    return tuple(dp[-1])


def cycle_dp(rows, root=0):
    if len(rows) < 3:
        return 0
    ends = endpoint_dp(rows, root)
    return sum(ways for v, ways in enumerate(ends) if rows[v] >> root & 1)


def total_paths(rows):
    return sum(sum(endpoint_dp(rows, root)) for root in range(len(rows)))


def permutation_cycles(rows, root):
    # Independent full permutation universe, with rotation fixed at root.
    count = 0
    words = []
    for tail in itertools.permutations(v for v in range(len(rows)) if v != root):
        word = (root,) + tail
        if all((rows[a] >> b) & 1 for a, b in zip(word, word[1:] + (root,))):
            count += 1
            words.append(word)
    return count, set(words)


def formula_cycles(h, m):
    z, x = m + 2, m + 3
    words = set()
    for mask in range(1 << m):
        if mask and not ((h - 1) >> (mask.bit_length() - 1)) & 1:
            continue
        before = tuple(j for j in range(2, m + 2) if mask >> (j - 2) & 1)
        after = tuple(j for j in range(2, m + 2) if not mask >> (j - 2) & 1)
        words.add((x, 1) + before + (z, 0) + after)
    return words


def strong(rows):
    for r in range(len(rows)):
        visited, todo = 0, 1 << r
        while todo:
            bit = todo & -todo
            todo -= bit
            if visited & bit:
                continue
            visited |= bit
            todo |= rows[bit.bit_length() - 1] & ~visited
        if visited != (1 << len(rows)) - 1:
            return False
    return True


def gates(rows):
    n = len(rows)
    return [(a, b) for a in range(n) for b in range(n)
            if rows[a] >> b & 1 and rows[a].bit_count() == 1
            and n - 1 - rows[b].bit_count() == 1]


def run():
    checks = []
    permutation_objects = 0
    endpoint_objects = 0
    for m in range(9):
        for h in range(1, (1 << m) + 1):
            base = build(h, m)
            forward = build(h, m, True)
            for rows in (base, forward):
                check_tournament(rows)
                require(strong(rows), ("strong", m, h, rows))
            z, x = m + 2, m + 3
            core = tuple(row & ((1 << x) - 1) for row in base[:-1])
            require(sum(endpoint_dp(core, 1)) == h, ("endpoint", m, h))
            require(strong(core) == (m == 0 or h > (1 << (m - 1))),
                    ("strong core criterion", m, h))
            require(endpoint_dp(base, 1)[x] == h, ("fixed endpoints", m, h))
            require(cycle_dp(base) == h, ("cycle", m, h))
            expected_forward = h + (m > 0)
            require(cycle_dp(forward) == expected_forward, ("forward", m, h))
            if m:
                require(not gates(base), ("unexpected base gate", m, h))
                require(gates(forward) == ([(z, 0)] if h == (1 << m) else []),
                        ("forward gate classification", m, h))
            require(len(formula_cycles(h, m)) == h, ("word family", m, h))
            if m <= 4:
                direct, words = permutation_cycles(base, x)
                require(direct == h and words == formula_cycles(h, m),
                        ("permutation comparison", m, h))
                direct2, words2 = permutation_cycles(forward, x)
                require(direct2 == expected_forward, ("forward permutations", m, h))
                if m:
                    extra = (x, z) + tuple(range(m + 2))
                    require(words2 == words | {extra}, ("forward bijection", m, h))
                else:
                    require(words2 == words, ("m=0 boundary",))
                permutation_objects += 2
            endpoint_objects += 1
            checks.append([m, h, list(base), h, list(forward), expected_forward])
    extra_checks = []
    for h in [511, 512, 513, 1023, 1024]:
        rows = build(h)
        require(cycle_dp(rows) == h, ("large control", h))
        extra_checks.append({"h": h, "order": len(rows), "HC": h})
    witnesses = []
    expected_H = {1: 5, 2: 13, 7: 91, 21: 587}
    for h in (1, 2, 7, 21):
        rows = build(h)
        hp = total_paths(rows)
        require(hp == expected_H[h], ("H witness", h, hp))
        witnesses.append({"h": h, "order": len(rows), "rows": list(rows),
                          "HC": cycle_dp(rows), "H": hp,
                          "gate_arcs": gates(rows)})
    require(cycle_dp((2, 0)) == 0 and sum(endpoint_dp((0,), 0)) == 1,
            "one-vertex closure boundary")
    require(cycle_dp((6, 4, 0)) == 0 and not strong((6, 4, 0)), "zero control")
    payload = json.dumps(checks, separators=(",", ":"), sort_keys=True).encode()
    return {"status": "PASS", "binary_widths": [0, 8],
            "endpoint_objects": endpoint_objects, "cycle_objects": 2 * endpoint_objects,
            "full_permutation_cycle_objects": permutation_objects,
            "finite_universe": "Every binary word of every width 0..8; both outer-arc directions",
            "all_case_sha256": hashlib.sha256(payload).hexdigest(),
            "large_controls": extra_checks, "witnesses": witnesses,
            "boundaries": {"n=1_endpoint_closure": "must be excluded",
                           "m=0_forward_increment": 0,
                           "m>=1_forward_increment": 1,
                           "nongated_base": "all m>=1",
                           "strong_core_iff": "m=0 or h>2^(m-1)",
                           "fixed_endpoints": "start1 to x counts h",
                           "forward_gate_iff": "h=2^m for m>=1"}}


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, sort_keys=True))
