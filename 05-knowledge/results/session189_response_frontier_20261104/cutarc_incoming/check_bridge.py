#!/usr/bin/env python3
"""Small independent exact check of the THM-4868 -> sigma-one ready bridge.

No imports from incoming code, no downloads, no third-party packages.
This is a bridge check through order 8, not a replay of incoming large censuses.
"""
import collections
import json
import math


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def induced(a, vertices):
    return tuple(sum(((a[u] >> v) & 1) << j for j, v in enumerate(vertices))
                 for u in vertices)


def strong(a):
    n = len(a)
    if n < 3:
        return False
    full = (1 << n) - 1
    for start in range(n):
        reached = 1 << start
        while True:
            nxt = reached
            for u in range(n):
                if reached >> u & 1:
                    nxt |= a[u]
            if nxt == reached:
                break
            reached = nxt
        if reached != full:
            return False
    return True


def endpoint_table(a, fixed_start=None):
    n = len(a)
    dp = [[0] * n for _ in range(1 << n)]
    for u in range(n):
        if fixed_start is None or u == fixed_start:
            dp[1 << u][u] = 1
    for mask in range(1 << n):
        for u, count in enumerate(dp[mask]):
            if count:
                left = a[u] & ~mask
                while left:
                    bit = left & -left
                    left -= bit
                    dp[mask | bit][bit.bit_length() - 1] += count
    return dp


def endpoint_counts(a, fixed_start=None):
    return endpoint_table(a, fixed_start)[-1]


def two_outneighbor_path_law(a, x, p):
    n = len(a)
    full = (1 << n) - 1
    end = endpoint_table(a)
    reverse = tuple(sum((1 << v) for v in range(n) if a[v] >> u & 1)
                    for u in range(n))
    start = endpoint_table(reverse)
    h = [sum(row) for row in end]
    h[0] = 1

    def f(z):
        return sum(h[mask] * start[full ^ mask][z]
                   for mask in range(1 << n) if not (mask >> z & 1))

    def g(q, z):
        return sum(end[mask][q] * start[full ^ mask][z]
                   for mask in range(1 << n)
                   if mask >> q & 1 and not (mask >> z & 1))

    require(f(x) == h[full] + 1, "digit-one cover identity")
    return h[full] - end[full][x] - end[full][p] + f(x) + f(p) - g(x, p) - g(p, x)


def hc(a):
    e = endpoint_counts(a, 0)
    return sum(e[u] for u in range(len(a)) if a[u] & 1)


def phi(a, exit_vertex):
    n = len(a)
    return tuple(sum((int(not (a[u] >> v & 1)) if
                      u != exit_vertex and v != exit_vertex else
                      (a[u] >> v & 1)) << v
                     for v in range(n) if u != v) for u in range(n))


def cut_count(a):
    answer = 0
    for u in range(len(a)):
        for v in range(len(a)):
            if a[u] >> v & 1:
                b = list(a)
                b[u] ^= 1 << v
                answer += not strong(tuple(b))
    return answer


def canonical_cycle_code(a, cycle):
    return min(induced(a, cycle[i:] + cycle[:i]) for i in range(len(a)))


def choose(n, k):
    return math.comb(n, k) if 0 <= k <= n else 0


def formula(n, k):
    return sum(choose(n - 4, k - 3 + 2*d) * choose(k - 2 + d, d)
               for d in range(n + 1))


def describe(a, cycle, word):
    n = len(a)
    nc = [u for u in range(n) if strong(induced(a, [v for v in range(n) if v != u]))]
    entries = [u for u in range(n) if n - 1 - a[u].bit_count() == 1]
    exits = [u for u in range(n) if a[u].bit_count() == 1]
    y = [u for u in entries if u in nc]
    x = [u for u in exits if u in nc]
    require(len(y) == len(x) == 1, "canonical endpoints")
    y, x = y[0], x[0]
    rdy = phi(a, x)
    remain = [u for u in range(n) if u != x]
    core = induced(rdy, remain)
    root = a[x].bit_length() - 1
    sigma = endpoint_counts(core)[remain.index(root)]
    k = cut_count(a)
    require(hc(a) == 1 and strong(rdy) and sigma == 1, "Phi image")
    require(k + len(nc) == n + 1, "run-cut dictionary")
    return {"n": n, "word": word, "uh_rows": a, "ready_rows": rdy,
            "cycle": cycle, "entry": y, "exit": x, "root": root,
            "run": len(nc), "cut_arcs": k, "entries": len(entries),
            "exits": len(exits), "H_UH": sum(endpoint_counts(a)),
            "H_ready": sum(endpoint_counts(rdy)), "sigma": sigma}


def main():
    levels = [((6, 12, 8, 1), (0, 1, 2, 3), "")]
    rows = []
    receipts = []
    earliest_same_state = None
    total_edges = 0
    for n in range(4, 9):
        data = [describe(*item) for item in levels]
        codes = {canonical_cycle_code(a, c) for a, c, _ in levels}
        require(len(codes) == len(levels), "no duplicated UH classes")
        hist = collections.Counter(d["cut_arcs"] for d in data)
        require(all(hist[k] == formula(n, k) for k in range(2, n)), "distribution")
        same_state = {}
        for d in data:
            state = (d["run"], d["entries"], d["exits"])
            if state in same_state and same_state[state]["H_ready"] != d["H_ready"]:
                if earliest_same_state is None:
                    earliest_same_state = [same_state[state], d]
            else:
                same_state[state] = d
        receipts.append({"n": n, "classes": len(levels), "cut_histogram": dict(sorted(hist.items())),
                         "ready_H_min": min(d["H_ready"] for d in data),
                         "ready_H_max": max(d["H_ready"] for d in data)})
        rows.extend(data)
        if n == 8:
            break
        next_level = []
        for (a, cycle, word), d in zip(levels, data):
            yi = cycle.index(d["entry"])
            parents = [cycle[yi], cycle[(yi - 1) % n]]
            if n - 1 - a[parents[1]].bit_count() == 1:
                parents.append(cycle[(yi - 2) % n])
            exact_admissible = [u for u, e in enumerate(endpoint_counts(a)) if e == 1]
            require(set(parents) == set(exact_admissible), "admissible parents")
            for j, p in enumerate(parents):
                require(p != d["exit"], "parent is not the canonical exit")
                child = list(a)
                child[p] |= 1 << n
                child.append(((1 << n) - 1) ^ (1 << p))
                child = tuple(child)
                pi = cycle.index(p)
                cc = cycle[:pi + 1] + (n,) + cycle[pi + 1:]
                cd = describe(child, cc, word + str(j))
                require(cd["exit"] == d["exit"], "exit persists")
                require(cd["run"] == d["run"] + j, "run transition")
                require(cd["ready_rows"][n] == (1 << p) | (1 << d["exit"]), "two-outneighbor child")
                require(induced(cd["ready_rows"], list(range(n))) == d["ready_rows"], "ready induced parent")
                require(cd["H_ready"] == two_outneighbor_path_law(d["ready_rows"], d["exit"], p),
                        "two-outneighbor path count update")
                next_level.append((child, cc, word + str(j)))
                total_edges += 1
        levels = next_level
    output = {"scope": "Independent bridge check, generated UH classes at orders 4-8 only",
              "classes_checked": len(rows), "insertion_edges_checked": total_edges,
              "receipts": receipts, "same_finite_state_different_ready_H": earliest_same_state}
    witness_inputs = [((6, 60, 8, 1, 13, 29), (0, 1, 5, 4, 2, 3),
                       (60, 5, 8, 50, 38, 6), 33),
                      ((22, 12, 8, 1, 46, 15), (0, 4, 5, 1, 2, 3),
                       (44, 53, 8, 50, 5, 20), 29)]
    for uh, cycle, ready, expected_h in witness_inputs:
        d = describe(uh, cycle, "explicit witness")
        require(d["ready_rows"] == ready and d["H_ready"] == expected_h,
                "explicit ready witness")
        require((d["n"], d["run"], d["cut_arcs"], d["entries"], d["exits"], d["H_UH"],
                 d["exit"], d["root"], d["sigma"]) == (6, 3, 4, 2, 2, 15, 2, 3, 1),
                "explicit common witness state")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
