"""Exact 223/233 local tournament and path-contraction probes.

Run with python or python -O; no external dependencies and no disabled asserts.
"""
from collections import Counter
from itertools import combinations, permutations, product
import json


def check(value, message):
    if not value:
        raise RuntimeError(message)


def emit(label, obj):
    print(label + "=" + json.dumps(obj, sort_keys=True))


def hp(adj):
    n = len(adj)
    if n == 0:
        return 1
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1, 1 << n):
        for i, count in enumerate(dp[mask]):
            if not count:
                continue
            rem = adj[i] & ~mask
            while rem:
                bit = rem & -rem
                rem -= bit
                dp[mask | bit][bit.bit_length() - 1] += count
    return sum(dp[-1])


def flip(adj, edges):
    out = list(adj)
    for i, j in edges:
        out[i] ^= 1 << j
        out[j] ^= 1 << i
    return out


def staircase(k):
    n = 2 * k
    ranks = [i // 2 + (k if i % 2 else 0) for i in range(n)]
    out = [0] * n
    for i, j in combinations(range(n), 2):
        u, v = (j, i) if i // 2 == j // 2 else ((i, j) if ranks[i] < ranks[j] else (j, i))
        out[u] |= 1 << v
    return out


def contract(adj, block):
    """The block becomes vertex zero: incoming from its first, outgoing from last."""
    outside = [i for i in range(len(adj)) if i not in block]
    out = [0] * (len(outside) + 1)
    for j, v in enumerate(outside, 1):
        if adj[block[-1]] & (1 << v):
            out[0] |= 1 << j
        if adj[v] & (1 << block[0]):
            out[j] |= 1
        for k, w in enumerate(outside, 1):
            if adj[v] & (1 << w):
                out[j] |= 1 << k
    return out


def brute_paths(adj):
    return [p for p in permutations(range(len(adj)))
            if all(adj[u] & (1 << v) for u, v in zip(p, p[1:]))]


def count_block(paths, block):
    return sum(any(p[i:i + len(block)] == tuple(block)
                   for i in range(len(p) - len(block) + 1)) for p in paths)


def odd_cycles(adj):
    found = []
    for start in range(len(adj)):
        def visit(path, mask):
            last = path[-1]
            if len(path) >= 3 and len(path) % 2 and adj[last] & (1 << start):
                found.append((tuple(path), mask))
            for v in range(start + 1, len(adj)):
                if not mask & (1 << v) and adj[last] & (1 << v):
                    visit(path + [v], mask | (1 << v))
        visit([start], 1 << start)
    return found


def packet(adj):
    cycles = odd_cycles(adj)
    lengths = Counter(len(c) for c, _ in cycles)
    pairs = Counter(tuple(sorted((len(c), len(d))))
                    for (c, x), (d, y) in combinations(cycles, 2) if not x & y)
    # Used only for <=8 nontrivial vertices, so three disjoint cycles are impossible.
    check(hp(adj) == 1 + 2 * len(cycles) + 4 * sum(pairs.values()), "independent OCF packet")
    return {"odd_cycle_counts": dict(sorted(lengths.items())),
            "disjoint_odd_pairs": {str(k): v for k, v in sorted(pairs.items())}}


def ordered_padding(adj, before=1, after=1):
    n, m = len(adj), before + len(adj) + after
    result = [0] * m
    for i in range(m):
        for j in range(i + 1, m):
            if before <= i < before + n and before <= j < before + n:
                u, v = (i, j) if adj[i - before] & (1 << (j - before)) else (j, i)
            else:
                u, v = i, j
            result[u] |= 1 << v
    return result


def scc(adj):
    reach = [row | (1 << i) for i, row in enumerate(adj)]
    for k in range(len(adj)):
        for i in range(len(adj)):
            if reach[i] & (1 << k):
                reach[i] |= reach[k]
    unused, blocks = set(range(len(adj))), []
    while unused:
        i = min(unused)
        block = [j for j in sorted(unused) if reach[i] & (1 << j) and reach[j] & (1 << i)]
        blocks.append(block)
        unused.difference_update(block)
    return blocks


def mixed_response_audit():
    rows = []
    for n in range(3, 6):
        edges = list(combinations(range(n), 2))
        cache = []
        graphs = []
        for mask in range(1 << len(edges)):
            adj = [0] * n
            for bit, (i, j) in enumerate(edges):
                u, v = (j, i) if mask & (1 << bit) else (i, j)
                adj[u] |= 1 << v
            graphs.append(adj)
            cache.append(hp(adj))
        index = {edge: i for i, edge in enumerate(edges)}
        tests = 0
        for mask, adj in enumerate(graphs):
            for v in range(n):
                for u, w in combinations([i for i in range(n) if i != v], 2):
                    # Choose the orientation in which this is a directed two-edge path.
                    if adj[u] & (1 << v) and adj[v] & (1 << w):
                        first, last = u, w
                    elif adj[w] & (1 << v) and adj[v] & (1 << u):
                        first, last = w, u
                    else:
                        continue
                    e = 1 << index[tuple(sorted((u, v)))]
                    f = 1 << index[tuple(sorted((v, w)))]
                    actual = cache[mask ^ e ^ f] - cache[mask ^ e] - cache[mask ^ f] + cache[mask]
                    predicted = hp(contract(adj, (first, v, last))) + hp(contract(adj, (last, v, first)))
                    check(actual == predicted and actual >= 0, "conditional two-edge path curvature")
                    tests += 1
        rows.append({"n": n, "all_labelled_tournaments": len(graphs), "directed_wedge_tests": tests})
    emit("mixed_response_exhaustive", rows)


def main():
    base = staircase(4)
    check(base == [252, 169, 242, 164, 202, 144, 42, 64] and hp(base) == 233, "THM316 source")
    edges = list(combinations(range(8), 2))
    one = [(list(e), hp(flip(base, [e]))) for e in edges]
    two = [(list(e), list(f), flip(base, [e, f])) for e, f in combinations(edges, 2)
           if hp(flip(base, [e, f])) == 223]
    check(not any(h == 223 for _, h in one), "one-flip 223 hostile")
    check([(e, f) for e, f, _ in two] == [([0, 3], [2, 3]), ([4, 5], [4, 7])], "complete distance-two witnesses")
    emit("one_flip_complete_28", one)
    emit("two_flip_complete_378_matches", [{"edges": [e, f], "adjacency": a} for e, f, a in two])

    e, f = (0, 3), (2, 3)
    corners = [base, flip(base, [e]), flip(base, [f]), flip(base, [e, f])]
    corner_counts = [hp(a) for a in corners]
    check(corner_counts == [233, 291, 123, 223], "four corner counts")
    corner_paths = [brute_paths(a) for a in corners]
    check([len(p) for p in corner_paths] == corner_counts, "permutation audit")
    block_counts = [count_block(corner_paths[0], (0, 3, 2)), count_block(corner_paths[3], (2, 3, 0))]
    check(block_counts == [31, 11], "literal contiguous block counts")
    check([hp(contract(base, (0, 3, 2))), hp(contract(base, (2, 3, 0)))] == block_counts,
          "six-vertex contraction counts")
    check(corner_counts[3] - corner_counts[1] - corner_counts[2] + corner_counts[0] == sum(block_counts),
          "mixed curvature42")
    emit("four_corner_packets", [{"flips": flips, "H": h, **packet(a)}
                                for flips, h, a in zip([[], [e], [f], [e, f]], corner_counts, corners)])
    emit("mixed_response_certificate", {"curvature": 42, "forward_block": [0, 3, 2], "forward_count": 31,
                                        "reverse_block": [2, 3, 0], "reverse_count": 11})

    # Changing any subset of the middle vertex's outside arcs cannot change this mixed response.
    middle_edges = [(3, j) for j in range(8) if j not in (0, 2, 3)]
    for bits in range(1 << len(middle_edges)):
        varied = flip(base, [edge for i, edge in enumerate(middle_edges) if bits & (1 << i)])
        values = [hp(varied), hp(flip(varied, [e])), hp(flip(varied, [f])), hp(flip(varied, [e, f]))]
        check(values[3] - values[1] - values[2] + values[0] == 42, "middle-neighborhood blindness")
    emit("middle_neighborhood_controls", {"orientations": 32, "same_curvature": 42})

    padded = []
    for a in (base, corners[-1]):
        ten = ordered_padding(a)
        check(hp(ten) == hp(a), "ordered singleton padding")
        check(list(map(len, scc(ten))) == [1, 8, 1], "padding is not strong")
        padded.append({"H": hp(ten), "adjacency": ten, "scc_sizes": list(map(len, scc(ten)))})
    emit("literal_ten_vertex_lifts", padded)
    family = []
    for k in range(2, 8):
        a = staircase(k)
        b = flip(a, [e, f])
        family.append({"pairs": k, "H_before": hp(a), "H_after": hp(b), "difference": hp(a) - hp(b)})
    emit("same_local_move_other_orders", family)
    mixed_response_audit()
    print("PASS: exact local graph identities; no arithmetic or Collatz dynamics transfer asserted")


if __name__ == "__main__":
    main()
