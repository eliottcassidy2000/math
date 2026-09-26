"""Exact graph/height probes for crossroads233_20260926_graph.md.

No third-party dependencies. All checks use require(), so python -O retains
them. The SCC implementation uses reachability, independently of cut coverage.
"""
from fractions import Fraction as F
from itertools import combinations, permutations, product
import json


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def emit(label, value):
    print(label + "=" + json.dumps(value, sort_keys=True))


def spine_tournament(n, backward):
    backward = set(backward)
    adj = [0] * n
    for i, j in combinations(range(n), 2):
        u, v = (j, i) if (i, j) in backward else (i, j)
        adj[u] |= 1 << v
    return adj


def components(adj):
    n = len(adj)
    reach = [adj[i] | (1 << i) for i in range(n)]
    for k in range(n):
        for i in range(n):
            if reach[i] & (1 << k):
                reach[i] |= reach[k]
    unseen = set(range(n))
    blocks = []
    while unseen:
        i = min(unseen)
        block = [j for j in sorted(unseen)
                 if reach[i] & (1 << j) and reach[j] & (1 << i)]
        blocks.append(block)
        unseen.difference_update(block)
    return blocks


def cut_blocks(n, backward):
    cuts = [r for r in range(1, n)
            if not any(i < r <= j for i, j in backward)]
    ends = [0] + cuts + [n]
    return [list(range(i, j)) for i, j in zip(ends, ends[1:])]


def cycle_triangle(adj, triple):
    i, j, k = triple
    a, b, c = bool(adj[i] & (1 << j)), bool(adj[j] & (1 << k)), bool(adj[i] & (1 << k))
    return a == b and a != c


def hamilton_count(adj):
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1, 1 << n):
        for i in range(n):
            if not dp[mask][i]:
                continue
            remaining = adj[i] & ~mask
            while remaining:
                bit = remaining & -remaining
                remaining -= bit
                dp[mask | bit][bit.bit_length() - 1] += dp[mask][i]
    return sum(dp[-1]), dp[-1]


def graph_census():
    rows = []
    for n in range(2, 8):
        chords = [(i, j) for i in range(n) for j in range(i + 2, n)]
        reduced = [(a, b + 1) for a, b in combinations(range(n - 1), 2)]
        require(chords == reduced, "edge-index bijection")
        for mask in range(1 << len(chords)):
            back = [edge for p, edge in enumerate(chords) if mask & (1 << p)]
            adj = spine_tournament(n, back)
            require(components(adj) == cut_blocks(n, back), "SCC cut theorem")
            triangle = any(cycle_triangle(adj, (i, i + 1, j))
                           for i in range(n - 2) for j in range(i + 2, n))
            require(triangle == bool(back), "spine triangle test family")
            if n <= 6:
                h, _ = hamilton_count(adj)
                require((h == 1) == (not back), "unique Hamilton path equivalence")
        rows.append({"n": n, "free_chords": len(chords), "tournaments": 1 << len(chords)})
    emit("fixed_spine_exhaustive_census", rows)
    adj = spine_tournament(4, [(0, 3)])
    require(not cycle_triangle(adj, (0, 1, 2)), "local triangle control")
    require(not cycle_triangle(adj, (1, 2, 3)), "local triangle control")
    require(cycle_triangle(adj, (0, 1, 3)), "long triangle hostile")
    emit("consecutive_triangle_hostile", {"n": 4, "only_backward_chord": [0, 3]})


def metric_controls():
    equilateral = [[int(i != j) for j in range(3)] for i in range(3)]
    c4 = [[min((i - j) % 4, (j - i) % 4) for j in range(4)] for i in range(4)]
    for d in (equilateral, c4):
        require(all(d[i][k] <= d[i][j] + d[j][k]
                    for i, j, k in product(range(len(d)), repeat=3)), "metric triangle inequality")
    require(all(sum(sorted(c4[i][j] for i, j in combinations(t, 2))[:2]) ==
                max(c4[i][j] for i, j in combinations(t, 2))
                for t in combinations(range(4), 3)), "C4 degenerate triangles")
    def is_line_order(p):
        return all(c4[p[i]][p[k]] == c4[p[i]][p[j]] + c4[p[j]][p[k]]
                   for i, j, k in combinations(range(4), 3))
    require(not any(is_line_order(p) for p in permutations(range(4))), "C4 non-line metric")
    emit("metric_controls", {"equilateral_K3": "triangle inequality with graph cycle",
                              "C4_metric": "all triangles degenerate, no compatible line order"})


def staircase(k):
    n = 2 * k
    rank = [i // 2 + (k if i % 2 else 0) for i in range(n)]
    adj = [0] * n
    for i, j in combinations(range(n), 2):
        u, v = (j, i) if i // 2 == j // 2 else ((i, j) if rank[i] < rank[j] else (j, i))
        adj[u] |= 1 << v
    return adj


def affine_word(word):
    a, b = [F(1)], [F(0)]
    for k in word:
        a.append(3 * a[-1] / 2 ** k)
        b.append((3 * b[-1] + 1) / 2 ** k)
    return a, b


def thresholds(a, b):
    return {(i, j): None if a[j] > a[i] else (b[j] - b[i]) / (a[i] - a[j])
            for i in range(len(a)) for j in range(i + 2, len(a))}


def cut_thresholds(n, theta):
    result = []
    for r in range(1, n):
        values = [t for (i, j), t in theta.items() if i < r <= j]
        result.append(None if None in values else max([F(0)] + values))
    return result


def direct_backward(a, b, x):
    heights = [u * x + v for u, v in zip(a, b)]
    return {(i, j) for i in range(len(a)) for j in range(i + 2, len(a))
            if heights[j] > heights[i]}


def affine_census():
    words, chambers, max_splits = 0, 0, 0
    for length in range(1, 6):
        for word in product(range(1, 5), repeat=length):
            words += 1
            a, b = affine_word(word)
            require(all(b[i] / a[i] < b[i + 1] / a[i + 1] for i in range(length)), "normalized carry order")
            theta = thresholds(a, b)
            positive = sorted({t for t in theta.values() if t is not None and t > 0})
            endpoints = [F(0)] + positive
            samples = [(u + v) / 2 for u, v in zip(endpoints, endpoints[1:])]
            samples.append(endpoints[-1] + 1)
            previous_back, previous_blocks = None, None
            all_blocks = []
            cut_theta = cut_thresholds(len(a), theta)
            for x in samples:
                chambers += 1
                back = direct_backward(a, b, x)
                predicted = {edge for edge, t in theta.items() if t is None or x < t}
                require(back == predicted, "affine threshold prediction")
                blocks = components(spine_tournament(len(a), back))
                cuts = [r for r, t in enumerate(cut_theta, 1) if t is not None and x > t]
                require(cuts == [block[0] for block in blocks[1:]], "exact cut thresholds")
                if previous_back is not None:
                    require(back <= previous_back, "backward set monotone")
                    require(all(any(set(block) <= set(old) for old in previous_blocks)
                                for block in blocks), "SCCs only split")
                previous_back, previous_blocks = back, blocks
                all_blocks.append(blocks)
            distinct_blocks = len({tuple(tuple(block) for block in blocks) for blocks in all_blocks})
            require(distinct_blocks <= length + 1, "SCC chamber bound")
            max_splits = max(max_splits, distinct_blocks - 1)
    emit("affine_word_census", {"lengths": [1, 5], "exponents": [1, 4], "words": words,
                               "positive_real_chambers_checked": chambers, "max_component_change_events": max_splits})


def odd_orbit(n, steps):
    nodes, word = [n], []
    for _ in range(steps):
        u = 3 * n + 1
        k = (u & -u).bit_length() - 1
        n = u >> k
        nodes.append(n)
        word.append(k)
    return nodes, word


def witness_233():
    word = [1, 1, 2, 1, 1, 2, 6, 1, 1, 1, 1, 2, 2, 1, 2, 1, 1]
    a, b = affine_word(word)
    n = len(a)
    theta = thresholds(a, b)
    ct = cut_thresholds(n, theta)
    contact = F(1441352971, 5077565)
    require(sum(word) == 27 and 3 ** 17 == 129140163, "witness clock")
    require(b[-1] * 2 ** 27 == 1441352971, "witness carry")
    require(theta[0, 17] == contact, "witness contact")
    require([r for r, t in enumerate(ct, 1) if t is not None] == [7], "only finite cut")
    require(ct[6] == contact, "cut7 contact")
    all_contacts = [(b[j] - b[i]) / (a[i] - a[j])
                    for i in range(n) for j in range(i + 1, n) if a[i] > a[j]]
    require(max(all_contacts) == contact, "all positive contacts below the next cylinder point")
    records = []
    for t in range(4):
        x = 231 + (2 ** 28) * t
        nodes, actual_word = odd_orbit(x, 17)
        require(actual_word == word, "exact valuation cylinder")
        require(len(set(nodes)) == len(nodes), "injective cylinder control")
        require(nodes == [int(u * x + v) for u, v in zip(a, b)], "direct affine agreement")
        blocks = components(spine_tournament(n, direct_backward(a, b, F(x))))
        require(blocks == ([list(range(18))] if t == 0 else [list(range(7)), list(range(7, 18))]), "233 split")
        records.append({"t": t, "source": x, "endpoint": nodes[-1], "scc_sizes": list(map(len, blocks))})
    emit("recovered_231_to_233", {"word": word, "nodes": odd_orbit(231, 17)[0],
                                  "carry": 1441352971, "contact": str(contact),
                                  "only_finite_cut": 7, "cylinder_modulus": 2 ** 28,
                                  "sampled_realizations": records})


if __name__ == "__main__":
    metric_controls()
    graph_census()
    staircase_rows = []
    for k in range(1, 6):
        adj = staircase(k)
        h, ends = hamilton_count(adj)
        require(all(bool(adj[i] & (1 << j)) == bool(adj[2*k-1-j] & (1 << (2*k-1-i)))
                    for i in range(2*k) for j in range(2*k) if i != j), "staircase anti-automorphism")
        staircase_rows.append({"k": k, "n": 2*k, "Hamilton_paths": h})
        if k == 4:
            require(h == 233 and list(reversed(ends)) == [80, 29, 39, 24, 21, 18, 11, 11], "THM316 recovery")
    emit("THM316_recovery", staircase_rows)
    affine_census()
    witness_233()
    print("STATUS=PASS; elementary proofs are in the companion note; no Collatz closure")
