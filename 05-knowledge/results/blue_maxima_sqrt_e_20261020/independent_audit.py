#!/usr/bin/env python3
"""Independent checks not importing the producer's implementation."""
import json
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, permutations, product
from math import factorial
from pathlib import Path

HERE = Path(__file__).resolve().parent


def rows_from_bits(k, a, b, matching=None):
    rows = [0] * (2 * k + (matching is not None))
    pair_bit = {pair: t for t, pair in enumerate(combinations(range(k), 2))}
    for u in range(2*k):
        for v in range(2*k):
            i, s = divmod(u, 2)
            j, t = divmod(v, 2)
            if i == j:
                edge = matching is not None and s != t and ((matching >> i & 1) == (s == 0))
            elif s != t:
                edge = (b >> pair_bit[tuple(sorted((i, j)))] & 1) == (s == 0)
            else:
                edge = (a >> pair_bit[tuple(sorted((i, j)))] & 1) == ((i < j) != bool(s))
            if edge:
                rows[u] |= 1 << v
    if matching is not None:
        for u in range(2*k):
            if u % 2:
                rows[2*k] |= 1 << u
            else:
                rows[u] |= 1 << (2*k)
    return rows


def h_by_recursion(k, a, b):
    rows = rows_from_bits(k, a, b)
    # Backward recurrence: choose a legal predecessor for the current first vertex.
    from functools import lru_cache
    @lru_cache(None)
    def count(remaining, first):
        if remaining == 0:
            return 1
        result = 0
        for i in range(k):
            if remaining >> i & 1:
                for v in (2*i, 2*i+1):
                    if rows[v] >> first & 1:
                        result += count(remaining ^ (1 << i), v)
        return result
    return sum(count(((1 << k)-1) ^ (1 << i), 2*i) for i in range(k))


def stable_colors(rows):
    """Canonical equitable vertex coloring, independently for each tournament."""
    n = len(rows)
    colors = [rows[u].bit_count() for u in range(n)]
    while True:
        palette = sorted(set(colors))
        signature = [
            (colors[u], tuple(sum((rows[u] >> v & 1) for v in range(n) if colors[v] == c)
                              for c in palette))
            for u in range(n)
        ]
        lookup = {x: i for i, x in enumerate(sorted(set(signature)))}
        new_colors = [lookup[x] for x in signature]
        if len(set(new_colors)) == len(set(colors)):
            return new_colors, sorted(signature)
        colors = new_colors


def isomorphic(left, right):
    """Exact tournament isomorphism with degree/refinement pruning."""
    lc, ls = stable_colors(left)
    rc, rs = stable_colors(right)
    if ls != rs:
        return False
    n = len(left)
    candidates = [[v for v in range(n) if lc[u] == rc[v]] for u in range(n)]
    mapping = {}
    used = set()
    def search():
        if len(mapping) == n:
            return True
        best_u = None
        best_c = None
        for u in range(n):
            if u in mapping:
                continue
            choices = [v for v in candidates[u] if v not in used and
                       all((left[u] >> x & 1) == (right[v] >> y & 1)
                           for x, y in mapping.items())]
            if not choices:
                return False
            if best_c is None or len(choices) < len(best_c):
                best_u, best_c = u, choices
        for v in best_c:
            mapping[best_u] = v
            used.add(v)
            if search():
                return True
            used.remove(v)
            del mapping[best_u]
        return False
    return search()


def triangle_count(rows):
    return sum(
        ((rows[a] >> b & 1) and (rows[b] >> c & 1) and (rows[c] >> a & 1))
        or ((rows[b] >> a & 1) and (rows[a] >> c & 1) and (rows[c] >> b & 1))
        for a, b, c in combinations(range(len(rows)), 3)
    )


def completion_audit(k, a, b, maximum):
    rows = rows_from_bits(k, a, b)
    automorphisms = []
    for p in permutations(range(k)):
        if all(((rows[2*i+s] >> (2*j+t)) & 1) ==
               ((rows[2*p[i]+s] >> (2*p[j]+t)) & 1)
               for i, j in permutations(range(k), 2)
               for s, t in product(range(2), repeat=2)):
            automorphisms.append(p)
    assert len(automorphisms) == 5
    orbit_reps = set()
    full_rows = []
    counts = Counter()
    for matching in range(1 << k):
        orbit = {sum(((matching >> i) & 1) << p[i] for i in range(k))
                 for p in automorphisms}
        rep = min(orbit)
        if rep in orbit_reps:
            continue
        orbit_reps.add(rep)
        completed = rows_from_bits(k, a, b, rep)
        stabilizer = len(automorphisms) // len(orbit)
        mass = maximum // stabilizer
        counts[(mass, triangle_count(completed) % 2)] += 1
        full_rows.append((rep, completed))
    comparisons = 0
    for (a1, rows1), (a2, rows2) in combinations(full_rows, 2):
        assert not isomorphic(rows1, rows2), (k, a1, a2, "completion collision")
        comparisons += 1
    return {
        "k": k,
        "h_by_independent_backward_recursion": h_by_recursion(k, a, b),
        "frame_automorphism_order": len(automorphisms),
        "completion_orbits": len(orbit_reps),
        "pairwise_full_tournament_nonisomorphism_checks": comparisons,
        "mass_and_triangle_parity": [
            {"mass": m, "parity": p, "classes": count}
            for (m, p), count in sorted(counts.items())
        ],
    }


def local_audit():
    base = [56, 200, 35, 224, 14, 194, 21, 21]
    assert base == rows_from_bits(4, 2, 11)
    assert all(r.bit_count() == 3 for r in base)
    frames = []
    for p in permutations(range(4)):
        for s in product(range(2), repeat=4):
            lift = [2*p[i]+(t ^ s[i]) for i in range(4) for t in range(2)]
            rows = [sum(((base[lift[u]] >> lift[v]) & 1) << v for v in range(8))
                    for u in range(8)]
            assert all(((rows[u] >> v) & 1) == ((rows[v ^ 1] >> (u ^ 1)) & 1)
                       for u, v in permutations(range(8), 2))
            assert all(not (rows[2*i] >> (2*i+1) & 1) and
                       not (rows[2*i+1] >> (2*i) & 1) for i in range(4))
            frames.append(rows)
    patterns = {
        "edge": ([(0, 2)], Fraction(1, 2)),
        "wedge": ([(0, 2), (2, 4)], Fraction(5, 16)),
        "disjoint_pair": ([(0, 2), (4, 6)], Fraction(1, 4)),
        "three_edge_path": ([(0, 2), (2, 4), (4, 6)], Fraction(3, 16)),
    }
    output = {}
    for label, (edges, expected) in patterns.items():
        probability = Fraction(sum(all(rows[u] >> v & 1 for u, v in edges)
                                   for rows in frames), len(frames))
        assert probability == expected
        output[label] = str(probability)
    excluded = []
    for i, j in combinations(range(4), 2):
        count = 0
        for u, v, w in permutations([2*i, 2*i+1, 2*j, 2*j+1], 3):
            if u // 2 == w // 2 and u // 2 != v // 2:
                count += (base[u] >> v & 1) * (base[v] >> w & 1)
        assert count == 2
        excluded.append(count)
    output["excluded_wedges_per_pair_of_parts"] = excluded
    output["valid_signed_permutations"] = len(frames)
    return output


output = {
    "local_distribution": local_audit(),
    "completion_checks": [
        completion_audit(6, 314, 3618, 1185),
        completion_audit(7, 268, 39375, 8295),
    ],
    "search_witnesses": [
        {"k": 8, "independent_h": h_by_recursion(8, 27827700, 242982930)},
        {"k": 9, "independent_h": h_by_recursion(9, 23122628159, 22262400520)},
        {"k": 10, "independent_h": h_by_recursion(10, 2996670331387, 17299192215552)},
    ],
}
assert output["completion_checks"][0]["h_by_independent_backward_recursion"] == 1185
assert output["completion_checks"][1]["h_by_independent_backward_recursion"] == 8295
assert output["search_witnesses"][0]["independent_h"] == 66123
assert output["search_witnesses"][1]["independent_h"] == 597249
assert output["search_witnesses"][2]["independent_h"] == 5975877
(HERE / "independent_audit.json").write_text(json.dumps(output, indent=2) + "\n")
print(json.dumps(output, indent=2))
