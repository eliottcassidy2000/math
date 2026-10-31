"""Finite checks of small-gain and triangle-packing corollaries."""

from collections import Counter
from functools import lru_cache
import hashlib
from itertools import combinations
import json
from pathlib import Path
import random

from cut_response import compile_response, from_code
from verify import direct_h


def induced(out, mask):
    vertices = [u for u in range(len(out)) if mask >> u & 1]
    return tuple(sum(1 << j for j, v in enumerate(vertices) if out[u] >> v & 1)
                 for u in vertices)


def triangles(out):
    result = []
    for triple in combinations(range(len(out)), 3):
        mask = sum(1 << v for v in triple)
        if all((out[v] & mask).bit_count() == 1 for v in triple):
            result.append(mask)
    return result


def triangle_packing(out, triangle_masks):
    @lru_cache(None)
    def maximum(mask):
        if not mask:
            return 0
        bit = mask & -mask
        best = maximum(mask ^ bit)
        for triangle in triangle_masks:
            if triangle & bit and triangle & mask == triangle:
                best = max(best, 1 + maximum(mask ^ triangle))
        return best
    return maximum((1 << len(out)) - 1)


def ordinal_triangles(number):
    out = [0] * (3 * number)
    for block in range(number):
        for offset in range(3):
            u = 3 * block + offset
            out[u] |= 1 << (3 * block + (offset + 1) % 3)
            out[u] |= sum(1 << v for v in range(3 * (block + 1), len(out)))
    return tuple(out)


def main():
    checks = Counter()
    digest = hashlib.sha256()

    def check(ok, label, data=None):
        if not ok:
            raise RuntimeError((label, data))
        checks[label] += 1

    rng = random.Random(124032)
    parents = [from_code(n, code) for n in range(6)
               for code in range(1 << (n * (n - 1) // 2))]
    parents += [from_code(n, rng.randrange(1 << (n * (n - 1) // 2)))
                for n in range(6, 10) for _ in range(12)]
    parents += [ordinal_triangles(3), ordinal_triangles(4)]
    controls = []
    for out in parents:
        n = len(out)
        full = (1 << n) - 1
        response = compile_response(out)
        check(response.h == direct_h(out), "parent_independent_dp")
        triangle_masks = triangles(out)
        nu = triangle_packing(out, triangle_masks)
        check(response.h >= 3 ** nu, "packing_path_lower_bound")
        check((response.h == 1) == (not triangle_masks), "count_one_iff_transitive")
        pair_h = {}
        for u in range(n):
            for v in range(u + 1, n):
                rest = full ^ (1 << u) ^ (1 << v)
                h = direct_h(induced(out, rest))
                pair_h[u, v] = pair_h[v, u] = h
                check(h >= 3 ** max(nu - 2, 0), "pair_deletion_packing_bound")
        values = response.all_values()
        positive_gains = []
        for cut, value in enumerate(values):
            gain = value - response.h
            if not gain:
                checks["zero_gain_cuts_excluded"] += 1
                continue
            positive_gains.append(gain)
            check(gain >= 2 * 3 ** max(nu - 2, 0), "positive_gain_packing_bound")
            k = gain // 2
            crossings = [(u, v) for u in range(n) if cut >> u & 1
                         for v in range(n) if out[u] >> v & 1 and not cut >> v & 1]
            check(bool(crossings), "positive_gain_has_crossing_arc")
            for u, v in crossings:
                check(pair_h[u, v] <= k, "every_crossing_pair_small_deletion")
                if gain <= 4:
                    check(pair_h[u, v] == 1, "gain_two_or_four_transitive_deletion")
        digest.update(json.dumps([out, response.h, nu,
                                  min(positive_gains, default=None)]).encode() + b"\n")
        if n in (9, 12) and out == ordinal_triangles(n // 3):
            controls.append({"out": out, "h": response.h, "nu": nu,
                             "minimum_positive_gain": min(positive_gains)})
    print(json.dumps({"status": "PASS", "parents": len(parents),
                      "universe": "all labelled n0..5; 12 seeded random per n6..9; ordinal sums of 3 and 4 triangles",
                      "seed": 124032, "checks": dict(sorted(checks.items())),
                      "structured_controls": controls, "data_sha256": digest.hexdigest(),
                      "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
                     sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
