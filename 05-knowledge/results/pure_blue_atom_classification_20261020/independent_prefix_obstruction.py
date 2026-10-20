#!/usr/bin/env python3
"""Independent constructive audit of the all-order prefix obstruction.

No producer imports or producer data. Default exhaustive universe: every
normalized skew connection set (1 in S) on each odd order 9 through 31.
The finite order-seven words are literal proof data, tested independently.
"""
import argparse
import itertools
import json
import math
from pathlib import Path


def permutation_order(p):
    seen = set()
    result = 1
    for start in range(len(p)):
        if start in seen:
            continue
        cur, length = start, 0
        while cur not in seen:
            seen.add(cur)
            length += 1
            cur = p[cur]
        assert cur == start
        result = math.lcm(result, length)
    return result


def skew_set(m, mask):
    n = 2 * m + 1
    return {1} | {d if mask & (1 << (d - 2)) else n - d
                  for d in range(2, m + 1)}


def prefix(n, S):
    """Find one of the three configurations, otherwise return residual type."""
    m = n // 2
    if 2 not in S and 3 in S:
        return [2, 0, 3, 1] + list(range(4, n)), "F1"
    for k in range(4, m + 1):
        if k in S and k - 2 not in S:
            q = list(range(n))
            q[0], q[k - 1] = q[k - 1], q[0]
            return q, "F2"
    neg = next((d for d in range(2, m + 1) if d not in S), None)
    if neg is not None and 3 <= neg <= m - 1:
        k = neg
        if k % 2:
            q = [1, k, 0] + list(range(2, k)) + list(range(k + 1, n))
        else:
            q = [1] + list(range(3, k + 1)) + [0, 2] + list(range(k + 1, n))
        return q, "F3odd" if k % 2 else "F3even"
    low = S & set(range(1, m + 1))
    if low == {1}:
        return None, "A"
    if low == set(range(1, m)):
        return None, "B"
    if low == set(range(1, m + 1)):
        return None, "C"
    raise AssertionError(("uncovered sign pattern", n, sorted(S)))


def fail_pair(path, S):
    n = len(path)
    reversal = [0] * n
    for i, x in enumerate(path):
        reversal[x] = path[n - 1 - i]
    return next(((x, y) for x in range(n) for y in range(n)
                 if (y - x) % n in S
                 and (reversal[x] - reversal[y]) % n not in S), None)


def check_path(path, S):
    n = len(path)
    assert sorted(path) == list(range(n))
    assert all((b - a) % n in S for a, b in zip(path, path[1:]))


def construct(n, original):
    m = n // 2
    S = set(original)
    multiplier = 1
    normalizations = []
    for depth in range(3):
        q, family = prefix(n, S)
        if q is not None:
            break
        a = -2 if family == "A" else 2
        assert a % n in S and math.gcd(a, n) == 1
        S = {d for d in range(1, n) if a * d % n in S}
        multiplier = multiplier * a % n
        normalizations.append(a)
    else:
        raise AssertionError(("normalization did not stop", n, original))
    assert len(normalizations) <= 2
    assert all(q[i] == i for i in range(m, n))
    assert permutation_order(q) % 2 == 0
    qi = [0] * n
    for i, v in enumerate(q):
        qi[v] = i
    r = lambda i: n - 1 - i
    mirror_inverse = [r(qi[r(i)]) for i in range(n)]
    assert all(mirror_inverse[i] == i for i in range(m + 1))
    assert not ({i for i in range(n) if q[i] != i}
                & {i for i in range(n) if mirror_inverse[i] != i})
    g = [q[mirror_inverse[i]] for i in range(n)]
    assert permutation_order(g) == permutation_order(q)
    check_path(q, S)
    assert fail_pair(q, S) is not None
    # Check the proof's direct, prescribed failed arc, not merely existence.
    if family == "F2":
        pair = (0, q[0]) if q[0] in S else (q[0], 0)
    else:
        pair = (0, 1)
    u, v = pair
    reversal = [q[r(qi[i])] for i in range(n)]
    assert (v - u) % n in S
    assert (reversal[u] - reversal[v]) % n not in S
    path = [multiplier * v % n for v in q]
    check_path(path, original)
    bad = fail_pair(path, original)
    assert bad is not None
    return family, len(normalizations), path, bad


def small_checks():
    n7 = [
        ({1, 4, 5}, [0, 1, 2, 6, 4, 5, 3], (0, 5)),
        ({1, 2, 4}, [0, 1, 2, 4, 5, 6, 3], (0, 2)),
        ({1, 3, 5}, [0, 1, 2, 5, 3, 6, 4], (2, 0)),
        ({1, 2, 3}, [0, 1, 3, 4, 6, 2, 5], (0, 2)),
    ]
    controls = []
    for S, path, pair in n7:
        check_path(path, S)
        rev = {v: path[6 - i] for i, v in enumerate(path)}
        u, v = pair
        assert (v - u) % 7 in S
        assert (rev[u] - rev[v]) % 7 not in S
        assert fail_pair(path, S) is not None
    for n in (3, 5):
        for mask in range(1 << (n // 2 - 1)):
            S = skew_set(n // 2, mask)
            hp = [p for p in itertools.permutations(range(n))
                  if all((b - a) % n in S for a, b in zip(p, p[1:]))]
            assert all(fail_pair(p, S) is None for p in hp)
            assert len(hp) == (3 if n == 3 else 15)
            controls.append({"n": n, "S": sorted(S), "HP": len(hp)})
    return {"n7_literal_certificates": 4, "positive_controls": controls}


def center_five_checks():
    """Directly distinguish all three free-flip partners of the C5 center."""
    rows = []
    for path in ([0, 1, 2, 3, 4], [0, 2, 3, 4, 1], [0, 2, 4, 1, 3]):
        retained = {frozenset((a, b)) for a, b in zip(path, path[1:])}
        adj = [[u != v and (((v - u) % 5 in {1, 2})
                            if frozenset((u, v)) in retained
                            else ((u - v) % 5 in {1, 2}))
                for v in range(5)] for u in range(5)]
        hp = sum(all(adj[a][b] for a, b in zip(p, p[1:]))
                 for p in itertools.permutations(range(5)))
        aut = sum(all(adj[u][v] == adj[p[u]][p[v]]
                      for u in range(5) for v in range(5))
                  for p in itertools.permutations(range(5)))
        rows.append({"path": path, "HP": hp, "Aut": aut})
    assert [(x["HP"], x["Aut"]) for x in rows] == [(11, 1), (15, 3), (13, 1)]
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-n", type=int, default=31)
    args = ap.parse_args()
    result = {"small": small_checks(), "center_five": center_five_checks(),
              "exhaustive": [],
              "large_structured_residuals": []}
    total = 0
    for n in range(9, args.max_n + 1, 2):
        m = n // 2
        families, depths = {}, {}
        for mask in range(1 << (m - 1)):
            fam, depth, _, _ = construct(n, skew_set(m, mask))
            families[fam] = families.get(fam, 0) + 1
            depths[depth] = depths.get(depth, 0) + 1
        row = {"n": n, "normalized_presentations": 1 << (m - 1),
               "families": families, "normalization_counts": depths}
        result["exhaustive"].append(row)
        total += row["normalized_presentations"]
        print(json.dumps(row), flush=True)
    for m in range(4, 201):
        n = 2 * m + 1
        for low in ({1}, set(range(1, m)), set(range(1, m + 1))):
            S = {d if d in low else n - d for d in range(1, m + 1)}
            construct(n, S)
    result["large_structured_residuals"] = {
        "m_min": 4, "m_max": 200, "cases_per_m": 3, "passed": 591}
    result["total_normalized_presentations"] = total
    result["failures"] = 0
    Path(__file__).with_suffix(".json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"total": total, "failures": 0, "positive_and_n7": "passed",
                      "large_residual_cases": 591}), flush=True)


if __name__ == "__main__":
    main()
