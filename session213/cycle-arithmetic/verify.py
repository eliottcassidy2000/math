#!/usr/bin/env python3
"""Exact independent control of the binary compiler correction. Standard library."""
import itertools
import json


def require(ok, label):
    if not ok:
        raise RuntimeError(label)


def apex(bits):
    m = len(bits)
    rows = [0] * (m + 1)
    for i in range(m):
        for j in range(i + 1, m):
            rows[i] |= 1 << j
        if bits[i]:
            rows[i] |= 1 << m
        else:
            rows[m] |= 1 << i
    return rows


def endpoint_core(q, m):
    return apex([0, 1] + [(q >> i) & 1 for i in range(m)])


def sole_exit(rows, root):
    n = len(rows)
    return [r | ((1 << n) if i != root else 0)
            for i, r in enumerate(rows)] + [1 << root]


def count(rows, start=None, cycles=False):
    n = len(rows)
    if not n:
        return 1
    if cycles and n < 3:
        return 0
    if cycles:
        start = 0
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        if start is None or v == start:
            dp[1 << v][v] = 1
    full = (1 << n) - 1
    for mask in range(1, full + 1):
        for v, value in enumerate(dp[mask]):
            if value:
                remaining = rows[v] & (full ^ mask)
                while remaining:
                    bit = remaining & -remaining
                    remaining -= bit
                    dp[mask | bit][bit.bit_length() - 1] += value
    return sum(value for v, value in enumerate(dp[full])
               if not cycles or rows[v] & 1)


def brute(rows, start=None, cycles=False):
    n = len(rows)
    total = 0
    for p in itertools.permutations(range(n)):
        if start is not None and p[0] != start:
            continue
        if cycles and p[0] != 0:
            continue
        if all((rows[p[i]] >> p[i + 1]) & 1 for i in range(n - 1)):
            if not cycles or (rows[p[-1]] >> p[0]) & 1:
                total += 1
    return total


def weighted_pairs(bits):
    return sum(1 if j - i == 1 else 1 << (j - i - 2)
               for j in range(len(bits)) if bits[j]
               for i in range(j) if not bits[i])


def D(q):
    if not q:
        return 1
    if q & 1:
        return D(q // 2)
    half = q // 2
    return D(half) + half + (half & 1)


def induced(rows, keep):
    return [sum(((rows[v] >> w) & 1) << j for j, w in enumerate(keep))
            for v in keep]


def all_tournaments(n):
    pairs = list(itertools.combinations(range(n), 2))
    for code in range(1 << len(pairs)):
        rows = [0] * n
        for i, (v, w) in enumerate(pairs):
            if (code >> i) & 1:
                rows[v] |= 1 << w
            else:
                rows[w] |= 1 << v
        yield rows


def main():
    report = {"status": "PASS", "scope": "FINITE-EXACT controls; proofs in PROOF.md",
              "endpoint_rows": [], "tournament_rows": []}
    family_total = 0
    permutation_total = 0
    for m in range(9):
        for q in range(1 << m):
            family_total += 1
            bits = [(q >> i) & 1 for i in range(m)]
            a = apex(bits)
            t = endpoint_core(q, m)
            h = q + 1
            d = D(q)
            k = count(t)
            require(count(a) == 1 + 2 * weighted_pairs(bits) == d,
                    ("one-apex", m, q))
            require(count(t, start=1) == h, ("endpoint", m, q))
            require(k == 2 * h + d, ("correction", m, q))
            require((k == 2 * h + 1) == (h & (h - 1) == 0),
                    ("equality locus", m, q))
            b = sole_exit(t, 1)
            require(count(b, cycles=True) == h, ("cycle closure", m, q))
            require(count(b) >= len(b) * h, ("closure floor", m, q))
            for newq in (2 * q, 2 * q + 1):
                newd = D(newq)
                newh = newq + 1
                if newq & 1:
                    require(newd == d and 2 * newh + newd == 2 * k - d,
                            ("one branch", m, q))
                else:
                    require(newd == d + h - 1 + ((h - 1) & 1),
                            ("zero correction", m, q))
                    require(2 * newh + newd == 2 * k - d + h - 3 + ((h - 1) & 1),
                            ("zero branch", m, q))
            if m <= 4:
                require(brute(t) == k and brute(t, start=1) == h,
                        ("core permutations", m, q))
                require(brute(b, cycles=True) == h, ("closure permutations", m, q))
                permutation_total += 1
        report["endpoint_rows"].append({"m": m, "all_codes": 1 << m,
                                        "core_order": m + 3})
    for n in range(1, 6):
        checked = 0
        deletions = 0
        for rows in all_tournaments(n):
            checked += 1
            h = count(rows)
            c = count(rows, cycles=True)
            require(h % 2 == 1 and h >= n * c, ("parity/closure", n, rows))
            for removed in range(n):
                old = count(induced(rows, [v for v in range(n) if v != removed]))
                require(h >= old, ("monotonicity", n, rows, removed))
                require(not(old >= 5 and h - old == 2), ("increment-two", n, rows, removed))
                deletions += 1
        report["tournament_rows"].append({"n": n, "labeled_tournaments": checked,
                                         "deletions": deletions})
    controls = []
    for h in (1, 2, 3, 4, 10, 11, 16, 23, 32, 100):
        q = h - 1
        controls.append({"h": h, "D": D(q), "H_core": 2 * h + D(q),
                         "proposed_target": 2 * h + 1})
    report.update({"family_cases": family_total, "permutation_family_cases": permutation_total,
                   "controls": controls})
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
