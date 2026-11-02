#!/usr/bin/env python3
"""Small proof controls; no repository imports and no assert-dependent gates.

Coverage of orders 5,6 uses Busch 2006 Lemma 1, not a claim that 272
labelled tournaments are the entire labelled universe. Cross arcs range
over all orientations after the two intersecting strong pieces are fixed.
"""
from itertools import combinations, permutations
from math import isqrt
import json


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def decode(n, code):
    rows = [0] * n
    for k, (i, j) in enumerate(combinations(range(n), 2)):
        u, v = (i, j) if (code >> k) & 1 else (j, i)
        rows[u] |= 1 << v
    return tuple(rows)


def encode(rows):
    return sum(((rows[i] >> j) & 1) << k
               for k, (i, j) in enumerate(combinations(range(len(rows)), 2)))


def strong(rows):
    full = (1 << len(rows)) - 1
    for s in range(len(rows)):
        reached = 1 << s
        while True:
            nxt = reached
            for i, row in enumerate(rows):
                if (reached >> i) & 1:
                    nxt |= row
            if nxt == reached:
                break
            reached = nxt
        if reached != full:
            return False
    return True


def count_dp(rows):
    n = len(rows)
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, full + 1):
        for v in range(n):
            count = dp[mask][v]
            if not count:
                continue
            rest = rows[v] & (full ^ mask)
            while rest:
                bit = rest & -rest
                rest -= bit
                dp[mask | bit][bit.bit_length() - 1] += count
    return sum(dp[full])


def count_literal(rows):
    return sum(all((rows[u] >> v) & 1 for u, v in zip(p, p[1:]))
               for p in permutations(range(len(rows))))


def relabel(rows, labels):
    ans = [0] * len(rows)
    for u, row in enumerate(rows):
        for v in range(len(rows)):
            if (row >> v) & 1:
                ans[labels[u]] |= 1 << labels[v]
    return tuple(ans)


def near(n):
    rows = list(decode(n, (1 << (n * (n - 1) // 2)) - 1))
    rows[0] ^= 1 << (n - 1)
    rows[n - 1] |= 1
    return tuple(rows)


def overlap(second, root, mask):
    """A directed triangle on 0,1,2, glued at 0 to a rooted strong B."""
    m = len(second)
    n = m + 2
    labels = {root: 0}
    labels.update(zip((i for i in range(m) if i != root), range(3, n)))
    ans = [0] * n
    ans[0] |= 1 << 1
    ans[1] |= 1 << 2
    ans[2] |= 1 << 0
    for u in range(m):
        for v in range(m):
            if (second[u] >> v) & 1:
                ans[labels[u]] |= 1 << labels[v]
    for k, (u, v) in enumerate((u, v) for u in (1, 2) for v in range(3, n)):
        a, b = (u, v) if (mask >> k) & 1 else (v, u)
        ans[a] |= 1 << b
    return tuple(ans)


def factor(n):
    fs = []
    p = 3
    while p * p <= n:
        while n % p == 0:
            fs.append(p)
            n //= p
        p += 2
    if n != 1:
        fs.append(n)
    return fs


def allowed(n):
    return n > 0 and n % 2 == 1 and n not in (7, 21)


def atom(n):
    fs = factor(n)
    return (len(fs) == 1 and n != 7) or n in (63, 343) or (
        len(fs) == 2 and 7 in fs and n != 21)


def arithmetic_max_length(n):
    fs = factor(n)
    a = fs.count(7)
    if a == 1 and all(p in (3, 7) for p in fs):
        return len(fs) - 2
    return len(fs) - (a + 1) // 2


def main():
    ledger = {}
    strong4 = {encode(decode(4, c)) for c in range(64) if strong(decode(4, c))}
    canonical = decode(4, min(strong4))
    orbit = {encode(relabel(canonical, p)) for p in permutations(range(4))}
    need(len(strong4) == 24 and strong4 == orbit, "strong4 orbit cover")
    need({count_literal(decode(4, c)) for c in strong4} == {5}, "strong4 spectrum")
    ledger["strong4"] = {"canonical_code": min(strong4), "labelled_count": 24,
                         "single_relabelling_orbit": True, "spectrum": [5]}
    expected = {5: [9, 11, 13, 15],
                6: [15, 17, 19, 23, 25, 27, 29, 31, 33, 37, 41, 43, 45]}
    for n, second, roots in ((5, (2, 4, 1), [0]), (6, canonical, range(4))):
        values = set()
        checks = 0
        for root in roots:
            for mask in range(1 << (2 * (n - 3))):
                rows = overlap(second, root, mask)
                need(strong(rows), "overlap strongness")
                d = count_dp(rows)
                need(d == count_literal(rows), "DP/literal overlap disagreement")
                need(d % 2 == 1, "parity")
                values.add(d)
                checks += 1
        nt = near(n)
        need(strong(nt), "near strongness")
        hnt = count_literal(nt)
        need(hnt == count_dp(nt) == 2 ** (n - 2) + 1, "near count")
        values.add(hnt)
        need(sorted(values) == expected[n], "small spectrum")
        need(7 not in values and 21 not in values, "forbidden value")
        ledger[f"order{n}_structural_cover"] = {"overlap_cases": checks,
            "nearly_transitive_count": hnt, "strong_spectrum": sorted(values)}
    # Independent explicit carry witnesses. The 63 rows come from the
    # repository's h63_counterexample_audit_s8.py; 49 and 343 were found
    # deterministically with random.Random(142), then frozen here.
    carries = {
        49: (110, 44, 32, 116, 7, 16, 54),
        63: (162, 252, 249, 17, 1, 152, 185, 24),
        343: (236, 161, 74, 130, 143, 220, 26, 68),
    }
    ledger["strong_carries"] = {}
    for h, rows in carries.items():
        need(decode(len(rows), encode(rows)) == rows, "carry tournament validity")
        need(strong(rows), "carry strongness")
        need(count_dp(rows) == count_literal(rows) == h, "carry count")
        ledger["strong_carries"][h] = {"rows": rows, "code": encode(rows),
                                         "order": len(rows)}
    # Exact target-monoid controls, independent multiplicative dynamic program.
    limit = 10001
    best = {1: 0}
    atoms = []
    for n in range(3, limit + 1, 2):
        if not allowed(n):
            continue
        best[n] = 1
        composite = False
        for d in range(3, isqrt(n) + 1, 2):
            if n % d == 0 and d in best and n // d in best:
                composite = True
                best[n] = max(best[n], best[d] + best[n // d])
        need(atom(n) == (not composite), "target atom classification")
        need(best[n] == arithmetic_max_length(n), "target factor length")
        if not composite:
            atoms.append(n)
    ledger["target_monoid"] = {"bound": limit, "allowed_values_checked": len(best),
        "atoms_checked": len(atoms), "formula_matches": True,
        "not_a_tournament_existence_check": True}
    print(json.dumps(ledger, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
