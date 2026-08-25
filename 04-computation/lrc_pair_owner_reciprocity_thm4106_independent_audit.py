#!/usr/bin/env python3
"""Independent exact audit of the proposed THM-4106 pair/row formulas."""

from fractions import Fraction
from itertools import combinations
from math import gcd


def req(flag, message):
    if not flag:
        raise RuntimeError(message)


def vtwo(x):
    ans = 0
    while x % 2 == 0:
        ans += 1
        x //= 2
    return ans


def norm(x):
    x %= 1
    return min(x, 1 - x)


def owner_by_cells(u, v):
    cuts = sorted(
        {Fraction(j, v - u) for j in range(v - u)}
        | {Fraction(j, v + u) for j in range(v + u)}
    )
    mass = Fraction(0)
    for i, left in enumerate(cuts):
        right = cuts[i + 1] if i + 1 < len(cuts) else Fraction(1)
        mid = (left + right) / 2
        if norm(u * mid) < norm(v * mid):
            mass += right - left
    return mass, len(cuts)


def brute_max(u, v):
    points = {Fraction(0)}
    for k in (u, v):
        points |= {Fraction(j, 2 * k) for j in range(2 * k)}
    for k in (v - u, v + u):
        points |= {Fraction(j, k) for j in range(k)}
    values = {t: min(norm(u * t), norm(v * t)) for t in points}
    best = max(values.values())
    return best, tuple(sorted(t for t, value in values.items() if value == best))


def formula(u, v):
    g = gcd(u, v)
    m, n = u // g, v // g
    if m % 2 and n % 2:
        maximum = Fraction(1, 2)
        owner = Fraction(1, 2)
        ties = 2 * v - 2 * g
        core = (Fraction(1, 2),)
    else:
        q = m + n
        maximum = Fraction(q - 1, 2 * q)
        owner = Fraction(1, 2) + Fraction(1, 2 * (n * n - m * m))
        ties = 2 * v - g
        k = ((q - 1) // 2) * pow(m, -1, q) % q
        core = tuple(sorted({Fraction(k, q), Fraction((-k) % q, q)}))
    lifted = tuple(sorted({(t + j) / g for t in core for j in range(g)}))
    return maximum, owner, ties, lifted


def directed_summary(left, right):
    u, v = sorted((left, right))
    maximum, owner, _, _ = formula(u, v)
    return maximum, owner if left == u else 1 - owner


def decode(summary):
    maximum, owner = summary
    defect = Fraction(1, 2) - maximum
    imbalance = owner - Fraction(1, 2)
    req(defect > 0 and imbalance, "blind edge passed to decoder")
    q = Fraction(1, 2 * defect)
    r = defect / abs(imbalance)
    m, n = (q - r) / 2, (q + r) / 2
    req(m.denominator == n.denominator == 1, "nonintegral decoded core")
    return m / n if imbalance > 0 else n / m


def gcd_many(row):
    ans = 0
    for x in row:
        ans = gcd(ans, x)
    return ans


def recover(row):
    # Deliberately choose the last vertex as root to stress arbitrary label order.
    root = len(row) - 1
    levels = tuple(vtwo(x) for x in row)
    bridge = next(i for i in range(len(row)) if levels[i] != levels[root])
    edges = []
    for i in range(len(row)):
        if i == root:
            continue
        edges.append((root, i) if levels[i] != levels[root] else (bridge, i))
    adj = [[] for _ in row]
    for i, j in edges:
        ratio = decode(directed_summary(row[i], row[j]))
        adj[i].append((j, ratio))
        adj[j].append((i, 1 / ratio))
    coord = {root: Fraction(1)}
    stack = [root]
    while stack:
        i = stack.pop()
        for j, i_over_j in adj[i]:
            if j not in coord:
                coord[j] = coord[i] / i_over_j
                stack.append(j)
    lcm = 1
    for x in coord.values():
        lcm = lcm * x.denominator // gcd(lcm, x.denominator)
    ans = tuple(int(coord[i] * lcm) for i in range(len(row)))
    common = gcd_many(ans)
    return tuple(x // common for x in ans), tuple(edges)


def main():
    pair_count = 0
    for u in range(1, 101):
        for v in range(u + 1, 101):
            pair_count += 1
            actual_mass, actual_ties = owner_by_cells(u, v)
            actual_max, actual_points = brute_max(u, v)
            expected_max, expected_mass, expected_ties, expected_points = formula(u, v)
            req((actual_mass, actual_ties) == (expected_mass, expected_ties),
                f"owner/ties {(u, v)}")
            req((actual_max, actual_points) == (expected_max, expected_points),
                f"maximum/maximizers {(u, v)}")
            if vtwo(u) != vtwo(v):
                req(decode(directed_summary(u, v)) == Fraction(u // gcd(u, v), v // gcd(u, v)),
                    f"forward label decode {(u, v)}")
                req(decode(directed_summary(v, u)) == Fraction(v // gcd(u, v), u // gcd(u, v)),
                    f"reverse label decode {(u, v)}")

    row_count = 0
    for chosen in combinations(range(1, 19), 5):
        if gcd_many(chosen) != 1 or len({vtwo(x) for x in chosen}) == 1:
            continue
        # Three deterministic label orders, including a root in a nonminimal v2 class.
        for row in (chosen, chosen[::-1], chosen[2:] + chosen[:2]):
            got, edges = recover(row)
            req(got == row, f"row decode {row} -> {got}")
            req(len(edges) == len(row) - 1 and len(set(edges)) == len(edges),
                f"edge duplication {row}: {edges}")
            req(all(vtwo(row[i]) != vtwo(row[j]) for i, j in edges),
                f"blind edge {row}: {edges}")
            row_count += 1

    # Exact common-grid correlation, independent of Fourier summation.
    corr_count = 0
    for r in range(1, 61):
        for s in range(1, 61):
            h = gcd(r, s)
            grid = r * s // h
            total = 0
            for cell in range(grid):
                t = Fraction(2 * cell + 1, 2 * grid)
                er = 1 if (r * t.numerator // t.denominator) % 2 == 0 else -1
                es = 1 if (s * t.numerator // t.denominator) % 2 == 0 else -1
                total += er * es
            actual = Fraction(total, grid)
            expected = Fraction(h * h, r * s) if (r // h) % 2 and (s // h) % 2 else Fraction(0)
            req(actual == expected, f"correlation {(r, s)}: {actual} != {expected}")
            corr_count += 1

    print(f"pairs={pair_count} labelled-row-orders={row_count} correlations={corr_count}: PASS")


if __name__ == "__main__":
    main()
