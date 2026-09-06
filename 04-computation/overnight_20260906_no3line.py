"""Exact no-three-in-line probes in row/column-saturated grid ensembles.

Universe: every labelled binary n by n matrix with row and column sum two,
2 <= n <= 6. No symmetry or collinearity filter precedes the census.
Two independent grid-line counts and a literal matching path audit controls.
All truth checks remain active under python -O. Standard library only.
"""
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations
from math import comb, factorial, gcd

CHECKS = 0


def require(test, message):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(message)


def line_key(p, q):
    dx, dy = q[0] - p[0], q[1] - p[1]
    d = gcd(dx, dy)
    a, b = dy // d, -dx // d
    if a < 0 or (a == 0 and b < 0):
        a, b = -a, -b
    return a, b, a * p[0] + b * p[1]


def line_sizes(points):
    pairs = Counter(line_key(p, q) for p, q in combinations(points, 2))
    counts = Counter()
    for (a, b, _), c in pairs.items():
        k = 2
        while comb(k, 2) < c:
            k += 1
        require(comb(k, 2) == c, "line pair multiplicity")
        counts[a, b, k] += 1
    return counts


def tuple_count(lines, k, nonaxis=False):
    return sum(count * comb(size, k) for (a, b, size), count in lines.items()
               if size >= k and (not nonaxis or a * b != 0))


def primitive_count(n, k, nonaxis=False):
    total = 0
    directions = [(0, 1)] + [(a, b) for a in range(1, n)
                             for b in range(1 - n, n) if gcd(a, b) == 1]
    for a, b in directions:
        if nonaxis and a * b == 0:
            continue
        for span in range(k - 1, n):
            x, y = n - span * abs(a), n - span * abs(b)
            if min(x, y) <= 0:
                break
            total += comb(span - 1, k - 2) * x * y
    return total


def boards(n):
    counts = [0] * n
    rows = []

    def visit(i):
        if i == n:
            require(all(v == 2 for v in counts), "terminal column degrees")
            yield tuple(rows)
            return
        for a, b in combinations(range(n), 2):
            if counts[a] == 2 or counts[b] == 2:
                continue
            counts[a] += 1
            counts[b] += 1
            rows.append((a, b))
            # Necessary remaining-capacity test, not a geometry filter.
            if all(v + n - i - 1 >= 2 for v in counts):
                yield from visit(i + 1)
            rows.pop()
            counts[a] -= 1
            counts[b] -= 1

    yield from visit(0)


def cycle_profile(rows):
    n = len(rows)
    adj = [[] for _ in range(2 * n)]
    for i, cols in enumerate(rows):
        for j in cols:
            adj[i].append(n + j)
            adj[n + j].append(i)
    seen, parts = set(), []
    for v in range(2 * n):
        if v in seen:
            continue
        stack, component = [v], set()
        while stack:
            u = stack.pop()
            if u in component:
                continue
            component.add(u)
            stack.extend(adj[u])
        seen.update(component)
        require(len(component) % 2 == 0, "bipartite cycle size")
        parts.append(len(component) // 2)
    return tuple(sorted(parts))


def matching_count(points, k):
    return sum(len({x for x, _ in t}) == k == len({y for _, y in t})
               for t in combinations(points, k))


def m3(n):
    return 2 * n * (n - 2) * (2 * n - 5) // 3


def m4(n, four_cycles):
    return n * (n - 3) * (2 * n - 7) * (2 * n - 5) // 6 + four_cycles


def all_boards_formula(n):
    # n!^2 [z^n] exp(-z/2) (1-z)^(-1/2).
    coefficient = sum(F((-1) ** j, 2 ** j * factorial(j)) *
                      F(comb(2 * (n - j), n - j), 4 ** (n - j))
                      for j in range(n + 1))
    value = factorial(n) ** 2 * coefficient
    require(value.denominator == 1, "integer ensemble cardinality")
    return int(value)


def main():
    print("Universe: all labelled row/column-degree-two boards, n=2,...,6")
    print("No geometry/symmetry filters. Counts are not relabelling or D4 orbits.")
    for n in range(2, 13):
        grid = [(x, y) for x in range(n) for y in range(n)]
        lines = line_sizes(grid)
        for k in (3, 4):
            for nonaxis in (False, True):
                require(tuple_count(lines, k, nonaxis) == primitive_count(n, k, nonaxis),
                        ("independent full-grid count", n, k, nonaxis))
    print("Independent pair-line and primitive-endpoint grid counts: n<=12, k=3,4")
    for n in range(2, 7):
        # profile -> [number, sum triples, sum quadruples, sum triples squared, zero]
        stats = defaultdict(lambda: [0] * 5)
        first_zero = {}
        square_orbits = set()
        for rows in boards(n):
            points = [(i, j) for i, cols in enumerate(rows) for j in cols]
            profile = cycle_profile(rows)
            lines = line_sizes(points)
            x, y = tuple_count(lines, 3), tuple_count(lines, 4)
            stat = stats[profile]
            for j, value in enumerate((1, x, y, x * x, int(x == 0))):
                stat[j] += value
            if x == 0:
                first_zero.setdefault(profile, rows)
                images = []
                for swap in (False, True):
                    for rx in (False, True):
                        for ry in (False, True):
                            images.append(tuple(sorted((n-1-u if rx else u, n-1-v if ry else v)
                                                       for u,v in ((y,x) if swap else (x,y)
                                                                   for x,y in points))))
                square_orbits.add(min(images))
            if n <= 5:
                literal = sum((b[0]-a[0])*(c[1]-a[1]) ==
                              (b[1]-a[1])*(c[0]-a[0])
                              for a,b,c in combinations(points, 3))
                require(x == literal, "independent collinearity determinant")
                require(matching_count(points, 3) == m3(n), "three matching law")
                require(matching_count(points, 4) == m4(n, profile.count(2)),
                        "four matching cycle correction")
        total = sum(s[0] for s in stats.values())
        require(total == all_boards_formula(n), "cycle-fugacity ensemble count")
        require(len(square_orbits) == {2:1,3:1,4:4,5:5,6:11}[n],
                "historical square-symmetry positive control")
        t3, t4 = primitive_count(n, 3, True), primitive_count(n, 4, True)
        print("n=", n, "boards=", total, "no3line=", sum(s[4] for s in stats.values()),
              "square_orbits=", len(square_orbits),
              "grid_nonaxis_triples=", t3, "quadruples=", t4)
        for profile, (count, sx, sy, sx2, zeros) in sorted(stats.items()):
            ex, ey = F(sx, count), F(sy, count)
            predicted3 = F(t3 * m3(n), comb(n, 3)**2 * factorial(3)) if n >= 3 else F(0)
            predicted4 = F(t4 * m4(n, profile.count(2)), comb(n, 4)**2 * factorial(4)) if n >= 4 else F(0)
            require(ex == predicted3, "cycle-blind triple mean")
            require(ey == predicted4, "four-cycle-sensitive quadruple mean")
            print("  cycle_row_sizes=", profile, "count=", count, "E[X3]=", str(ex),
                  "E[X4]=", str(ey), "Var[X3]=", str(F(sx2,count)-ex*ex), "zero=", zeros)
        if n == 4:
            require(stats[(2, 2)][4] == 9 and stats[(4,)][4] == 2,
                    "n4 same-mean unequal-zero-probability hostile")
            print("  n4 no3line examples=", first_zero)
            weighted_count = sum(count * 2**len(profile) for profile,(count,*_) in stats.items())
            weighted_zeros = sum(zeros * 2**len(profile) for profile,(*_,zeros) in stats.items())
            require(F(weighted_zeros,weighted_count) != F(sum(s[4] for s in stats.values()),total),
                    "ordered permutation pair sampling bias")
            print("  ordered-pair zero_probability=",str(F(weighted_zeros,weighted_count)),
                  "uniform-board zero_probability=",str(F(sum(s[4] for s in stats.values()),total)))
    print("Active exact checks:", CHECKS)
    print("PASS. No asymptotic independence or universal no-three-in-line bound claimed.")


if __name__ == '__main__':
    main()
