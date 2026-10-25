"""Literal odd actions, literal tournament classes, and signed interpolation."""
from collections import Counter
from fractions import Fraction
from itertools import permutations, product
from math import comb, factorial, gcd, prod
import json


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def shifted(P):
    n = len(P) - 1
    return [sum((-1)**(n + k) * comb(k, j) * P[k]
                for k in range(j, n + 1)) for j in range(n + 1)]


def perm_mask(mask, p):
    return sum(1 << p[i] for i in range(len(p)) if mask >> i & 1)


def subset_poly(group, n):
    unseen = set(range(1 << n))
    P = [0] * (n + 1)
    while unseen:
        x = min(unseen)
        orbit = {perm_mask(x, g) for g in group}
        unseen.difference_update(orbit)
        P[x.bit_count()] += 1
    return P


def action_cases():
    actions = [("empty", [()], 0), ("trivial6", [tuple(range(6))], 6)]
    for n in (3, 5, 9):
        actions.append((f"C{n}", [tuple((i + a) % n for i in range(n))
                                  for a in range(n)], n))
    points = list(product(range(3), repeat=2))
    actions.append(("C3_squared", [tuple(points.index(tuple((v[i] + a[i]) % 3
                        for i in range(2))) for v in points) for a in points], 9))
    actions.append(("F21", [tuple((a * x + b) % 7 for x in range(7))
                            for a in (1, 2, 4) for b in range(7)], 7))
    actions.append(("C15_two_fixed", [tuple([(i+a)%3 for i in range(3)]
        + [3+(i+b)%5 for i in range(5)] + [8, 9])
        for a in range(3) for b in range(5)], 10))
    rows = []
    for name, group, n in actions:
        pool = set(group)
        need(len(pool) % 2 == 1, "not odd group")
        need(all(tuple(g[h[i]] for i in range(n)) in pool for g in group for h in group),
             "action not closed")
        P = subset_poly(group, n)
        Q = shifted(P)
        difference = [comb(n, i) - P[i] for i in range(n+1)]
        need(all(x >= 0 and x % 2 == 0 for x in difference), "positive half correction")
        need([x % 2 for x in Q] == [0] * n + [1], "universal parity jet")
        need(all(x >= 0 for x in Q), "positive jet")
        orbit_rank = len({frozenset(g[i] for g in group) for i in range(n)})
        for F in range(1 << n):
            invariant = all(perm_mask(F, g) == F for g in group)
            if not invariant:
                g = next(g for g in group if perm_mask(F, g) != F)
                need((F ^ perm_mask(F, g)) != 0, "noninvariant toggle descent hostile")
        rows.append({"action": name, "gaps": n, "group_order": len(group),
                     "toggle_rank": orbit_rank, "P": P, "Q": Q})
    need(subset_poly([(0, 1), (1, 0)], 2) == [1, 1, 1], "even action hostile")
    # Enumerate both candidate images of a C5 generator on two points, then
    # all maps from its regular five-point orbit. Positivity alone would
    # supply two abstract points; this tests the actual quotient operation.
    surjections = 0
    for generator in permutations(range(2)):
        fifth = tuple(range(2))
        for _ in range(5):
            fifth = tuple(generator[fifth[i]] for i in range(2))
        if fifth != (0, 1):
            continue
        for f in product(range(2), repeat=5):
            if all(f[(i+1)%5] == generator[f[i]] for i in range(5)):
                surjections += len(set(f)) == 2
    need(surjections == 0, "C5 positive correction has no equivariant surjection")
    return rows


def odd_parts(n, least=1):
    if n == 0:
        yield ()
    for d in range(least, n+1, 2):
        for rest in odd_parts(n-d, d):
            yield (d,) + rest


def signed_rows(N):
    b = []
    for M in range(N+1):
        E = N-M
        value = Fraction(0)
        for parts in odd_parts(M):
            counts = Counter(parts)
            exponent = sum((d-1)//2 for d in parts)
            exponent += sum(gcd(a, b) for i, a in enumerate(parts) for b in parts[i+1:])
            denominator = prod(d**k * factorial(k) for d, k in counts.items())
            excess = [1] + [0] * E
            for d, k in counts.items():
                nxt = [0] * (E+1)
                for e in range(E//d+1):
                    beta = (-1)**e * comb(k+e-1, e)
                    for j in range(E-d*e+1):
                        nxt[j+d*e] += excess[j] * beta
                excess = nxt
            value += Fraction((1 << exponent) * excess[E], denominator)
        need(value.denominator == 1, ("signed coefficient not integral", N, M))
        b.append(value.numerator)
    return b


def interpolation(N):
    b = signed_rows(N)
    R = [sum(b[M] * comb(M, j) for M in range(j, N+1)) for j in range(N+1)]
    need(all(x >= 0 for x in R), ("interpolation positivity", N))
    expected = [0] * (N+1)
    if N == 0:
        expected[0] = 1
    elif N % 2:
        expected[0] = expected[1] = 1
    else:
        expected[1] = expected[2] = 1
    need([x % 2 for x in R] == expected, ("interpolation parity", N))
    need(R[0] % 2 == (N == 0 or N % 2 == 1), "pair-free parity")
    return R


def literal_tournaments(n):
    pairs = [(a, b) for a in range(n) for b in range(a+1, n)]
    index = {p: i for i, p in enumerate(pairs)}
    perms = list(permutations(range(n)))
    def arc(mask, a, b):
        return bool(mask >> index[min(a,b), max(a,b)] & 1) == (a < b)
    def transform(mask, perm):
        out = 0
        for a, b in pairs:
            x, y = perm[a], perm[b]
            if arc(mask, a, b) == (x < y):
                out |= 1 << index[min(x,y), max(x,y)]
        return out
    unseen = set(range(1 << len(pairs)))
    total = [0] * (n+1)
    classes = core_count = 0
    while unseen:
        mask = min(unseen)
        orbit = {transform(mask, p) for p in perms}
        unseen.difference_update(orbit)
        aut = [p for p in perms if transform(mask, p) == mask]
        need(len(aut) % 2 == 1, "tournament automorphisms not odd")
        gaps = [(a,b) for a,b in pairs
                if all(arc(mask,a,k) == arc(mask,b,k) for k in range(n) if k not in (a,b))]
        adj = {i: set() for i in range(n)}
        for a,b in gaps:
            adj[a].add(b); adj[b].add(a)
        remaining = set(range(n)); q = 0
        while remaining:
            stack = [remaining.pop()]; q += 1
            while stack:
                for v in adj[stack.pop()] & remaining:
                    remaining.remove(v); stack.append(v)
        need(len(gaps) == n-q, "homogeneous-pair forest not path forest")
        gap_index = {e:i for i,e in enumerate(gaps)}
        action = [tuple(gap_index[tuple(sorted((p[a],p[b])))] for a,b in gaps) for p in aut]
        P = subset_poly(action, len(gaps))
        term = mul([comb(q,j) for j in range(q+1)], shifted(P))
        for j,v in enumerate(term):
            total[j] += v
        classes += 1
        core_count += not gaps
    need(total == interpolation(n), ("literal class interpolation", n, total))
    return {"n": n, "classes": classes, "pair_free": core_count, "R": total}


def main():
    actions = action_cases()
    rows = [interpolation(n) for n in range(21)]
    literal = [literal_tournaments(n) for n in range(7)]
    need([x['classes'] for x in literal] == [1,1,1,2,4,12,56], "class prefix")
    need([x['pair_free'] for x in literal] == [1,1,0,1,2,5,26], "core prefix")
    need(1+3*literal[3]['pair_free']+literal[4]['pair_free'] == 6 != literal[4]['classes'],
         "integer inflation-transform hostile")
    print(json.dumps({"status":"PASS", "action_cases":actions,
        "signed_interpolation_range":[0,20], "literal_tournaments":literal,
        "pair_free_counts":[r[0] for r in rows],
        "parity_rows":[[v%2 for v in r] for r in rows]}, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
