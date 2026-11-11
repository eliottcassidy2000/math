#!/usr/bin/env python3
"""S92 (README sections 5, 6): small tournaments, their arc-usage profiles, and the Collatz maps on the lines of the
path spectrum.  Everything here is exact and elementary; nothing is a theorem about Chan-Pak or Collatz.

Part 1 (tournaments of order <= 5).  For every isomorphism class: H, number of 3-cycles, number of directed cycles,
Pi_2 (ordered pairs of possibly empty disjoint paths covering T), |Aut|, self-converse or not, strong or not, and
the arc-usage profile W = multiset over the arcs of the number of Hamiltonian paths through the arc (invariant
under converse; the "profile" of the pasted note and W(u, v) of Session265).  Checks: order 4 has profiles
{1,1,1}, {2,2,2,1,1,1} (twice, a converse pair), {4,3,3,3,1,1}; the numbers of directed cycles are 0, 1, 1, 3;
Pi_2 = 8 (2 + c_3) = 16, 24, 24, 32; order 5 has 12 classes, 6 strong; the generating function of the classes by
strong components, B(x) = 1/(1 - A(x)), A = x + x^3 + x^4 + 6 x^5, gives 1, 1, 2, 4, 12 through order 5.

Part 2 (lines).  For odd N > 1 put c_minus(N) = odd part of N - 1 (N is on the seed line 1 + c 2^s of THM-5126)
and c_plus(N) = odd part of N + 1 (N is on the anti-line c 2^s - 1 of THM-5130).  For every odd n:
        c_minus(3n) = T_minus(n) = odd part of 3n - 1,        c_plus(3n) = T_plus(n) = odd part of 3n + 1.
So the ordered join with a 3-cycle (H -> 3H) moves a path count from its line to the line labelled by the image of
the count under the Syracuse map of 3n - 1, and from its anti-line to the anti-line labelled by its image under
the Syracuse map of 3n + 1.  Checked for n < 10^6.  Printed: the three cycles of T_minus as cycles of lines, the
transfer points 3n, and where the two non-values 7 and 21 sit.
Exit status 1 on failure.
"""
import itertools
import sys
from collections import Counter, defaultdict


def tournaments(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for bits in range(1 << len(pairs)):
        adj = [[0] * n for _ in range(n)]
        for k, (i, j) in enumerate(pairs):
            if (bits >> k) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
        yield adj


def canon(adj, n):
    best = None
    for p in itertools.permutations(range(n)):
        key = tuple(adj[p[i]][p[j]] for i in range(n) for j in range(n))
        if best is None or key < best:
            best = key
    return best


def ham_path_list(adj, n):
    out = []
    for p in itertools.permutations(range(n)):
        if all(adj[p[i]][p[i + 1]] for i in range(n - 1)):
            out.append(p)
    return out


def count_cycles(adj, n):
    tot = 0
    c3 = 0
    for k in range(3, n + 1):
        for sub in itertools.combinations(range(n), k):
            first = sub[0]
            for rest in itertools.permutations(sub[1:]):
                cyc = (first,) + rest
                if all(adj[cyc[i]][cyc[(i + 1) % k]] for i in range(k)):
                    tot += 1
                    if k == 3:
                        c3 += 1
    return c3, tot


def ham(adj, verts):
    m = len(verts)
    if m == 0:
        return 1
    return sum(1 for p in itertools.permutations(verts) if all(adj[p[i]][p[i + 1]] for i in range(m - 1)))


def pi2(adj, n):
    tot = 0
    for mask in range(1 << n):
        a = [v for v in range(n) if (mask >> v) & 1]
        b = [v for v in range(n) if not (mask >> v) & 1]
        tot += ham(adj, a) * ham(adj, b)
    return tot


def strong(adj, n):
    def reach(src, fwd):
        seen = {src}
        stack = [src]
        while stack:
            v = stack.pop()
            for w in range(n):
                if w not in seen and (adj[v][w] if fwd else adj[w][v]):
                    seen.add(w)
                    stack.append(w)
        return len(seen) == n
    return reach(0, True) and reach(0, False)


def odd_part(x):
    while x % 2 == 0:
        x //= 2
    return x


def main():
    ok = True
    table = {}
    for n in (1, 2, 3, 4, 5):
        classes = {}
        for adj in tournaments(n):
            key = canon(adj, n)
            if key not in classes:
                classes[key] = [adj, 0]
            classes[key][1] += 1
        rows = []
        for key, (adj, labelled) in classes.items():
            paths = ham_path_list(adj, n)
            use = Counter()
            for p in paths:
                for i in range(n - 1):
                    use[frozenset((p[i], p[i + 1]))] += 1
            prof = tuple(sorted(use.values(), reverse=True))
            conv = [[adj[j][i] for j in range(n)] for i in range(n)]
            c3, cyc = count_cycles(adj, n)
            fact = 1
            for k in range(2, n + 1):
                fact *= k
            rows.append(dict(H=len(paths), c3=c3, cycles=cyc, pi2=pi2(adj, n), aut=fact // labelled,
                             selfconv=canon(conv, n) == key, strong=strong(adj, n) if n > 1 else True, profile=prof))
        rows.sort(key=lambda r: (r["H"], r["cycles"], r["profile"]))
        table[n] = rows
        print(f"order {n}: {len(rows)} classes, {sum(r['strong'] for r in rows)} strong, "
              f"{sum(r['selfconv'] for r in rows)} self-converse")
        for r in rows:
            print(f"   H = {r['H']:2d}  3-cycles {r['c3']}  directed cycles {r['cycles']:2d}  Pi_2 = {r['pi2']:3d}  |Aut| = {r['aut']}"
                  f"  {'strong' if r['strong'] else 'not strong'}  {'self-converse' if r['selfconv'] else 'converse pair'}"
                  f"  arc usage {list(r['profile'])}")
    r4 = table[4]
    ok &= [r["H"] for r in r4] == [1, 3, 3, 5]
    ok &= [list(r["profile"]) for r in r4] == [[1, 1, 1], [2, 2, 2, 1, 1, 1], [2, 2, 2, 1, 1, 1], [4, 3, 3, 3, 1, 1]]
    ok &= [r["cycles"] for r in r4] == [0, 1, 1, 3]
    ok &= [r["pi2"] for r in r4] == [16, 24, 24, 32] and all(r["pi2"] == 8 * (2 + r["c3"]) for r in r4)
    ok &= len(table[5]) == 12 and sum(r["strong"] for r in table[5]) == 6
    # generating function 1/(1 - A)
    A = [0, 1, 0, 1, 1, 6]
    B = [1] + [0] * 5
    for n in range(1, 6):
        B[n] = sum(A[k] * B[n - k] for k in range(1, n + 1))
    print(f"classes by strong components: 1/(1 - (x + x^3 + x^4 + 6x^5)) = {B} (orders 0..5); counted: "
          f"{[1] + [len(table[n]) for n in range(1, 6)]}")
    ok &= B == [1] + [len(table[n]) for n in range(1, 6)]
    s5 = [r for r in table[5] if r["strong"]]
    print(f"the six strong classes of order 5: H = {[r['H'] for r in s5]}, self-converse: "
          f"{[r['selfconv'] for r in s5]}, |Aut| = {[r['aut'] for r in s5]}")

    # Part 2
    def tm(n):
        return odd_part(3 * n - 1)

    def tp(n):
        return odd_part(3 * n + 1)
    good = all(odd_part(3 * n - 1) == tm(n) and odd_part(3 * n + 1) == tp(n) for n in range(1, 10 ** 6, 2))
    print(f"c_minus(3n) = T_minus(n) and c_plus(3n) = T_plus(n) for odd n < 10^6: {good}")
    ok &= good
    for start in (1, 5, 17):
        cyc = [start]
        while tm(cyc[-1]) != start:
            cyc.append(tm(cyc[-1]))
        steps = []
        for n in cyc:
            N = 3 * n
            s = 0
            x = N - 1
            while x % 2 == 0:
                x //= 2
                s += 1
            steps.append(f"3*{n} = {N} = 1 + {x}*2^{s}")
        print(f"cycle of 3n - 1 through {start}: lines L_c, c = {cyc}; transfers: " + "; ".join(steps))
    print("non-values: 7 = 1 + 3*2 (line L_3, whose label has no preimage under T_minus: no point of L_3 is a multiple of 3)"
          " = 1*2^3 - 1 (anti-line of the fixed point 1 of T_plus);")
    print("            21 = 3*7 = 1 + 5*2^2 (line L_5, 5 = T_minus(7): the transfer point of the 2-cycle {5, 7}) "
          "= 11*2 - 1 (anti-line L'_11, 11 = T_plus(7))")
    ok &= tm(7) == 5 and tm(5) == 7 and tp(7) == 11 and odd_part(7 - 1) == 3 and odd_part(21 - 1) == 5
    # basins of the three cycles among the odd numbers
    cyc_of = {}
    for start in (1, 5, 17):
        n = start
        while True:
            cyc_of[n] = start
            n = tm(n)
            if n == start:
                break
    cnt = Counter()
    for n in range(1, 200001, 2):
        x = n
        seen = []
        while x not in cyc_of:
            seen.append(x)
            x = tm(x)
        cnt[cyc_of[x]] += 1
    print(f"basins of the cycles of 3n - 1 among the odd n < 200000 (NUMERICAL): "
          f"{ {k: round(v / 100000, 4) for k, v in sorted(cnt.items())} }")
    print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
