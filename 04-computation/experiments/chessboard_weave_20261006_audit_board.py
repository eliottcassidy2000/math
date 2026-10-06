#!/usr/bin/env python3
"""Independent audit (written from scratch) of Section 1 of
05-knowledge/results/chessboard_weave_20261006.md:
Prop 1.1 (rings 2-factor, diagonal components, odd/even line partitions, D4 orbits),
Prop 1.2 (bishop/queen/rook reach, affine identity with lambda),
Prop 1.3 (knight ring law on 2m x 2m boards; within-ring structure; square/diamond naming).
Exact integer / Fraction arithmetic only.
"""
from fractions import Fraction as Fr
from math import comb
from collections import Counter, defaultdict

KN = [(1, 2), (2, 1), (-1, 2), (-2, 1), (1, -2), (2, -1), (-1, -2), (-2, -1)]
ORTH = [(1, 0), (-1, 0), (0, 1), (0, -1)]
DIAG = [(1, 1), (1, -1), (-1, 1), (-1, -1)]
fails = []


def check(cond, msg):
    if not cond:
        fails.append(msg)
        print("FAIL:", msg)


def ring(i, j, N):
    r2 = max(abs(2 * i - (N - 1)), abs(2 * j - (N - 1)))  # = 2k+1 for N even
    assert r2 % 2 == 1
    return (r2 - 1) // 2


def edges(N, vecs):
    E = set()
    for i in range(N):
        for j in range(N):
            for dx, dy in vecs:
                a, b = i + dx, j + dy
                if 0 <= a < N and 0 <= b < N:
                    E.add(frozenset([(i, j), (a, b)]))
    return E


def components(V, E):
    adj = defaultdict(set)
    for e in E:
        u, v = tuple(e)
        adj[u].add(v)
        adj[v].add(u)
    seen, comps = set(), []
    for v in V:
        if v in seen:
            continue
        stack, comp = [v], set([v])
        seen.add(v)
        while stack:
            x = stack.pop()
            for y in adj[x]:
                if y not in seen:
                    seen.add(y)
                    comp.add(y)
                    stack.append(y)
        comps.append(comp)
    return comps, adj


def is_single_cycle(Vs, E):
    # induced subgraph on Vs
    Ein = [e for e in E if all(v in Vs for v in e)]
    deg = Counter()
    for e in Ein:
        for v in e:
            deg[v] += 1
    if any(deg[v] != 2 for v in Vs):
        return False, len(Ein)
    comps, _ = components(list(Vs), Ein)
    return len(comps) == 1 and len(Ein) == len(Vs), len(Ein)


# ---------------- Prop 1.1 -----------------
print("=== Prop 1.1 ===")
for m in range(1, 21):
    N = 2 * m
    V = [(i, j) for i in range(N) for j in range(N)]
    R = {v: ring(*v, N) for v in V}
    sizes = Counter(R.values())
    check(all(sizes[k] == 4 * (2 * k + 1) for k in range(m)) and len(sizes) == m, f"ring sizes N={N}")
    G = edges(N, ORTH)
    check(len(G) == 2 * N * (N - 1), f"grid edges N={N}")
    ring_edges = 0
    for k in range(m):
        Vs = set(v for v in V if R[v] == k)
        ok, ne = is_single_cycle(Vs, G)
        check(ok and ne == 4 * (2 * k + 1), f"ring {k} is a cycle C_{4*(2*k+1)} N={N}")
        ring_edges += ne
    spokes = Counter()
    for e in G:
        u, v = tuple(e)
        if R[u] != R[v]:
            spokes[tuple(sorted((R[u], R[v])))] += 1
    check(all(abs(a - b) == 1 for (a, b) in spokes), f"spokes only adjacent rings N={N}")
    check(all(spokes[(k, k + 1)] == 8 * (k + 1) for k in range(m - 1)), f"spoke counts 8(k+1) N={N}")
    check(ring_edges + sum(spokes.values()) == len(G), f"edge split N={N}")
    if N == 8:
        print(f"N=8: grid edges {len(G)}, ring edges {ring_edges}, spokes {sum(spokes.values())}, by pair {dict(spokes)}")
    # (b) diagonal-adjacency (ferz) components = colour classes
    D = edges(N, DIAG)
    comps, _ = components(V, D)
    if N >= 2:
        ok = len(comps) == 2 and all(len(set((i + j) % 2 for (i, j) in c)) == 1 for c in comps)
        check(ok, f"ferz components = colour classes N={N}")
    # line lengths per colour class
    for col in (0, 1):
        anti = Counter()  # i+j = s
        dia = Counter()   # i-j = d
        for (i, j) in V:
            if (i + j) % 2 == col:
                anti[i + j] += 1
                dia[i - j] += 1
        anti_l = [anti[s] for s in sorted(anti)]
        dia_l = [dia[d] for d in sorted(dia)]
        odd_seq = list(range(1, 2 * m, 2)) + list(range(2 * m - 1, 0, -2))
        even_seq = list(range(2, 2 * m + 1, 2)) + list(range(2 * m - 2, 0, -2))
        if col == 0:
            check(anti_l == odd_seq and dia_l == even_seq, f"even class lines N={N}")
        else:
            check(anti_l == even_seq and dia_l == odd_seq, f"odd class lines N={N}")
        if N == 8:
            print(f"N=8 colour {col}: anti-diagonal lengths {anti_l}, diagonal lengths {dia_l}")
    # (c) odd-length lines partition; even-length lines partition
    lines = []
    for s in range(0, 2 * N - 1):
        sq = [(i, s - i) for i in range(N) if 0 <= s - i < N]
        lines.append(('anti', s, sq))
    for d in range(-(N - 1), N):
        sq = [(i, i - d) for i in range(N) if 0 <= i - d < N]
        lines.append(('diag', d, sq))
    odd_lines = [L for L in lines if len(L[2]) % 2 == 1]
    even_lines = [L for L in lines if len(L[2]) % 2 == 0]
    for name, LL in (("odd", odd_lines), ("even", even_lines)):
        cover = Counter()
        for L in LL:
            for v in L[2]:
                cover[v] += 1
        check(all(cover[v] == 1 for v in V), f"{name}-length lines partition N={N}")
    # odd lines: anti-diagonals of colour 0 and diagonals of colour 1
    ok = all(((L[0] == 'anti' and L[1] % 2 == 0) or (L[0] == 'diag' and L[1] % 2 != 0)) for L in odd_lines)
    check(ok, f"odd lines = even-class anti-diagonals + odd-class diagonals N={N}")
    oc = Counter(len(L[2]) for L in odd_lines)
    check(len(odd_lines) == 4 * m and all(oc[2 * k + 1] == 4 for k in range(m)), f"odd line census N={N}")
    check(len(even_lines) == 4 * m - 2, f"even line count N={N}")
    if N == 8:
        print(f"N=8: odd lines {len(odd_lines)} census {dict(oc)}; even lines {len(even_lines)} census {dict(Counter(len(L[2]) for L in even_lines))}")
        len1 = sorted(v for L in odd_lines if len(L[2]) == 1 for v in L[2])
        print("N=8: length-1 lines =", len1, "rings", [R[v] for v in len1])
    # (d) D4 orbits
    def d4(v):
        i, j = v
        n1 = N - 1
        imgs = [(i, j), (j, n1 - i), (n1 - i, n1 - j), (n1 - j, i), (j, i), (n1 - i, j), (n1 - j, n1 - i), (i, n1 - j)]
        return frozenset(imgs)
    orbits = set(d4(v) for v in V)
    check(len(orbits) == m * (m + 1) // 2, f"D4 orbit count N={N}")
    per_ring = defaultdict(list)
    for o in orbits:
        rs = set(R[v] for v in o)
        assert len(rs) == 1
        per_ring[rs.pop()].append(len(o))
    check(all(sorted(per_ring[k]) == [4] + [8] * k for k in range(m)), f"D4 orbits per ring N={N}")
print("Prop 1.1 checked for N = 2..40 (even)")

# ---------------- Prop 1.2 -----------------
print("=== Prop 1.2 ===")
for N in range(2, 41, 2):
    for i in range(N):
        for j in range(N):
            k = ring(i, j, N)
            b = 0
            for dx, dy in DIAG:
                a, c = i + dx, j + dy
                while 0 <= a < N and 0 <= c < N:
                    b += 1
                    a += dx
                    c += dy
            r = 0
            for dx, dy in ORTH:
                a, c = i + dx, j + dy
                while 0 <= a < N and 0 <= c < N:
                    r += 1
                    a += dx
                    c += dy
            check(b == 2 * N - 3 - 2 * k, f"bishop reach N={N} {(i, j)}")
            check(r == 2 * (N - 1), f"rook reach N={N}")
            check(b + r == 4 * N - 5 - 2 * k, f"queen reach N={N}")
            x, y = Fr(2 * i + 1, 2 * N), Fr(2 * j + 1, 2 * N)
            lam = min(min(x, 1 - x), min(y, 1 - y))
            check(lam == Fr(N - 1 - 2 * k, 2 * N), f"lambda centre N={N}")
            check(b == (N - 2) + 2 * N * lam, f"affine identity N={N}")
print("Prop 1.2 + affine identity checked for N = 2..40 (even), all squares")
print("8x8 bishop reach by ring:", sorted(set((ring(i, j, 8), 13 - 2 * ring(i, j, 8)) for i in range(8) for j in range(8))))

# ---------------- Prop 1.3 -----------------
print("=== Prop 1.3 ===")


def knight_profile(N):
    E = edges(N, KN)
    R = {}
    pair = Counter()
    for e in E:
        u, v = tuple(e)
        ru, rv = ring(*u, N), ring(*v, N)
        pair[tuple(sorted((ru, rv)))] += 1
    return E, pair


for m in range(1, 21):
    N = 2 * m
    E, pair = knight_profile(N)
    check(len(E) == 4 * (N - 1) * (N - 2) if N >= 2 else True, f"knight edge total N={N}")
    within = {k: pair.get((k, k), 0) for k in range(m)}
    check(within.get(0, 0) == 0 and all(within[k] == 8 for k in range(1, m)), f"within-ring 8 per ring N={N}")
    adj = sum(c for (a, b), c in pair.items() if b - a == 1)
    two = sum(c for (a, b), c in pair.items() if b - a == 2)
    big = sum(c for (a, b), c in pair.items() if b - a >= 3)
    check(sum(within.values()) == 8 * (m - 1), f"within total N={N}")
    check(adj == 16 * comb(m, 2), f"adjacent total N={N}")
    check(two == 16 * comb(m - 1, 2), f"two-ring total N={N}")
    check(big == 0, f"no 3-ring jumps N={N}")
    check(all(((u[0] + u[1]) - (v[0] + v[1])) % 2 == 1 for e in E for (u, v) in [tuple(e)]), f"colour alternation N={N}")
    if N == 8:
        print("N=8 knight pair profile:", dict(sorted(pair.items())), "total", len(E))
        claimed = {(0, 1): 16, (0, 2): 16, (1, 1): 8, (1, 2): 32, (1, 3): 32, (2, 2): 8, (2, 3): 48, (3, 3): 8}
        check(dict(pair) == claimed, "N=8 pair profile matches note")
print("Prop 1.3 totals checked for N = 2..40 (even)")

# within-ring structure
print("--- within-ring knight graph structure ---")
for m in range(2, 13):
    N = 2 * m
    E, _ = knight_profile(N)
    for k in range(1, m):
        Ek = [e for e in E if all(ring(*v, N) == k for v in e)]
        deg = Counter(v for e in Ek for v in e)
        comps, _ = components(list(deg), Ek)
        if k == 1:
            ok = len(Ek) == 8 and len(deg) == 8 and all(d == 2 for d in deg.values()) and sorted(len(c) for c in comps) == [4, 4]
            check(ok, f"ring 1 = two 4-cycles N={N}")
        else:
            ok = len(Ek) == 8 and len(deg) == 16 and all(d == 1 for d in deg.values())
            check(ok, f"ring {k} = perfect matching on 16 squares N={N}")

# square vs diamond classification of the two ring-1 4-cycles (central 4x4 of 8x8)
N = 8
E, _ = knight_profile(N)
E1 = [tuple(e) for e in E if all(ring(*v, N) == 1 for v in e)]
comps, adj = components(sorted(set(v for e in E1 for v in e)), [frozenset(e) for e in E1])


def cyc_order(comp, adj):
    start = min(comp)
    order = [start]
    prev = None
    cur = start
    while True:
        nxt = [y for y in adj[cur] if y in comp and y != prev]
        nxt = sorted(nxt)[0] if prev is None else [y for y in adj[cur] if y in comp and y != prev][0]
        if nxt == start:
            break
        order.append(nxt)
        prev, cur = cur, nxt
        if len(order) > 10:
            break
    return order


def shape(order):
    vecs = [(order[(t + 1) % 4][0] - order[t][0], order[(t + 1) % 4][1] - order[t][1]) for t in range(4)]
    dots = [vecs[t][0] * vecs[(t + 1) % 4][0] + vecs[t][1] * vecs[(t + 1) % 4][1] for t in range(4)]
    return "square" if all(d == 0 for d in dots) else "diamond(rhombus)"


for comp in comps:
    order = cyc_order(comp, adj)
    loc = [(i - 2, j - 2) for (i, j) in order]
    print("ring-1 within-ring 4-cycle (central-4x4 local coords):", loc, "->", shape(order))
# classical decomposition of a 4x4 block into 4 knight 4-cycles
blk = [(i, j) for i in range(4) for j in range(4)]
E4 = [tuple(e) for e in edges(4, KN)]
cycles4 = []
import itertools
for quad in itertools.combinations(blk, 4):
    qs = set(quad)
    sub = [e for e in E4 if e[0] in qs and e[1] in qs]
    deg = Counter(v for e in sub for v in e)
    if len(sub) == 4 and all(deg[v] == 2 for v in qs):
        c, a = components(list(qs), [frozenset(e) for e in sub])
        if len(c) == 1:
            cycles4.append((quad, shape(cyc_order(qs, a))))
print("all knight 4-cycles of a 4x4 block:")
for q, s in cycles4:
    print("   ", q, s, " corners/centre involved:", any(v in [(0, 0), (0, 3), (3, 0), (3, 3)] for v in q))
print()
print("TOTAL FAILURES:", len(fails))
