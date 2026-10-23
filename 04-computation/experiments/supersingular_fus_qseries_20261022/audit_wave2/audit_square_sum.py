#!/usr/bin/env python3
"""audit_wave2: independent re-implementation of the square-sum claims of
05-knowledge/results/monster_levels_square_sum_7_21_20261022.md, section 4 and rows D, E of section 5.

Different code paths from square_sum_structure.py: BFS components (not union-find), a window scan for
s(n) and its jumps (not the isqrt formula), the gain/loss coincidences found by intersecting the two event
sets up to 10^7 (not by the Pell test), an exhaustive Hamiltonian-path COUNT (every path, both directions
identified) for n <= 30 instead of a witness search, an exhaustive Hamiltonian-cycle search for n <= 40,
a Held-Karp DP count for H(QR_7), and Aut(QR_7) found by brute force over all 7! permutations.
Run: python3 -B audit_square_sum.py > audit_square_sum.out
"""
from __future__ import annotations

import itertools
import math
import sys
from collections import deque

FAIL = []


def check(cond, msg):
    if not cond:
        FAIL.append(msg)
        print("   !! CHECK FAILED:", msg)


def sq(x):
    r = math.isqrt(x)
    return x >= 0 and r * r == x


def adj(n, S=sq):
    A = {v: [] for v in range(1, n + 1)}
    for a in range(1, n + 1):
        for b in range(a + 1, n + 1):
            if S(a + b):
                A[a].append(b)
                A[b].append(a)
    return A


def comps(n, S=sq):
    A, seen, out = adj(n, S), set(), []
    for v in range(1, n + 1):
        if v in seen:
            continue
        C, dq = [], deque([v])
        seen.add(v)
        while dq:
            u = dq.popleft()
            C.append(u)
            for w in A[u]:
                if w not in seen:
                    seen.add(w)
                    dq.append(w)
        out.append(sorted(C))
    return out


# ---------------------------------------------------------------- 1. components
print("== 1. components ==")
cnt = [len(comps(n)) for n in range(1, 41)]
print("  #components n = 1..40:", cnt)
check(cnt[:14] == [1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 1] and all(c == 1 for c in cnt[13:]), "component counts")
c12 = comps(12)
print("  Q_12:", c12, "sizes", sorted(map(len, c12)))
E12 = [(a, b) for a in range(1, 13) for b in range(a + 1, 13) if sq(a + b)]
print("  squares used in Q_12:", sorted({a + b for a, b in E12}))
for v in (13, 14):
    print(f"  vertex {v}: partners {[(a, a + v) for a in range(1, v) if sq(a + v)]}")
# connectivity for all n >= 14 up to 20000 (incremental: n attaches to k^2 - n)
ok = all(any(sq(a + n) for a in range(1, n)) for n in range(14, 20001))
print("  every n in [14, 20000] has a partner below it (so Q_n connected for n >= 14 by induction):", ok)
check(ok, "connectivity")

# ---------------------------------------------------------------- 2. s(n)
print("\n== 2. s(n), runs, gains/losses, coincidences ==")
s = [0] * 200001
# window scan: s(n) = #{k : n+1 <= k^2 <= 2n-1}
for n in range(1, 200001):
    lo, hi = n + 1, 2 * n - 1
    k0 = math.isqrt(lo - 1) + 1
    c = 0
    k = k0
    while k * k <= hi:
        c += 1
        k += 1
    s[n] = c
brute_ok = all(s[n] == sum(1 for a in range(1, n) if sq(a + n)) for n in range(1, 3001))
form_ok = all(s[n] == math.isqrt(2 * n - 1) - math.isqrt(n) for n in range(1, 200001))
print("  window count = brute force (n <= 3000):", brute_ok, "; = floor(sqrt(2n-1)) - floor(sqrt(n)) (n <= 2*10^5):", form_ok)
check(brute_ok and form_ok, "s(n)")
runs, st = [], 1
for n in range(2, 61):
    if s[n] != s[n - 1]:
        runs.append((st, n - 1, s[n - 1]))
        st = n
runs.append((st, 60, s[60]))
print("  runs to 60:", runs)
gains = {(j * j + 2) // 2 for j in range(2, 4500)}            # ceil((j^2+1)/2)
losses = {k * k for k in range(2, 3200)}
law = all(s[n] - s[n - 1] == (n in gains) - (n in losses) for n in range(2, 200001))
print("  jump law s(n)-s(n-1) = [gain] - [loss] for 2 <= n <= 2*10^5:", law)
check(law, "jump law")
# coincidences up to 10^7 by intersecting the two event sets
G = {(j * j + 2) // 2 for j in range(2, 4473)}               # all gains <= 10^7
Lset = {k * k for k in range(1, 3163)}
co = sorted(x for x in G & Lset if x <= 10 ** 7)
print("  gain & loss coincide (j >= 2) at n <= 10^7:", co, "  k =", [math.isqrt(x) for x in co])
# Pell side: k^2 - 2 t^2 = 1 (t >= 1) or 2 k^2 - m^2 = 1
pk = sorted({k for k in range(1, 3163) for t in [math.isqrt((k * k - 1) // 2)] if t >= 1 and k * k - 2 * t * t == 1}
            | {k for k in range(1, 3163) if sq(2 * k * k - 1)})
print("  Pell k (k^2-2t^2=1, t>=1, or 2k^2-1 square):", pk)
check(co == [k * k for k in pk if k > 1], "Pell characterisation (n = 1 aside)")
pell = [0, 1]
while len(pell) < 14:
    pell.append(2 * pell[-1] + pell[-2])
hcp = [1, 1]
while len(hcp) < 14:
    hcp.append(2 * hcp[-1] + hcp[-2])
print("  Pell numbers P_n:", pell[:12], "; half-companion Pell H_n:", hcp[:12])
print("  coincidence k that are Pell numbers:", [k for k in pk if k in pell], "; that are not:", [k for k in pk if k not in pell])
print("  n = 1 is degenerate: the square 1 never lies in a window [n+1, 2n-1]; gains/losses start at j, k = 2")

# ---------------------------------------------------------------- 3. Hamiltonian paths, exhaustive counts
print("\n== 3. Hamiltonian paths of Q_n, exhaustive count (undirected), n <= 30 ==")


def connected(mask, nb):
    if mask == 0:
        return True
    seen = fr = mask & -mask
    while fr:
        nx, m = 0, fr
        while m:
            b = m & -m
            nx |= nb[b.bit_length() - 1]
            m ^= b
        nx &= mask & ~seen
        seen |= nx
        fr = nx
    return seen == mask


def ham_paths(n, limit=None, collect=False):
    A = adj(n)
    nb = [0] * (n + 1)
    for v in range(1, n + 1):
        for w in A[v]:
            nb[v] |= 1 << w
    full = sum(1 << v for v in range(1, n + 1))
    count = 0
    found = []

    def rec(v, mask, path):
        nonlocal count
        if mask == full:
            count += 1
            if collect:
                found.append(list(path))
            return
        rest = full & ~mask
        # sound prunes: every unvisited vertex needs a neighbour in rest | {v}; a vertex with exactly one
        # such neighbour must be the final end, so at most one of them (when |rest| >= 2); rest connected
        m, ones, multi = rest, 0, rest & (rest - 1)
        while m:
            b = m & -m
            u = b.bit_length() - 1
            c = bin(nb[u] & (rest | (1 << v))).count("1")
            if c == 0:
                return
            if c == 1 and multi:
                ones += 1
                if ones > 1:
                    return
            m ^= b
        if not connected(rest, nb):
            return
        m = nb[v] & rest
        while m:
            b = m & -m
            u = b.bit_length() - 1
            path.append(u)
            rec(u, mask | b, path)
            path.pop()
            m ^= b

    for st in range(1, n + 1):
        rec(st, 1 << st, [st])
    # each undirected path counted twice (once per direction) unless n == 1
    return (count if n == 1 else count // 2), found


H = {}
for n in range(1, 31):
    c, f = ham_paths(n, collect=(n in (15, 16, 17)))
    H[n] = c
    if n in (15, 16, 17):
        uniq = sorted({tuple(p) if p[0] < p[-1] else tuple(reversed(p)) for p in f})
        print(f"  n = {n}: {c} path(s): {[list(u) for u in uniq]}")
print("  counts n = 1..30:", [H[n] for n in range(1, 31)])
has = [n for n in range(1, 31) if H[n] > 0]
print("  n <= 30 with a Hamiltonian path:", has)
check(has == [1, 15, 16, 17, 23] + list(range(25, 31)), "path set")
check([H[n] for n in (15, 16, 17, 23, 25, 26, 27, 28, 29, 30)] == [1, 1, 1, 3, 10, 12, 35, 52, 19, 20], "OEIS A071983 values")
onsets = [n for n in range(2, 31) if H[n] and not H[n - 1]]
lossn = [n for n in range(2, 31) if not H[n] and H[n - 1]]
print("  onsets:", onsets, " s:", [s[n] for n in onsets], "; losses:", lossn)
print("  18's partners below 31:", [a for a in range(1, 31) if a != 18 and sq(a + 18)],
      "; leaves of Q_18:", [v for v, w in adj(18).items() if len(w) == 1])

# ---------------------------------------------------------------- 3b. Hamiltonian cycles n <= 40
print("\n== 3b. Hamiltonian cycles of Q_n, n <= 40 (exhaustive existence search) ==")


def has_cycle(n):
    if n < 3:
        return False
    A = adj(n)
    if any(len(A[v]) < 2 for v in A):
        return False
    nb = [0] * (n + 1)
    for v in range(1, n + 1):
        for w in A[v]:
            nb[v] |= 1 << w
    full = sum(1 << v for v in range(1, n + 1))
    start = min(range(1, n + 1), key=lambda v: len(A[v]))

    def rec(v, mask):
        if mask == full:
            return bool(nb[v] & (1 << start))
        rest = full & ~mask
        m = rest
        while m:
            b = m & -m
            u = b.bit_length() - 1
            # u needs two neighbours among rest | {v, start} (one if it is the only remaining vertex)
            need = 1 if rest == b else 2
            if bin(nb[u] & (rest | (1 << v) | (1 << start))).count("1") < need:
                return False
            m ^= b
        m = nb[v] & rest
        while m:
            b = m & -m
            if rec(b.bit_length() - 1, mask | b):
                return True
            m ^= b
        return False

    return rec(start, 1 << start)


cyc = [n for n in range(1, 41) if has_cycle(n)]
print("  n <= 40 with a Hamiltonian cycle:", cyc)
check(cyc == list(range(32, 41)), "cycles from 32")
mind = {n: min(len(w) for w in adj(n).values()) for n in range(2, 41)}
print("  n <= 40 whose Q_n has a vertex of degree <= 1:", [n for n in range(2, 41) if mind[n] <= 1])

# ---------------------------------------------------------------- 4. induced edges at 25
print("\n== 4. induced edges in Q_25 ==")
A25 = adj(25)
for S in ([10, 6, 3, 1, 8, 17], [12, 4, 5, 11], [2, 7, 18, 9, 16]):
    Ei = sorted((a, b) for a in S for b in A25[a] if b in S and a < b)
    deg = {v: sum(1 for e in Ei if v in e) for v in S}
    print(f"  {sorted(S)}: {len(Ei)} edges {Ei}; degrees {deg}")
check(sorted(len([1 for a in S for b in A25[a] if b in S and a < b]) for S in ([10, 6, 3, 1, 8, 17], [12, 4, 5, 11], [2, 7, 18, 9, 16])) == [3, 4, 5], "edge counts")

# ---------------------------------------------------------------- 5. hostile probe over sum graphs
print("\n== 5. hostile probe: component sizes of sum graphs a + b in S, 3 <= n < 120 ==")


def icbrt(x):
    r = round(x ** (1 / 3))
    while r ** 3 > x:
        r -= 1
    while (r + 1) ** 3 <= x:
        r += 1
    return r


fam = {
    "squares": sq,
    "squares+1": lambda x: sq(x - 1),
    "squares-1": lambda x: sq(x + 1),
    "twice squares": lambda x: x % 2 == 0 and sq(x // 2),
    "triangular": lambda x: sq(8 * x + 1),
    "pentagonal k(3k-1)/2": lambda x: any(k * (3 * k - 1) // 2 == x for k in range(1, 40)),
    "gen. pentagonal": lambda x: any(k * (3 * k - 1) // 2 == x for k in range(-40, 40) if k),
    "cubes": lambda x: icbrt(x) ** 3 == x,
    "primes": lambda x: x > 1 and all(x % d for d in range(2, math.isqrt(x) + 1)),
    "powers of 2": lambda x: x & (x - 1) == 0,
    "odd squares": lambda x: sq(x) and x % 2 == 1,
    "pronic k(k+1)": lambda x: any(k * (k + 1) == x for k in range(1, 20)),
}
pyth = {(a, b, c) for c in range(1, 200) for b in range(1, c) for a in range(1, b + 1) if a * a + b * b == c * c}
for name, S in fam.items():
    exact3, contain = [], []
    for n in range(3, 120):
        sz = sorted(map(len, comps(n, S)))
        if len(sz) == 3 and tuple(sz) in pyth:
            exact3.append(n)
        if any(t in pyth for t in itertools.combinations(sz, 3)):
            contain.append((n, sz if len(sz) <= 6 else f"{len(sz)} comps"))
    print(f"  {name:22s}: exactly three components of Pythagorean sizes at n = {exact3}; "
          f"sizes containing a Pythagorean triple: {contain[:3]}{' ...' if len(contain) > 3 else ''} ({len(contain)})")

# ---------------------------------------------------------------- 6. QR_7
print("\n== 6. Paley tournament QR_7 ==")
R = {1, 2, 4}
arc = [[(v - u) % 7 in R for v in range(7)] for u in range(7)]
# Held-Karp count of Hamiltonian paths
dp = [[0] * 7 for _ in range(1 << 7)]
for v in range(7):
    dp[1 << v][v] = 1
for mask in range(1 << 7):
    for v in range(7):
        if dp[mask][v]:
            for w in range(7):
                if not mask >> w & 1 and arc[v][w]:
                    dp[mask | 1 << w][w] += dp[mask][v]
HQ = sum(dp[(1 << 7) - 1])
print("  H(QR_7) by Held-Karp DP:", HQ)
check(HQ == 189, "H = 189")
aut = [p for p in itertools.permutations(range(7)) if all(arc[u][v] == arc[p[u]][p[v]] for u in range(7) for v in range(7) if u != v)]
affine = sorted(tuple((a * x + b) % 7 for x in range(7)) for a in (1, 2, 4) for b in range(7))
print("  |Aut(QR_7)| by brute force over 7! permutations:", len(aut), "; equals {x -> ax+b : a in {1,2,4}}:", sorted(aut) == affine)
check(len(aut) == 21 and sorted(aut) == affine, "Aut")
paths = [p for p in itertools.permutations(range(7)) if all(arc[p[i]][p[i + 1]] for i in range(6))]
orb = {min(tuple(g[x] for x in p) for g in aut) for p in paths}
free = all(tuple(g[x] for x in p) != p for p in paths for g in aut if list(g) != list(range(7)))
print("  paths:", len(paths), "; free action:", free, "; orbits:", len(orb))
check(free and len(orb) == 9, "free, 9 orbits")
print("  (freeness is automatic: an automorphism fixing a directed Hamiltonian path fixes its vertices in order;")
print("   so |Aut(T)| divides H(T) for every tournament T)")
# check that general fact on random tournaments with nontrivial automorphisms: circulant tournaments on 5, 7, 9 vertices
import random
random.seed(1)
gen_ok = True
for m in (5, 7, 9):
    for Sset in itertools.combinations(range(1, m), (m - 1) // 2):
        if any((m - x) % m in Sset for x in Sset):
            continue
        a2 = [[(v - u) % m in Sset for v in range(m)] for u in range(m)]
        dpm = [[0] * m for _ in range(1 << m)]
        for v in range(m):
            dpm[1 << v][v] = 1
        for mask in range(1 << m):
            for v in range(m):
                if dpm[mask][v]:
                    for w in range(m):
                        if not mask >> w & 1 and a2[v][w]:
                            dpm[mask | 1 << w][w] += dpm[mask][v]
        h = sum(dpm[(1 << m) - 1])
        gen_ok &= h % m == 0
print("  every circulant tournament on 5, 7, 9 vertices has m | H (rotations act freely):", gen_ok)

print("\nfailures:", len(FAIL), FAIL)
if not FAIL:
    print("ALL AUDIT CHECKS PASSED")
print("python", sys.version.split()[0])
