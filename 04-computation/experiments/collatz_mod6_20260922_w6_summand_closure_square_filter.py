#!/usr/bin/env python3
"""collatz_mod6_20260922_w6_summand_closure_square_filter.py

Lane summand_closure_square_filter (wave 6, 2026-09-22, session collatz-mod6-20260917).

Blocks
  A  the summand closure of THM-2422: holes {1,4,6}, asynchronous vs synchronous
     closure, M_t = 27*2^(t-4)+1, the adjoin-1/4/6 corrections (finite re-check).
  B  the square-sum graph Q_n: exact components for n <= 13, the three chains as
     a linear forest, merge squares at 13 and 14, connectivity 14..2000.
  C  incremental build: degree of the newest vertex, founders (born isolated),
     born-leaf set, leaf law, leaf sets L(n), no leaves iff n >= 31.
  D  self-similarity: the 2m-square leaf window (m = 2j^2), finite scar j <= 3;
     the square count in (2j^2, 4j^2) grows; no dyadic law.
  E  Hamiltonian paths for n <= 32 (existence + counts vs OEIS A090460),
     leaf obstruction at 18, forced edges at 15, cycles at 31/32.
  F  analogy audit numbers (Collatz odd arrow as summand arrow).

Every numeric claim of the note is printed here.  python3 -O gives identical
output modulo the timing line.
"""
import math
import sys
import time

T0 = time.time()
sys.setrecursionlimit(10000)


def out(*a):
    print(*a)
    sys.stdout.flush()


def isqrt(x):
    return math.isqrt(x)


def is_square(x):
    return x >= 0 and isqrt(x) ** 2 == x


# ---------------------------------------------------------------- block A
out("=== A. THM-2422 summand closure, finite re-check (cited by path, not re-proved) ===")


def closure_strict(seeds, bound):
    S = set(seeds)
    changed = True
    while changed:
        changed = False
        L = sorted(S)
        for i in range(len(L)):
            for j in range(i + 1, len(L)):
                z = L[i] + L[j]
                if z <= bound and z not in S:
                    S.add(z)
                    changed = True
    return S


B = 200
cl = closure_strict({2, 3}, B)
holes = sorted(set(range(1, B + 1)) - cl)
out("A1 Cl^<({2,3}) cap [1,%d] misses exactly %s" % (B, holes))
if holes != [1, 4, 6]:
    raise RuntimeError("closure holes changed")
for seeds, expect in (({2, 3, 1}, []), ({2, 3, 4}, [1]), ({2, 3, 6}, [1, 4])):
    h = sorted(set(range(1, B + 1)) - closure_strict(seeds, B))
    out("A2 seeds %s -> holes %s" % (sorted(seeds), h))
    if h != expect:
        raise RuntimeError("adjoin correction failed")

# synchronous closure S_t (19)-(26)
S = {2, 3}
out("A3 synchronous S_t: t, |S_t|, max, gaps in [7,max], M_t formula check")
for t in range(0, 10):
    mx = max(S)
    gaps = [z for z in range(7, mx + 1) if z not in S]
    formula = 27 * 2 ** (t - 4) + 1 if t >= 4 else None
    ok = (formula == mx) if t >= 4 else True
    out("   t=%d |S_t|=%d max=%d gaps=%s M_t=%s ok=%s" % (t, len(S), mx, gaps[:6], formula, ok))
    if t >= 4 and (not ok or gaps or S != {2, 3, 5} | set(range(7, mx + 1))):
        raise RuntimeError("synchronous law failed at t=%d" % t)
    S = S | {a + b for a in S for b in S if a < b}

# ---------------------------------------------------------------- block B
out()
out("=== B. Square-sum graph Q_n: components, chains, merges, connectivity ===")


def nbrs(m, n):
    """neighbours of m in Q_n (x~y iff x+y square, x != y)."""
    res = []
    k = isqrt(m + 1)
    if k * k < m + 1:
        k += 1
    while k * k <= m + n:
        y = k * k - m
        if 1 <= y <= n and y != m:
            res.append(y)
        k += 1
    return res


def components(n):
    seen = [False] * (n + 1)
    comps = []
    for s in range(1, n + 1):
        if seen[s]:
            continue
        stack = [s]
        seen[s] = True
        comp = []
        while stack:
            v = stack.pop()
            comp.append(v)
            for w in nbrs(v, n):
                if not seen[w]:
                    seen[w] = True
                    stack.append(w)
        comps.append(sorted(comp))
    return sorted(comps)


def edges(n):
    return sorted((x, y) for x in range(1, n + 1) for y in nbrs(x, n) if x < y)


out("B1 exact components of Q_n, n = 1..14 (with max degree and edge list)")
ccount = {}
for n in range(1, 15):
    cs = components(n)
    ccount[n] = len(cs)
    md = max(len(nbrs(m, n)) for m in range(1, n + 1))
    out("   n=%2d c=%d maxdeg=%d comps=%s" % (n, len(cs), md, cs))
    if n == 14:
        out("      edges of Q_14 with their squares: %s" % [(x, y, x + y) for x, y in edges(14)])
    if n in (12, 13):
        out("      edges of Q_%d with their squares: %s" % (n, [(x, y, x + y) for x, y in edges(n)]))

expect_three = {n: 3 for n in range(4, 13)}
for n, c in expect_three.items():
    if ccount[n] != c:
        raise RuntimeError("component count at n=%d is %d" % (n, ccount[n]))
if ccount[13] != 2 or ccount[14] != 1:
    raise RuntimeError("merge events wrong")

A12 = components(12)
out("B2 the three chains of Q_12 (named by founder = least element): %s" % A12)
out("   founders (least elements): %s" % [c[0] for c in A12])
# linear forest test: max degree <= 2 iff n <= 12
lf = [n for n in range(1, 41) if max(len(nbrs(m, n)) for m in range(1, n + 1)) <= 2]
out("B3 n <= 40 with max degree <= 2 (Q_n a linear forest = union of chains): %s" % lf)
if lf != list(range(1, 13)):
    raise RuntimeError("linear forest range wrong")
out("   degree-3 vertices of Q_13: %s" % [(m, nbrs(m, 13)) for m in range(1, 14) if len(nbrs(m, 13)) == 3])
out("B4 merge events: vertex 13 attaches to %s (squares %s); vertex 14 attaches to %s (squares %s)"
    % (nbrs(13, 13), [13 + y for y in nbrs(13, 13)], nbrs(14, 14), [14 + y for y in nbrs(14, 14)]))
# component containing 2 over n
for n in (12, 13, 14):
    comp2 = [c for c in components(n) if 2 in c][0]
    out("   component of 2 in Q_%d: %s" % (n, comp2))

# connectivity 14..2000 (union-find incremental)
NMAX = 2000
parent = list(range(NMAX + 1))


def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x


ncomp = 0
cseq = []
merges = []
for n in range(1, NMAX + 1):
    ncomp += 1
    for y in nbrs(n, n):
        a, b = find(n), find(y)
        if a != b:
            parent[a] = b
            ncomp -= 1
            merges.append((n, y, n + y))
    cseq.append(ncomp)
out("B5 c(Q_n) for n=1..20: %s" % cseq[:20])
bad = [n for n in range(14, NMAX + 1) if cseq[n - 1] != 1]
out("   n in [14,%d] with c(Q_n) != 1: %s  (FINITE-EXACT connectivity)" % (NMAX, bad))
if bad:
    raise RuntimeError("disconnected Q_n found")
out("   union-find merge events n<=14 (new vertex, old vertex, square): %s" % [m for m in merges if m[0] <= 14])
# increases of c: births in isolation
iso_birth = [n for n in range(1, NMAX + 1) if (cseq[n - 1] - (cseq[n - 2] if n >= 2 else 0)) == 1]
out("B6 vertices born isolated (c increases): %s" % iso_birth)
if iso_birth != [1, 2, 4]:
    raise RuntimeError("founders are not {1,2,4}")

# ---------------------------------------------------------------- block C
out()
out("=== C. Incremental build: degrees, founders, leaves ===")


def deg_new(n):
    """degree of the newest vertex n in Q_n = #squares in (n, 2n) = #{k: n < k^2 <= 2n-1}."""
    return isqrt(2 * n - 1) - isqrt(n)


def deg_in(m, n):
    """degree of m in Q_n: squares in (m, m+n] minus the diagonal 2m."""
    return isqrt(m + n) - isqrt(m) - (1 if is_square(2 * m) else 0)


# verify formulas against brute force
for n in range(1, 400):
    for m in range(1, n + 1):
        if deg_in(m, n) != len(nbrs(m, n)):
            raise RuntimeError("deg formula fails at m=%d n=%d" % (m, n))
    if deg_new(n) != len(nbrs(n, n)):
        raise RuntimeError("deg_new fails at n=%d" % n)
out("C1 formulas deg_n(m) = floor(sqrt(m+n)) - floor(sqrt(m)) - [2m square], deg_n(n) = floor(sqrt(2n-1)) - floor(sqrt(n)): verified for all m <= n < 400")
out("C2 deg_n(n) for n=1..40: %s" % [deg_new(n) for n in range(1, 41)])
born_iso = [n for n in range(1, NMAX + 1) if deg_new(n) == 0]
born_leaf = [n for n in range(1, NMAX + 1) if deg_new(n) == 1]
born_two = [n for n in range(1, NMAX + 1) if deg_new(n) == 2]
out("C3 born isolated (deg_n(n)=0): %s" % born_iso)
out("C4 born as a leaf (deg_n(n)=1): %s (count %d, max %d)" % (born_leaf, len(born_leaf), max(born_leaf)))
out("C5 born with degree 2: count %d, max %d, first %s ... last %s" % (len(born_two), max(born_two), born_two[:6], born_two[-4:]))
# bound check: (sqrt(n)+1)^2 < 2n iff n >= 6
first_always = min(n for n in range(1, 100) if all((isqrt(k) + 1) ** 2 < 2 * k for k in range(n, 100)))
out("C6 smallest n0 with (floor(sqrt(k))+1)^2 < 2k for all n0 <= k < 100: %d" % first_always)
# leaves
leafsets = {}
for n in range(1, NMAX + 1):
    leafsets[n] = [m for m in range(1, n + 1) if deg_in(m, n) == 1]
out("C7 leaf sets L(n) for n=1..32:")
for n in range(1, 33):
    iso = [m for m in range(1, n + 1) if deg_in(m, n) == 0]
    out("   n=%2d leaves=%s isolated=%s" % (n, leafsets[n], iso))
with_leaves = [n for n in range(1, NMAX + 1) if leafsets[n]]
out("C8 n <= %d having at least one leaf: %s" % (NMAX, with_leaves))
if with_leaves != [n for n in range(3, 31)]:
    raise RuntimeError("leaf range not [3,30]")
out("   => Q_n has a leaf iff 3 <= n <= 30 (FINITE-EXACT to %d; PROVED for all n >= 31 in the note)" % NMAX)
# leaf windows per vertex m: for each m, the n-interval [m, ...] where m is a leaf
out("C9 leaf windows: m -> set of n >= m with m a leaf of Q_n (n <= 60)")
for m in range(1, 31):
    w = [n for n in range(m, 61) if deg_in(m, n) == 1]
    if w:
        out("   m=%2d leaf for n in [%d,%d] (%d values), 2m square=%s, squares>m: %s" %
            (m, w[0], w[-1], len(w), is_square(2 * m), [k * k for k in range(isqrt(m) + 1, isqrt(m) + 5)]))
# session lead probe: degree<=1 vertices n=18..32
out("C10 session lead probe: degree<=1 vertices for n=18..32: %s" %
    [(n, [m for m in range(1, n + 1) if deg_in(m, n) <= 1]) for n in range(18, 33)])
# proof helper for n >= 31: n in [24,52] and m=2i^2 <= n direct check
out("C11 proof helper: for 24 <= n <= 60 and every m <= n with 2m square, deg_n(m) = %s" %
    sorted(set(deg_in(m, n) for n in range(31, 61) for m in range(1, n + 1) if is_square(2 * m))))
out("   for 24 <= n <= 60 and every m <= n with 2m NOT square, min deg_n(m) = %d" %
    min(deg_in(m, n) for n in range(24, 61) for m in range(1, n + 1) if not is_square(2 * m)))
out("   n0 with 4*sqrt(n)+4 <= n for all n >= n0: %d" % min(n for n in range(1, 200) if all(4 * math.sqrt(k) + 4 <= k for k in range(n, 200))))
out("   n0 with 6*sqrt(n)+9 <= n for all n >= n0: %d" % min(n for n in range(1, 200) if all(6 * math.sqrt(k) + 9 <= k for k in range(n, 200))))
out("   m=18: squares in (18, 18+n] for n=30: %s; n=31: %s" %
    ([k * k for k in range(5, 8) if k * k <= 48], [k * k for k in range(5, 8) if k * k <= 49]))

# ---------------------------------------------------------------- block D
out()
out("=== D. Self-similarity: the 2m-square scar vs THM-2422's dyadic law ===")
out("D1 THM-2422 frontier M_t = 27*2^(t-4)+1 for t=4..12: %s" % [27 * 2 ** (t - 4) + 1 for t in range(4, 13)])
out("D2 m = 2j^2: squares strictly between 2j^2 and 4j^2, and the leaf window of m")
for j in range(1, 13):
    m = 2 * j * j
    between = [k * k for k in range(isqrt(m) + 1, 2 * j) if m < k * k < 4 * j * j]
    w = [n for n in range(m, 6 * m + 60) if deg_in(m, n) == 1]
    out("   j=%2d m=%4d #squares in (2j^2,4j^2)=%d %s leaf window=%s" %
        (j, m, len(between), between[:4], ([w[0], w[-1]] if w else None)))
# the count of squares in (2j^2, 4j^2) = 2j-1-floor(j*sqrt2)
out("D3 #squares in (2j^2,4j^2) = 2j-1-floor(j*sqrt(2)); j=1..12: %s" %
    [2 * j - 1 - isqrt(2 * j * j) for j in range(1, 13)])
# degrees at n=2000 grow like (sqrt2-1)sqrt(n)
out("D4 deg_n(n)/sqrt(n) at n=100,400,900,1600,2000: %s" %
    ["%.4f" % (deg_new(n) / math.sqrt(n)) for n in (100, 400, 900, 1600, 2000)])
out("   sqrt(2)-1 = %.4f" % (math.sqrt(2) - 1))
# composition of two reflections x -> k^2-x, x -> l^2-x is translation by l^2-k^2
diffs = sorted(set((l * l - k * k) for k in range(1, 40) for l in range(k + 1, 40) if l * l - k * k <= 60))
out("D5 two-step translations l^2-k^2 <= 60: %s" % diffs)
out("   residues mod 4 of these: %s ; missing residue: %s" % (sorted(set(d % 4 for d in diffs)), sorted({0, 1, 2, 3} - set(d % 4 for d in diffs))))
# no dyadic repetition: the component count is 1 for all n >= 14 and the born-leaf set is finite
out("D6 born-leaf set is finite: max %d; born-isolated set {1,2,4}; c(Q_n)=1 for 14<=n<=%d: no self-similar repetition of the three-chain phase" % (max(born_leaf), NMAX))

# ---------------------------------------------------------------- block E
out()
out("=== E. Hamiltonian paths, leaf obstruction, forced edges (n <= 32) ===")


def ham_paths_count(n, count_all=True, budget=240.0):
    """number of Hamiltonian paths of Q_n up to reversal (count_all) or existence (first found)."""
    adj = [0] * (n + 1)
    for m in range(1, n + 1):
        for y in nbrs(m, n):
            adj[m] |= 1 << y
    full = ((1 << (n + 1)) - 2)
    t0 = time.time()
    found = [0]
    first = [None]

    def dfs(v, visited, path):
        if visited == full:
            found[0] += 1
            if first[0] is None:
                first[0] = list(path)
            return not count_all
        if time.time() - t0 > budget:
            raise TimeoutError
        rem = full & ~visited
        # prune: every unvisited vertex needs an unvisited-or-current neighbour;
        # at most one unvisited vertex may have exactly one such neighbour (the last vertex)
        ones = 0
        r = rem
        while r:
            lb = r & -r
            u = lb.bit_length() - 1
            r ^= lb
            d = bin(adj[u] & (rem | (1 << v))).count("1")
            if d == 0:
                return False
            if d == 1:
                ones += 1
                if ones > 1 and rem != lb:
                    return False
        cand = adj[v] & rem
        while cand:
            lb = cand & -cand
            w = lb.bit_length() - 1
            cand ^= lb
            path.append(w)
            if dfs(w, visited | lb, path):
                return True
            path.pop()
        return False

    for s in range(1, n + 1):
        if dfs(s, 1 << s, [s]):
            break
    return found[0] // 2 if count_all else found[0], first[0]


oeis_A090460 = {15: 1, 16: 1, 17: 1, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 3, 24: 0, 25: 10, 26: 12,
                27: 35, 28: 52, 29: 19, 30: 20, 31: 349, 32: 361}
out("E1 Hamiltonian path existence, n=1..32 (backtracking)")
exist = []
for n in range(1, 33):
    c, p = ham_paths_count(n, count_all=False, budget=60.0)
    if c:
        exist.append(n)
    if n in (15, 23):
        out("   n=%d witness path: %s" % (n, p))
out("   n with a Hamiltonian path: %s" % exist)
out("   n without: %s" % [n for n in range(1, 33) if n not in exist])
if exist != [1, 15, 16, 17, 23] + list(range(25, 33)):
    raise RuntimeError("Hamiltonian existence differs from session lead / A090461")
out("E2 Hamiltonian path counts up to reversal, n=15..30, vs OEIS A090460 (CITED, fetched 2026-09-22)")
for n in range(15, 31):
    c, _ = ham_paths_count(n, count_all=True, budget=240.0)
    out("   n=%d paths=%d A090460=%d %s" % (n, c, oeis_A090460[n], "ok" if c == oeis_A090460[n] else "MISMATCH"))
    if c != oeis_A090460[n]:
        raise RuntimeError("path count mismatch at n=%d" % n)
# leaf obstruction
out("E3 leaf obstruction: a Hamiltonian path has at most 2 leaves; leaf counts n=15..32: %s" %
    [(n, len(leafsets[n])) for n in range(15, 33)])
out("   n=18 leaves %s with neighbours %s -> no Hamiltonian path (PROVED by leaf count)" %
    (leafsets[18], [(m, nbrs(m, 18)) for m in leafsets[18]]))
out("   n=19 leaves %s neighbours %s: both endpoints forced" % (leafsets[19], [(m, nbrs(m, 19)) for m in leafsets[19]]))
# forced edges at n=15: degree-2 vertices
d2 = [(m, nbrs(m, 15)) for m in range(1, 16) if deg_in(m, 15) == 2]
out("E4 Q_15 degree-2 vertices (both edges forced in any Hamiltonian path): %s" % d2)
out("   Q_15 degree of 4: %d neighbours %s; degree list: %s" % (deg_in(4, 15), nbrs(4, 15), [deg_in(m, 15) for m in range(1, 16)]))
out("   Q_15 edge count %d, path uses %d edges, unused edges: %s" %
    (len(edges(15)), 14, [e for e in edges(15) if not any({e[0], e[1]} == {a, b} for a, b in zip([8, 1, 15, 10, 6, 3, 13, 12, 4, 5, 11, 14, 2, 7, 9], [1, 15, 10, 6, 3, 13, 12, 4, 5, 11, 14, 2, 7, 9]))]))
# session lead's paths verified
for n, path in ((15, [8, 1, 15, 10, 6, 3, 13, 12, 4, 5, 11, 14, 2, 7, 9]),
                (23, [18, 7, 9, 16, 20, 5, 11, 14, 2, 23, 13, 12, 4, 21, 15, 10, 6, 19, 17, 8, 1, 3, 22])):
    okp = sorted(path) == list(range(1, n + 1)) and all(is_square(a + b) for a, b in zip(path, path[1:]))
    out("E5 session lead path n=%d valid=%s sums=%s" % (n, okp, [a + b for a, b in zip(path, path[1:])]))
    if not okp:
        raise RuntimeError("session lead path invalid")


def ham_cycle_exists(n, budget=120.0):
    adj = [0] * (n + 1)
    for m in range(1, n + 1):
        for y in nbrs(m, n):
            adj[m] |= 1 << y
    full = ((1 << (n + 1)) - 2)
    t0 = time.time()
    s = 1
    res = [None]

    def dfs(v, visited, path):
        if visited == full:
            if adj[v] >> s & 1:
                res[0] = list(path)
                return True
            return False
        if time.time() - t0 > budget:
            raise TimeoutError
        rem = full & ~visited
        r = rem
        while r:
            lb = r & -r
            u = lb.bit_length() - 1
            r ^= lb
            if bin(adj[u] & (rem | (1 << v) | (1 << s))).count("1") < 2:
                return False
        cand = adj[v] & rem
        while cand:
            lb = cand & -cand
            w = lb.bit_length() - 1
            cand ^= lb
            path.append(w)
            if dfs(w, visited | lb, path):
                return True
            path.pop()
        return False

    dfs(s, 1 << s, [s])
    return res[0]


for n in (30, 31, 32):
    cyc = ham_cycle_exists(n)
    out("E6 Hamiltonian cycle in Q_%d: %s %s" % (n, cyc is not None, cyc if cyc else ""))
    if (cyc is not None) != (n == 32):
        raise RuntimeError("cycle existence differs from A090460 comment")

# ---------------------------------------------------------------- block F
out()
out("=== F. Analogy audit numbers ===")
# Collatz odd arrow n -> n + (n+1)/2 is a summand arrow; which targets are squares?
sq_targets = [n for n in range(1, 2001, 2) if is_square(n + (n + 1) // 2)]
out("F1 odd n <= 2000 whose odd-arrow target (3n+1)/2 is a square: %s (count %d)" % (sq_targets[:12], len(sq_targets)))
out("   PROVED: (3n+1)/2 = k^2 needs 2k^2 = 1 mod 3, but k^2 mod 3 in %s; so the odd-arrow target is never a square" % sorted(set(k * k % 3 for k in range(30))))
sq31 = [n for n in range(1, 2001, 2) if is_square(3 * n + 1)]
out("   contrast: odd n <= 2000 with 3n+1 itself a square: %s (count %d), k = %s" % (sq31[:8], len(sq31), [isqrt(3 * n + 1) for n in sq31[:8]]))
# Q_n out-degree vs Collatz functional graph
out("F2 Q_n max degree at n=2000: %d ; Collatz odd-arrow out-degree: 1 (functional graph)" % max(deg_in(m, 2000) for m in range(1, 2001)))
# pasted claims
out("F3 pasted 'prime ladder 3,7,11,17' differences: %s (not an AP)" % [7 - 3, 11 - 7, 17 - 11])
out("F4 pasted Lean square_sum_graph on Fin n with values x+1: value 2 has 2+2=4 square -> Adj 2 2 true, not loopless; vertex values with 2v square, v<=32: %s" %
    [v for v in range(1, 33) if is_square(2 * v)])
out("F5 pasted '196 = 1+3+5+7+11+13+17+19+23+29+31+37' check: %d" % sum([1, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]))
out("F6 pasted 'N=14 first unification (3 components join)': c(Q_12)=%d c(Q_13)=%d c(Q_14)=%d -> two merges, one at 13 and one at 14" % (ccount[12], ccount[13], ccount[14]))
out("F7 pasted 'N>=25 conjectured connected': connected already for all n >= 14 (B5, and PROVED in note S8)")

out()
out("=== G. Threshold constants used in the proofs ===")
out("G1 (1+sqrt2)^2 = %.3f : (sqrt n +1)^2 < 2n iff n > this, so every n >= 6 has a square strictly inside (n,2n)" % ((1 + math.sqrt(2)) ** 2))
out("G2 (2+sqrt8)^2 = %.3f : 4 sqrt n + 4 <= n iff n >= 24 (two squares in (m,m+n] for all m <= n)" % ((2 + math.sqrt(8)) ** 2))
out("G3 (3+sqrt18)^2 = %.3f : 6 sqrt n + 9 <= n iff n >= 53 (three squares in (m,m+n] for all m <= n)" % ((3 + math.sqrt(18)) ** 2))
out("G4 born-leaf finiteness: (s+2)^2 <= 2n-1 follows from n - 4 sqrt n - 5 >= 0 i.e. n >= %d; direct check deg_n(n) >= 2 for 19 <= n <= 24: %s" % (25, [deg_new(n) for n in range(19, 25)]))
out("G5 squares <= 2*12-1 = 23: %s (three squares, hence three reflections x -> k^2 - x for n <= 12)" % [k * k for k in range(1, 5) if k * k <= 23])
out("G6 closure holes {1,4,6} vs Q founders {1,2,4}: intersection %s, symmetric difference %s; excluded diagonals 2+2=%d (isolates 2 in Q_2; is the missing route to 4 in Cl^<), 3+3=%d (missing route to 6 in Cl^<)" % (sorted({1, 4, 6} & {1, 2, 4}), sorted({1, 4, 6} ^ {1, 2, 4}), 2 + 2, 3 + 3))
out("G7 Q_15: vertices of degree 3: %s; leaves: %s; forced edges = all %d edges except (1,3): the unique path" % ([m for m in range(1, 16) if deg_in(m, 15) == 3], leafsets[15], len(edges(15))))
out("G8 c(Q_n) is non-increasing for n >= 5 (new vertex has deg >= 1): check over n <= 2000: %s" % all(cseq[n - 1] <= cseq[n - 2] for n in range(6, NMAX + 1)))
out("G9 closure of Q_13 component sizes: %s; T_3 = 1+2+3 = %d lies in the chain founded by 1" % ([len(c) for c in components(13)], 6))
out()
out("elapsed %.1f s" % (time.time() - T0))
