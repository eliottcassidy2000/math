#!/usr/bin/env python3
"""collatz_mod6_20260922_w6_summand_closure_square_filter_audit.py

Adversarial audit of lane summand_closure_square_filter (wave 6, 2026-09-22).
Independent recomputation: adjacency by a plain double loop over x+y square
(no isqrt formula), BFS components, brute-force degrees, a separately written
Hamiltonian path counter (start < end, so each undirected path is counted
once, no halving) and cycle search, exact integer thresholds (no floats),
and the counterexample to the explorer's "even-length walk" wording (S16).
python3 -O gives identical output modulo the timing line.
"""
import math
import sys
import time

T0 = time.time()
sys.setrecursionlimit(10000)


def out(*a):
    print(*a)
    sys.stdout.flush()


def is_sq(x):
    return x >= 0 and math.isqrt(x) ** 2 == x


def adjacency(n):
    """plain double loop; x ~ y iff x != y and x + y is a square."""
    adj = {v: [] for v in range(1, n + 1)}
    for x in range(1, n + 1):
        for y in range(x + 1, n + 1):
            if is_sq(x + y):
                adj[x].append(y)
                adj[y].append(x)
    return adj


def comps_bfs(adj):
    seen = set()
    cs = []
    for s in adj:
        if s in seen:
            continue
        q = [s]
        seen.add(s)
        c = []
        while q:
            v = q.pop()
            c.append(v)
            for w in adj[v]:
                if w not in seen:
                    seen.add(w)
                    q.append(w)
        cs.append(sorted(c))
    return sorted(cs)


out("=== AUDIT 1: components, chains, merges (BFS on double-loop adjacency) ===")
cc = {}
for n in range(1, 41):
    adj = adjacency(n)
    cs = comps_bfs(adj)
    cc[n] = len(cs)
    if n <= 14:
        out("   n=%2d c=%d comps=%s maxdeg=%d" % (n, len(cs), cs, max(len(v) for v in adj.values())))
out("A1 c(Q_n), n=1..40: %s" % [cc[n] for n in range(1, 41)])
if [cc[n] for n in range(4, 13)] != [3] * 9 or cc[13] != 2 or cc[14] != 1 or cc[1] != 1 or cc[2] != 2 or cc[3] != 2:
    raise RuntimeError("component counts differ from note S5/S6/S8")
if comps_bfs(adjacency(12)) != [[1, 3, 6, 8, 10], [2, 7, 9], [4, 5, 11, 12]]:
    raise RuntimeError("chains of Q_12 differ")
adj13 = adjacency(13)
out("A2 Q_13 neighbours of 13: %s ; Q_14 neighbours of 14: %s" % (adj13[13], adjacency(14)[14]))
if adj13[13] != [3, 12] or adjacency(14)[14] != [2, 11]:
    raise RuntimeError("merge vertices differ")
lin = [n for n in range(1, 41) if max(len(v) for v in adjacency(n).values()) <= 2]
out("A3 n<=40 with max degree <= 2: %s" % lin)
if lin != list(range(1, 13)):
    raise RuntimeError("linear forest range differs")

out()
out("=== AUDIT 2: connectivity to 2000 (newest vertex has a neighbour; BFS spot checks) ===")


def newest_deg(n):
    return sum(1 for y in range(1, n) if is_sq(n + y))


nd = [newest_deg(n) for n in range(1, 2001)]
out("B1 deg_n(n), n=1..40: %s" % nd[:40])
zero = [n for n in range(1, 2001) if nd[n - 1] == 0]
out("B2 newest vertex isolated (founders) n<=2000: %s" % zero)
if zero != [1, 2, 4]:
    raise RuntimeError("founders differ")
one = [n for n in range(1, 2001) if nd[n - 1] == 1]
out("B3 newest vertex a leaf n<=2000: %s (count %d)" % (one, len(one)))
if one != [3, 5, 6, 7, 8, 9, 10, 11, 12, 16, 17, 18]:
    raise RuntimeError("born-leaf set differs")
for n in (14, 100, 500, 2000):
    c = len(comps_bfs(adjacency(n)))
    out("   BFS: c(Q_%d) = %d" % (n, c))
    if c != 1:
        raise RuntimeError("disconnected")
out("B4 c(Q_14)=1 and deg_n(n)>=1 for all 15<=n<=2000 => connected 14..2000 (min newest degree over 15..2000 = %d)" % min(nd[14:]))
# exact integer threshold for (s+1)^2 < 2n with s = isqrt(n)
thr = [n for n in range(1, 200) if (math.isqrt(n) + 1) ** 2 >= 2 * n]
out("B5 n<200 with (isqrt(n)+1)^2 >= 2n (newest vertex may lack the (s+1)^2 neighbour): %s" % thr)
# (sqrt n + 1)^2 < 2n  <=>  n - 2 sqrt n - 1 > 0  <=>  (n-1)^2 > 4n and n > 1  (exact)
thr2 = [n for n in range(2, 200) if not ((n - 1) ** 2 > 4 * n)]
out("   n<200 failing the sufficient bound (n-1)^2 > 4n (i.e. n > (1+sqrt2)^2): %s" % thr2)

out()
out("=== AUDIT 3: degree formula, leaves, leaf windows (brute force to 400, formula beyond) ===")


def deg_bf(m, n):
    return sum(1 for y in range(1, n + 1) if y != m and is_sq(m + y))


def deg_formula(m, n):
    return math.isqrt(m + n) - math.isqrt(m) - (1 if is_sq(2 * m) else 0)


bad = [(m, n) for n in range(1, 401) for m in range(1, n + 1) if deg_bf(m, n) != deg_formula(m, n)]
out("C1 degree formula mismatches for m<=n<=400: %s" % bad)
if bad:
    raise RuntimeError("degree formula wrong")
leaves = {}
for n in range(1, 2001):
    f = deg_bf if n <= 400 else deg_formula
    leaves[n] = [m for m in range(1, n + 1) if f(m, n) == 1]
out("C2 leaf sets n=1..32: %s" % [(n, leaves[n]) for n in range(1, 33)])
if leaves[18] != [16, 17, 18] or leaves[19] != [16, 18] or any(leaves[n] != [18] for n in range(20, 31)) or leaves[31] or leaves[32]:
    raise RuntimeError("session lead leaf lists differ")
hasleaf = [n for n in range(1, 2001) if leaves[n]]
out("C3 n<=2000 with a leaf: min %d max %d, count %d (contiguous=%s)" % (min(hasleaf), max(hasleaf), len(hasleaf), hasleaf == list(range(3, 31))))
if hasleaf != list(range(3, 31)):
    raise RuntimeError("leaf range differs")
win = {}
for m in range(1, 31):
    w = [n for n in range(m, 200) if deg_formula(m, n) == 1]
    if w:
        win[m] = (w[0], w[-1], len(w))
out("C4 leaf windows m -> (first n, last n, length): %s" % win)
longest = sorted(win.items(), key=lambda kv: -kv[1][2])[:6]
out("C5 windows sorted by length (top 6): %s" % longest)
sevens = sorted(m for m, v in win.items() if v[2] == 7)
out("   windows of length exactly 7: m in %s  => the explorer's 'three longest windows are 2, 8, 18' is WRONG (4 and 9 tie with 2); the two longest (9 and 13) are 8 and 18" % sevens)
if sevens != [2, 4, 9] or win[8][2] != 9 or win[18][2] != 13:
    raise RuntimeError("window lengths differ")
# S13 thresholds, exact integers: 4 sqrt n + 4 <= n  <=>  n >= 4 and 16 n <= (n-4)^2
t24 = min(n for n in range(1, 400) if all(k >= 4 and 16 * k <= (k - 4) ** 2 for k in range(n, 400)))
t53 = min(n for n in range(1, 400) if all(k >= 9 and 36 * k <= (k - 9) ** 2 for k in range(n, 400)))
t25 = min(n for n in range(1, 400) if all(k >= 5 and 16 * k <= (k - 5) ** 2 for k in range(n, 400)))
out("C6 exact thresholds: 4sqrt(n)+4<=n from n=%d; 6sqrt(n)+9<=n from n=%d; 4sqrt(n)+5<=n from n=%d" % (t24, t53, t25))
if (t24, t53, t25) != (24, 53, 25):
    raise RuntimeError("thresholds differ")
sqm = sorted(set(deg_formula(m, n) for n in range(31, 61) for m in range(1, n + 1) if is_sq(2 * m)))
sqm24 = sorted(set(deg_formula(m, n) for n in range(24, 61) for m in range(1, n + 1) if is_sq(2 * m)))
out("C7 deg_n(m) over 31<=n<=60, 2m square: %s ; over 24<=n<=60: %s (explorer's C11 label said 24 but computed 31: label bug, fixed)" % (sqm, sqm24))
if sqm != [2, 3, 4, 5] or 1 not in sqm24:
    raise RuntimeError("C11 audit differs")
out("C8 deg_n(18) for n=29,30,31,32: %s" % [deg_bf(18, n) for n in (29, 30, 31, 32)])

out()
out("=== AUDIT 4: self-similarity claims ===")
cnt = [sum(1 for k in range(1, 2 * j) if 2 * j * j < k * k < 4 * j * j) for j in range(1, 51)]
frm = [2 * j - 1 - math.isqrt(2 * j * j) for j in range(1, 51)]
out("D1 #squares in (2j^2,4j^2) j=1..12 brute: %s formula ok for j<=50: %s" % (cnt[:12], cnt == frm))
if cnt != frm:
    raise RuntimeError("count formula differs")
jl = [j for j in range(1, 51) if cnt[j - 1] <= 1]
out("D2 j<=50 with at most one square in (2j^2,4j^2): %s" % jl)
for j in (1, 2, 3, 4, 5, 6):
    m = 2 * j * j
    w = [n for n in range(m, 8 * m + 100) if deg_formula(m, n) == 1]
    out("   m=%d leaf window %s" % (m, [w[0], w[-1]] if w else None))
# S16: differences of squares l^2-k^2 (k>=1) are never 2 mod 4; but an even-length walk CAN move by 2 mod 4
diffs = sorted(set(l * l - k * k for k in range(1, 40) for l in range(k + 1, 40)))
out("D3 residues mod 4 of l^2-k^2 (1<=k<l<40): %s" % sorted(set(d % 4 for d in diffs)))
walk = [1, 3, 6, 10, 15]
out("D4 COUNTEREXAMPLE to 'every even-length walk moves by an integer != 2 mod 4': walk %s, sums %s, length %d, displacement %d = %d mod 4"
    % (walk, [a + b for a, b in zip(walk, walk[1:])], len(walk) - 1, walk[-1] - walk[0], (walk[-1] - walk[0]) % 4))
if not all(is_sq(a + b) for a, b in zip(walk, walk[1:])) or (walk[-1] - walk[0]) % 4 != 2:
    raise RuntimeError("counterexample invalid")
out("   correct statement: a walk of length exactly 2 moves by l^2-k^2, never 2 mod 4; length-2 walks in Q_15 with displacement 2 mod 4: %d" %
    sum(1 for x in range(1, 16) for y in adjacency(15)[x] for z in adjacency(15)[y] if z != x and (z - x) % 4 == 2))

out()
out("=== AUDIT 5: Hamiltonian paths/cycles, independent counter (start < end, each path once) ===")


def ham_count(n, existence_only=False, budget=300.0):
    adj = adjacency(n)
    bits = {v: sum(1 << w for w in adj[v]) for v in adj}
    full = sum(1 << v for v in range(1, n + 1))
    t0 = time.time()
    cnt = [0]
    wit = [None]

    def rec(v, vis, path, start):
        if vis == full:
            if v > start:
                cnt[0] += 1
                if wit[0] is None:
                    wit[0] = list(path)
                return existence_only
            return False
        if time.time() - t0 > budget:
            raise TimeoutError("n=%d" % n)
        rem = full & ~vis
        # each unvisited vertex needs a neighbour in rem or v; at most one unvisited
        # vertex (the future end) may have a single such neighbour
        low = 0
        r = rem
        while r:
            b = r & -r
            u = b.bit_length() - 1
            r ^= b
            d = bin(bits[u] & (rem | (1 << v))).count("1")
            if d == 0 or (d == 1 and (low := low + 1) > 1 and rem != b):
                return False
        c = bits[v] & rem
        while c:
            b = c & -c
            w = b.bit_length() - 1
            c ^= b
            path.append(w)
            if rec(w, vis | b, path, start):
                return True
            path.pop()
        return False

    if n == 1:
        return 1, [1]
    for s in range(1, n + 1):
        if rec(s, 1 << s, [s], s):
            break
    return cnt[0], wit[0]


A090460 = [1, 1, 1, 0, 0, 0, 0, 0, 3, 0, 10, 12, 35, 52, 19, 20, 349, 361]
ex = []
for n in range(1, 33):
    c, w = ham_count(n, existence_only=True, budget=60.0)
    if c:
        ex.append(n)
out("E1 Hamiltonian path exists for n<=32: %s" % ex)
if ex != [1, 15, 16, 17, 23] + list(range(25, 33)):
    raise RuntimeError("existence differs")
counts = []
for n in range(15, 31):
    c, w = ham_count(n, budget=300.0)
    counts.append(c)
    if n == 15:
        out("   n=15 unique path: %s" % w)
out("E2 counts n=15..30: %s ; A090460 (fetched 2026-09-22, offset 15): %s ; equal=%s" % (counts, A090460[:16], counts == A090460[:16]))
if counts != A090460[:16]:
    raise RuntimeError("counts differ")
# leaves in Q_15..Q_32 and forced edges at 15
adj15 = adjacency(15)
degs = [len(adj15[v]) for v in range(1, 16)]
out("E3 Q_15 degrees: %s edges=%d leaves=%s deg3=%s" % (degs, sum(degs) // 2, [v for v in range(1, 16) if degs[v - 1] == 1], [v for v in range(1, 16) if degs[v - 1] == 3]))
forced = set()
for v in range(1, 16):
    if degs[v - 1] <= 2:
        for w in adj15[v]:
            forced.add((min(v, w), max(v, w)))
alle = set((v, w) for v in range(1, 16) for w in adj15[v] if v < w)
out("   forced edges %d, unforced %s" % (len(forced), sorted(alle - forced)))
if sorted(alle - forced) != [(1, 3)]:
    raise RuntimeError("forced-edge argument differs")
for n, p in ((15, [8, 1, 15, 10, 6, 3, 13, 12, 4, 5, 11, 14, 2, 7, 9]),
             (23, [18, 7, 9, 16, 20, 5, 11, 14, 2, 23, 13, 12, 4, 21, 15, 10, 6, 19, 17, 8, 1, 3, 22]),
             (23, [2, 23, 13, 12, 4, 21, 15, 10, 6, 19, 17, 8, 1, 3, 22, 14, 11, 5, 20, 16, 9, 7, 18])):
    ok = sorted(p) == list(range(1, n + 1)) and all(is_sq(a + b) for a, b in zip(p, p[1:]))
    out("E4 path n=%d valid=%s" % (n, ok))
    if not ok:
        raise RuntimeError("path invalid")


def ham_cycle(n, budget=200.0):
    adj = adjacency(n)
    bits = {v: sum(1 << w for w in adj[v]) for v in adj}
    full = sum(1 << v for v in range(1, n + 1))
    t0 = time.time()
    res = [None]
    s = 1

    def rec(v, vis, path):
        if vis == full:
            if bits[v] >> s & 1:
                res[0] = list(path)
                return True
            return False
        if time.time() - t0 > budget:
            raise TimeoutError
        rem = full & ~vis
        r = rem
        while r:
            b = r & -r
            u = b.bit_length() - 1
            r ^= b
            if bin(bits[u] & (rem | (1 << v) | (1 << s))).count("1") < 2:
                return False
        c = bits[v] & rem
        while c:
            b = c & -c
            w = b.bit_length() - 1
            c ^= b
            path.append(w)
            if rec(w, vis | b, path):
                return True
            path.pop()
        return False

    rec(s, 1 << s, [s])
    return res[0]


for n in (30, 31, 32):
    cyc = ham_cycle(n)
    ok = cyc is not None and sorted(cyc) == list(range(1, n + 1)) and all(is_sq(a + b) for a, b in zip(cyc, cyc[1:] + cyc[:1]))
    out("E5 Hamiltonian cycle Q_%d: %s" % (n, ok))
    if ok != (n == 32):
        raise RuntimeError("cycle existence differs")
c32 = [1, 8, 28, 21, 4, 32, 17, 19, 30, 6, 3, 13, 12, 24, 25, 11, 5, 31, 18, 7, 29, 20, 16, 9, 27, 22, 14, 2, 23, 26, 10, 15]
out("   explorer's Q_32 cycle valid: %s" % (sorted(c32) == list(range(1, 33)) and all(is_sq(a + b) for a, b in zip(c32, c32[1:] + c32[:1]))))

out()
out("=== AUDIT 6: analogy numbers, target-filter founders (S22, S25), pasted arithmetic ===")
bad = [n for n in range(1, 100001, 2) if is_sq((3 * n + 1) // 2)]
out("F1 odd n<=100000 with (3n+1)/2 square: %s ; k^2 mod 3 values: %s" % (bad, sorted(set(k * k % 3 for k in range(100)))))
if bad:
    raise RuntimeError("S22 false")
sq31 = [n for n in range(1, 2001, 2) if is_sq(3 * n + 1)]
out("   odd n<=2000 with 3n+1 square: first %s count %d" % (sq31[:8], len(sq31)))


def founders(T, upto):
    return [n for n in range(1, upto + 1) if not any((n + y) in T for y in range(1, n))]


hits = 0
for mask in range(1 << 14):  # T subset {2..15}; founders among n <= 8 depend only on T cap [2,15]
    T = {v + 2 for v in range(14) if mask >> v & 1}
    if founders(T, 8) == [1, 4, 6]:
        hits += 1
out("G1 T subset {2..15} with founders among n<=8 exactly {1,4,6}: %d" % hits)
if hits:
    raise RuntimeError("S25 false")
out("   proof check: windows [n+1,2n-1] for n=4,5,6: %s (window of 5 lies inside the union of those of 4 and 6: %s)" %
    ([(n, n + 1, 2 * n - 1) for n in (4, 5, 6)], set(range(6, 10)) <= set(range(5, 8)) | set(range(7, 12))))
out("   square filter founders n<=20: %s" % founders({k * k for k in range(1, 8)}, 20))
out("F3 3,7,11,17 differences %s ; 196 check %d ; 2v square v<=32: %s" %
    ([4, 4, 6], sum([1, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]), [v for v in range(1, 33) if is_sq(2 * v)]))

out()
out("=== AUDIT 7: THM-2422 re-check with an independent closure (asynchronous, bound 300) ===")


def closure(seeds, bound):
    S = set(seeds)
    frontier = True
    while frontier:
        new = {a + b for a in S for b in S if a < b and a + b <= bound} - S
        frontier = bool(new)
        S |= new
    return S


for seeds, exp in (({2, 3}, [1, 4, 6]), ({1, 2, 3}, []), ({2, 3, 4}, [1]), ({2, 3, 6}, [1, 4])):
    h = sorted(set(range(1, 301)) - closure(seeds, 300))
    out("H1 seeds %s holes in [1,300]: %s" % (sorted(seeds), h))
    if h != exp:
        raise RuntimeError("closure differs")
S = {2, 3}
Ms = []
for t in range(0, 8):
    if t >= 4:
        Ms.append(max(S))
        if S != {2, 3, 5} | set(range(7, max(S) + 1)):
            raise RuntimeError("synchronous law")
    S = S | {a + b for a in S for b in S if a < b}
out("H2 synchronous M_t t=4..7: %s = 27*2^(t-4)+1: %s" % (Ms, [27 * 2 ** (t - 4) + 1 for t in range(4, 8)]))
out()
out("AUDIT VERDICT: all explorer numbers reproduced; two wording defects found (S12 'three longest windows', S16 'even-length walk') and one .out label defect (C11 range), all corrected in the note/script.")
out("elapsed %.1f s" % (time.time() - T0))
