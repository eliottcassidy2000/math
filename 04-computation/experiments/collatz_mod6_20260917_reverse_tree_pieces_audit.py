#!/usr/bin/env python3
"""
collatz_mod6_20260917_reverse_tree_pieces_audit.py
Independent adversarial recomputation for lane reverse_tree_pieces (session collatz-mod6-20260917, mac-mini),
written 2026-09-21 by the single adversarial verify-and-fix pass.  Everything here is recomputed from scratch
with different code paths from the lane script (standard Collatz map instead of the accelerated T where possible,
forward simulation instead of inverse closure, numpy eigenvalues instead of power iteration, and so on).
All checks are explicit `raise` (survive python -O).

A1  row law, base table, rotation (brute force)                       A7  pruned tree: 4/3, (4/3)^k, root-8 counts
A2  tree to depth 6, j<=8 (independent generator)                     A8  difference inequality: numeric verification + exponent + base case
A3  Theorem 1.2 modulus (row needs node mod 3 only; sharpness)        A9  Perron roots with numpy on a grid: all roots of rho(g)=1 per split
A4  depth census (standard Collatz map, odd-step count)               A10 section 4: G vs child, guarded/TO-1-in-E/FROM-1-in-E closures, SCC
A5  fibre uniformity (all internal u < 2000)                          A11 section 5: image law, orbits, survival, parity independence
A6  truncated fibre count 416666                                      A12 boundary cases and near misses (u=1,2,4,9; 17; 11)
"""
import sys, math, time
from fractions import Fraction
from collections import Counter, defaultdict
import numpy as np

def fail(msg):
    raise RuntimeError("AUDIT CHECK FAILED: " + msg)

def check(cond, msg):
    if not cond:
        fail(msg)

def banner(t):
    print("=" * 78); print(t); print("=" * 78)

def T(n):
    m = 3 * n + 1
    while m % 2 == 0:
        m //= 2
    return m

def h0(u):
    r = u % 3
    if r == 0:
        fail("h0 of multiple of 3")
    return 1 if r == 1 else 0

def child(u, j):
    return ((1 << (h0(u) + 1 + 2 * j)) * u - 1) // 3

# ------------------------------------------------------------------ A1
banner("A1  row law and base table")
for n in range(1, 20001, 2):
    r = n % 6
    check(r == {4: 1, 1: 3, 7: 5}[(3 * n + 1) % 9], "rho at %d" % n)
    check(r == {1: 1, 0: 3, 2: 5}[n % 3], "row is a function of n mod 3 at %d" % n)
base = {}
for u in range(1, 2000, 2):
    if u % 3 == 0:
        continue
    rows = [child(u, j) % 6 for j in range(9)]
    for j in range(8):
        check(rows[j + 1] == (rows[j] + 4) % 6, "rotation +4 at u=%d" % u)
    base.setdefault(u % 18, set()).add(rows[0])
check(base == {1: {1}, 11: {1}, 13: {5}, 17: {5}, 5: {3}, 7: {3}}, "base table %s" % base)
print("row(n) = rho(3n+1 mod 9) and row(n) = n mod 3 lifted (1->1, 0->3, 2->5) for all odd n < 20001: OK")
print("base row by u mod 18: %s ; rotation +4 mod 6 along j for all internal u < 2000: OK" % {k: sorted(v) for k, v in sorted(base.items())})
print("NOTE: the row of an odd node is a function of the node mod 3 (not mod 9); the lane note's Theorem 1.2 proof")
print("      sentence 'the row of a node needs the node mod 9' is wrong as written; with mod 3 the modulus 3^(D+1) is right.")

# ------------------------------------------------------------------ A2
banner("A2  tree to depth 6 with j<=8 (independent BFS)")
levels = [[1]]
seen = {1}
rowc = []
for d in range(1, 7):
    nxt = []
    rc = Counter()
    for u in levels[-1]:
        if u % 3 == 0:
            continue
        for j in range(9):
            n = child(u, j)
            check(T(n) == u, "T(child)")
            if n == 1:
                check(u == 1 and j == 0, "unexpected 1")
                continue
            check(n not in seen, "duplicate %d" % n)
            seen.add(n)
            nxt.append(n)
            rc[n % 6] += 1
    levels.append(nxt)
    rowc.append(dict(sorted(rc.items())))
sizes = [len(l) for l in levels]
check(sizes == [1, 8, 45, 270, 1620, 9720, 58320], "level sizes %s" % sizes)
check(sum(sizes) == 69984, "total")
check(rowc[0] == {1: 2, 3: 3, 5: 3}, "level-1 rows")
for d in range(1, 6):
    c = sizes[d + 1] // 3
    check(rowc[d] == {1: c, 3: c, 5: c}, "level %d rows" % (d + 1))
print("levels %s total %d, row counts %s: OK" % (sizes, sum(sizes), rowc))

# ------------------------------------------------------------------ A3
banner("A3  Theorem 1.2: rows along a fixed index path are a function of u mod 3^(D+1), sharp")
def pattern(u, path):
    out = []
    x = u
    for j in path:
        if x % 3 == 0:
            out.append('L'); break
        x = child(x, j)
        out.append(x % 6)
    return tuple(out)
for D in range(1, 6):
    for path in [(0,) * D, (2, 1, 0, 1, 2)[:D], (1,) * D]:
        good = {}
        for u in range(1, 4 * 3 ** (D + 2), 2):
            if u % 3 == 0:
                continue
            k = u % 3 ** (D + 1)
            p = pattern(u, path)
            check(good.setdefault(k, p) == p, "not a function of u mod 3^(D+1), D=%d" % D)
        # sharpness: every class mod 3^D (coprime to 3) splits
        bad = defaultdict(set)
        for u in range(1, 4 * 3 ** (D + 2), 2):
            if u % 3 == 0:
                continue
            bad[u % 3 ** D].add(pattern(u, path))
        splitting = sorted(k for k, s in bad.items() if len(s) > 1)
        check(len(splitting) > 0, "no class mod 3^D splits, D=%d path=%s" % (D, path))
        if path == (0,) * D:
            check(1 % 3 ** D in splitting, "witness class 1 mod 3^D on the zero path, D=%d" % D)
        print("D=%d path=%s: function of u mod %d; classes mod %d that split: %d of %d (a class whose path reaches a leaf cannot split)" %
              (D, path, 3 ** (D + 1), 3 ** D, len(splitting), len(bad)))
print("Theorem 1.2 modulus 3^(D+1) confirmed for D<=5 on three paths; sharp (some class mod 3^D splits; class 1 splits on the zero path).")

# ------------------------------------------------------------------ A4
banner("A4  depth census of odd u <= 10^6 via the standard Collatz map (odd-step count)")
X = 10 ** 6
t0 = time.time()
depth = {1: 0}
for u in range(3, X + 1, 2):
    x, s = u, 0
    while x >= u:
        x = 3 * x + 1
        s += 1
        while x % 2 == 0:
            x //= 2
    depth[u] = s + depth[x]
maxd = max(depth.values())
arg = [u for u, d in depth.items() if d == maxd]
check(maxd == 195 and arg == [837799], "max depth %d at %s" % (maxd, arg))
# standard total stopping time of 837799 is 524 (classical record below 10^6): recompute
x, tot = 837799, 0
while x != 1:
    x = x // 2 if x % 2 == 0 else 3 * x + 1
    tot += 1
check(tot == 524, "total stopping time of 837799 = %d" % tot)
print("odd u <= 10^6: max odd-step depth 195 at 837799 (its full Collatz step count is %d): OK, %.1fs" % (tot, time.time() - t0))
cnt = Counter(depth.values())
by9 = defaultdict(Counter)
for u, d in depth.items():
    by9[d][u % 9] += 1
expected_counts = {0: 1, 1: 9, 2: 34, 3: 78, 4: 176, 5: 282, 6: 508, 7: 871, 8: 1056, 9: 1657, 10: 1860,
                   11: 2593, 12: 3415, 13: 3514, 14: 4882, 15: 4295}
for d, c in expected_counts.items():
    check(cnt[d] == c, "depth %d count %d vs %d" % (d, cnt[d], c))
check([by9[2][r] for r in range(9)] == [4, 4, 3, 6, 5, 5, 2, 2, 3], "depth-2 residues")
check([by9[15][r] for r in range(9)] == [502, 483, 471, 480, 488, 457, 468, 478, 468], "depth-15 residues")
print("depth counts 0..15 and the mod-9 rows for depths 2 and 15 match the lane .out: OK")
# depth 2 by hand: children <= X of the internal depth-1 nodes <= X give 33; the 34th has parent 1398101 > X
d2 = sorted(u for u, d in depth.items() if d == 2)
parents = Counter(T(u) for u in d2)
check(len(d2) == 34 and parents[1398101] == 1 and 932067 in d2, "depth-2 anatomy %s" % parents)
print("depth-2 anatomy: parents %s; the node 932067 has parent 1398101 > X, so the truncated level is NOT the" % dict(parents))
print("   union of truncated fibres of nodes <= X (children can be smaller than parents when u = 2 mod 3).")
buckets = {}
for lo in range(0, maxd + 1, 20):
    c = sum(cnt[d] for d in range(lo, lo + 20))
    m9 = [sum(by9[d][r] for d in range(lo, lo + 20)) for r in range(9)]
    buckets[lo] = (c, m9)
check(buckets[0][0] == 52108 and buckets[20][0] == 168230 and buckets[180][0] == 4, "bucket counts")
check(abs((buckets[20][1][0] + buckets[20][1][3] + buckets[20][1][6]) / 168230 - 0.3333) < 5e-5, "bucket 20-39 frac 0 mod 3")
print("bucket counts 0-19: %d, 20-39: %d, 180-199: %d: OK" % (buckets[0][0], buckets[20][0], buckets[180][0]))
tot9 = [sum(by9[d][r] for d in by9) for r in range(9)]
check(tot9 == [55556, 55556, 55555, 55556, 55555, 55556, 55555, 55556, 55555], "pooled mod 9 %s" % tot9)
n_int = sum(1 for u in depth if u % 3)
n_int1 = sum(1 for u in depth if u % 3 == 1)
n_int_d1 = sum(1 for u in depth if u % 3 and u != 1)
n_int1_d1 = n_int1 - 1
print("internal nodes: %d, of which 1 mod 3: %d (fraction %.6f); excluding depth 0: %d / %d = %.6f" %
      (n_int, n_int1, n_int1 / n_int, n_int1_d1, n_int_d1, n_int1_d1 / n_int_d1))
check(n_int == 333333 and n_int1 == 166667, "internal split")
print("   this is 166667/333333 = the count of odd u = 1 mod 6 vs 5 mod 6 up to X: a TRIVIAL fact about [1,X], not a tree statement.")
worst = max((max(abs(by9[d][r] / cnt[d] - 1 / 9) for r in range(9)), d) for d in cnt if cnt[d] >= 1000)
check(abs(worst[0] - 0.0194) < 5e-5 and worst[1] == 91, "worst deviation %s" % (worst,))
print("largest mod-9 deviation over depths with >= 1000 nodes: %.4f at depth %d: OK" % worst)

# ------------------------------------------------------------------ A5
banner("A5  fibre uniformity mod 9 for all internal u < 2000; ord_27(4)")
o = 1; x = 4
while x != 1:
    x = x * 4 % 27; o += 1
check(o == 9, "ord_27(4)")
for u in range(1, 2000, 2):
    if u % 3 == 0:
        continue
    check(sorted(child(u, j) % 9 for j in range(9)) == list(range(9)), "fibre uniformity u=%d" % u)
    check(sorted(child(u, j) % 27 for j in range(27)) == list(range(27)), "fibre uniformity mod 27 u=%d" % u)
print("ord_27(4) = 9; every internal u < 2000 has children j=0..8 covering Z/9 and j=0..26 covering Z/27: OK")
print("   (period 27 in j mod 27 follows from ord_81(4) = 27, the same argument one digit deeper)")

# ------------------------------------------------------------------ A6
banner("A6  truncated fibre count")
cross = sum(1 for n in range(1, X + 1, 2) if T(n) <= X)
check(cross == 416666, "cross %d" % cross)
print("#{odd n <= 10^6 : T(n) <= 10^6} = %d = 333333 * 1.25 (1.2500 mean truncated fibre): OK" % cross)

# ------------------------------------------------------------------ A7
banner("A7  pruned T-tree: children counts, (4/3)^k exact, root-8 depth counts")
def pch(y):
    out = [2 * y]
    if y % 3 == 2 and ((2 * y - 1) // 3) % 3 != 0:
        out.append((2 * y - 1) // 3)
    return out
cnts = {r: len(pch(r + 90)) for r in (1, 2, 4, 5, 7, 8)}
check(cnts == {1: 1, 2: 2, 4: 1, 5: 1, 7: 1, 8: 2}, "children by residue %s" % cnts)
for k in range(1, 8):
    M = 3 ** (k + 1)
    total = 0
    for a in range(1, M):
        if a % 3 == 0:
            continue
        fr = [a + M]          # a different representative than the lane script (a + M instead of a)
        for _ in range(k):
            fr = [c for y in fr for c in pch(y)]
        total += len(fr)
    check(Fraction(total, 2 * 3 ** k) == Fraction(4, 3) ** k, "(4/3)^k at k=%d" % k)
print("sum over root classes of #depth-k nodes = 8*4^(k-1), average (4/3)^k for k<=7 (representatives a+3^(k+1)): OK")
fr = [8]; counts = []
for k in range(33):
    counts.append(len(fr))
    fr = [c for y in fr for c in pch(y)]
exp_counts = [1, 2, 2, 2, 3, 4, 6, 10, 13, 14, 18, 25, 33, 46, 61, 77, 107, 144, 189, 253, 331, 441, 591, 802, 1066, 1412, 1876, 2492, 3345, 4453, 5936, 7925, 10563]
check(counts == exp_counts, "root-8 counts %s" % counts)
check(abs(counts[32] ** (1 / 32) - 1.3358) < 5e-5, "32nd root")
print("root-8 pruned tree depth counts 0..32 match; 10563^(1/32) = %.4f: OK" % counts[32] ** (1 / 32))

# ------------------------------------------------------------------ A8
banner("A8  modulus-9 difference inequality: numeric check on N restricted to [1,10^6], exponent, base case")
# all n <= 10^6 reach 1 under T (A4 covers odd n; even n halve to an odd n < n); so N cap [1,10^6] = non-multiples of 3
LIM = 10 ** 6
pref = [[0] * (LIM + 1) for _ in range(9)]
for r in (1, 2, 4, 5, 7, 8):
    p = pref[r]
    for n in range(1, LIM + 1):
        p[n] = p[n - 1] + (1 if n % 9 == r else 0)
def f(r, x):
    x = int(math.floor(x))
    return pref[r][x] if x >= 0 else 0
def A(x):
    return f(1, x) + f(4, x) + f(7, x)
def B(x):
    return f(2, x) + f(5, x) + f(8, x)
viol = 0
for x in list(range(1, 3000)) + list(range(3000, 600000, 997)):
    if not (A(x) >= B(x / 2) + f(2, 3 * x / 2)): viol += 1
    if not (B(x) >= A(x / 2) + f(8, 3 * x / 2)): viol += 1
    if not (f(2, x) >= max(f(1, x / 2), f(7, x / 8), f(4, x / 32))): viol += 1
    if not (f(8, x) >= max(f(4, x / 2), f(1, x / 8), f(7, x / 32))): viol += 1
    if not (A(x) >= A(x / 4) + A(3 * x / 128) / 3 + A(3 * x / 64) / 3): viol += 1
check(viol == 0, "%d violations" % viol)
print("inequalities (i), (ii), (iii) hold at every integer x < 3000 and a sparse grid to 6*10^5 (universe: all n <= 10^6, which reach 1): OK")
# these are lower bounds for the full N; inside [1,10^6] N is all non-multiples of 3, so the checks are on the true N truncated.
def phi(g):
    return 4 ** (-g) + (3 / 128) ** g / 3 + (3 / 64) ** g / 3 - 1
lo, hi = 0.0, 1.0
for _ in range(100):
    mid = (lo + hi) / 2
    if phi(mid) > 0: lo = mid
    else: hi = mid
gam = (lo + hi) / 2
check(abs(gam - 0.246227) < 2e-6, "gamma %.6f" % gam)
# uniqueness of the root on (0, 2]: phi is strictly decreasing (each term is decreasing in g)
check(all(phi(g / 100) > phi((g + 1) / 100) for g in range(0, 200)), "phi monotone")
# base case: c = 256^(-gamma); check c*y^gamma <= 1 <= A(y) on [1,256) and the induction step numerically to 10^6
c = 256 ** (-gam)
check(all(c * y ** gam <= 1 <= A(y) for y in range(1, 256)), "base case")
check(all(A(x) >= c * x ** gam for x in range(1, LIM + 1, 13)), "A(x) >= c x^gamma numerically")
print("gamma = %.6f (unique root, phi strictly decreasing); base case c = 256^-gamma = %.4f; A(x) >= c x^gamma on a grid to 10^6: OK" % (gam, c))

# ------------------------------------------------------------------ A9
banner("A9  Perron roots of the constant-split system with numpy: all roots of rho(g) = 1 on [0, 3]")
states = [1, 2, 4, 5, 7, 8]
idx = {r: i for i, r in enumerate(states)}
def rho(g, s2, s8):
    M = np.zeros((6, 6))
    for r in states:
        M[idx[2 * r % 9], idx[r]] += 2.0 ** (-g)
    for s, a in s2.items():
        M[idx[s], idx[2]] += a * 1.5 ** g
    for s, a in s8.items():
        M[idx[s], idx[8]] += a * 1.5 ** g
    return max(abs(np.linalg.eigvals(M)))
def roots(s2, s8):
    gs = [i / 1000 for i in range(0, 3001)]
    vals = [rho(g, s2, s8) - 1 for g in gs]
    out = []
    for i in range(len(gs) - 1):
        if abs(vals[i]) < 1e-12:
            out.append(round(gs[i], 6)); continue
        if vals[i] * vals[i + 1] < 0 and abs(vals[i + 1]) >= 1e-12:
            a, b = gs[i], gs[i + 1]
            fa = vals[i]
            for _ in range(60):
                m = (a + b) / 2
                fm = rho(m, s2, s8) - 1
                if (fm > 0) == (fa > 0): a, fa = m, fm
                else: b = m
            out.append(round((a + b) / 2, 6))
    return out, min(vals) + 1
u3 = {1: 1 / 3, 4: 1 / 3, 7: 1 / 3}; v3 = {2: 1 / 3, 5: 1 / 3, 8: 1 / 3}
r_u, min_u = roots(u3, v3)
print("uniform split: roots of rho(g)=1 on [0,3]: %s  (rho(g) = 2^-g + 3^(g-1) 2^-g exactly, since the all-ones vector is an eigenvector)" % r_u)
check(r_u == [1.0, 2.0], "uniform roots %s" % r_u)
check(all(abs(rho(g, u3, v3) - (2 ** (-g) + 1.5 ** g / 3)) < 1e-9 for g in (0.3, 1.0, 1.7, 2.5)), "uniform closed form")
table = {}
for s2 in (1, 4, 7):
    for s8 in (2, 5, 8):
        rr, mn = roots({s2: 1.0}, {s8: 1.0})
        table[(s2, s8)] = rr
        print("vertex split 2->%d, 8->%d: roots %s ; min rho on [0,3] = %.4f" % (s2, s8, rr, mn))
lane_vals = {(1, 2): 1.244017, (1, 5): 0.774576, (4, 5): 0.576032, (7, 2): 0.602605, (7, 5): 0.436588}
for k, v in lane_vals.items():
    check(table[k] and abs(table[k][0] - v) < 2e-6, "smallest root at %s: %s vs lane %f" % (k, table[k], v))
for s2 in (1, 4, 7):
    check(table[(s2, 8)] == [], "8->8 has no root: %s" % table[(s2, 8)])
    check(min(rho(g, {s2: 1.0}, {8: 1.0}) for g in (0, 0.5, 1, 1.5, 2, 3)) > 1, "8->8 rho > 1")
check(table[(4, 2)] == [], "2->4, 8->2 has no root: %s" % table[(4, 2)])
# the lane's 0.867659 for 2->4, 8->2 is a power-iteration artifact: the support graph has the 3-cycle 2->4->8->2 and the
# 6-cycle, period gcd = 3, so the one-step ratio of the power iteration oscillates with period 3 and hits 1 at one phase
Mx = np.zeros((6, 6)); g = 0.867659
for r in states:
    Mx[idx[2 * r % 9], idx[r]] += 2.0 ** (-g)
Mx[idx[4], idx[2]] += 1.5 ** g; Mx[idx[2], idx[8]] += 1.5 ** g
v = np.ones(6)
ratios = [float((np.linalg.matrix_power(Mx, k + 1) @ v).sum() / (np.linalg.matrix_power(Mx, k) @ v).sum()) for k in range(300, 306)]
rad = max(abs(np.linalg.eigvals(Mx)))
check(abs(rad - 1.168535) < 1e-5 and abs(ratios[2] - 1.0) < 1e-6 and abs(ratios[0] - 1.3233) < 1e-3, "artifact anatomy %s %s" % (rad, ratios))
print("split 2->4, 8->2 at the lane's g = 0.867659: true spectral radius %.6f; power-iteration one-step ratios cycle %s" %
      (rad, [round(r, 4) for r in ratios]))
smallest = min(v[0] for v in table.values() if v)
check(abs(smallest - 0.436588) < 2e-6, "vertex minimum")
print("REFUTED (four table entries): the lane's 'g = 1.500000' for the three 8->8 splits is the bisection cap, not a root")
print("   (rho(g) >= 1.5^g > 1 for g > 0 and rho(g) >= 2^-g > 1 for g < 0), and 'g = 0.867659' for 2->4, 8->2 is a period-3")
print("   power-iteration artifact (min rho on [0,3] is 1.1208 > 1: no root).  The five remaining vertex values and the minimum")
print("   0.436588 at 2->7, 8->5 are the SMALLEST roots and are confirmed; rho(g) is not monotone (uniform split has roots 1 and 2),")
print("   so the lane's bisection on [0, 1.5] found the smaller root only by luck of the bracket.")

# ------------------------------------------------------------------ A10
banner("A10  section 4: G vs minimal child, guarded closure, TO-1-in-E and FROM-1-in-E closures, giant SCC")
def G(m):
    k = 0
    while ((m << k) % 9) not in (4, 7):
        k += 1
    return ((m << k) - 1) // 3, k
for u in range(1, 100001, 2):
    if u % 3 == 0:
        continue
    g, k = G(u)
    check((g == child(u, 0)) == (u % 9 in (1, 2, 8)), "G vs child at %d" % u)
    if u % 9 in (4, 7):
        check(k == 0 and g % 2 == 0, "k=0 even image at %d" % u)
    if u % 9 == 5:
        check(k == 3 and g == child(u, 1), "u=5 mod 9 at %d" % u)
print("G(u) = child(u,0) iff u = 1,2,8 mod 9 for all odd u < 10^5; u = 4,7 mod 9 -> k=0 even image; u = 5 mod 9 -> G = child(u,1): OK")
check(G(17)[0] == child(17, 0) == 11, "u=17 near miss")
check(child(11, 0) == 7 and G(11)[0] == 7 and 7 % 9 not in (1, 2, 8), "m(11)=G(11)=7 leaves the coincidence locus")
XB = 10 ** 5
# guarded closure of 1 in [1,XB] == {y : the forward Collatz orbit of y stays <= XB}: forward simulation, no closure code
guarded = set()
for y in range(1, XB + 1):
    x = y; ok = True
    while x != 1:
        x = x // 2 if x % 2 == 0 else 3 * x + 1
        if x > XB:
            ok = False; break
    if ok:
        guarded.add(y)
check(len(guarded) == 39706 and sum(1 for y in guarded if y % 3) == 26478, "guarded %d" % len(guarded))
miss_g = [y for y in range(1, XB + 1) if y % 3 and y not in guarded][:8]
check(miss_g == [703, 871, 937, 1055, 1249, 1307, 1406, 1471], "guarded misses %s" % miss_g)
print("guarded closure of 1 in [1,10^5] (= forward Collatz orbit stays in the box): %d nodes, %d non-multiples of 3; first misses %s: OK" %
      (len(guarded), 26478, miss_g))
# TO-1-in-E inside [1,XB]: BFS backwards from 1 along E-predecessors 2x and (x-1)/3 (x = 1 mod 3), non-multiples of 3 only
to1 = {1}; st = [1]
while st:
    x = st.pop()
    for y in ([2 * x] + ([(x - 1) // 3] if x % 3 == 1 and x > 1 else [])):
        if y <= XB and y % 3 and y not in to1:
            to1.add(y); st.append(y)
# FROM-1-in-E inside [1,XB]: BFS forwards from 1 along E-arrows x -> 3x+1 (all x) and x -> x/2 (even x)
from1 = {1}; st = [1]
while st:
    x = st.pop()
    for y in ([3 * x + 1] + ([x // 2] if x % 2 == 0 else [])):
        if y <= XB and y not in from1:
            from1.add(y); st.append(y)
check(all(y % 3 for y in from1), "FROM-1 set avoids multiples of 3 (3Z is forward-closed in E)")
check(len(to1) == 27472, "TO-1 size %d" % len(to1))
check(all(y in to1 for y in guarded if y % 3), "guarded subset of TO-1")
miss_to1 = [y for y in range(1, XB + 1) if y % 3 and y not in to1][:8]
check(miss_to1 == [1535, 2047, 2207, 2287, 2303, 2527, 2687, 2815], "TO-1 misses %s" % miss_to1)
miss_from1 = [y for y in range(1, XB + 1) if y % 3 and y not in from1][:8]
giant = to1 & from1
print("TO-1-in-E closure (predecessor closure of 1, the lane's 'unguarded closure'): %d of 66667; first misses %s" % (len(to1), miss_to1))
print("FROM-1-in-E closure (successor closure of 1 = the Q2 set inside the box): %d of 66667; first misses %s" % (len(from1), miss_from1))
print("intersection (giant SCC of E|[1,10^5] containing 1): %d  (E-graph lane .out S3: giant size 17077)" % len(giant))
check(len(giant) == 17077, "giant SCC size %d" % len(giant))
for y in miss_to1[:5]:
    check(y in from1, "E-lane outsider %d should be FROM-1 reachable (its greedy inverse peak < N)" % y)
print("the five smallest TO-1 misses 1535, 2047, 2207, 2287, 2303 ARE in the FROM-1 set (E-lane: their binding constraint is the")
print("   forward Collatz peak).  So the lane's 'unguarded closure' illustrates Q1 in the E-graph lane's sense (every n reaches 1 in E),")
print("   NOT Q2 (1 reaches every m in E).  Direction error in the lane note's finite illustration: REFUTED as labelled, numbers correct.")
# no E-arrow enters 3Z from outside: 3n+1 is never 0 mod 3, and n/2 = 0 mod 3 forces n = 0 mod 3 (3Z is backward-closed,
# not forward-closed: 3 -> 10 leaves it).  Hence a path from 1 never meets a multiple of 3 and the clause is automatic.
for n in range(1, 20000):
    check((3 * n + 1) % 3 != 0, "3n+1 at %d" % n)
    if n % 2 == 0 and (n // 2) % 3 == 0:
        check(n % 3 == 0, "halving into 3Z at %d" % n)
check((3 * 3 + 1) == 10 and 10 % 3 == 1, "3Z is not forward-closed (3 -> 10)")
print("no E-arrow enters 3Z from outside (3n+1 != 0 mod 3; n/2 = 0 mod 3 only if n = 0 mod 3), while 3 -> 10 leaves 3Z;")
print("   so 'through non-multiples of 3' in Q2 is automatic (the E-lane's transient-singleton statement), not an extra hypothesis: OK")
p7 = [7]
while p7[-1] != 1:
    p7.append(p7[-1] // 2 if p7[-1] % 2 == 0 else 3 * p7[-1] + 1)
check(p7 == [7, 22, 11, 34, 17, 52, 26, 13, 40, 20, 10, 5, 16, 8, 4, 2, 1], "path of 7")
check(G(7) == (2, 0) and G(2) == (1, 1), "G-orbit of 7")
print("m = 7 example paths: OK")

# ------------------------------------------------------------------ A11
banner("A11  section 5: image law, orbits to 10^6, survival classes, parity independence")
m = lambda u: child(u, 0)
img = set(m(u) for u in range(1, 400001, 2) if u % 3)
for n in range(1, 200001, 2):
    check((n in img) == (n % 8 != 5), "image law at %d" % n)
    if n % 8 == 5:
        check((n - 1) // 4 % 2 == 1, "5 mod 8 = 4*(odd)+1")
print("image(m) = odd n != 5 mod 8 for odd n < 2*10^5 (domain odd u < 4*10^5): OK; 5 mod 8 = 4*odd + 1: OK")
check(m(1) == 1 and all(m(u) != u for u in range(3, 10001, 2) if u % 3), "fixed points")
for u in range(1, 3000, 2):
    if u % 3 == 0:
        continue
    check(T(m(u)) == u, "T o m")
    if u % 3 == 1 and u > 1:
        check(m(u) > u, "m(u) > u")
    if u % 3 == 2:
        check(m(u) < u, "m(u) < u")
t0 = time.time()
lengths = Counter(); best = (-1, None); bestpk = (0, None)
for u in range(1, X + 1, 2):
    if u % 3 == 0:
        continue
    if u == 1:
        lengths[0] += 1; continue
    x, L, pk = u, 0, u
    while x % 3:
        x = m(x); L += 1
        if x > pk: pk = x
        if L > 10000: fail("cap")
    lengths[L] += 1
    if L > best[0]: best = (L, u)
    if pk / u > bestpk[0]: bestpk = (pk / u, u)
check(best == (29, 774883), "longest %s" % (best,))
check(abs(bestpk[0] - 33.26) < 5e-3 and bestpk[1] == 86023, "peak %s" % (bestpk,))
exp_len = {1: 111112, 2: 74074, 3: 49383, 4: 32924, 5: 21949, 6: 14624, 7: 9754, 8: 6498, 9: 4336, 10: 2889, 11: 1931, 12: 1286}
for L, c in exp_len.items():
    check(lengths[L] == c, "L=%d count %d" % (L, lengths[L]))
check(sum(c for L, c in lengths.items() if L >= 13) == 2572, "tail")
print("m-orbits of 333333 odd non-multiples of 3 <= 10^6: longest L=29 at 774883, peak/u max %.2f at %d, length table matches: OK (%.1fs)" %
      (bestpk[0], bestpk[1], time.time() - t0))
# exact geometric law is a class count: #{u <= X in the L-classes} -- compare with the residue-class prediction (1/3)(2/3)^(L-1)
for L in range(1, 9):
    M = 3 ** (L + 1)
    surv_odd = 0; surv_even = 0; tot = 0
    for a in range(1, M):
        if a % 3 == 0:
            continue
        tot += 1
        for rep, acc in ((a if a % 2 else a + M), 'odd'), ((a + M if a % 2 else a), 'even'):
            x = rep; ok = True
            for _ in range(L):
                x = ((1 << (h0(x) + 1)) * x - 1) // 3
                if x % 3 == 0:
                    ok = False; break
            if acc == 'odd': surv_odd += ok
            else: surv_even += ok
    check(surv_odd == surv_even == 2 ** (L + 1) and tot == 2 * 3 ** L, "survival L=%d: %d %d" % (L, surv_odd, surv_even))
print("survival classes mod 3^(L+1): exactly 2^(L+1) of 2*3^L survive L steps for L<=8, identical for odd and even representatives: OK")
# lifting lemma directly: m(a + t 3^s) = m(a) + 2^(h0+1) t 3^(s-1)
for s in range(1, 6):
    for a in range(1, 3 ** s):
        if a % 3 == 0: continue
        for t in range(3):
            lhs = ((1 << (h0(a) + 1)) * (a + t * 3 ** s) - 1) // 3
            rhs = ((1 << (h0(a) + 1)) * a - 1) // 3 + (1 << (h0(a) + 1)) * t * 3 ** (s - 1)
            check(lhs == rhs, "lifting lemma")
print("lifting lemma m(a + t 3^s) = m(a) + 2^(h0+1) t 3^(s-1) (so the three lifts of a class mod 3^s map onto the three lifts of m(a) mod 3^(s-1)): OK")

# ------------------------------------------------------------------ A12
banner("A12  boundary cases")
check(child(1, 0) == 1 and G(1) == (1, 2) and m(1) == 1 and T(1) == 1, "u=1")
check(T(9) == 7 and 9 % 3 == 0, "n=9 is a leaf with parent 7")
check(T(3) == 5 and T(5) == 1 and T(21) == 1, "small leaves")
# even 'u' is outside the domain: the formula still returns an integer but T does not invert it
check(((1 << 1) * 2 - 1) // 3 == 1 and T(1) != 2, "u=2 is not a valid parent: (2*2-1)/3 = 1 but T(1) = 1")
check(((1 << 2) * 4 - 1) // 3 == 5 and T(5) == 1 != 4, "u=4: (4*4-1)/3 = 5 but T(5) = 1 (the E-lane index j=-1 phenomenon)")
check(child(5, 0) == 3 and child(7, 0) == 9, "u=5,7 mod 9 minimal child is a leaf")
print("u=1 loop, leaves 3, 9, 21, even u outside the domain (child formula not inverted by T), minimal leaf children of 5 and 7: OK")
print()
print("ALL AUDIT CHECKS PASSED")
