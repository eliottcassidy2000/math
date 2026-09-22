#!/usr/bin/env python3
"""collatz_mod6_20260922_w6_sign_specific_probes_audit.py

Independent recomputation of the key numbers of lane sign_specific_probes
(collatz_mod6_20260922_w6_sign_specific_probes.py) with different code paths:
  A0  square-sum graph Q_n: components, leaves, Hamiltonian existence (plain
      backtracking, no Warnsdorff), ALL Hamiltonian paths of Q_15 (uniqueness and
      the squares used), the pasted-Lean loops
  A1  root asymmetry: the n = 7 witness against "1 -> 4 is the only unhalved
      arrow with companion a multiple of 3"; the (a,b) fixed-point table
  A2  (3,b,k) edge counts by enumeration over TARGETS y (not sources), C_edge,
      per-k differences, the B3 mixing counts
  A3  T_- basins to 10^6 (C-form iteration, dict memo), greedy Q2_- histogram,
      the exact one-step law 444444, the E_- paths
  A4  carry B_L against direct iteration; EXACT equality of the |B_10| mod 4 and
      mod 3 counts on complete residue classes (K_10 <= 18, odd n < 2^19)
  A5  greedy mean k_1 as an exact rational from residue counts mod 9
      (77779/66667 vs 77778/66667); 2-cycles of the inverse relation for ALL
      exponents (plus: divisor argument; minus: a + c <= 3)
  A6  residue-of-minimum law on the cycles; last-odd-before-1 census to 10^5
All checks use explicit raise.
"""
import math
import sys
from fractions import Fraction


def check(c, m):
    if not c:
        raise RuntimeError("AUDIT CHECK FAILED: " + m)


def v2(n):
    return (n & -n).bit_length() - 1


def Tb(n, b):
    m = 3 * n + b
    k = v2(m)
    return m >> k, k


print("A0  square-sum graph Q_n (x ~ y iff x + y square, x != y), n <= 32")


def adj_of(n):
    A = {x: set() for x in range(1, n + 1)}
    for x in range(1, n + 1):
        for y in range(x + 1, n + 1):
            r = math.isqrt(x + y)
            if r * r == x + y:
                A[x].add(y)
                A[y].add(x)
    return A


def ncomp(A):
    seen = set()
    c = 0
    for s in A:
        if s in seen:
            continue
        c += 1
        st = [s]
        seen.add(s)
        while st:
            u = st.pop()
            for w in A[u]:
                if w not in seen:
                    seen.add(w)
                    st.append(w)
    return c


def ham_count(A, n, stop_at=None, budget=20_000_000):
    """count directed Hamiltonian paths by plain backtracking (no heuristics)."""
    cnt = [0]
    nodes = [0]
    vis = [False] * (n + 1)

    def rec(u, depth):
        nodes[0] += 1
        if nodes[0] > budget:
            raise RuntimeError("budget")
        if depth == n:
            cnt[0] += 1
            return stop_at is not None and cnt[0] >= stop_at
        for w in A[u]:
            if not vis[w]:
                vis[w] = True
                if rec(w, depth + 1):
                    return True
                vis[w] = False
        return False

    for s in range(1, n + 1):
        vis[s] = True
        if rec(s, 1):
            break
        vis[s] = False
    return cnt[0]


comps = {}
leaves = {}
ham = {}
for n in range(1, 33):
    A = adj_of(n)
    comps[n] = ncomp(A)
    leaves[n] = sorted(x for x in A if len(A[x]) <= 1)
    if n > 1 and (any(len(A[x]) == 0 for x in A) or len(leaves[n]) >= 3):
        ham[n] = False
    else:
        ham[n] = ham_count(A, n, stop_at=1) >= 1
print("components:", [(n, comps[n]) for n in range(1, 33)])
print("leaves n>=14:", [(n, leaves[n]) for n in range(14, 33)])
print("Hamiltonian path exists:", [n for n in range(1, 33) if ham[n]])
check([n for n in range(1, 33) if ham[n]] == [1, 15, 16, 17, 23] + list(range(25, 33)), "ham list")
check(all(comps[n] == 3 for n in range(4, 13)) and comps[13] == 2 and all(comps[n] == 1 for n in range(14, 33)), "comps")
check(leaves[18] == [16, 17, 18] and leaves[19] == [16, 18] and all(leaves[n] == [18] for n in range(20, 31))
      and leaves[31] == [] and leaves[32] == [], "leaves")
check(leaves[14] == [8, 9, 10] and leaves[15] == [8, 9] and leaves[16] == [8, 16] and leaves[17] == [16, 17], "leaves 14-17")
A15 = adj_of(15)
tot15 = ham_count(A15, 15)
print("directed Hamiltonian paths of Q_15:", tot15, "(= 2 means one undirected path)")
check(tot15 == 2, "Q_15 unique path")
# squares available in Q_15 and the squares on the unique path
lead15 = [8, 1, 15, 10, 6, 3, 13, 12, 4, 5, 11, 14, 2, 7, 9]
sq15 = sorted(set(lead15[i] + lead15[i + 1] for i in range(14)))
edges15 = sum(len(A15[x]) for x in A15) // 2
sq_avail = sorted(set(x + y for x in A15 for y in A15[x]))
print("Q_15: %d edges, squares available %s, squares on the unique path %s (edge 1-3 with sum 4 unused)" % (edges15, sq_avail, sq15))
check(sq_avail == [4, 9, 16, 25] and sq15 == [9, 16, 25], "squares Q_15")
A23 = adj_of(23)
tot23 = ham_count(A23, 23, stop_at=3)
print("Q_23 has at least %d directed Hamiltonian paths (search stopped at 3)" % tot23)
check(tot23 >= 2, "Q_23 path exists")
loops = [v for v in range(1, 33) if math.isqrt(2 * v) ** 2 == 2 * v]
print("pasted Lean loops (2v square), v <= 32:", loops)
check(loops == [2, 8, 18, 32], "loops")
print("3,7,11,17 gaps:", [4, 4, 6])

print()
print("A1  root asymmetry")
for n in (1, 3, 5, 7):
    print("  n=%d: halved c+=%d c-=%d unhalved c+=%d c-=%d" % (n, (n + 1) // 2, (n - 1) // 2, 2 * n + 1, 2 * n - 1))
odd_mult3 = [n for n in range(1, 30, 2) if (2 * n + 1) % 3 == 0]
print("  odd n < 30 whose unhalved plus companion 2n+1 is a multiple of 3 (n = 1 mod 3):", odd_mult3)
check(odd_mult3 == [1, 7, 13, 19, 25], "companion multiple of 3 is NOT unique to n = 1: witness 7 -> 22 = 7 + 15")
fam = {}
for a, b in [(3, 1), (3, -1), (1, 1), (1, -1), (5, 1), (5, -1), (7, 1), (7, -1)]:
    s = a + b
    fam[(a, b)] = s > 0 and s & (s - 1) == 0
print("  1 fixed by (a n + b)/2^v2:", fam)
check(fam[(1, 1)] and fam[(3, 1)] and fam[(3, -1)] and not fam[(1, -1)] and not fam[(5, 1)] and fam[(5, -1)], "family")

print()
print("A2  (3,b,k) edges x -> y, 3x + b = 2^k y, x, y odd, x != y, (x^2 + y^2)/2 <= X; enumeration over targets y")


def edges_by_target(X, b):
    out = {}
    ymax = math.isqrt(2 * X) + 1
    for y in range(1, ymax + 1, 2):
        k = 1
        while True:
            m = (y << k) - b
            if m > 3 * ymax + 3 and y * y + (m // 3) ** 2 > 2 * X:
                break
            if m > 0 and m % 3 == 0:
                x = m // 3
                if x % 2 == 1 and x != y and x * x + y * y <= 2 * X:
                    out.setdefault(k, []).append((x, y))
            k += 1
            if k > 60:
                break
    return out


C_edge = sum(math.sqrt(2.0 / (1.0 + 9.0 / 4 ** k)) / 2 ** (k + 1) for k in range(1, 80))
print("  C_edge = %.6f" % C_edge)
check(abs(C_edge - 0.507819) < 1e-6, "C_edge")
for X in (10 ** 6, 10 ** 7):
    E = {b: edges_by_target(X, b) for b in (1, -1)}
    tot = {b: sum(len(v) for v in E[b].values()) for b in (1, -1)}
    tri = {b: len(set((max(x, y), min(x, y)) for v in E[b].values() for x, y in v)) for b in (1, -1)}
    ks = sorted(set(E[1]) | set(E[-1]))
    perk = [(k, len(E[1].get(k, [])), len(E[-1].get(k, []))) for k in ks]
    print("  X=%d: plus %d (tri %d) minus %d (tri %d); per k (k, plus, minus): %s" % (X, tot[1], tri[1], tot[-1], tri[-1], perk))
    check(all(abs(p - m) <= 1 for _, p, m in perk), "per-k diff in {-1,0,1}")
    check(abs(tot[1] - tot[-1]) <= 2 * math.log2(3 * math.sqrt(2 * X) + 1), "note's O(log X) bound")
    if X == 10 ** 6:
        check((tot[1], tri[1], tot[-1], tri[-1]) == (507, 507, 506, 505), "507/507/506/505")
        union = set((max(x, y), min(x, y)) for b in (1, -1) for v in E[b].values() for x, y in v)
        check(len(union) == 1012, "union 1012")
        print("  union of triangles: %d" % len(union))
        check([p for _, p, _ in perk] == [196, 141, 82, 44, 22, 11, 5, 3, 1, 1, 0, 1], "plus per-k at 10^6")
    else:
        check((tot[1], tri[1], tot[-1], tri[-1]) == (1605, 1605, 1605, 1604), "1605/1605/1605/1604")
        check([p for _, p, _ in perk] == [620, 447, 261, 138, 69, 35, 17, 9, 4, 3, 1, 1, 0], "plus per-k at 10^7")
        check([m for _, _, m in perk] == [619, 447, 262, 137, 70, 35, 18, 8, 5, 2, 1, 0, 1], "minus per-k at 10^7")
    # B3 mixing: k=1 edges of sheet b whose child (2x+b) -> y (k=2, sheet -b) is also below X
    for b in (1, -1):
        k1 = E[b].get(1, [])
        inside = 0
        for x, y in k1:
            u = 2 * x + b
            check((3 * u - b) == 4 * y, "B3 identity")
            if u * u + y * y <= 2 * X:
                inside += 1
        k2other = set(E[-b].get(2, []))
        check(all(((2 * x + b, y) in k2other) == (u2 * u2 + y * y <= 2 * X) for x, y in k1 for u2 in [2 * x + b]), "child membership")
        print("  X=%d sheet %+d: %d k=1 edges, %d children below X" % (X, b, len(k1), inside))
        if X == 10 ** 6:
            check((len(k1), inside) == ((196, 141) if b == 1 else (195, 141)), "mixing counts 10^6")
        else:
            check((len(k1), inside) == ((620, 447) if b == 1 else (619, 447)), "mixing counts 10^7")
    # every k=2 edge of sheet -b is a B3 child of a k=1 edge of sheet b (converse, on the window)
    for b in (1, -1):
        for u, y in E[-b].get(2, []):
            x = (u - b) // 2
            check((u - b) % 2 == 0 and x % 2 == 1 and 3 * x + b == 2 * y, "converse B3 (control S6)")
    print("  converse B3 (every k=2 edge is a child of a k=1 edge of the other sheet) verified on the window")
print("  7 -> 11 (plus, k=1): triangle (%d, %d, %d)" % (7 * 11, (121 - 49) // 2, (121 + 49) // 2))
check((7 * 11, (121 - 49) // 2, (121 + 49) // 2) == (77, 36, 85), "paste PPT")

print()
print("A3  T_- basins, greedy Q2_-, E_- paths")
X3 = 10 ** 6
cyc_min = {1: 1, 5: 5, 7: 5, 17: 17, 25: 17, 37: 17, 55: 17, 41: 17, 61: 17, 91: 17}
first_below = {}
basin = {}
maxdrop = 0
for n in range(1, X3 + 1, 2):
    if n in cyc_min:
        basin[n] = cyc_min[n]
        continue
    # C-form: 3n-1 then halvings, count accelerated steps until < n or a cycle member
    m = n
    st = 0
    while m >= n and m not in cyc_min:
        m = 3 * m - 1
        while m % 2 == 0:
            m //= 2
        st += 1
    maxdrop = max(maxdrop, st)
    basin[n] = cyc_min[m] if m in cyc_min else basin[m]
cnt = {1: 0, 5: 0, 17: 0}
for n in range(1, X3 + 1, 2):
    cnt[basin[n]] += 1
print("  basins odd n <= 10^6: %s, max accelerated steps to first drop %d" % (cnt, maxdrop))
check(cnt == {1: 163486, 5: 162122, 17: 174392} and maxdrop == 125, "basins")


def Gm(m):
    # minus-sheet greedy: (2^k m + 1)/3 integer, not 0 mod 3, k minimal
    k = 0
    while True:
        t = (m << k) + 1
        if t % 3 == 0 and (t // 3) % 3 != 0:
            return t // 3, k
        k += 1


hist = {}
cyc = []
for m in range(2, X3 + 1):
    if m % 3 == 0:
        continue
    x = m
    st = 0
    seen = set()
    while x >= m:
        if x in seen:
            cyc.append(m)
            break
        seen.add(x)
        x, _ = Gm(x)
        st += 1
    else:
        hist[st] = hist.get(st, 0) + 1
print("  greedy Q2_-: certified %d, cycled %s, histogram %s max %d" % (sum(hist.values()), cyc, sorted(hist.items())[:6], max(hist)))
check(sum(hist.values()) == 666665 and cyc == [4] and hist[1] == 444444 and hist[2] == 148147 and hist[3] == 24692
      and hist[4] == 24691 and hist[5] == 5487 and hist[6] == 8230 and max(hist) == 35, "greedy histogram")
one_step = sum(1 for m in range(2, X3 + 1) if m % 9 in (1, 2, 5, 7))
print("  exact one-step law: #m in [2,10^6] with m mod 9 in {1,2,5,7} = %d" % one_step)
check(one_step == 444444, "444444")
esc5 = [5, 14, 7, 20, 10, 29, 86, 43, 128, 64, 32, 16, 8, 4, 2, 1]
esc17 = [17, 50, 25, 74, 37, 110, 55, 164, 82, 41, 122, 61, 182, 91, 272, 136, 68, 203, 608, 304, 911, 2732, 1366, 683,
         2048, 1024, 512, 256, 128, 64, 32, 16, 8, 4, 2, 1]
p4 = [1, 2, 5, 14, 7, 20, 10, 29, 86, 43, 128, 64, 32, 16, 8, 4]
for path in (esc5, esc17, p4):
    for i in range(len(path) - 1):
        u, w = path[i], path[i + 1]
        check((u % 2 == 0 and w == u // 2) or w == 3 * u - 1, "E_- arrow %d -> %d" % (u, w))
    ev = [(path[i], path[i + 1]) for i in range(len(path) - 1) if path[i] % 2 == 0 and path[i + 1] == 3 * path[i] - 1]
    print("  path %d -> %d: %d arrows, even->3n-1 arrows %s" % (path[0], path[-1], len(path) - 1, ev))
check(len(esc5) - 1 == 15 and len(esc17) - 1 == 35 and len(p4) - 1 == 15, "path lengths")
print("  SCC note section 6 already states Q1_- and Q2_- to 10^6 (inherited; this lane reconfirms)")

print()
print("A4  carry B_L")


def carry_direct(n, b, L):
    """2^K T^L(n) - 3^L n by direct iteration, plus the word."""
    x = n
    K = 0
    w = []
    for _ in range(L):
        x, k = Tb(x, b)
        K += k
        w.append(k)
    return (x << K) - 3 ** L * n, K, tuple(w)


def carry_formula(w, b):
    L = len(w)
    B = 0
    K = 0
    for i, k in enumerate(w):
        B += 3 ** (L - 1 - i) << K
        K += k
    return b * B, K


rows = []
for b, n0, L in [(1, 1, 1), (-1, 1, 1), (-1, 5, 2), (-1, 17, 7), (-1, 17, 14), (-1, 17, 21), (1, 27, 8), (1, 837799, 8)]:
    Bd, K, w = carry_direct(n0, b, L)
    Bf, Kf = carry_formula(w, b)
    check((Bd, K) == (Bf, Kf), "carry formula vs direct")
    check(Bd % 2 == 1, "B odd")
    rows.append((b, n0, L, Bd, K))
print("  formula = direct on", rows)
check(rows[3][3] == -2363 and rows[4][3] == -10007305 and rows[5][3] == -31797116387, "seven-cycle carries")
for (b, n0, L, B, K) in rows[:6]:
    check(Fraction(B, 2 ** K - 3 ** L) == n0, "gate")
# B mod 3 = b * 2^(K_(L-1)) mod 3 : parity of K_(L-1)
for n in range(1, 2001, 2):
    for b in (1, -1):
        B, K, w = carry_direct(n, b, 6)
        check(B % 3 == (b * 2 ** (K - w[-1])) % 3, "B mod 3 law")
print("  B_L mod 3 = b 2^(K_(L-1)) mod 3 (parity of K_(L-1)) checked, n <= 2000, L = 6")
# exact equality of the |B_10| mod 4 / mod 3 counts on COMPLETE residue classes
J = 19
LL = 10
Kcap = J - 1
cnt4 = {}
cnt3 = {}
for b in (1, -1):
    c4 = {}
    c3 = {}
    for n in range(1, 1 << J, 2):
        B, K, w = carry_direct(n, b, LL)
        if K <= Kcap:
            a = abs(B)
            c4[a % 4] = c4.get(a % 4, 0) + 1
            c3[a % 3] = c3.get(a % 3, 0) + 1
    cnt4[b] = sorted(c4.items())
    cnt3[b] = sorted(c3.items())
print("  odd n < 2^%d with K_10 <= %d (complete classes mod 2^(K+1)): |B| mod 4 plus %s minus %s; |B| mod 3 plus %s minus %s"
      % (J, Kcap, cnt4[1], cnt4[-1], cnt3[1], cnt3[-1]))
check(cnt4[1] == cnt4[-1] and cnt3[1] == cnt3[-1], "exact mirror on complete classes")

print()
print("A5  greedy means as exact rationals; 2-cycles for all exponents")
kp = {1: 2, 2: 1, 4: 0, 5: 3, 7: 0, 8: 1}
km = {1: 1, 2: 0, 4: 3, 5: 0, 7: 1, 8: 2}
rc = {r: sum(1 for m in range(1, 100001) if m % 9 == r) for r in kp}
print("  residue counts mod 9 in [1,10^5]:", rc)
mp = Fraction(sum(kp[r] * rc[r] for r in kp), sum(rc.values()))
mm = Fraction(sum(km[r] * rc[r] for r in km), sum(rc.values()))
print("  mean k_1: plus %s = %.6f, minus %s = %.6f, difference %s" % (mp, float(mp), mm, float(mm), mp - mm))
check(mp == Fraction(77779, 66667) and mm == Fraction(77778, 66667), "means")
# 2-cycles m(2^(a+c) - 9) = b(2^c + 3): plus, all a, c >= 0
plus = []
for c in range(0, 6):
    for a in range(0, 10):
        s = a + c
        den = 2 ** s - 9
        if den > 0 and (2 ** c + 3) % den == 0:
            plus.append(((2 ** c + 3) // den, a, c))
print("  plus 2-cycles found with c <= 5, a <= 9:", plus)
print("  PROVED complete: c >= 4 needs 2^c (2^a - 1) <= 12 (a >= 1) or 2^c - 9 | 12 (a = 0), impossible; c <= 3 is the census")
check(plus == [(1, 2, 2)], "plus 2-cycle unique")
minus = []
for s in range(0, 4):
    for c in range(0, s + 1):
        a = s - c
        den = 2 ** s - 9
        num = -(2 ** c + 3)
        if num % den == 0 and num // den > 0:
            m = num // den
            m2 = 2 ** a * m + 1
            if m2 % 3 == 0 and (2 ** c * (m2 // 3) + 1) == 3 * m:
                minus.append((m, m2 // 3, a, c))
print("  minus 2-cycles (a + c <= 3 is exhaustive since 2^(a+c) < 9):", minus)
check(sorted(set((min(x, y), max(x, y)) for x, y, _, _ in minus)) == [(1, 1), (4, 11), (5, 7)], "minus 2-cycles")

print()
print("A6  residue-of-minimum law; last odd before 1")
for cyc_ in [(5, 7), (17, 25, 37, 55, 41, 61, 91)]:
    n0, M = min(cyc_), max(cyc_)
    _, k0 = Tb(n0, -1)
    _, kM = Tb(M, -1)
    check(n0 % 4 == 1 and k0 == 1 and M % 4 == 3 and kM >= 2, "law on minus cycles")
    print("  cycle %s: min %d = %d mod 4 k=%d, max %d = %d mod 4 k=%d" % (cyc_, n0, n0 % 4, k0, M, M % 4, kM))
# fixed points: n(2^k - 3) = b has n = 1 only
fp = [(b, k) for b in (1, -1) for k in range(1, 5) if (b % (2 ** k - 3) == 0 if 2 ** k != 3 else False) and b // (2 ** k - 3) == 1]
print("  fixed points n = b/(2^k - 3) = 1:", fp)
check(fp == [(1, 2), (-1, 1)], "fixed points")
h = {}
for n in range(3, 10 ** 5 + 1, 2):
    x = n
    prev = x
    while x != 1:
        prev = x
        x, _ = Tb(x, 1)
    h[prev % 4] = h.get(prev % 4, 0) + 1
print("  last odd before 1, odd 3 <= n <= 10^5:", sorted(h.items()))
check(h == {1: 49999}, "last odd 1 mod 4")
print("  (classical: a 3n+1 cycle minimum is 3 mod 4 and its maximum 1 mod 4; elementary, UNCITED-RECOLLECTION as folklore)")

print()
print("ALL AUDIT CHECKS PASSED")
