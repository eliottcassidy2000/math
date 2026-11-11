"""aud_defg.py - audit tasks 2(d), 2(e), 2(f), 2(g) with the auditor's own counters (aud_core, aud_pairstate).

D   D4: order-8 pair (code 14705281, O = {1,3,4}): state (-1, 1024, 217, 64), out-lists of the theorem, closed
    tournaments with 2^(10+t) - 1 paths (t = 1, 2, 3); the tower states of D4 by the level map and, for n = 1 and 7,
    D2 / D3 literal spot checks.
E   H3: the five 9-grid pairs (levels 7..11), closed counts for r = 1, 2; the 43-grid carrier (root code
    45174751748, root 2, port 8), the four low sources, the head 175; H1 examples (5,9), (17,18), (33,36).
F   H2: domination ray on all labelled pairs of order <= 4 (t = 1, 2, 3) by brute-force dominating-set counting;
    M >= delta_dom on all labelled pairs of order <= 5; the inequality H(T) - H(T - z) >= d(G_z) - 1 on every
    (T, z), T labelled of order <= 5, and on the closed tournaments.
G   H4: H(X) >= H(X[S]) H(X - S); n + 1 <= Pi_2(X)/H(X) <= 2^n, equality on the right iff transitive, on all
    labelled tournaments of order <= 5 (permutation enumeration for H of every induced sub-tournament);
    multiplicativity over X => Y; merge map injective (explicit merge, all pairs of paths) for order <= 5.
"""
import itertools, random, subprocess, os, sys, time
from aud_core import *

HERE = os.path.dirname(os.path.abspath(__file__))
EXE = os.path.join(HERE, "aud_pairstate.exe")
allok = True


def report(name, ok, extra=""):
    global allok
    allok &= bool(ok)
    print(("PASS  " if ok else "FAIL  ") + name + ("  " + extra if extra else ""))
    sys.stdout.flush()


def rec(mode, A, O):
    n = len(A)
    return " ".join([str(mode), str(n)] + [str(sum(1 << j for j in range(n) if A[i][j])) for i in range(n)]
                    + [str(sum(1 << o for o in O))])


def cstate(A, O, mode=0):
    out = subprocess.run([EXE], input=rec(mode, A, O) + "\n", capture_output=True, text=True).stdout.strip()
    if mode == 1:
        a, b = out.split("|")
        a = tuple(map(int, a.split()))[1:]
        b = tuple(map(int, b.split()))
        assert a == b, (a, b)
        return a
    return tuple(map(int, out.split()))[1:]


def cH(A):
    out = subprocess.run([EXE], input=rec(2, A, []) + "\n", capture_output=True, text=True).stdout.strip()
    return int(out.split()[1])


def odd_part(x):
    v = 0
    while x % 2 == 0:
        x //= 2
        v += 1
    return x, v


# ====================================================================== D
A8 = adj_from_code(8, 14705281)
outl = {0: [1], 1: [2, 4], 2: [0, 3, 4], 3: [0, 1, 7], 4: [0, 3, 5, 6], 5: [0, 1, 2, 3], 6: [0, 1, 2, 3, 5],
        7: [0, 1, 2, 4, 5, 6]}
ok = all(sorted(j for j in range(8) if A8[i][j]) == outl[i] for i in range(8))
report("D0 out-lists of D4 in THM-5130 = code 14705281", ok)
O8 = [1, 3, 4]
s_b = brute_state(A8, O8)
s_d = dp_state(A8, O8)
s_c = cstate(A8, O8, 1)
report("D1 state of the order-8 pair = (-1, 1024, 217, 64) (permutations, Python recursion, C recursion + literal)",
       s_b == s_d == s_c == (-1, 1024, 217, 64), str(s_b))
ok = True
vals = []
for t in (1, 2, 3):
    T = close(A8, O8, transitive(t))
    h1 = cH(T)
    h2 = dp_H(T)
    vals.append(h1)
    ok &= (h1 == h2 == 2 ** (10 + t) - 1 and len(T) == 9 + t)
hb = brute_H(close(A8, O8, transitive(1)))
ok &= (hb == 2047)
report("D2 closed tournaments on 9 + t vertices have 2^(10+t) - 1 paths, t = 1, 2, 3 (t = 1 also by 10! permutations)",
       ok, str(vals))


def level(s, delta, i, m):
    r, M, h, E = s
    a = i + delta + m
    al, ga = 2 ** delta + 1, 3 ** delta + 2 ** delta
    return (-r, al * 2 ** m * (h + r) + 2 * ga * 2 ** (i - 1) * 3 ** m * M - (E if i == 1 else 0), r + 2 ** a * M,
            al * 2 ** (i - 1) * M + ((r - h) if m == 0 else 0))


ok = True
for n in range(1, 8):
    s = (1, 2 ** (n - 1), 1, 1 if n == 1 else 0)
    order = n
    trace = [s]
    for (d, i, m) in ((1, 8 - n, 0), (1, 4, 0), (2, 2, 1)):
        s = level(s, d, i, m)
        order += i + d + m + 1
        trace.append(s)
    ok &= (s[0] == -1 and s[1] == 2 ** 23 and order == 22 and trace[1][1] == 646)
report("D3 D4 towers over (TT_n, {1}), n = 1..7, levels (1,8-n,0),(1,4,0),(2,2,1): offset -1, slope 2^23, base order 22, M_1 = 646 (level map)", ok)
# literal check of the first two levels of those towers (orders 10 and 16) for n = 1, 4, 7
ok = True
for n in (1, 4, 7):
    A = transitive(n)
    O = [0]
    s = (1, 2 ** (n - 1), 1, 1 if n == 1 else 0)
    for (d, i, m) in ((1, 8 - n, 0), (1, 4, 0)):
        a = i + d + m
        nn = len(A)
        A = close(A, O, transitive(a))
        O = [nn + i - 1, nn + i - 1 + d]
        s = level(s, d, i, m)
        ok &= (cstate(A, O) == s)
report("D4 first two levels of the D4 towers recounted literally (orders 10 and 16; n = 1, 4, 7)", ok)

# D2 / D3 families, literal
rng = random.Random(77)
ok = True
cnt = 0
for _ in range(60):
    n = rng.randint(1, 4)
    p = rng.randint(1, n)
    d = rng.randint(1, 2)
    i = rng.randint(1, 3)
    m = rng.randint(0, 2)
    t = rng.randint(1, 2)
    a = i + d + m
    B1 = close(transitive(n), [p - 1], transitive(a))
    T = close(B1, [n + i - 1, n + i - 1 + d], transitive(t))
    if len(T) > 14:
        continue
    M1 = (2 ** d + 1) * 2 ** (m + 1) + 2 * (3 ** d + 2 ** d) * 2 ** (i - 1 + n - p) * 3 ** m - (1 if (i == 1 and p == n) else 0)
    cnt += 1
    ok &= (cH(T) == 2 ** t * M1 - 1)
report("D5 family D2 (two closers over a fan), literal H", ok, "%d tournaments" % cnt)
ok = True
cnt = 0
for _ in range(60):
    d1 = rng.randint(1, 2)
    i1 = rng.randint(1, 2)
    m1 = rng.randint(0, 2)
    d2 = rng.randint(1, 2)
    i2 = rng.randint(1, 3)
    m2 = rng.randint(0, 2)
    t = rng.randint(1, 2)
    a1 = i1 + d1 + m1
    a2 = i2 + d2 + m2
    B1 = close(transitive(a1), [i1 - 1, i1 - 1 + d1], transitive(a2))
    T = close(B1, [a1 + i2 - 1, a1 + i2 - 1 + d2], transitive(t))
    if len(T) > 15:
        continue
    q = (2 ** d1 + 1) * (3 ** d2 + 2 ** d2)
    M2 = 2 * q * 2 ** (m1 + i2 - 1) * 3 ** m2 - (1 if (i2 == 1 and m1 == 0) else 0)
    cnt += 1
    ok &= (cH(T) == 1 + 2 ** t * M2)
report("D6 family D3 (two closers over a transitive pair), literal H", ok, "%d tournaments" % cnt)
ok = True
for d in (1, 2):
    for i in (1, 2, 3):
        for m in (0, 1, 2):
            for t in (1, 2, 3):
                T = close(transitive(i + d + m), [i - 1, i - 1 + d], transitive(t))
                ok &= (dp_H(T) == (2 ** d + 1) * 2 ** (m + t) - 1)
report("D7 family D1 (one closer), literal H", ok)

# ====================================================================== E
print("--- 9-grid sources (levels 7..11)")
src = {7: (8, 11026436, [2], 129, 288, 129), 8: (8, 6562209, [2], 257, 288, 257), 9: (9, 442533377, [0], 513, 576, 513),
       10: (9, 6793626886, [1, 2], 1025, 4608, 1565), 11: (9, 194216770, [0], 2049, 2304, 2049)}
ok = True
for s, (m, code, O, off, slope, hB) in src.items():
    A = adj_from_code(m, code)
    st = cstate(A, O, 1)
    d, v = odd_part(st[1])
    good = (st[0] == off == 1 + 2 ** s and st[1] == slope and st[2] == hB and d == 9 and v <= s)
    if m == 8:
        good &= (brute_state(A, O) == st)
    closed = []
    for rr in (1, 2):
        t = rr + s - v
        T = close(A, O, transitive(t))
        hT = cH(T)
        closed.append((rr, t, len(T), hT))
        good &= (hT == 1 + 2 ** s + 9 * 2 ** (rr + s))
    ok &= good
    print("   s=%d: order %d code %d O=%s state (r,M,h,E)=%s  v=%d  closed (r, t, order, H)=%s  %s"
          % (s, m, code, O, st, v, closed, "ok" if good else "BAD"))
report("E1 five 9-grid pairs: offset 1 + 2^s, slope 9 * 2^v (v <= s), H(B) as tabulated; N(9, r, s) for r = 1, 2", ok)

print("--- H1 examples")
ok = True
pt = [[0]]
st = dp_state(close(pt, [0], transitive(2)), [1])           # (point, a = 2, p = 1)
ok &= (st[0] == 5 and st[1] == 9)
st2 = dp_state(close(pt, [0], transitive(4)), [1 + 3])       # (point, a = 4, p = 4)
ok &= (st2[0] == 17 and st2[1] == 18)
A4 = adj_from_code(4, 34)
s4 = brute_state(A4, [2, 3])
st3 = dp_state(close(A4, [2, 3], transitive(1)), [4])
ok &= (s4[:3] == (1, 16, 5) and st3[0] == 33 and st3[1] == 36)
print("   (point; a=2,p=1) ->", st, " (point; a=4,p=4) ->", st2, " seed (4, code 34, O={2,3}) =", s4, "; a=p=1 ->", st3)
report("E2 H1: the pairs (5, 9), (17, 18), (33, 36) are one-vertex gadgets on seeds of slope 2^k", ok)
# H1 formula on random bases
rng = random.Random(5)
ok = True
for _ in range(120):
    n = rng.randint(1, 4)
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    O = [i for i in range(n) if rng.random() < 0.5]
    r, M, h, E = brute_state(A, O)
    a = rng.randint(1, 4)
    p = rng.randint(1, a)
    lit = dp_state(close(A, O, transitive(a)), [n + p - 1])
    pred = (r + 2 ** a * M, 2 ** (a - p) * (r + h) + 2 ** p * 3 ** (a - p) * M - (E if p == 1 else 0), r + 2 ** a * M,
            ((r - h) if p == a else 0) + 2 ** (p - 1) * M)
    ok &= (lit == pred)
report("E3 H1 formula (one-vertex gadget (TT_a, {p})) on random bases, all four coordinates", ok)

print("--- 43-grid")
m, code, root, port = 9, 45174751748, 2, 8
L = adj_from_code(m, code)
indeg_port = sum(L[x][port] for x in range(m))
inn = [x for x in range(m) if L[x][port]]
ok = (indeg_port == 1 and L[root][port] == 0 and L[port][root] == 1)
print("   port in-degree in the rooted tournament: %d (in-neighbour %s); root beats port: %d; out-set of the root: %s"
      % (indeg_port, inn, L[root][port], [x for x in range(m) if L[root][x]]))
report("E4 the port (vertex 8) has in-degree 1 and is not beaten by the root (vertex 2)", ok)
ok = True
Ms = []
for a in (1, 2, 3, 4):
    verts = []
    for x in range(m):
        if x == root:
            continue
        if x == port:
            verts += [(x, k) for k in range(a)]
        else:
            verts.append((x, 0))
    n = len(verts)
    B = [[0] * n for _ in range(n)]
    for i, (x, xi) in enumerate(verts):
        for j, (y, yj) in enumerate(verts):
            if i == j:
                continue
            if x == y:
                B[i][j] = 1 if xi < yj else 0
            else:
                B[i][j] = L[x][y]
    assert is_tournament(B)
    O = [i for i, (x, _) in enumerate(verts) if L[root][x]]
    st = cstate(B, O, 1 if a <= 3 else 0)
    Ms.append(st)
    ok &= (st[0] == 1 and st[1] == 16 + 172 * 2 ** a)
    if a >= 3:   # closed count for (r, s) = (a - 2, 5): t = 1
        T = close(B, O, transitive(1))
        hT = cH(T)
        rr = a - 2
        ok &= (hT == 1 + 2 ** 5 + 43 * 2 ** (rr + 5))
        print("   a=%d: closed tournament (t=1) on %d vertices has H = %d = N(43, %d, 5)" % (a, len(T), hT, rr))
print("   states (r, M, h, E) for a = 1, 2, 3, 4:", Ms)
report("E5 carrier: offset 1 and slope 16 + 172 * 2^a for a = 1, 2, 3 (and 4)", ok)
low = {1: (7, 532488, 10, 3, 172), 2: (7, 313346, 10, 5, 172), 3: (6, 16392, 33, 9, 86), 4: (6, 17034, 9, 17, 86)}
ok = True
for s, (mm, code2, Om, r0, M0) in low.items():
    A = adj_from_code(mm, code2)
    O = [i for i in range(mm) if Om >> i & 1]
    st = brute_state(A, O)
    d, v = odd_part(st[1])
    good = (st[0] == r0 == 1 + 2 ** s and st[1] == M0 and d == 43)
    rs = (2, 3) if v > s else (1, 2)
    closed = []
    for rr in rs:
        t = rr + s - v
        hT = dp_H(close(A, O, transitive(t)))
        closed.append((rr, t, hT))
        good &= (hT == 1 + 2 ** s + 43 * 2 ** (rr + s))
    ok &= good
    print("   level %d: order %d code %d O=%s state=%s v=%d (v <= s: %s) closed (r,t,H)=%s" % (s, mm, code2, O, st, v, v <= s, closed))
report("E6 the four low sources of the 43-grid: offsets 3, 5, 9, 17, slopes 172, 172, 86, 86", ok)
h175 = brute_H(adj_from_code(7, 17976))
report("E7 H(order 7, code 17976) = 175 = N(43, 1, 1)", h175 == 175 == 1 + 2 + 43 * 4, str(h175))
# logical completeness of the 43-grid argument on a box: which mechanism covers (r, s)
miss = []
for s in range(1, 41):
    for rr in range(1, 41):
        cov = False
        if s >= 5:                               # carrier: N = 1 + 2^t (16 + 172 2^a), s = t + 4, r = a - 2, t >= 1, a >= 1
            t, a = s - 4, rr + 2
            cov = (t >= 1 and a >= 1 and 1 + 2 ** t * (16 + 172 * 2 ** a) == 1 + 2 ** s + 43 * 2 ** (rr + s))
        elif s in low:
            off, sl = low[s][3], low[s][4]
            d, v = odd_part(sl)
            t = rr + s - v
            cov = (t >= 1 and off + 2 ** t * sl == 1 + 2 ** s + 43 * 2 ** (rr + s))
        if not cov and (rr, s) == (1, 1):
            cov = True                            # the head 175
        if not cov:
            miss.append((rr, s))
report("E8 43-grid bookkeeping: carrier (s >= 5, r = a - 2 >= 1), sources (s = 2, 3, 4 all r; s = 1, r >= 2), head (1, 1): nothing left out for r, s <= 40",
       not miss, str(miss[:10]))

# ====================================================================== F
print("--- domination ray")


def dom_count(left, right, edges):
    """number of dominating sets of the bipartite graph (closed neighbourhoods), brute force over all subsets."""
    nodes = left + right
    idx = {x: k for k, x in enumerate(nodes)}
    N = len(nodes)
    nb = [1 << k for k in range(N)]
    for (a, b) in edges:
        nb[idx[a]] |= 1 << idx[b]
        nb[idx[b]] |= 1 << idx[a]
    full = (1 << N) - 1
    d = 0
    for S in range(1 << N):
        cov = 0
        T = S
        while T:
            low_ = T & -T
            cov |= nb[low_.bit_length() - 1]
            T ^= low_
        d += (cov == full)
    return d


def G_of(T, z):
    """the graph at z: left = out-neighbours of z, right = in-neighbours of z, o ~ x iff o -> x."""
    n = len(T)
    left = [x for x in range(n) if T[z][x]]
    right = [x for x in range(n) if T[x][z]]
    edges = [(o, x) for o in left for x in right if T[o][x]]
    return left, right, edges


t0 = time.time()
bad_ray = bad_slope = bad_ineq = eq_nontrans = trans_noneq = pairs = 0
bad_single = 0
for m in range(1, 5):
    for code in range(1 << (m * (m - 1) // 2)):
        A = adj_from_code(m, code)
        hB = brute_H(A)
        for Om in range(1, 1 << m):
            O = [i for i in range(m) if Om >> i & 1]
            I = [i for i in range(m) if i not in O]
            r, M, h, E = brute_state(A, O)
            delta = 0
            for k in range(1, len(O) + 1):
                for S in itertools.combinations(O, k):
                    delta += 2 ** sum(1 for x in I if any(A[s][x] for s in S))
            # G_0: no tail
            edges0 = [(o, x) for o in O for x in I if A[o][x]]
            d0 = dom_count(O, I, edges0)
            iota = 1 if all(any(A[o][x] for x in I) for o in O) else 0
            rho = d0 - iota - delta + 1
            pairs += 1
            bad_slope += (M < delta)
            if M == delta and hB > 1:
                eq_nontrans += 1
            if hB == 1 and M != delta:
                trans_noneq += 1
            if len(O) == 1:
                bad_single += (delta != 2 ** sum(A[O[0]]))
            for t in (1, 2, 3):
                T = close(A, O, transitive(t))
                z = len(T) - 1
                left, right, edges = G_of(T, z)
                assert sorted(left) == sorted(O)
                d = dom_count(left, right, edges)
                bad_ray += (d != rho + 2 ** t * delta)
                HT = dp_H(T)
                assert HT == r + 2 ** t * M
                bad_ineq += (HT - hB < d - 1)
                bad_ineq += ((r - h - rho + 1) + 2 ** t * (M - delta) < 0)
report("F1 d(G_t) = rho_dom + 2^t delta_dom (t = 1, 2, 3), rho_dom = d(G_0) - iota - delta_dom + 1; M >= delta_dom; "
       "equality iff B transitive; U_B(o) >= 2^outdeg(o); Session261 inequality on the closed tournaments; all labelled pairs of order <= 4",
       bad_ray == 0 and bad_slope == 0 and eq_nontrans == 0 and trans_noneq == 0 and bad_ineq == 0 and bad_single == 0,
       "%d pairs, %.0fs" % (pairs, time.time() - t0))
t0 = time.time()
bad_slope = eq_nontrans = trans_noneq = pairs = 0
for m in (5,):
    for code in range(1 << (m * (m - 1) // 2)):
        A = adj_from_code(m, code)
        endat, startat, tot = dp_tables(A)
        full = (1 << m) - 1
        hB = tot[full]
        Uo = [sum(tot[S] * startat[full ^ S][o] for S in range(full) if not (S >> o) & 1) for o in range(m)]
        outm = [sum(1 << j for j in range(m) if A[i][j]) for i in range(m)]
        for Om in range(1, 1 << m):
            O = [i for i in range(m) if Om >> i & 1]
            M = sum(Uo[o] for o in O)
            delta = 0
            sub = Om
            while sub:
                un = 0
                for o in O:
                    if sub >> o & 1:
                        un |= outm[o]
                delta += 2 ** bin(un & ~Om).count("1")
                sub = (sub - 1) & Om
            pairs += 1
            bad_slope += (M < delta)
            if M == delta and hB > 1:
                eq_nontrans += 1
            if hB == 1 and M != delta:
                trans_noneq += 1
report("F2 M >= delta_dom, equality iff B transitive: all labelled pairs of order 5", bad_slope == 0 and eq_nontrans == 0 and trans_noneq == 0,
       "%d pairs, %.0fs" % (pairs, time.time() - t0))
# the cited inequality itself on every (T, z), T labelled of order <= 5 (and equality iff T - z transitive for H(T) >= d(G_z))
t0 = time.time()
bad1 = bad2 = bad3 = cnt = 0
for m in range(1, 6):
    for code in range(1 << (m * (m - 1) // 2)):
        T = adj_from_code(m, code)
        HT = dp_H(T)
        for z in range(m):
            Tz = [[T[i][j] for j in range(m) if j != z] for i in range(m) if i != z]
            Hz = dp_H(Tz)
            left, right, edges = G_of(T, z)
            d = dom_count(left, right, edges)
            cnt += 1
            bad1 += (HT - Hz < d - 1)
            bad2 += (HT < d)
            bad3 += ((HT == d) != (Hz == 1))
report("F3 cited inequality H(T) - H(T - z) >= d(G_z) - 1 and H(T) >= d(G_z) with equality iff T - z transitive, all (T, z) of order <= 5",
       bad1 == 0 and bad2 == 0 and bad3 == 0, "%d rooted tournaments, %.0fs" % (cnt, time.time() - t0))

# ====================================================================== G
print("--- merge inequality / abundancy")


def sub(A, S):
    return [[A[i][j] for j in S] for i in S]


def paths_of(A, S):
    return [p for p in itertools.permutations(S) if all(A[p[k]][p[k + 1]] for k in range(len(p) - 1))]


def merge(A, P, Q):
    """the merge of two disjoint paths keeping both orders (as in merge sort)."""
    out = []
    i = j = 0
    while i < len(P) and j < len(Q):
        if A[P[i]][Q[j]]:
            out.append(P[i]); i += 1
        else:
            out.append(Q[j]); j += 1
    out += list(P[i:]) + list(Q[j:])
    return tuple(out)


t0 = time.time()
bad_merge = bad_inj = bad_lo = bad_hi = bad_eq = cnt = 0
bad_pi = 0
for n in range(1, 6):
    for code in range(1 << (n * (n - 1) // 2)):
        A = adj_from_code(n, code)
        V = list(range(n))
        Hs = {}
        for k in range(0, n + 1):
            for S in itertools.combinations(V, k):
                Hs[S] = len(paths_of(A, S))            # H of every induced sub-tournament by permutations (H(empty) = 1)
        HX = Hs[tuple(V)]
        pi2 = 0
        allpaths = set(paths_of(A, tuple(V)))
        for k in range(0, n + 1):
            for S in itertools.combinations(V, k):
                R = tuple(v for v in V if v not in S)
                pi2 += Hs[S] * Hs[R]
                bad_merge += (Hs[S] * Hs[R] > HX)
                if n <= 5:
                    imgs = set()
                    for P in paths_of(A, S):
                        for Q in paths_of(A, R):
                            mg = merge(A, P, Q)
                            if mg not in allpaths:
                                bad_inj += 1
                            if tuple(v for v in mg if v in S) != P or tuple(v for v in mg if v in R) != Q:
                                bad_inj += 1
                            imgs.add(mg)
                    bad_inj += (len(imgs) != Hs[S] * Hs[R])
        g = brute_gadget(A, [])
        bad_pi += (g["Pi"] != pi2 or g["H"] != HX)
        cnt += 1
        bad_lo += (pi2 < (n + 1) * HX)
        bad_hi += (pi2 > 2 ** n * HX)
        bad_eq += ((pi2 == 2 ** n * HX) != (HX == 1))
report("G1 H(X) >= H(X[S]) H(X - S) for every S; the merge keeps both orders, lands on a Hamiltonian path and is injective; "
       "(n + 1) H <= Pi_2 <= 2^n H with equality on the right iff transitive; all labelled tournaments of order <= 5",
       bad_merge == 0 and bad_inj == 0 and bad_lo == 0 and bad_hi == 0 and bad_eq == 0 and bad_pi == 0,
       "%d tournaments, %.0fs" % (cnt, time.time() - t0))
rng = random.Random(11)
ok = True
for _ in range(150):
    n = rng.randint(1, 4)
    k = rng.randint(1, 4)
    X = [[0] * n for _ in range(n)]
    Y = [[0] * k for _ in range(k)]
    for M_ in (X, Y):
        for i in range(len(M_)):
            for j in range(i + 1, len(M_)):
                if rng.random() < 0.5:
                    M_[i][j] = 1
                else:
                    M_[j][i] = 1
    J = [[0] * (n + k) for _ in range(n + k)]
    for i in range(n):
        for j in range(n):
            J[i][j] = X[i][j]
        for j in range(k):
            J[i][n + j] = 1
    for i in range(k):
        for j in range(k):
            J[n + i][n + j] = Y[i][j]
    ok &= (dp_pi2(J) * dp_H(X) * dp_H(Y) == dp_pi2(X) * dp_pi2(Y) * dp_H(J) and dp_H(J) == dp_H(X) * dp_H(Y)
           and dp_pi2(J) == dp_pi2(X) * dp_pi2(Y))
report("G2 a(X => Y) = a(X) a(Y) (Pi_2 and H both multiplicative), random X, Y of order <= 4", ok)
print("ALL PASS" if allok else "SOME FAILED")
