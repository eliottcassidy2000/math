"""sources_check.py - checks for SOURCES_AND_DOMINATION.md (kind-pasteur-2026-11-10-S22, continuation of 2026-11-11).
Pure Python, no input files; every count is a literal subset dynamic programme (towers.py).

S1  one-vertex gadget on a seed of slope 2^k: D = 1 + 2^(a+k), U = 2^(a-p)(1+h) + 2^(p+k) 3^(a-p) - [p=1] E
S2  the sources (5, 9), (17, 18), (33, 36) of the 9-grid (Session261) are instances of S1
S3  new sources of the 9-grid at levels s = 7, 8, 9, 10, 11: states and the first two grid points
S4  a carrier for the 43-grid: slope family 16 + 172 * 2^a at one port of an order-9 unit root
S5  domination ray: d(G_t) = rho + 2^t delta_dom, delta_dom = sum over nonempty S in O of 2^|N+(S) - O|,
    M >= delta_dom with equality iff B transitive (all labelled tournaments of order <= 5, every O)
S6  merge inequality H(X) >= H(X[S]) H(X - S); abundancy a(X) = Pi_2(X)/H(X) in [n + 1, 2^n]; a(X => Y) = a(X) a(Y)
S7  transitive interfaces: s(TT_n, O) = (1 - w(O), sum 2^(n-o), 1, [n in O]) with w = sum over pairs of max(2, 2^(dist-1))"""
import itertools, random
from fractions import Fraction
from towers import tables, H, state, gadget, apply, level_up, TT, rand_t
rng = random.Random(20261111)
allok = True


def report(name, ok, extra=""):
    global allok
    allok &= bool(ok)
    print(("PASS  " if ok else "FAIL  ") + name + ("  " + extra if extra else ""))


def adj(m, code):
    A = [[0] * m for _ in range(m)]
    b = 0
    for i in range(m):
        for j in range(i + 1, m):
            if code >> b & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            b += 1
    return A


def odd(x):
    v = 0
    while x % 2 == 0:
        x //= 2
        v += 1
    return x, v


def bits(mask, m):
    return [i for i in range(m) if mask >> i & 1]


# S1: one-vertex gadget formula against literal states, on seeds with a power of two as slope
seeds = [((1, 0, 1), "transitive point"), ((2, 0, 1), "TT_2, source"), ((3, 0, 1), "TT_3, source"),
         ((4, 34, 12), "h = 5, M = 16"), ((5, 642, 10), "h = 9, M = 32"), ((6, 17410, 18), "h = 15, M = 64")]
bad = 0
cnt = 0
for (m, code, Om), _ in seeds:
    A = adj(m, code)
    O = bits(Om, m)
    r, M, h, E = state(A, O)
    assert r == 1 and M & (M - 1) == 0
    k = M.bit_length() - 1
    for a in range(1, 4):
        for p in range(1, a + 1):
            D = 1 + 2 ** (a + k)
            U = 2 ** (a - p) * (1 + h) + 2 ** (p + k) * 3 ** (a - p) - (E if p == 1 else 0)
            st = state(level_up(A, O, TT(a)), [m + p - 1])
            cnt += 1
            bad += not (st[0] == D and st[1] == U and st[2] == D)
report("S1 one-vertex gadget on a seed of slope 2^k: D = 1 + 2^(a+k), U = 2^(a-p)(1+h) + 2^(p+k) 3^(a-p) - [p=1]E", bad == 0,
       "%d literal states" % cnt)

# S2: three sources of Session261 as gadget instances
ok = True
for (m, code, Om), (a, p), (D, U) in [((1, 0, 1), (2, 1), (5, 9)), ((1, 0, 1), (4, 4), (17, 18)), ((4, 34, 12), (1, 1), (33, 36))]:
    A = adj(m, code)
    st = state(level_up(A, bits(Om, m), TT(a)), [m + p - 1])
    ok &= (st[0] == D and st[1] == U)
report("S2 (5, 9), (17, 18), (33, 36) are one-vertex gadgets on seeds of slope a power of two", ok)

# S3: new 9-grid sources, levels 7..11
src = {7: (8, 11026436, 4), 8: (8, 6562209, 4), 9: (9, 442533377, 1), 10: (9, 6793626886, 6), 11: (9, 194216770, 1)}
ok = True
lines = []
for s, (m, code, Om) in src.items():
    A = adj(m, code)
    O = bits(Om, m)
    r, M, h, E = state(A, O)
    d, v = odd(M)
    good = (r == 1 + 2 ** s and d == 9 and v <= s)
    for rr in (1, 2):
        t = rr + s - v
        good &= (H(level_up(A, O, TT(t))) == 1 + 2 ** s + 9 * 2 ** (rr + s))
    ok &= good
    lines.append("s=%d: order %d, code %d, O=%s, (offset, slope, H(B)) = (%d, %d, %d)" % (s, m, code, O, r, M, h))
report("S3 9-grid sources at levels 7..11 (offset 1 + 2^s, slope 9 * 2^v, v <= s); N(9, r, s) recounted for r = 1, 2", ok)
for l in lines:
    print("        " + l)

# S4: carrier for the 43-grid: order-9 unit root, code 45174751748, root 2, port 8; module TT_a at the port
m, code, root, port = 9, 45174751748, 2, 8
L = adj(m, code)
ok = True
Ms = []
for a in (1, 2, 3):
    orig = []
    for x in range(m):
        orig += [(x, i) for i in range(a if x == port else 1)]
    n = len(orig)
    T = [[0] * n for _ in range(n)]
    for i, (x, xi) in enumerate(orig):
        for j, (y, yj) in enumerate(orig):
            if i != j:
                T[i][j] = (1 if xi < yj else 0) if x == y else L[x][y]
    rv = [i for i, (x, _) in enumerate(orig) if x == root][0]
    B = [[T[i][j] for j in range(n) if j != rv] for i in range(n) if i != rv]
    keep = [i for i in range(n) if i != rv]
    O = [k for k, i in enumerate(keep) if T[rv][i]]
    r, M, h, E = state(B, O)
    ok &= (r == 1 and M == 16 + 172 * 2 ** a)
    Ms.append(M)
report("S4 carrier for d = 43: unit offset for module sizes a = 1, 2, 3 at the port, slope 16 + 172 * 2^a", ok, str(Ms))

# S5: domination ray on all labelled tournaments of order <= 5
bad_slope = 0
eq_nontrans = 0
eq_trans_fail = 0
pairs = 0
bad_ray = 0
for m in range(1, 6):
    for code in range(1 << (m * (m - 1) // 2)):
        A = adj(m, code)
        hB = H(A)
        for Om in range(1, 1 << m):
            O = bits(Om, m)
            I = [i for i in range(m) if i not in O]
            r, M, h, E = state(A, O)
            dl = 0
            for k in range(1, len(O) + 1):
                for S in itertools.combinations(O, k):
                    dl += 2 ** sum(1 for i in I if any(A[s][i] for s in S))
            pairs += 1
            bad_slope += (M < dl)
            if M == dl and hB > 1:
                eq_nontrans += 1
            if hB == 1 and M != dl:
                eq_trans_fail += 1
            if m <= 3:
                # literal dominating sets of G_t, t = 1, 2: left O, right I + tail (tail adjacent to all of O)
                for t in (1, 2):
                    nodes = len(O) + len(I) + t
                    nb = [1 << v for v in range(nodes)]
                    for a_, o in enumerate(O):
                        for b_, i in enumerate(I):
                            if A[o][i]:
                                nb[a_] |= 1 << (len(O) + b_)
                                nb[len(O) + b_] |= 1 << a_
                        for z in range(t):
                            nb[a_] |= 1 << (len(O) + len(I) + z)
                            nb[len(O) + len(I) + z] |= 1 << a_
                    d = 0
                    for Sm in range(1 << nodes):
                        cov = 0
                        for v in range(nodes):
                            if Sm >> v & 1:
                                cov |= nb[v]
                        d += (cov == (1 << nodes) - 1)
                    if t == 1:
                        d1 = d
                    else:
                        bad_ray += (d - d1 != 2 * dl)
                    bad_ray += ((r + 2 ** t * M) - d < hB - 1)
report("S5 M >= delta_dom for every pair; equality iff B transitive; d(G_2) - d(G_1) = 2 delta_dom; H(T_t) - d(G_t) >= H(B) - 1",
       bad_slope == 0 and eq_nontrans == 0 and eq_trans_fail == 0 and bad_ray == 0, "%d pairs" % pairs)

# S6: merge inequality and abundancy
bad = 0
lo = hi = True
for m in range(1, 6):
    for code in range(1 << (m * (m - 1) // 2)):
        A = adj(m, code)
        full = (1 << m) - 1
        en, st, hs = tables(A)
        pi2 = sum(hs[X] * hs[full ^ X] for X in range(full + 1))
        bad += any(hs[X] * hs[full ^ X] > hs[full] for X in range(full + 1))
        a = Fraction(pi2, hs[full])
        lo &= (a >= m + 1)
        hi &= (a <= 2 ** m and ((a == 2 ** m) == (hs[full] == 1)))
mult = True
for _ in range(100):
    X = rand_t(rng.randint(1, 4), rng)
    Y = rand_t(rng.randint(1, 4), rng)
    n, k = len(X), len(Y)
    J = [[0] * (n + k) for _ in range(n + k)]
    for i in range(n):
        for j in range(n):
            J[i][j] = X[i][j]
        for j in range(k):
            J[i][n + j] = 1
    for i in range(k):
        for j in range(k):
            J[n + i][n + j] = Y[i][j]

    def ab(T):
        full = (1 << len(T)) - 1
        en, st, hs = tables(T)
        return Fraction(sum(hs[S] * hs[full ^ S] for S in range(full + 1)), hs[full])
    mult &= (ab(J) == ab(X) * ab(Y))
report("S6 H(X) >= H(X[S]) H(X - S); n + 1 <= Pi_2(X)/H(X) <= 2^n with equality on the right iff X transitive; multiplicative",
       bad == 0 and lo and hi and mult)

# S7: transitive interfaces
ok = True
for n in range(1, 8):
    for Om in range(1, 1 << n):
        O = bits(Om, n)
        w = sum(max(2, 2 ** (b - a - 1)) for a, b in itertools.combinations(O, 2))
        ok &= (state(TT(n), O) == (1 - w, sum(2 ** (n - 1 - o) for o in O), 1, 1 if (n - 1) in O else 0))
report("S7 s(TT_n, O) = (1 - w(O), sum_{o in O} 2^(n-o), 1, [n in O]), w = sum over pairs of max(2, 2^(distance-1)), n <= 7", ok)
# S8: the 43-grid N(43, r, s) = 1 + 2^s + 43 * 2^(r+s), all r, s >= 1
#   carrier (S4): the port has in-degree 1 in the rooted tournament and is not selected, so H, W, U are affine in 2^a
#   (the module is visited at most twice); unit offset at a = 1, 2 then holds for every a: s >= 5, r >= 1.
ok = (sum(L[x][port] for x in range(len(L))) == 1 and not L[root][port])
low = {1: (7, 532488, 10, 3, 172), 2: (7, 313346, 10, 5, 172), 3: (6, 16392, 33, 9, 86), 4: (6, 17034, 9, 17, 86)}
for s, (mm, code, Om, r0, M0) in low.items():
    A = adj(mm, code)
    O = bits(Om, mm)
    r, M, h, E = state(A, O)
    ok &= (r == r0 == 1 + 2 ** s and M == M0)
    d, v = odd(M)
    for rr in ((2, 3) if s == 1 else (1, 2)):          # level 1 starts at r = 2 (v2 = 2); the head N(43,1,1) = 175 is below
        ok &= (H(level_up(A, O, TT(rr + s - v))) == 1 + 2 ** s + 43 * 2 ** (rr + s))
ok &= (H(adj(7, 17976)) == 175 == 1 + 2 + 43 * 4)
report("S8 43-grid: carrier (s >= 5), sources at s = 1 (r >= 2), 2, 3, 4, and the head N(43, 1, 1) = 175 (order 7)", ok)
print("ALL PASS" if allok else "SOME FAILED")
