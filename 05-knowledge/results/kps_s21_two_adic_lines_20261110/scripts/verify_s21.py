"""verify_s21.py - kind-pasteur-2026-11-10-S21. Exact checks for THM-5126 (C)-(E) and (B1): the extension-mean
identity, doubling vertices, the normal forms of one-apex counts, simultaneous doubling. Pure Python, no dependencies.
usage: python verify_s21.py            (about 20 seconds)
Every check prints PASS/FAIL with its universe; the script exits 1 on any failure.
A deliberately wrong statement is tested at the end (negative control) and must FAIL.
"""
import itertools, math, random, sys
from functools import lru_cache

FAILS = []
def report(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + ("  | " + detail if detail else ""))
    if not ok: FAILS.append(name)

# ---------------------------------------------------------------- basic tournament tools
def all_tournaments(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for bits in range(1 << len(pairs)):
        A = [[0] * n for _ in range(n)]
        for k, (i, j) in enumerate(pairs):
            if bits >> k & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

def rand_tournament(n, rng):
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

def path_tables(A):
    """end[mask][v] = # Hamiltonian paths of T[mask] ending at v; start[mask][v] likewise starting at v."""
    n = len(A); full = 1 << n
    out = [sum(1 << u for u in range(n) if A[v][u]) for v in range(n)]
    inn = [sum(1 << u for u in range(n) if A[u][v]) for v in range(n)]
    end = [[0] * n for _ in range(full)]; start = [[0] * n for _ in range(full)]
    for v in range(n): end[1 << v][v] = 1; start[1 << v][v] = 1
    for m in range(1, full):
        for v in range(n):
            c = end[m][v]
            if c:
                nx = out[v] & ~m
                while nx:
                    u = (nx & -nx).bit_length() - 1; nx &= nx - 1
                    end[m | 1 << u][u] += c
            c = start[m][v]
            if c:
                nx = inn[v] & ~m
                while nx:
                    u = (nx & -nx).bit_length() - 1; nx &= nx - 1
                    start[m | 1 << u][u] += c
    return end, start

def H(A):
    n = len(A)
    if n == 0: return 1
    end, _ = path_tables(A)
    return sum(end[(1 << n) - 1])

def blowup(A, reps, kind):
    """Replace vertex v by reps[v] copies. kind[v] = 'T' (transitive module) or 'I' (pairwise non-adjacent copies)."""
    n = len(A); orig = []; idx = []
    for v in range(n):
        for i in range(reps[v]): orig.append(v); idx.append(i)
    N = len(orig); B = [[0] * N for _ in range(N)]
    for p in range(N):
        for q in range(N):
            if p == q: continue
            if orig[p] == orig[q]:
                B[p][q] = 1 if (kind[orig[p]] == 'T' and idx[p] < idx[q]) else 0
            else: B[p][q] = A[orig[p]][orig[q]]
    return B

def strong_components(A):
    n = len(A); reach = [[A[i][j] or i == j for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            if reach[i][k]:
                for j in range(n):
                    if reach[k][j]: reach[i][j] = 1
    comp = [None] * n
    for i in range(n):
        if comp[i] is None:
            for j in range(n):
                if reach[i][j] and reach[j][i]: comp[j] = i
    return comp

def on_cycle(A):
    comp = strong_components(A); n = len(A)
    size = {}
    for c in comp: size[c] = size.get(c, 0) + 1
    return [size[comp[v]] >= 3 for v in range(n)]

def stirling2(n, k):
    @lru_cache(None)
    def S(n, k):
        if n == k: return 1
        if k == 0 or k > n: return 0
        return k * S(n - 1, k) + S(n - 1, k - 1)
    return S(n, k)

def ordered_covers(A):
    """P[j] = number of ordered j-tuples of vertex-disjoint paths covering T (j = 1..n)."""
    n = len(A); full = (1 << n) - 1
    end, _ = path_tables(A)
    hsub = [sum(end[m]) for m in range(full + 1)]; hsub[0] = 0
    cur = [1 if m == 0 else 0 for m in range(full + 1)]   # covers of mask by j paths
    P = [0] * (n + 1)
    for j in range(1, n + 1):
        nxt = [0] * (full + 1)
        for m in range(1, full + 1):
            s = m; tot = 0
            while s:
                if hsub[s] and cur[m ^ s]: tot += hsub[s] * cur[m ^ s]
                s = (s - 1) & m
            nxt[m] = tot
        cur = nxt; P[j] = cur[full]
    return P

# ---------------------------------------------------------------- (A) strict submodularity: Q(a,b)+Q(b,a) >= 2, even
def Qmatrix(A):
    n = len(A); full = (1 << n) - 1
    end, start = path_tables(A)
    Q = [[0] * n for _ in range(n)]
    for X in range(1, full):
        Y = full ^ X
        for a in range(n):
            if X >> a & 1 and end[X][a]:
                for b in range(n):
                    if Y >> b & 1 and start[Y][b]: Q[a][b] += end[X][a] * start[Y][b]
    return Q, end, start

def check_A():
    ok = True; mn = None; cnt = 0; eq2 = 0
    for n in range(2, 7):
        for A in all_tournaments(n):
            Q, _, _ = Qmatrix(A)
            for a in range(n):
                for b in range(a + 1, n):
                    M = Q[a][b] + Q[b][a]; cnt += 1
                    if M < 2 or M % 2: ok = False
                    if M == 2: eq2 += 1
                    mn = M if mn is None else min(mn, M)
    rng = random.Random(1)
    for n in (7, 8, 9):
        for _ in range(40):
            Q, _, _ = Qmatrix(rand_tournament(n, rng))
            for a in range(n):
                for b in range(a + 1, n):
                    M = Q[a][b] + Q[b][a]; cnt += 1
                    if M < 2 or M % 2: ok = False
    report("(A) M_ab = Q(a,b)+Q(b,a) is even and >= 2", ok,
           "all labelled tournaments n=2..6 + 120 random n=7..9; %d pairs; minimum %d attained %d times (n<=6)" % (cnt, mn, eq2))

# ---------------------------------------------------------------- (A') cut identity and (C) mean identities
def check_C():
    ok_cut = True; ok_mean = True; ok_lb = True; ok_k2 = True; ok_pi = True; tested = 0
    for n in range(1, 6):
        for A in all_tournaments(n):
            Q, end, start = Qmatrix(A); full = (1 << n) - 1; h = sum(end[full])
            P = ordered_covers(A); tested += 1
            # P_j >= j! S(n,j), equality for all j iff transitive
            trans = (h == 1)
            eq_all = all(P[j] == math.factorial(j) * stirling2(n, j) for j in range(1, n + 1))
            if any(P[j] < math.factorial(j) * stirling2(n, j) for j in range(1, n + 1)) or (eq_all != trans): ok_lb = False
            if P[1] != h or (n >= 2 and P[2] != sum(Q[a][b] for a in range(n) for b in range(n) if a != b)): ok_cut = False
            # pi_T(2) = 2H + P_2 >= 2^n
            if 2 * h + (P[2] if n >= 2 else 0) < 2 ** n: ok_pi = False
            tot = 0
            for O in range(1 << n):   # v beats O
                B = [row[:] + [0] for row in A] + [[0] * (n + 1)]
                for u in range(n):
                    if O >> u & 1: B[n][u] = 1
                    else: B[u][n] = 1
                hv = H(B); tot += hv
                cut = sum(start[full][b] for b in range(n) if O >> b & 1) + sum(end[full][a] for a in range(n) if not O >> a & 1) \
                    + sum(Q[a][b] for a in range(n) for b in range(n) if a != b and not O >> a & 1 and O >> b & 1)
                if cut != hv: ok_cut = False
            if 4 * tot != 2 ** n * (4 * h + (P[2] if n >= 2 else 0)): ok_mean = False
    # k = 2 new vertices: sum over all 2^(2n+1) extensions = 2^(2n+1) * 2!*2^(1-2) * sum_j C(3,j) 2^-j P_j
    for n in range(1, 5):
        for A in all_tournaments(n):
            P = ordered_covers(A) + [0, 0, 0]
            tot = 0
            for O1 in range(1 << n):
                for O2 in range(1 << n):
                    for arc in (0, 1):
                        B = [row[:] + [0, 0] for row in A] + [[0] * (n + 2), [0] * (n + 2)]
                        for u in range(n):
                            if O1 >> u & 1: B[n][u] = 1
                            else: B[u][n] = 1
                            if O2 >> u & 1: B[n + 1][u] = 1
                            else: B[u][n + 1] = 1
                        if arc: B[n][n + 1] = 1
                        else: B[n + 1][n] = 1
                        tot += H(B)
            # E = 2!*2^(1-2)*(3 P1/2 + 3 P2/4 + P3/8) = (12 P1 + 6 P2 + P3)/8
            if 8 * tot != 2 ** (2 * n + 1) * (12 * P[1] + 6 * P[2] + P[3]): ok_k2 = False
    report("(A') H(T+v_O) = capacity of the directed s-t cut (THM-4114 (32)), P_2 = sum Q", ok_cut, "all labelled n=1..5, all 2^n out-sets")
    report("(C1) sum_O H(T+v_O) = 2^(n-2) (4H + P_2)", ok_mean, "all labelled n=1..5 (%d tournaments)" % tested)
    report("(C2) P_j(T) >= j! S(n,j) for all j, equality for all j iff T transitive; 2H+P_2 >= 2^n", ok_lb and ok_pi, "all labelled n=1..5")
    report("(C3) two new vertices: sum over all extensions = 2^(2n+1) (12 P_1 + 6 P_2 + P_3)/8", ok_k2, "all labelled n=1..4")

# ---------------------------------------------------------------- (B) simultaneous doubling
def check_B():
    ok_pair = True; ok_all = True; npairs = 0; nall = 0
    for n in range(3, 6):
        for A in all_tournaments(n):
            cyc = on_cycle(A)
            for x in range(n):
                for y in range(x + 1, n):
                    reps = [1] * n; reps[x] = 2; reps[y] = 2
                    c22 = H(blowup(A, reps, ['I'] * n)) // 4
                    npairs += 1
                    if (c22 > 0) != (cyc[x] and cyc[y]): ok_pair = False
            if n <= 4 and all(cyc):   # strong: every vertex twice
                nall += 1
                if H(blowup(A, [2] * n, ['I'] * n)) == 0: ok_all = False
    rng = random.Random(2)
    for _ in range(60):
        n = 6; A = rand_tournament(n, rng); cyc = on_cycle(A)
        R = [v for v in range(n) if cyc[v]]
        if len(R) >= 3:
            R = R[:3]; reps = [2 if v in R else 1 for v in range(n)]
            nall += 1
            if H(blowup(A, reps, ['I'] * n)) == 0: ok_all = False
    report("(B1) some Hamiltonian walk visits x and y twice each  iff  x and y both lie on cycles", ok_pair, "all labelled n=3..5, %d pairs" % npairs)
    report("(B2) simultaneous doubling of a set of cycle vertices", ok_all, "%d cases (all strong n<=4 with R=V; random n=6 with |R|=3)" % nall)

# ---------------------------------------------------------------- (D) normal forms of one-apex counts
@lru_cache(None)
def g(v):
    if v == 0: return 1
    if v & 1: return g(v >> 1)
    return g(v >> 1) + (v >> 1) + ((v >> 1) & 1)
@lru_cache(None)
def pi(b):      # pair penalty of the one-bits of b
    if b == 0: return 0
    if b & 1: return pi(b >> 1) + (b >> 1) + ((b >> 1) & 1)
    return pi(b >> 1)
def G(x): return g(x - 1)
def one_apex(word):   # chain 0..L-1 increasing, apex L beats chain vertex i iff word[i] == 1
    L = len(word); A = [[0] * (L + 1) for _ in range(L + 1)]
    for i in range(L):
        for j in range(i + 1, L): A[i][j] = 1
        if word[i]: A[L][i] = 1
        else: A[i][L] = 1
    return A
def check_D():
    ok_g = True; ok_pi = True; ok_traj = True; ok_sym = True
    for L in range(1, 11):
        for w in itertools.product((0, 1), repeat=L):
            v = int("".join(map(str, w)), 2)
            if L <= 9 and H(one_apex(w)) != g(v): ok_g = False
    for b in range(0, 1 << 16):
        # direct pair penalty
        ones = [i for i in range(16) if b >> i & 1]
        direct = sum((1 << (q - p - 1)) + (1 if q - p == 1 else 0) for i, p in enumerate(ones) for q in ones[i + 1:])
        if direct != pi(b) or g(b) != 1 + 2 * (b // 2) - pi(b): ok_pi = False
    for m in range(1, 1 << 15, 2):
        o = m; tot = 0
        while o != 1:
            tot += 2 * ((o + 1) // 4)
            o = o + 1
            while o % 2 == 0: o //= 2
        if G(m) - 1 != tot or G(2 * m) != G(m): ok_traj = False
    # symmetric product form for words u 1^k v: N = (1+m)(1+y) - pi(m) - pi(y), y = A(u) even, m = 2^k C(v)
    for Lu in range(0, 6):
        for u in itertools.product((0, 1), repeat=Lu):
            if Lu and u[-1] != 0: continue
            y = int("".join(map(str, u)), 2) if Lu else 0
            for k in range(1, 4):
                for Lv in range(1, 5):
                    for vv in itertools.product((0, 1), repeat=Lv):
                        C = sum((1 - vv[j]) << j for j in range(Lv)); m = (1 << k) * C
                        w = u + (1,) * k + vv; val = g(int("".join(map(str, w)), 2))
                        if val != (1 + m) * (1 + y) - pi(m) - pi(y): ok_sym = False
    report("(D1) one-apex count H(w) = g(A(w)), g(2v+1)=g(v), g(2v)=g(v)+v+(v mod 2)", ok_g, "brute force, all words of length <= 9")
    report("(D2) g(b) = 1 + 2 floor(b/2) - pi(b), pi = pair penalty of the one-bits", ok_pi, "all b < 2^16")
    report("(D3) G(m)-1 = sum over the (x+1)-trajectory of m of 2 floor((o+1)/4); G(2m)=G(m)", ok_traj, "all odd m < 2^15")
    report("(D4) H(u 1^k v) = (1+m)(1+y) - pi(m) - pi(y)", ok_sym, "u of length <= 5 ending in 0 (or empty), k <= 3, v of length <= 4")

# ---------------------------------------------------------------- (E) doubling vertices
def check_E():
    ok_tower = True; ok_sub = True; nd = 0
    C3 = [[0, 1, 0], [0, 0, 1], [1, 0, 0]]
    T4 = [[0, 1, 1, 0], [0, 0, 1, 1], [0, 0, 0, 1], [1, 0, 0, 0]]
    mods = [C3, T4, [[0, 1], [0, 0]], [[0]]]
    def pi2(M): P = ordered_covers(M) + [0, 0]; return 2 * P[1] + (P[2] if len(M) >= 2 else 0)
    def substitute(A, x, M):
        n = len(A); k = len(M); verts = [(v, 0) for v in range(n) if v != x] + [(x, i) for i in range(k)]
        N = len(verts); B = [[0] * N for _ in range(N)]
        for p in range(N):
            for q in range(N):
                if p == q: continue
                (v, i), (u, j) = verts[p], verts[q]
                B[p][q] = M[i][j] if (v == x and u == x) else A[v][u]
        return B
    for n in range(3, 6):
        for A in all_tournaments(n):
            h = H(A)
            for x in range(n):
                reps = [1] * n; kind = ['T'] * n
                vals = []
                for a in (2, 3, 4):
                    reps[x] = a; vals.append(H(blowup(A, reps, kind)))
                if vals == [2 * h - 1, 4 * h - 3, 8 * h - 7]:
                    nd += 1
                    reps[x] = 5
                    if H(blowup(A, reps, kind)) != 16 * h - 15: ok_tower = False
                    for M in mods:
                        if len(M) + n - 1 <= 8 and H(substitute(A, x, M)) != H(M) + (h - 1) * pi2(M) // 2: ok_sub = False
    report("(E1) doubling vertex: H(K[TT_a at x]) = 1 + (H-1) 2^(a-1) persists at a = 5", ok_tower, "%d doubling vertices, all labelled n=3..5" % nd)
    report("(E2) doubling vertex: H(K[M at x]) = H(M) + (H(K)-1) pi_M(2)/2 for M in {C3, T4, TT2, point}", ok_sub)

# ---------------------------------------------------------------- negative control
def negative_control():
    # false claim: every vertex on a cycle is a doubling vertex (fails at the strong 4-tournament)
    T4 = [[0, 1, 1, 0], [0, 0, 1, 1], [0, 0, 0, 1], [1, 0, 0, 0]]
    h = H(T4); bad = False
    for x in range(4):
        reps = [1] * 4; reps[x] = 2
        if H(blowup(T4, reps, ['T'] * 4)) != 2 * h - 1: bad = True
    print(("PASS " if bad else "FAIL ") + "negative control: 'every cycle vertex is doubling' is correctly rejected")
    if not bad: FAILS.append("negative control")

if __name__ == "__main__":
    check_A(); check_C(); check_B(); check_D(); check_E(); negative_control()
    print("ALL PASS" if not FAILS else "FAILURES: %s" % FAILS)
    sys.exit(1 if FAILS else 0)
