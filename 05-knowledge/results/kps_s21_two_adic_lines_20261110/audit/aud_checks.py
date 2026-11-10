"""aud_checks.py -- independent audit checks for THM-5126 (auditor's own code; imports nothing from the packet).
Everything here is brute force from the DEFINITIONS (permutations and split points), plus two independent
Hamiltonian-path counters (subset DP and depth-first enumeration of the paths themselves).

  T1  tail law (A) on random rooted tournaments and words, state (H,d,E,U) from the definitions
  T2  offset formula c0 = d - delta', M = H - E + delta' and its cut form (1'):
      delta' = sum over pairs of O of W_ef, M = sum over O of U_B(y); c0 odd; c0 <= H(B), equality iff |O| <= 1;
      W_ef even and >= 2 H(B - e - f)
  T3  seeds of line_seeds.txt: c = 3, 5, 7, 13, ten random other lines, the four module seeds and ALL eleven
      order-ten seeds: (H,d,E,U), c0, M, and literal counts of the first ray points (DP and DFS)
  T4  extension-mean law (C) for k = 1, 2, 3 with P_j from the definition; equality case of P_j >= j! S(n,j)
  T5  trajectory formula (D) against literal one-apex counts
  T6  remark (E): the Camion construction is carried out literally and the resulting walk is checked
  T7  README section 6.2: mean of the path-cover polynomial = sum over odd-order permutations (closed forms, n <= 12;
      brute force over all labelled tournaments, n <= 4)
usage: python aud_checks.py [path to line_seeds.txt]
"""
import itertools, math, os, random, sys
from fractions import Fraction
sys.setrecursionlimit(10000)
FAILS = []
def report(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + ("  | " + detail if detail else "")); sys.stdout.flush()
    if not ok: FAILS.append(name)

# ------------------------------------------------------------------ tournaments
def dec(code, m):
    A = [[0] * m for _ in range(m)]; b = 0
    for i in range(m):
        for j in range(i + 1, m):
            if code >> b & 1: A[i][j] = 1
            else: A[j][i] = 1
            b += 1
    return A
def rand_t(n, rng):
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A
def all_t(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for bits in range(1 << len(pairs)):
        A = [[0] * n for _ in range(n)]
        for k, (i, j) in enumerate(pairs):
            if bits >> k & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A
def sub(A, keep): return [[A[i][j] for j in keep] for i in keep]
def is_path(A, seq): return all(A[seq[i]][seq[i + 1]] for i in range(len(seq) - 1))

def H_dp(A):
    n = len(A)
    if n == 0: return 1
    out = [sum(1 << u for u in range(n) if A[v][u]) for v in range(n)]
    f = [[0] * n for _ in range(1 << n)]
    for v in range(n): f[1 << v][v] = 1
    for s in range(1, 1 << n):
        fs = f[s]
        for v in range(n):
            c = fs[v]
            if c:
                nx = out[v] & ~s
                while nx:
                    lb = nx & -nx; nx ^= lb
                    f[s | lb][lb.bit_length() - 1] += c
    return sum(f[(1 << n) - 1])
def H_dfs(A):
    """literal enumeration of the Hamiltonian paths (no dynamic programming)"""
    n = len(A)
    if n == 0: return 1
    outl = [[u for u in range(n) if A[v][u]] for v in range(n)]; full = (1 << n) - 1; cnt = 0
    stack = [(v, 1 << v) for v in range(n)]
    while stack:
        v, mask = stack.pop()
        if mask == full: cnt += 1; continue
        for u in outl[v]:
            if not mask >> u & 1: stack.append((u, mask | 1 << u))
    return cnt
def H_perm(A): return sum(1 for p in itertools.permutations(range(len(A))) if is_path(A, p))

def state_bf(A, v):
    """(H, d, E, U, delta') straight from the definitions of THM-5126 (Setting)"""
    n = len(A); others = [u for u in range(n) if u != v]
    Hh = H_perm(A); d = 0; E = 0; U = 0; dp = 0
    for p in itertools.permutations(others):
        if is_path(A, p): d += 1
        if is_path(A, p + (v,)): E += 1
        for k in range(len(others) + 1):
            P1 = p[:k]; P2 = (v,) + p[k:]
            if is_path(A, P1) and is_path(A, P2):
                U += 1
                if 1 <= k < len(others) and A[v][P1[-1]]: dp += 1      # crossed cover
    return Hh, d, E, U, dp
def cut_bf(B):
    """Q_B(e,f) and U_B(y) from the definitions"""
    m = len(B); Q = [[0] * m for _ in range(m)]; UB = [0] * m
    for p in itertools.permutations(range(m)):
        for k in range(0, m):
            P1 = p[:k]; P2 = p[k:]
            if is_path(B, P1) and is_path(B, P2):
                UB[P2[0]] += 1
                if k >= 1: Q[P1[-1]][P2[0]] += 1
    return Q, UB
def state_dp(A, v):
    """(H, d, E, U) for larger orders: U by summing over the vertex set X of P2 (own DP)"""
    n = len(A); full = (1 << n) - 1
    out = [sum(1 << u for u in range(n) if A[x][u]) for x in range(n)]
    endc = [[0] * n for _ in range(1 << n)]
    for x in range(n): endc[1 << x][x] = 1
    for s in range(1, 1 << n):
        row = endc[s]
        for x in range(n):
            c = row[x]
            if c:
                nx = out[x] & ~s
                while nx:
                    lb = nx & -nx; nx ^= lb
                    endc[s | lb][lb.bit_length() - 1] += c
    hsub = [sum(r) for r in endc]; hsub[0] = 1
    # paths starting at v: DP forward from v only
    fromv = [[0] * n for _ in range(1 << n)]; fromv[1 << v][v] = 1
    for s in range(1, 1 << n):
        if not s >> v & 1: continue
        row = fromv[s]
        for x in range(n):
            c = row[x]
            if c:
                nx = out[x] & ~s
                while nx:
                    lb = nx & -nx; nx ^= lb
                    fromv[s | lb][lb.bit_length() - 1] += c
    U = sum(sum(fromv[X]) * hsub[full ^ X] for X in range(1, 1 << n) if X >> v & 1)
    return hsub[full], hsub[full ^ (1 << v)], endc[full][v], U
def tail(A, v, y):
    n = len(A); t = len(y); N = n + t; B = [[0] * N for _ in range(N)]
    for i in range(n):
        for j in range(n): B[i][j] = A[i][j]
    for a in range(t):
        z = n + a
        for b in range(a + 1, t): B[z][n + b] = 1
        for u in range(n):
            if u != v: B[u][z] = 1
        if y[a]: B[v][z] = 1
        else: B[z][v] = 1
    return B
def blow(A, reps):
    """replace vertex v by a transitive module of reps[v] vertices"""
    n = len(A); orig = []; idx = []
    for v in range(n):
        for i in range(reps[v]): orig.append(v); idx.append(i)
    N = len(orig); B = [[0] * N for _ in range(N)]
    for p in range(N):
        for q in range(N):
            if p != q: B[p][q] = (1 if idx[p] < idx[q] else 0) if orig[p] == orig[q] else A[orig[p]][orig[q]]
    return B

# ------------------------------------------------------------------ T1 tail law
def T1():
    rng = random.Random(20261110); ok = True; ok_dfs = True; cnt = 0; cnt_single = 0
    pool = [[[0]]] + [A for n in (2, 3) for A in all_t(n)] + [rand_t(n, rng) for n in (4, 5, 6) for _ in range(25)]
    for A in pool:
        n = len(A)
        for v in range(n):
            Hh, d, E, U, dp = state_bf(A, v)
            words = [y for t in (1, 2, 3) for y in itertools.product((0, 1), repeat=t)]
            words += [tuple(rng.randint(0, 1) for _ in range(rng.randint(4, 6))) for _ in range(6)]
            for y in words:
                C = sum((1 - y[j]) << j for j in range(len(y))); eps = 1 - y[0]
                gy = H_dp(tail([[0]], 0, y))                       # g(y): literal count of {v} o y
                lhs = H_dp(tail(A, v, y)); cnt += 1; cnt_single += (n == 1)
                if lhs != Hh + (U - d) * C + d * (gy - 1) + (d - E) * eps: ok = False
                if n + len(y) <= 8 and H_dfs(tail(A, v, y)) != lhs: ok_dfs = False
            for t in range(1, 6):                                   # rays (1)
                if H_dp(tail(A, v, (0,) * t)) != (Hh - E + 2 * d - U) + (U - d) * 2 ** t: ok = False
    report("T1 tail law H(L o y) = H + (U-d)C(y) + d(g(y)-1) + (d-E)eps(y); rays c0 + M 2^t, t = 1..5", ok,
           "%d (L,v,y) triples (%d with L a single vertex); state from the definitions; all labelled n<=3, 75 random n=4..6" % (cnt, cnt_single))
    report("T1b DP counter agrees with literal DFS enumeration on every extension of order <= 8", ok_dfs)

# ------------------------------------------------------------------ T2 offset formula and cut form
def T2():
    rng = random.Random(7); ok_a = True; ok_cut = True; ok_par = True; ok_eq = True; ok_W = True; cnt = 0
    pool = [A for n in (1, 2, 3, 4) for A in all_t(n)] + [rand_t(n, rng) for n in (5, 6) for _ in range(40)] + [rand_t(7, rng) for _ in range(12)]
    for A in pool:
        n = len(A)
        for v in range(n):
            Hh, d, E, U, dp = state_bf(A, v); c0 = Hh - E + 2 * d - U; M = U - d; cnt += 1
            if c0 != d - dp or M != Hh - E + dp: ok_a = False
            keep = [u for u in range(n) if u != v]; B = sub(A, keep); O = [i for i, u in enumerate(keep) if A[v][u]]
            if n >= 2:
                Q, UB = cut_bf(B); m = n - 1
                W = sum(Q[e][f] + Q[f][e] for i, e in enumerate(O) for f in O[i + 1:])
                if dp != W or M != sum(UB[y] for y in O) or c0 != d - W: ok_cut = False
                for e in range(m):
                    for f in range(e + 1, m):
                        w = Q[e][f] + Q[f][e]
                        if w % 2 or w < 2 * H_dp(sub(B, [x for x in range(m) if x not in (e, f)])): ok_W = False
            else:
                if dp != 0 or M != 0: ok_cut = False
            if c0 % 2 == 0: ok_par = False
            if c0 > d or ((c0 == d) != (len(O) <= 1)): ok_eq = False
    report("T2a c0 = d - delta', M = H - E + delta' (delta' = crossed covers counted from the definition)", ok_a, "%d rooted tournaments: all labelled n<=4, 80 random n=5,6, 12 random n=7, every root" % cnt)
    report("T2b cut form (1'): delta' = sum_{pairs of O} W_ef, M = sum_{y in O} U_B(y), c0 = H(B) - W(O)", ok_cut)
    report("T2c c0 is odd", ok_par)
    report("T2d c0 <= H(B), with equality iff |O| <= 1", ok_eq)
    report("T2e W_ef is even and W_ef >= 2 H(B - e - f) > 0 for every pair", ok_W)

# ------------------------------------------------------------------ T3 seeds
def T3(fn):
    rows = []
    for line in open(fn):
        if line.startswith("#"): continue
        p = line.split(); rows.append((int(p[0]), int(p[1]), p[2], int(p[3]), int(p[4]), int(p[5])) + tuple(int(x) for x in p[6:10]))
    byc = {r[0]: r for r in rows}
    rng = random.Random(5126)
    pick = [3, 5, 7, 13] + sorted(rng.sample([r[0] for r in rows if r[0] not in (3, 5, 7, 13) and r[3] <= 9 and r[2] == "tail"], 10))
    mods = [r[0] for r in rows if r[2] == "module"]
    ten = [19, 33, 71, 79, 83, 97, 101, 103, 107, 113, 115]
    ok = True; lit = 0; lines = []
    for c in pick + mods + ten:
        c_, s0, mech, m, code, v, Hh, d, E, U = byc[c]; A = dec(code, m)
        if mech == "tail":
            st = state_dp(A, v)
            if m <= 7:
                bf = state_bf(A, v)
                if bf[:4] != st or bf[4] != st[1] - 1: ok = False
            good = (st == (Hh, d, E, U)) and st[0] - st[2] + 2 * st[1] - st[3] == 1 and st[3] - st[1] == c * 2 ** (s0 - 1)
            vals = []
            for t in ((1, 2) if m == 10 else (1, 2, 3)):
                T = tail(A, v, (0,) * t); h1 = H_dp(T); h2 = H_dfs(T) if m + t <= 12 else h1; lit += 1
                vals.append(h1)
                if h1 != 1 + c * 2 ** (s0 - 1 + t) or h2 != h1: good = False
            if m == 10 and sum(A[v]) not in (2, 3, 4): good = False
            lines.append("  c=%d s0=%d tail order %d root %d outdeg %d (H,d,E,U)=%s c0=%d M=%d ray points %s %s" % (c, s0, m, v, sum(A[v]), st, st[0] - st[2] + 2 * st[1] - st[3], st[3] - st[1], vals, "ok" if good else "BAD"))
        else:
            h = H_dp(A); good = (h == Hh) and (h - 1 == c * 2 ** s0); vals = []
            for a in (2, 3, 4):
                reps = [1] * m; reps[v] = a; T = blow(A, reps); h1 = H_dp(T); h2 = H_dfs(T) if m + a - 1 <= 11 else h1; lit += 1; vals.append(h1)
                if h1 != 1 + (h - 1) * 2 ** (a - 1) or h1 != h2: good = False
            lines.append("  c=%d s0=%d module order %d x=%d H=%d K[TT_a at x], a=2,3,4: %s %s" % (c, s0, m, v, h, vals, "ok" if good else "BAD"))
        ok &= good
    for l in lines: print(l)
    report("T3 seeds c = 3, 5, 7, 13, ten random lines %s, the module seeds %s and the eleven order-ten seeds: states, c0 = 1, M = c 2^(s0-1), literal ray points (DP and DFS)" % (pick[4:], mods), ok, "%d literal counts" % lit)

# ------------------------------------------------------------------ T4 extension-mean law
def Pj_bf(A):
    n = len(A); P = [0] * (n + 2)
    for p in itertools.permutations(range(n)):
        b = sum(1 for i in range(n - 1) if not A[p[i]][p[i + 1]])
        for j in range(b + 1, n + 1): P[j] += math.comb(n - 1 - b, j - 1 - b)
    return P
def stirling2(n, k):
    if n == k: return 1
    if k == 0 or k > n: return 0
    return k * stirling2(n - 1, k) + stirling2(n - 1, k - 1)
def ext_sum(A, k):
    n = len(A); N = n + k; pairs = [(n + a, u) for a in range(k) for u in range(n)] + [(n + a, n + b) for a in range(k) for b in range(a + 1, k)]
    tot = 0
    for bits in range(1 << len(pairs)):
        B = [[0] * N for _ in range(N)]
        for i in range(n):
            for j in range(n): B[i][j] = A[i][j]
        for q, (a, b) in enumerate(pairs):
            if bits >> q & 1: B[a][b] = 1
            else: B[b][a] = 1
        tot += H_dp(B)
    return tot, len(pairs)
def T4():
    rng = random.Random(3); ok = True; ok_b = True; cnt = 0
    uni = [(1, A) for n in (1, 2, 3, 4) for A in all_t(n)] + [(1, rand_t(n, rng)) for n in (5, 6) for _ in range(8)]
    uni += [(2, A) for n in (1, 2, 3) for A in all_t(n)] + [(2, rand_t(4, rng)) for _ in range(6)]
    uni += [(3, A) for n in (1, 2) for A in all_t(n)] + [(3, rand_t(3, rng)) for _ in range(3)]
    for k, A in uni:
        n = len(A); P = Pj_bf(A) + [0] * 8; tot, nb = ext_sum(A, k); cnt += 1
        rhs = Fraction(math.factorial(k), 2 ** (k - 1)) * sum(Fraction(math.comb(k + 1, j), 2 ** j) * P[j] for j in range(1, k + 2))
        if Fraction(tot, 2 ** nb) != rhs: ok = False
        if P[1] != H_dp(A): ok = False
        alleq = all(P[j] == math.factorial(j) * stirling2(n, j) for j in range(1, n + 1))
        if any(P[j] < math.factorial(j) * stirling2(n, j) for j in range(1, n + 1)) or alleq != (P[1] == 1): ok_b = False
        if k == 1 and tot * 4 != 2 ** n * (4 * P[1] + P[2]): ok = False
        if k == 1 and tot < 2 ** (2 * n - 2) + 2 ** (n - 1): ok_b = False
    report("T4 (C) mean of H over all k-vertex extensions = k! 2^(1-k) sum_j C(k+1,j) 2^(-j) P_j, P_j from the definition", ok, "%d (k, T) cases: k=1 all n<=4 + 16 random n=5,6; k=2 all n<=3 + 6 random n=4; k=3 all n<=2 + 3 random n=3" % cnt)
    report("T4b P_j >= j! S(n,j), equality for all j iff transitive; k = 1 sum >= 2^(2n-2) + 2^(n-1)", ok_b)

# ------------------------------------------------------------------ T5 trajectory formula
def T5():
    def g(v):
        if v == 0: return 1
        return g(v >> 1) if v & 1 else g(v >> 1) + (v >> 1) + ((v >> 1) & 1)
    ok = True; cnt = 0
    for L in range(1, 9):
        for w in itertools.product((0, 1), repeat=L):
            val = int("".join(map(str, w)), 2); lit = H_dfs(tail([[0]], 0, w)); cnt += 1      # {v} o w: v -> z_j iff w_j = 1
            if lit != g(val): ok = False
            m = val + 1
            if m % 2:
                o = m; tot = 1
                while o != 1:
                    tot += 2 * ((o + 1) // 4); o += 1
                    while o % 2 == 0: o //= 2
                if tot != lit: ok = False
    report("T5 (D) literal one-apex count of {v} o w = g(value of w, first letter most significant) = trajectory formula at m = value + 1 odd", ok, "%d words of length <= 8, DFS counts" % cnt)

# ------------------------------------------------------------------ T6 remark (E): the Camion construction, literally
def components(A):
    n = len(A); reach = [[A[i][j] or i == j for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            if reach[i][k]:
                for j in range(n):
                    if reach[k][j]: reach[i][j] = 1
    comps = []; seen = set()
    for i in range(n):
        if i in seen: continue
        C = [j for j in range(n) if reach[i][j] and reach[j][i]]; seen.update(C); comps.append(C)
    comps.sort(key=lambda C: -sum(1 for j in range(n) if reach[C[0]][j]))     # sources first
    return comps
def ham_path(A, C): return next(p for p in itertools.permutations(C) if is_path(A, p))
def ham_cycle_from(A, C, x):
    rest = [u for u in C if u != x]
    return next((x,) + p for p in itertools.permutations(rest) if is_path(A, (x,) + p) and A[p[-1]][x])
def T6():
    rng = random.Random(11); ok = True; ok_iff = True; cnt = 0; cnt_no = 0
    pool = [A for n in (3, 4, 5) for A in all_t(n)] + [rand_t(n, rng) for n in (6, 7) for _ in range(40)]
    for A in pool:
        n = len(A); comps = components(A); where = {u: i for i, C in enumerate(comps) for u in C}
        for x in range(n):
            for y in range(n):
                if x == y: continue
                cyc = len(comps[where[x]]) >= 3 and len(comps[where[y]]) >= 3
                if n <= 5:      # literal count of walks through x twice and y twice: Hamiltonian paths of K[I_2 at x, I_2 at y] / 4
                    orig = [u for u in range(n)] + [x, y]; N = n + 2; B = [[0] * N for _ in range(N)]
                    for p in range(N):
                        for q in range(N):
                            if p != q and orig[p] != orig[q]: B[p][q] = A[orig[p]][orig[q]]
                    c22 = H_dfs(B) // 4
                    if (c22 > 0) != cyc: ok_iff = False
                if not cyc: cnt_no += 1; continue
                walk = []
                for i, C in enumerate(comps):
                    if i == where[x] and i == where[y]:
                        w = list(ham_cycle_from(A, C, x)) + [x]
                        w = [y] + w if A[y][x] else w + [y]
                    elif i == where[x]: w = list(ham_cycle_from(A, C, x)) + [x]
                    elif i == where[y]: w = list(ham_cycle_from(A, C, y)) + [y]
                    else: w = list(ham_path(A, C))
                    walk += w
                cnt += 1
                good = all(A[walk[i]][walk[i + 1]] for i in range(len(walk) - 1)) and walk.count(x) == 2 and walk.count(y) == 2 \
                    and all(walk.count(u) == 1 for u in range(n) if u not in (x, y)) and all(walk[i] != walk[i + 1] for i in range(len(walk) - 1))
                if not good: ok = False
    report("T6 (E) the Camion construction gives a walk through x twice, y twice, every other vertex once, copies non-adjacent", ok, "%d ordered pairs on cycles (all labelled n=3..5, 80 random n=6,7); %d pairs not both on cycles" % (cnt, cnt_no))
    report("T6b some walk visits x and y twice each iff both lie on cycles (literal count, all labelled n=3..5)", ok_iff)

# ------------------------------------------------------------------ T7 the odd-order mean identity (README 6.2)
def T7():
    ok = True
    # closed forms: E[P_j] = n! C(n-1,j-1) 2^(j-n) (linearity of expectation); both sides have EGF ((2+t)/(2-t))^x
    for n in range(1, 13):
        for x in range(1, 6):
            lhs = sum(Fraction(math.factorial(n) * math.comb(n - 1, j - 1), 2 ** (n - j)) * math.comb(x, j) for j in range(1, n + 1))
            # coefficient of t^n in ((1+t/2)/(1-t/2))^x times n!
            num = [Fraction(math.comb(x, i), 2 ** i) for i in range(n + 1)]                      # (1+t/2)^x
            den = [Fraction(math.comb(x + i - 1, i), 2 ** i) for i in range(n + 1)]              # (1-t/2)^(-x)
            rhs = math.factorial(n) * sum(num[i] * den[n - i] for i in range(n + 1))
            if lhs != rhs: ok = False
            if n <= 7:      # the sum over permutations of odd order, by brute force
                br = Fraction(0)
                for perm in itertools.permutations(range(n)):
                    seen = [0] * n; lens = []
                    for i in range(n):
                        if not seen[i]:
                            l = 0; j = i
                            while not seen[j]: seen[j] = 1; j = perm[j]; l += 1
                            lens.append(l)
                    if any(l % 2 == 0 for l in lens): continue
                    br += Fraction(2 ** sum(1 for l in lens if l > 1), 2 ** sum(l for l in lens if l > 1)) * x ** len(lens)
                if br != lhs: ok = False
    okb = True
    for n in range(1, 5):
        tot = [0] * (n + 2); N = 0
        for A in all_t(n):
            P = Pj_bf(A); N += 1
            for j in range(1, n + 1): tot[j] += P[j]
        for j in range(1, n + 1):
            if Fraction(tot[j], N) != Fraction(math.factorial(n) * math.comb(n - 1, j - 1), 2 ** (n - j)): okb = False
    report("T7 mean of pi_T(x) over labelled tournaments = n! [t^n] ((2+t)/(2-t))^x = sum over odd-order permutations of 2^(psi-moved) x^cyc", ok, "closed forms n <= 12, x = 1..5; permutation sum by brute force n <= 7")
    report("T7b E[P_j] = n! C(n-1,j-1) 2^(j-n) by brute force over all labelled tournaments, n <= 4", okb)

if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    fn = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "..", "line_seeds.txt")
    T1(); T2(); T3(fn); T4(); T5(); T6(); T7()
    print("ALL PASS" if not FAILS else "FAILURES: %s" % FAILS)
    sys.exit(1 if FAILS else 0)
