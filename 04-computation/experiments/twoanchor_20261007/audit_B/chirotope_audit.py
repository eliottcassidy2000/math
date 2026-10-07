#!/usr/bin/env python3
"""Audit B (independent) of THM-4602 (1)-(2): 4-Pfaffians, rank-2 chirotopes, local transitivity,
total positivity and SL2 friezes.  Written from scratch; shares no code with tournament_chirotope.py.

Conventions: B[i][j] = +1 iff i -> j, B[j][i] = -B[i][j].
  (a) realisable : sign det(v_i, v_j) = B[i][j] for some vectors v_i in R^2 (pairwise independent)
  (b) every 4-Pfaffian b_ij b_kl - b_ik b_jl + b_il b_jk has absolute value 1
  (c) locally transitive: every out- and in-neighbourhood induces a transitive tournament
Realisability is certified CONSTRUCTIVELY with integer vectors: switch T so that vertex 0 is a source;
if the result is transitive with source-first order r(.), put u_i = eps_i * (1, r(i)) and check signs exactly.
"""
import itertools, math, random, sys
from fractions import Fraction
import numpy as np

OUT = []
def say(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); OUT.append(s)

def perm_sign(p):
    s, seen = 1, [False]*len(p)
    for i in range(len(p)):
        if not seen[i]:
            j, L = i, 0
            while not seen[j]: seen[j] = True; j = p[j]; L += 1
            if L % 2 == 0: s = -s
    return s

def pf(B, i, j, k, l):
    return B[i][j]*B[k][l] - B[i][k]*B[j][l] + B[i][l]*B[j][k]

# ------------------------------------------------------------------ Part 1: all 64 labelled 4-tournaments
say("== Part 1: the 64 labelled 4-tournaments ==")
pairs4 = list(itertools.combinations(range(4), 2))
tally = {}
for bits in range(64):
    B = [[0]*4 for _ in range(4)]
    for idx, (i, j) in enumerate(pairs4):
        s = 1 if (bits >> idx) & 1 else -1
        B[i][j], B[j][i] = s, -s
    P = pf(B, 0, 1, 2, 3)
    sc = sorted(sum(1 for j in range(4) if B[i][j] == 1) for i in range(4))
    c3 = 4 - sum(s*(s-1)//2 for s in sc)            # Moon/Kendall-Babington Smith count of 3-cycles
    # brute-force 3-cycle count as a second, independent count
    c3b = sum(1 for (i, j, k) in itertools.combinations(range(4), 3)
              if (B[i][j] == B[j][k] == B[k][i]))
    assert c3 == c3b
    assert P in (-3, -1, 1, 3)
    assert (abs(P) == 3) == (c3 == 1)
    # |Pf| = 3 iff the three products are equal
    terms = (B[0][1]*B[2][3], -B[0][2]*B[1][3], B[0][3]*B[1][2])
    assert (abs(P) == 3) == (len(set(terms)) == 1)
    # covariance under relabelling: Pf(B^sigma) = sgn(sigma) Pf(B)  (so |Pf| is an isomorphism invariant)
    for p in itertools.permutations(range(4)):
        Bp = [[B[p[a]][p[b]] for b in range(4)] for a in range(4)]
        assert pf(Bp, 0, 1, 2, 3) == perm_sign(p) * P
    # switching a vertex negates Pf (each term uses each index exactly once)
    for v in range(4):
        Bs = [[B[a][b] * (-1 if (a == v) != (b == v) else 1) for b in range(4)] for a in range(4)]
        assert pf(Bs, 0, 1, 2, 3) == -P
    key = {(0, 1, 2, 3): "TT4", (1, 1, 1, 3): "3-cycle+source", (0, 2, 2, 2): "3-cycle+sink",
           (1, 1, 2, 2): "strong (2 three-cycles)"}[tuple(sc)]
    tally.setdefault(key, {}).setdefault(P, 0)
    tally[key][P] += 1
for k in tally: say(f"  {k:26s} c3={'1' if 'cycle+' in k else ('0' if k=='TT4' else '2')}  Pf distribution {dict(sorted(tally[k].items()))}")
say("  CONFIRMED: Pf in {+-1,+-3}; |Pf|=3 iff exactly one 3-cycle (3-cycle+source or 3-cycle+sink);"
    " Pf(B^sigma)=sgn(sigma)Pf(B); switching a vertex negates Pf")

# ------------------------------------------------------------------ helpers for (a), (b), (c)
def lt_by_definition(B, n):
    for v in range(n):
        for side in (1, -1):
            S = [w for w in range(n) if w != v and B[v][w] == side]
            sc = [sum(1 for x in S if x != u and B[u][x] == 1) for u in S]
            if len(set(sc)) != len(S):           # transitive iff scores distinct
                return False
    return True

def all_pf_pm1(B, n):
    return all(abs(pf(B, *q)) == 1 for q in itertools.combinations(range(n), 4))

def realise(B, n):
    """Constructive realisation by integer vectors, or None."""
    eps = [1] + [B[0][j] for j in range(1, n)]
    S = [[eps[i]*eps[j]*B[i][j] if i != j else 0 for j in range(n)] for i in range(n)]
    sc = [sum(1 for j in range(n) if S[i][j] == 1) for i in range(n)]
    if len(set(sc)) != n: return None
    r = [n - 1 - s for s in sc]
    U = [(eps[i], eps[i]*r[i]) for i in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                d = U[i][0]*U[j][1] - U[i][1]*U[j][0]
                if (d > 0) - (d < 0) != B[i][j]: raise AssertionError("certificate failed")
    return U

def pattern_from_vectors(V):
    n = len(V); B = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                d = V[i][0]*V[j][1] - V[i][1]*V[j][0]
                assert d != 0
                B[i][j] = 1 if d > 0 else -1
    return B

# ------------------------------------------------------------------ Part 2a: exhaustive n = 4..7 (numpy)
say("== Part 2a: exhaustive over ALL labelled tournaments, n = 4..7 ==")
for n in range(4, 8):
    pairs = list(itertools.combinations(range(n), 2)); m = len(pairs)
    codes = np.arange(1 << m, dtype=np.uint32)
    s = {}
    for idx, (i, j) in enumerate(pairs):
        a = (((codes >> idx) & 1).astype(np.int8) * 2 - 1)
        s[(i, j)] = a; s[(j, i)] = -a
    def S(i, j): return s[(i, j)]
    # (b)
    okb = np.ones(1 << m, dtype=bool)
    for (i, j, k, l) in itertools.combinations(range(n), 4):
        P = S(i, j)*S(k, l) - S(i, k)*S(j, l) + S(i, l)*S(j, k)
        okb &= (np.abs(P) == 1)
    # (c) by definition: no vertex v dominating / dominated by a 3-cycle {a,b,c}
    okc = np.ones(1 << m, dtype=bool)
    for v in range(n):
        rest = [w for w in range(n) if w != v]
        for (a, b, c) in itertools.combinations(rest, 3):
            cyc = (S(a, b) == S(b, c)) & (S(b, c) == S(c, a))
            dom = (S(v, a) == S(v, b)) & (S(v, b) == S(v, c))
            okc &= ~(cyc & dom)
    # (a) constructive: canonical switch making 0 a source, transitivity, integer vectors
    eps = {0: np.ones(1 << m, dtype=np.int8)}
    for j in range(1, n): eps[j] = S(0, j)
    score = {}
    for i in range(n):
        sc = np.zeros(1 << m, dtype=np.int8)
        for j in range(n):
            if j != i: sc += ((eps[i]*eps[j]*S(i, j)) == 1).astype(np.int8)
        score[i] = sc
    mask = np.zeros(1 << m, dtype=np.int64)
    for i in range(n): mask |= (np.int64(1) << score[i].astype(np.int64))
    oka = (mask == (1 << n) - 1)
    # verify the integer-vector certificate on every candidate: sign(eps_i eps_j (r_j - r_i)) == s_ij
    cert = np.ones(1 << m, dtype=bool)
    for i in range(n):
        for j in range(n):
            if i != j:
                ri = (n - 1 - score[i]).astype(np.int16); rj = (n - 1 - score[j]).astype(np.int16)
                d = eps[i].astype(np.int16)*eps[j].astype(np.int16)*(rj - ri)
                cert &= (np.sign(d).astype(np.int8) == S(i, j))
    assert np.all(cert[oka]), "integer certificate failed"
    assert np.array_equal(oka, okb) and np.array_equal(okb, okc)
    cnt = int(okb.sum()); expect = math.factorial(n - 1) * 2**(n - 1)
    say(f"  n={n}: {1 << m} tournaments; #(a realisable, certified) = #(b all |Pf|=1) = #(c LT) = {cnt}"
        f"  [(n-1)! 2^(n-1) = {expect}]  sets identical: True")
    assert cnt == expect

# ------------------------------------------------------------------ Part 2b: random tests n = 4..9
say("== Part 2b: random tests, n = 4..9 (seed 20261007) ==")
rnd = random.Random(20261007)
for n in range(4, 10):
    stats = dict(rand=0, rand_lt=0, planar=0, switchedTT=0, dfs=0)
    # (i) uniform random tournaments
    for _ in range(3000 if n <= 7 else 1500):
        B = [[0]*n for _ in range(n)]
        for i, j in itertools.combinations(range(n), 2):
            x = rnd.choice((1, -1)); B[i][j], B[j][i] = x, -x
        a = realise(B, n) is not None; b = all_pf_pm1(B, n); c = lt_by_definition(B, n)
        assert a == b == c
        stats['rand'] += 1; stats['rand_lt'] += int(b)
    # (ii) random real planar configurations (exact rationals: random integer vectors), incl. both half-planes
    for _ in range(600):
        while True:
            V = [(rnd.randint(-50, 50), rnd.randint(-50, 50)) for _ in range(n)]
            if all(V[i][0]*V[j][1] - V[i][1]*V[j][0] != 0 for i, j in itertools.combinations(range(n), 2)): break
        B = pattern_from_vectors(V)
        # GP identity (a)=>(b) exactly
        for (i, j, k, l) in itertools.combinations(range(n), 4):
            d = lambda x, y: V[x][0]*V[y][1] - V[x][1]*V[y][0]
            assert d(i, j)*d(k, l) - d(i, k)*d(j, l) + d(i, l)*d(j, k) == 0
        assert all_pf_pm1(B, n) and lt_by_definition(B, n) and realise(B, n) is not None
        stats['planar'] += 1
    # (iii) random switchings of random transitive tournaments (all must be LT / realisable)
    for _ in range(600):
        order = list(range(n)); rnd.shuffle(order); pos = {v: r for r, v in enumerate(order)}
        sw = [rnd.choice((1, -1)) for _ in range(n)]
        B = [[0 if i == j else sw[i]*sw[j]*(1 if pos[i] < pos[j] else -1) for j in range(n)] for i in range(n)]
        assert all_pf_pm1(B, n) and lt_by_definition(B, n) and realise(B, n) is not None
        stats['switchedTT'] += 1
    # (iv) random members of class (b) produced WITHOUT geometry: random-order DFS on arcs, prune |Pf|=3
    pairs = list(itertools.combinations(range(n), 2))
    for _ in range(300):
        B = [[0]*n for _ in range(n)]
        order = pairs[:]; rnd.shuffle(order)
        def ok_after(i, j):
            for k, l in itertools.combinations([x for x in range(n) if x not in (i, j)], 2):
                q = sorted((i, j, k, l))
                if all(B[x][y] != 0 for x, y in itertools.combinations(q, 2)) and abs(pf(B, *q)) == 3:
                    return False
            return True
        def dfs(t):
            if t == len(order): return True
            i, j = order[t]
            for x in rnd.sample((1, -1), 2):
                B[i][j], B[j][i] = x, -x
                if ok_after(i, j) and dfs(t + 1): return True
            B[i][j] = B[j][i] = 0
            return False
        assert dfs(0)
        assert all_pf_pm1(B, n) and lt_by_definition(B, n)
        assert realise(B, n) is not None          # (b) => (a), constructive
        stats['dfs'] += 1
    say(f"  n={n}: {stats}  -> (a)<=>(b)<=>(c) on every sample; every (b)-member realised by integer vectors")

# ------------------------------------------------------------------ Part 3: total positivity and SL2 friezes
say("== Part 3: totally positive part and SL2 friezes ==")
def det(u, v): return u[0]*v[1] - u[1]*v[0]
# 3a. Over R the standard affine chart (1, x_1 < ... < x_m), (0, 1) is totally positive and its pattern is the
#     transitive tournament with infinity a sink -- the real counterpart of THM-4602 (3)(a).
for m in range(2, 10):
    for _ in range(100):
        xs = sorted(rnd.sample(range(-1000, 1000), m))
        V = [(1, x) for x in xs] + [(0, 1)]
        B = pattern_from_vectors(V)
        assert all(B[i][j] == 1 for i in range(m + 1) for j in range(i + 1, m + 1))
say("  R: the chart (1,x_1<...<x_m),(0,1) has transitive pattern (infinity = sink), m = 2..9")

def frieze_from_quiddity(a):
    """v_1=(1,0), v_2=(0,1), v_{i+1} = a_i v_i - v_{i-1} (i = 2..n+1); returns vectors and closure flag."""
    n = len(a)
    V = [(1, 0), (0, 1)]
    for i in range(1, n + 1):          # uses a_2..a_{n}, a_1 (cyclic) -> vectors v_3..v_{n+2}
        ai = a[i % n]
        V.append((ai*V[-1][0] - V[-2][0], ai*V[-1][1] - V[-2][1]))
    return V

# 3b. Conway-Coxeter: positive INTEGER friezes <-> triangulations, counted by Catalan C_{n-2}
def cc_count(n):
    cnt = 0
    def rec(a, V):
        nonlocal cnt
        k = len(V)
        if k == n + 2:
            # closure: v_{n+1} = -v_1, v_{n+2} = -v_2 (antiperiodic, monodromy -I)
            if V[n] == (-V[0][0], -V[0][1]) and V[n + 1] == (-V[1][0], -V[1][1]):
                cnt += 1
            return
        for x in range(1, n):
            w = (x*V[-1][0] - V[-2][0], x*V[-1][1] - V[-2][1])
            # positivity of all entries det(v_i, w) for i in the current window (i > k-n+1)
            if k < n and any(det(V[i], w) <= 0 for i in range(0, k)):
                continue
            rec(a + [x], V + [w])
    rec([], [(1, 0), (0, 1)])
    return cnt
cat = lambda m: math.comb(2*m, m)//(m + 1)
row = []
for n in range(3, 10):
    c = cc_count(n); row.append((n, c, cat(n - 2))); assert c == cat(n - 2)
say(f"  Conway-Coxeter positive integer friezes with n points (width n-3): (n, #friezes, C_(n-2)) = {row}")

# 3c. positive REAL friezes = the slice {p_(i,i+1) = 1, p_(1n) = 1} of Gr+(2,n); torus normalisation
def torus_normalise(P, n):
    """Solve lam_i lam_(i+1) P[i][i+1] = 1 (i=0..n-2), lam_(n-1) lam_0 P[0][n-1] = 1 in positive reals (log-linear).
    Returns (lam, residual)."""
    A = np.zeros((n, n)); rhs = np.zeros(n)
    for i in range(n - 1): A[i, i] = A[i, i + 1] = 1; rhs[i] = -math.log(P[i][i + 1])
    A[n - 1, n - 1] = A[n - 1, 0] = 1; rhs[n - 1] = -math.log(P[0][n - 1])
    sol, res, rank, _ = np.linalg.lstsq(A, rhs, rcond=None)
    return np.exp(sol), float(np.max(np.abs(A @ sol - rhs))), rank
nodd_ok, neven_fail, neven_total = 0, 0, 0
for n in range(4, 12):
    for _ in range(200):
        th = sorted(rnd.uniform(0.01, math.pi - 0.01) for _ in range(n))
        rr = [rnd.uniform(0.3, 3) for _ in range(n)]
        V = [(rr[i]*math.cos(t), rr[i]*math.sin(t)) for i, t in enumerate(th)]
        P = [[det(V[i], V[j]) for j in range(n)] for i in range(n)]
        assert all(P[i][j] > 0 for i in range(n) for j in range(i + 1, n))       # TP <=> transitive pattern
        lam, resid, rank = torus_normalise(P, n)
        if n % 2 == 1:
            assert rank == n and resid < 1e-9
            W = [(lam[i]*V[i][0], lam[i]*V[i][1]) for i in range(n)]
            Wx = W + [(-W[0][0], -W[0][1]), (-W[1][0], -W[1][1])]          # antiperiodic closure
            assert all(abs(det(Wx[i], Wx[i + 1]) - 1) < 1e-9 for i in range(n + 1))
            quid = [det(Wx[i - 1], Wx[i + 1]) for i in range(1, n + 1)]
            # frieze recurrence closes: v_(i+1) = a_i v_i - v_(i-1)
            for i in range(1, n + 1):
                pred = (quid[i - 1]*Wx[i][0] - Wx[i - 1][0], quid[i - 1]*Wx[i][1] - Wx[i - 1][1])
                assert abs(pred[0] - Wx[i + 1][0]) < 1e-7 and abs(pred[1] - Wx[i + 1][1]) < 1e-7
            assert all(det(W[i], W[j]) > 0 for i in range(n) for j in range(i + 1, n))
            nodd_ok += 1
        else:
            neven_total += 1
            # solvable iff p12 p34 ... p_(n-1,n) = p23 p45 ... p_(n-2,n-1) p_(1n)
            lhs = math.prod(P[i][i + 1] for i in range(0, n - 1, 2))
            rhs = math.prod(P[i][i + 1] for i in range(1, n - 2, 2)) * P[0][n - 1]
            solvable = abs(math.log(lhs/rhs)) < 1e-9
            assert solvable == (resid < 1e-9)
            if not solvable: neven_fail += 1
say(f"  odd n (5..11): {nodd_ok} random TP points, each torus-normalises to a unique positive real frieze (closure checked)")
say(f"  even n (4..10): {neven_fail}/{neven_total} random TP points are NOT torus-equivalent to any frieze "
    "(alternating-product obstruction p12 p34.. = p23 p45.. p1n)")
# n = 4: every width-1 frieze (a, 2/a, a, 2/a) has cross-ratio p12 p34 / (p13 p24) = 1/2
for a in [Fraction(1), Fraction(2), Fraction(3, 7), Fraction(11, 5)]:
    q = [a, 2/a, a, 2/a]
    V = [(Fraction(1), Fraction(0)), (Fraction(0), Fraction(1))]
    for i in range(1, 4): V.append((q[i]*V[-1][0] - V[-2][0], q[i]*V[-1][1] - V[-2][1]))
    assert det(V[3], (-V[0][0], -V[0][1])) == 1 and all(det(V[i], V[i + 1]) == 1 for i in range(3))
    cr = det(V[0], V[1])*det(V[2], V[3]) / (det(V[0], V[2])*det(V[1], V[3]))
    assert cr == Fraction(1, 2)
say("  n=4: all width-1 friezes (a,2/a,a,2/a) have the same cross-ratio 1/2 -> the frieze slice is NOT Gr+(2,4)/T")
say("ALL CHIROTOPE CHECKS PASSED")
with open(__file__.replace('.py', '.out'), 'w') as f: f.write("\n".join(OUT) + "\n")
