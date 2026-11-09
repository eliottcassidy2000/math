"""AUDIT2 for THM-5041 (independent code; H by brute force over all vertex orderings, no DP).
Checks: (1a) bijection/bit count, (1c) class table and counts, (1d) transpose, complement, reversal relabelling
(U, L) -> (s L, s U) with s(X)_{ij} = X_{n+1-j, n+1-i}; (2a) loops irrelevant, H(U,L) = H(L,U), H = 2 hp on equal
pairs, THM-5071 (3) parity by explicit packings of defect cycles; (2b) tile-wise additivity at EVERY pair of random
digraphs; (2c) the flip corollary, parity and difference statements, and the signed expansion over reciprocal sets;
(3) Session234 examples and DG2 in (U, L) coordinates (exhaustive for semicomplete n <= 5).
"""
import itertools, random, sys

def H(A):
    n = len(A)
    if n == 0: return 1
    return sum(1 for p in itertools.permutations(range(n)) if all(A[p[i]][p[i + 1]] for i in range(n - 1)))

def paths(A):
    n = len(A)
    return [p for p in itertools.permutations(range(n)) if all(A[p[i]][p[i + 1]] for i in range(n - 1))]

def UL(A):
    n = len(A)
    U = {(i, j): A[i][j] for i in range(n) for j in range(i + 1, n)}
    L = {(i, j): A[j][i] for i in range(n) for j in range(i + 1, n)}
    return U, L

def s_refl(X, n):   # (sX)_{ij} = X_{n-1-j, n-1-i} (0-based form of n+1-j, n+1-i)
    return {(i, j): X[(n - 1 - j, n - 1 - i)] for (i, j) in X}

def rand_digraph(n, rng, loops=True):
    return [[(rng.randint(0, 1) if (i != j or loops) else 0) for j in range(n)] for i in range(n)]

def setpair(A, i, j, a, b):
    B = [r[:] for r in A]; B[i][j] = a; B[j][i] = b; return B

rng = random.Random(5041)
bad = {}
def flag(name, cond):
    bad.setdefault(name, [0, 0]); bad[name][1] += 1
    if not cond: bad[name][0] += 1

# (1a)/(1c) exhaustive at n = 3 (all 2^9 matrices)
n = 3
seen = set(); cls = {'tour': 0, 'semi': 0, 'orient': 0, 'sym': 0, 'loopless': 0}
for bits in range(2**(n * n)):
    A = [[bits >> (n * i + j) & 1 for j in range(n)] for i in range(n)]
    U, L = UL(A); diag = tuple(A[i][i] for i in range(n))
    key = (diag, tuple(sorted(U.items())), tuple(sorted(L.items())))
    seen.add(key)
    if any(diag): continue
    cls['loopless'] += 1
    tour = all(A[i][j] + A[j][i] == 1 for i in range(n) for j in range(i + 1, n))
    semi = all(A[i][j] + A[j][i] >= 1 for i in range(n) for j in range(i + 1, n))
    orient = all(A[i][j] + A[j][i] <= 1 for i in range(n) for j in range(i + 1, n))
    sym = all(A[i][j] == A[j][i] for i in range(n) for j in range(n))
    flag('class tournament <=> L = complement U', tour == all(L[p] == 1 - U[p] for p in U))
    flag('class semicomplete <=> U or L = 1 everywhere', semi == all(U[p] | L[p] for p in U))
    flag('class oriented <=> U and L = 0 everywhere', orient == all(not (U[p] & L[p]) for p in U))
    flag('class undirected <=> U = L', sym == all(U[p] == L[p] for p in U))
    cls['tour'] += tour; cls['semi'] += semi; cls['orient'] += orient; cls['sym'] += sym
flag('(1a) bijection at n=3', len(seen) == 2**9)
flag('(1c) counts 4^3, 3^3, 2^3 at n=3', (cls['loopless'], cls['semi'], cls['tour']) == (64, 27, 8))

# random digraphs: (1d), (2a), (2b)
for trial in range(400):
    n = rng.randint(2, 6)
    A = rand_digraph(n, rng)
    h = H(A)
    U, L = UL(A)
    # transpose
    At = [[A[j][i] for j in range(n)] for i in range(n)]
    Ut, Lt = UL(At)
    flag('(1d) transpose (U,L)->(L,U)', Ut == L and Lt == U)
    flag('(2a) H(transpose) = H', H(At) == h)
    # complement J - A
    Ac = [[1 - A[i][j] for j in range(n)] for i in range(n)]
    Uc, Lc = UL(Ac)
    flag('(1d) complement J-A', all(Uc[p] == 1 - U[p] and Lc[p] == 1 - L[p] for p in U) and all(Ac[i][i] == 1 - A[i][i] for i in range(n)))
    # reversal relabelling i -> n+1-i
    Ar = [[A[n - 1 - i][n - 1 - j] for j in range(n)] for i in range(n)]
    Ur, Lr = UL(Ar)
    flag('(1d) reversal (U,L)->(sL,sU)', Ur == s_refl(L, n) and Lr == s_refl(U, n))
    # loops
    Al = [r[:] for r in A]
    for i in range(n): Al[i][i] = rng.randint(0, 1)
    flag('(2a) loops irrelevant', H(Al) == h)
    # (2b) at every pair
    for i in range(n):
        for j in range(i + 1, n):
            v = {(a, b): H(setpair(A, i, j, a, b)) for a in (0, 1) for b in (0, 1)}
            flag('(2b) tile additivity (every pair)', v[(1, 0)] + v[(0, 1)] == v[(1, 1)] + v[(0, 0)])

# equal pairs: H = 2 hp
for trial in range(200):
    n = rng.randint(2, 7)
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            A[i][j] = A[j][i] = rng.randint(0, 1)
    hp = len({min(p, p[::-1]) for p in paths(A)})
    flag('(2a) H = 2 hp on equal pairs', H(A) == 2 * hp)

# THM-5071 (3) parity by explicit packings of defect cycles, on general random loopless digraphs
def cycles(A):
    """directed cycles of length >= 2, as (frozenset of vertices, tuple of arcs); each cycle once"""
    n = len(A); out = []
    for k in range(2, n + 1):
        for S in itertools.combinations(range(n), k):
            first = S[0]
            for rest in itertools.permutations(S[1:]):
                cyc = (first,) + rest
                arcs = [(cyc[t], cyc[(t + 1) % k]) for t in range(k)]
                if k == 2 and rest[0] < first: continue
                if all(A[a][b] for a, b in arcs):
                    out.append((frozenset(S), arcs))
    return out

def packings(cycs):
    """number of sets of pairwise vertex-disjoint cycles (empty set included)"""
    cycs = list(cycs)
    def rec(i, used):
        if i == len(cycs): return 1
        tot = rec(i + 1, used)
        if not (cycs[i][0] & used): tot += rec(i + 1, used | cycs[i][0])
        return tot
    return rec(0, frozenset())

for trial in range(300):
    n = rng.randint(2, 6)
    A = rand_digraph(n, rng, loops=False)
    Ab = [[(1 - A[i][j]) if i != j else 0 for j in range(n)] for i in range(n)]
    B = {frozenset((i, j)) for i in range(n) for j in range(i + 1, n) if A[i][j] and A[j][i]}
    Z = {frozenset((i, j)) for i in range(n) for j in range(i + 1, n) if not A[i][j] and not A[j][i]}
    defect = [c for c in cycles(A) if any(frozenset(a) in B for a in c[1])] + \
             [c for c in cycles(Ab) if any(frozenset(a) in Z for a in c[1])]
    flag('THM-5071 (3): H = #packings of defect cycles (mod 2), general digraphs', H(A) % 2 == packings(defect) % 2)

# (2c) flip corollary on tournaments, at every tile
for trial in range(150):
    n = rng.randint(2, 7)
    T = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if rng.randint(0, 1): T[i][j] = 1
            else: T[j][i] = 1
    hT = H(T)
    for i in range(n):
        for j in range(i + 1, n):
            Tp = setpair(T, i, j, T[j][i], T[i][j])
            hr = H(setpair(T, i, j, 1, 1)); hm = H(setpair(T, i, j, 0, 0)); hTp = H(Tp)
            flag('(2c) H(recip) + H(missing) = H(T) + H(T^+p)', hr + hm == hT + hTp)
            flag('(2c) one-defect totals have equal parity', (hr - hm) % 2 == 0)
            thru = sum(1 for P in paths(T) if any({P[t], P[t + 1]} == {i, j} for t in range(n - 1))) + \
                   sum(1 for P in paths(Tp) if any({P[t], P[t + 1]} == {i, j} for t in range(n - 1)))
            flag('(2c) H(recip) - H(missing) = HPs of T or T^+p through p', hr - hm == thru)

# (2c) signed expansion over reciprocal sets of random semicomplete digraphs
for trial in range(120):
    n = rng.randint(2, 6)
    D = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            st = rng.choice([(1, 0), (0, 1), (1, 1), (1, 1)])
            D[i][j], D[j][i] = st
    R = [(i, j) for i in range(n) for j in range(i + 1, n) if D[i][j] and D[j][i]]
    tot = 0
    for r in range(len(R) + 1):
        for Aset in itertools.combinations(R, r):
            rest = [p for p in R if p not in Aset]
            for orient in itertools.product((0, 1), repeat=len(rest)):
                X = [row[:] for row in D]
                for (i, j) in Aset: X[i][j] = X[j][i] = 0
                for (i, j), o in zip(rest, orient):
                    X[i][j], X[j][i] = (1, 0) if o else (0, 1)
                tot += (-1) ** r * H(X)
    flag('(2c) signed expansion over reciprocal sets', tot == H(D))

# (3) Session234 examples and DG2 in (U, L) coordinates
def parse(rows): return [[int(c) for c in r] for r in rows]
ex = {'7 (order 4)': parse(['0111', '0011', '1001', '1000']),
      '21 (order 5)': parse(['01111', '00110', '00011', '10001', '11000']),
      '21 (order 6, one reciprocal pair)': parse(['011110', '001111', '000111', '100011', '000001', '100000'])}
for name, A in ex.items():
    U, L = UL(A)
    st = [(U[p], L[p]) for p in U]
    print("Session234 example %s: H = %d, tiles (1,1) = %d, (0,0) = %d, diagonal = %s"
          % (name, H(A), st.count((1, 1)), st.count((0, 0)), [A[i][i] for i in range(len(A))]))
A6 = ex['21 (order 6, one reciprocal pair)']
print("   order-6 example: resolutions and missing state at {0,3}:", H(setpair(A6, 0, 3, 1, 0)), H(setpair(A6, 0, 3, 0, 1)), H(setpair(A6, 0, 3, 0, 0)))
# DG2 exhaustive on semicomplete digraphs n <= 5 (with H computed by brute force)
for n in range(2, 6):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    hist = {}
    for states in itertools.product(((1, 0), (0, 1), (1, 1)), repeat=len(pairs)):
        A = [[0] * n for _ in range(n)]
        for (i, j), (a, b) in zip(pairs, states): A[i][j], A[j][i] = a, b
        h = H(A)
        if h in (7, 21):
            b = states.count((1, 1)); hist.setdefault(h, set()).add(b)
    print("semicomplete n=%d: reciprocal-pair counts at H=7: %s, at H=21: %s" % (n, sorted(hist.get(7, [])), sorted(hist.get(21, []))))
    flag('DG2 (semicomplete n<=5): H=7 => >=2 (1,1) tiles', all(b >= 2 for b in hist.get(7, [])))
    flag('DG2 (semicomplete n<=5): H=21 => >=1 (1,1) tile', all(b >= 1 for b in hist.get(21, [])))

allok = True
for k, (nb, nt) in bad.items():
    print("%-75s tests %6d  failures %d" % (k, nt, nb))
    allok &= nb == 0
print("ALL:", "PASS" if allok else "FAIL")
