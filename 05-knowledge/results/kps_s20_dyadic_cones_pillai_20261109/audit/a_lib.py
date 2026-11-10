# a_lib.py -- S20 audit helpers (independent of the packet's tl.py).
# Tournaments are adjacency matrices A (list of lists), A[i][j] = 1 iff i -> j.
import itertools, subprocess, os, random

HERE = os.path.dirname(os.path.abspath(__file__))
EXE = os.path.join(HERE, 'a_paths.exe')

def encode(A):
    n = len(A)
    return '%d %s' % (n, ''.join('1' if A[i][j] else '0' for i in range(n) for j in range(i + 1, n))) if n >= 2 else '%d' % n

def run(tours):
    """batch call of a_paths.exe; returns list of dicts H, Q, e[], R[], s[]"""
    inp = '\n'.join(encode(A) for A in tours) + '\n'
    out = subprocess.run([EXE], input=inp, capture_output=True, text=True, check=True).stdout.strip().split('\n')
    res = []
    for line in out:
        a, e, R, s = line.split('|')
        H, Q = map(int, a.split())
        res.append(dict(H=H, Q=Q, e=list(map(int, e.split())), R=list(map(int, R.split())), s=list(map(int, s.split()))))
    assert len(res) == len(tours)
    return res

def triple(A, v, res=None):
    r = res if res is not None else run([A])[0]
    return (r['H'], r['e'][v], r['R'][v])

# ---------- naive (definition-level) counters, for validation on small n ----------
def ham_paths(A, verts):
    verts = list(verts)
    if not verts: return [()]
    out = []
    for perm in itertools.permutations(verts):
        if all(A[perm[i]][perm[i + 1]] for i in range(len(perm) - 1)): out.append(perm)
    return out

def naive_H(A, verts=None):
    if verts is None: verts = range(len(A))
    return len(ham_paths(A, verts))

def paths_from(A, v, avail):
    """all directed paths starting at v using vertices of avail (v in avail), as tuples"""
    res = []
    def dfs(path, used):
        res.append(tuple(path))
        for w in avail:
            if w not in used and A[path[-1]][w]:
                path.append(w); used.add(w); dfs(path, used); path.pop(); used.discard(w)
    dfs([v], {v})
    return res

def naive_stats(A):
    n = len(A); V = list(range(n))
    hp = ham_paths(A, V)
    H = len(hp)
    e = [sum(1 for p in hp if p[-1] == v) for v in V] if n else []
    s = [sum(1 for p in hp if p[0] == v) for v in V] if n else []
    R = []
    for v in V:
        tot = 0
        for P in paths_from(A, v, V):
            rest = [w for w in V if w not in P]
            tot += naive_H(A, rest)
        R.append(tot)
    Q = 0
    for mask in range(1 << n):
        S = [w for w in V if mask >> w & 1]; Sc = [w for w in V if not mask >> w & 1]
        Q += naive_H(A, S) * naive_H(A, Sc)
    return dict(H=H, Q=Q, e=e, R=R, s=s)

# ---------- constructions ----------
def empty(n):
    return [[0] * n for _ in range(n)]

def TT(m):
    """transitive: vertex i -> k for i < k (index 0 = g_1 = source, index m-1 = g_m = sink)"""
    A = empty(m)
    for i in range(m):
        for k in range(i + 1, m): A[i][k] = 1
    return A

def rand_tour(n, rng):
    A = empty(n)
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

def ordinal_sum(X, Y):
    """X => Y : every vertex of X beats every vertex of Y; X first, then Y (offset len(X))"""
    a, b = len(X), len(Y); A = empty(a + b)
    for i in range(a):
        for j in range(a): A[i][j] = X[i][j]
        for j in range(b): A[i][a + j] = 1
    for i in range(b):
        for j in range(b): A[a + i][a + j] = Y[i][j]
    return A

def one_arc_join(T, v, G, u):
    """direct convention: N = T => G with the single arc u -> v (u in G, v in T) instead of v -> u.
       G's vertices are offset by len(T)."""
    N = ordinal_sum(T, G); a = len(T)
    N[v][a + u] = 0; N[a + u][v] = 1
    return N

def from_code(code, n):
    """S19 class-file code: bits over pairs (i<j) in order i then j, bit set means i->j"""
    A = empty(n); k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if code >> k & 1: A[i][j] = 1
            else: A[j][i] = 1
            k += 1
    return A

def is_strong(A):
    n = len(A)
    if n == 0: return False
    def reach(fwd):
        seen = {0}; st = [0]
        while st:
            x = st.pop()
            for y in range(n):
                if y not in seen and (A[x][y] if fwd else A[y][x]): seen.add(y); st.append(y)
        return len(seen) == n
    return reach(True) and reach(False)

def ispow2(x): return x > 0 and x & (x - 1) == 0
def oddpart(x):
    while x % 2 == 0: x //= 2
    return x
def v2(x):
    k = 0
    while x % 2 == 0: x //= 2; k += 1
    return k
