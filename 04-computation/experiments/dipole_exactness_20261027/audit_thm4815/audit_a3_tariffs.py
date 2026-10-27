"""audit_a3_tariffs.py -- checks of THM-4815 A3 per-move laws on every labelled move at n = 4, 5, 6:
c3' - c3 = s(a) - s(b) - 1 exactly; Hamming distance to the (converse-closed) transitive bank, i.e. fas, changes by at
most n - 1; and the path skew matrix A_Q has operator norm 2cos(pi/(n+1)) and nuclear norm sum_k |2cos(k pi/(n+1))|.
Run: python3 -I audit_a3_tariffs.py"""
import itertools, math
try:
    import numpy as np
except Exception:
    np = None

def tours(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for code in range(1 << len(pairs)):
        A = [[0] * n for _ in range(n)]
        for k, (i, j) in enumerate(pairs):
            if code >> k & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

def hps(A):
    n = len(A); res = []
    def rec(p, used):
        if len(p) == n: res.append(tuple(p)); return
        for v in range(n):
            if not used >> v & 1 and A[p[-1]][v]: p.append(v); rec(p, used | 1 << v); p.pop()
    for s in range(n): rec([s], 1 << s)
    return res

def tau(A, Q):
    n = len(A); B = [[A[j][i] if i != j else 0 for j in range(n)] for i in range(n)]
    for i in range(n - 1):
        u, v = Q[i], Q[i + 1]; B[u][v] = 1; B[v][u] = 0
    return B

def c3(A):
    n = len(A); s = [sum(r) for r in A]
    return math.comb(n, 3) - sum(math.comb(x, 2) for x in s)

def fas(A):
    n = len(A); best = None
    for p in itertools.permutations(range(n)):
        pos = {v: i for i, v in enumerate(p)}
        b = sum(1 for u in range(n) for v in range(n) if A[u][v] and pos[u] > pos[v])
        best = b if best is None or b < best else best
    return best

for n in (4, 5, 6):
    bad_c3 = bad_fas = moves = 0
    fcache = {}
    for A in tours(n):
        s = [sum(r) for r in A]; c = c3(A)
        key = tuple(map(tuple, A))
        fA = fcache.setdefault(key, fas(A)) if n <= 5 else None
        for Q in hps(A):
            moves += 1
            B = tau(A, Q)
            if c3(B) - c != s[Q[0]] - s[Q[-1]] - 1: bad_c3 += 1
            if n <= 5:
                kB = tuple(map(tuple, B)); fB = fcache.setdefault(kB, fas(B))
                if abs(fB - fA) > n - 1: bad_fas += 1
    print(f"n={n}: labelled moves {moves}; c3-law failures {bad_c3}; fas tariff (<= n-1) failures {bad_fas if n <= 5 else 'not run'}")

if np is not None:
    for n in range(3, 13):
        S = np.zeros((n, n))
        for i in range(n - 1): S[i, i + 1] = 1; S[i + 1, i] = -1
        sv = np.linalg.svd(S, compute_uv=False)
        op, nuc = sv.max(), sv.sum()
        op_f = 2 * math.cos(math.pi / (n + 1)); nuc_f = sum(abs(2 * math.cos(k * math.pi / (n + 1))) for k in range(1, n + 1))
        print(f"n={n}: ||A_Q||_op={op:.12f} formula {op_f:.12f}; ||A_Q||_*={nuc:.12f} formula {nuc_f:.12f}")
else:
    print('numpy unavailable: norm check skipped')
