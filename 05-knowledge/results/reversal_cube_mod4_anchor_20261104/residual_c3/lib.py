"""Helpers: substitution tournaments, primality, exact coefficient calls (exploratory)."""
import itertools, subprocess, os, random
HERE = os.path.dirname(os.path.abspath(__file__))
EXE = os.path.join(HERE, 'fcoef.exe')

def subst(Q, children):
    """Q: m x m 0/1 matrix (Q[i][j]=1 iff i->j). children: list of k_i x k_i matrices.
    Returns (adj, blocks) with names assigned child by child."""
    blocks, n = [], 0
    for C in children:
        blocks.append(list(range(n, n + len(C)))); n += len(C)
    A = [[0] * n for _ in range(n)]
    for i, Bi in enumerate(blocks):
        for p, u in enumerate(Bi):
            for q, v in enumerate(Bi):
                A[u][v] = children[i][p][q]
        for j, Bj in enumerate(blocks):
            if i != j and Q[i][j]:
                for u in Bi:
                    for v in Bj:
                        A[u][v] = 1
    return A, blocks

def TT(k):
    return [[1 if p < q else 0 for q in range(k)] for p in range(k)]

def conv(M):
    k = len(M); return [[M[q][p] for q in range(k)] for p in range(k)]

def is_module(A, S):
    n = len(A); S = set(S)
    for v in range(n):
        if v in S: continue
        vals = {A[v][u] for u in S}
        if len(vals) > 1: return False
    return True

def is_prime(Q):
    m = len(Q)
    for r in range(2, m):
        for S in itertools.combinations(range(m), r):
            if is_module(Q, S): return False
    return True

def random_tournament(m, rng):
    Q = [[0] * m for _ in range(m)]
    for i in range(m):
        for j in range(i + 1, m):
            if rng.random() < 0.5: Q[i][j] = 1
            else: Q[j][i] = 1
    return Q

def coeffs(tours):
    """Exact coefficient vectors (dict (S,R)->int) for a list of adjacency matrices."""
    inp = [str(len(tours))]
    for A in tours:
        inp.append(str(len(A)))
        inp += [''.join(str(x) for x in row) for row in A]
    out = subprocess.run([EXE], input='\n'.join(inp) + '\n', capture_output=True, text=True, check=True).stdout.split('\n')
    res = []
    for k, A in enumerate(tours):
        n, d = map(int, out[2 * k].split())
        vals = list(map(int, out[2 * k + 1].split()))
        D = [()] + [(i,) for i in range(n)] + [p for p in itertools.combinations(range(n), 2)]
        res.append({(D[i], D[j]): vals[i * d + j] for i in range(d) for j in range(d)})
    return res
