# mac-mini-2026-11-06-microcosm. THM-4900 note, permutohedron proposition: the full descent polynomial of the gated tournament is a sum over ordered 3-partitions (faces of the core permutohedron) with junction weights.
import itertools, random
from functools import lru_cache
def desc_poly_sub(n, arc, S):
    """descent polynomial (list of coefficients) of the subtournament on vertex set S (tuple)."""
    m = len(S)
    if m == 0: return [1]
    a = [0] * m
    for p in itertools.permutations(S):
        k = sum(1 for i in range(m - 1) if not arc[p[i]][p[i + 1]])
        a[k] += 1
    return a
def padd(p, q):
    r = [0] * max(len(p), len(q))
    for i, c in enumerate(p): r[i] += c
    for i, c in enumerate(q): r[i] += c
    return r
def pmul(p, q):
    r = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q): r[i + j] += a * b
    return r
def shift(p, k): return [0] * k + p
def predicted_gate_poly(n, arc):
    V = range(n); tot = [0]
    for lab in itertools.product(range(3), repeat=n):          # ordered 3-partition (A, B, C)
        A = tuple(v for v in V if lab[v] == 0); B = tuple(v for v in V if lab[v] == 1); C = tuple(v for v in V if lab[v] == 2)
        base = pmul(pmul(desc_poly_sub(n, arc, A), desc_poly_sub(n, arc, B)), desc_poly_sub(n, arc, C))
        w1 = (1 if A else 0) + (1 if not B else 0) + (1 if C else 0)    # order alpha s beta t gamma
        w2 = 2 * (1 if B else 0)                                       # order alpha t beta s gamma
        tot = padd(tot, padd(shift(base, w1), shift(base, w2)))
    while len(tot) > 1 and tot[-1] == 0: tot.pop()
    return tot
def gate(n, arc):
    m = n + 2; s, t = n, n + 1
    A = [[False] * m for _ in range(m)]
    for i in range(n):
        for j in range(n): A[i][j] = arc[i][j]
        A[s][i] = True; A[i][t] = True
    A[t][s] = True
    return m, A
random.seed(5); ok = True
for n in range(1, 6):
    for trial in range(15):
        arc = [[False] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                if random.random() < 0.5: arc[i][j] = True
                else: arc[j][i] = True
        m, A = gate(n, arc)
        actual = desc_poly_sub(m, A, tuple(range(m)))
        while len(actual) > 1 and actual[-1] == 0: actual.pop()
        ok &= actual == predicted_gate_poly(n, arc)
print('gate transform of the full descent polynomial (3-fold subset convolution with junction weights) verified:', ok)
