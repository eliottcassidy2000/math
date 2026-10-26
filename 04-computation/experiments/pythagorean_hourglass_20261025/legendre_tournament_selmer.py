# "Quadratic reciprocity describes the 4-vertex tournaments": the precise version already in the repo
# (bsd_selmer_growth_20261013a.md section 3).  For primes p_i = 3 (mod 8), QR gives A_ij + A_ji = 1 (A_ij = [p_j/p_i],
# [a/p] = (1 - (a/p))/2), so the Legendre relation is a tournament; with A_ii = row parity, the pure 2-Selmer dimension
# of the congruent-number curve E_n (n = prod p_i) is nullity of the Monsky matrix = 2 dim_F4 ker(A^2 + A + I).
# Check: (i) the value for each of the four isomorphism classes of 4-tournaments; (ii) arithmetic witnesses with the
# full Monsky matrix (independent of the reduction); (iii) all four classes occur.  mac-mini-2026-10-25-pyth.
import itertools
import numpy as np
from sympy import primerange, legendre_symbol

def nullity_F2(M):
    M = [list(map(lambda v: v % 2, row)) for row in M]; n = len(M); m = len(M[0]); r = 0
    for c in range(m):
        piv = next((i for i in range(r, n) if M[i][c]), None)
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]
        for i in range(n):
            if i != r and M[i][c]: M[i] = [(a ^ b) for a, b in zip(M[i], M[r])]
        r += 1
    return m - r

def cls(adj):
    n = len(adj); sc = tuple(sorted(sum(adj[i][j] for j in range(n) if j != i) for i in range(n)))
    return {(0, 1, 2, 3): 'TT4', (1, 1, 1, 3): 'source+C3', (0, 2, 2, 2): 'sink+C3', (1, 1, 2, 2): 'strong C4'}[sc]

# (i) abstract: ker(A^2 + A + I) over F2 for each class (A = adjacency + row-parity diagonal)
seen = {}
for bits in itertools.product((0, 1), repeat=6):
    adj = [[0] * 4 for _ in range(4)]; k = 0
    for i in range(4):
        for j in range(i + 1, 4):
            adj[i][j] = bits[k]; adj[j][i] = 1 - bits[k]; k += 1
    A = [[adj[i][j] if i != j else sum(adj[i]) % 2 for j in range(4)] for i in range(4)]
    A2 = [[sum(A[i][t] * A[t][j] for t in range(4)) % 2 for j in range(4)] for i in range(4)]
    Q = [[(A2[i][j] + A[i][j] + (1 if i == j else 0)) % 2 for j in range(4)] for i in range(4)]
    seen.setdefault(cls(adj), set()).add(nullity_F2(Q))
print("(i) dim_F2 ker(A^2+A+I) by class:", {c: sorted(v) for c, v in seen.items()})
# (ii)/(iii) arithmetic witnesses: 4 primes = 3 mod 8, Monsky matrix nullity computed directly
def bracket(a, p): return (1 - legendre_symbol(a % p, p)) // 2
def monsky_nullity(ps):
    k = len(ps)
    A = [[bracket(ps[j], ps[i]) if i != j else 0 for j in range(k)] for i in range(k)]
    for i in range(k): A[i][i] = sum(A[i]) % 2
    D2 = [bracket(2, p) for p in ps]; Dm2 = [bracket(-2, p) for p in ps]
    M = [[0] * (2 * k) for _ in range(2 * k)]
    for i in range(k):
        for j in range(k):
            M[i][j] = (A[i][j] + (D2[i] if i == j else 0)) % 2
            M[i][k + j] = D2[i] if i == j else 0
            M[k + i][j] = D2[i] if i == j else 0
            M[k + i][k + j] = (A[i][j] + (Dm2[i] if i == j else 0)) % 2
    return nullity_F2(M)
P = [p for p in primerange(3, 400) if p % 8 == 3]
wit = {}
for quad in itertools.combinations(P, 4):
    adj = [[1 if (i != j and legendre_symbol(quad[j] % quad[i], quad[i]) == 1) else 0 for j in range(4)] for i in range(4)]
    assert all(adj[i][j] + adj[j][i] == 1 for i in range(4) for j in range(4) if i != j)   # QR: a tournament
    c = cls(adj)
    if c not in wit:
        wit[c] = (quad, monsky_nullity(list(quad)))
    if len(wit) == 4: break
print("(ii) first arithmetic witnesses (n = product, Monsky nullity = pure 2-Selmer dimension):")
for c, (quad, nu) in wit.items():
    print(f"     {c:10s}: primes {quad}, n = {np.prod(quad)}, nullity {nu}")
