# audit_H1: (i) d = 3 vs d = 4 normalised rate N^2 rho/log^2 N (HYP-9490 "d = 4 agrees to within 1%"), (ii) saddle-point reading N^2 rho ~ (log K + 2)^2 with K (log K + 2) = N, (iii) Gamma^k(K1) = transitive order with the k concentric chords (s_i, t_i) reversed.
import math
from audit_tournament_lib import gate

RMAX = 130
W = {1: 1}
Wl = [dict(W)]
for r in range(RMAX):
    new = {}
    for k, c in W.items():
        if k == 0:
            continue
        for kk, wt in ((k + 1, k), (k, k * k), (k - 1, -k)):
            new[kk] = new.get(kk, 0) + wt * c
    W = new
    Wl.append(W)
U = lambda d, r: sum(c * (k ** d if d else 1) for k, c in Wl[r].items())
print(" r   N3   d=3 norm   N4   d=4 norm   rel.diff   saddle K   N^2 rho/(log K+2)^2 (d=3)")
for r in (10, 20, 40, 60, 80, 100, 120):
    vals = []
    for d in (3, 4):
        N = d + 2 * r
        rho = U(d - 2, r) / U(d, r)
        vals.append((N, rho * N * N / math.log(N) ** 2, rho))
    N3 = vals[0][0]
    # solve K (log K + 2) = N3
    K = N3 / 5.0
    for _ in range(100):
        K = N3 / (math.log(K) + 2)
    sad = vals[0][2] * N3 * N3 / (math.log(K) + 2) ** 2
    print(f"{r:3d} {vals[0][0]:4d} {vals[0][1]:9.4f} {vals[1][0]:5d} {vals[1][1]:9.4f} {abs(vals[1][1]/vals[0][1]-1):9.2%} {K:9.2f} {sad:9.4f}")
# (iii)
for k in range(1, 9):
    X = [[0]]
    for _ in range(k):
        X = gate(X)
    n = 2 * k + 1
    # my gate() appends s = old n, t = old n+1; linear order: s_k..s_1, v, t_1..t_k
    order = [2 * i - 1 for i in range(k, 0, -1)] + [0] + [2 * i for i in range(1, k + 1)]
    pos = {v: p for p, v in enumerate(order)}
    rev = [(a, b) for a in range(n) for b in range(n) if X[a][b] and pos[a] > pos[b]]
    expected = {(2 * i, 2 * i - 1) for i in range(1, k + 1)}   # t_i -> s_i
    assert set(rev) == expected, (k, rev)
print("(iii) Gamma^k(K1) = TT_(2k+1) with exactly the k concentric chords (s_i, t_i) reversed: k = 1..8 OK")
