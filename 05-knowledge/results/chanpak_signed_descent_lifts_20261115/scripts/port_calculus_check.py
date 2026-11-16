# Exact checks of the port calculus on random tournaments (brute-force Hamiltonian path counts).
#   For a tournament A and a vertex u:
#     E_A(u) = # Hamiltonian paths ending at u,  S_A(u) = # starting at u,
#     U_A(u)  = # ordered covers (P1, u P2) of A by two paths, P1 possibly empty  (no junction condition),
#     U*_A(u) = # ordered covers (P1 u, P2), P2 possibly empty,
#     D_A(u) = U_A(u) - H(A),  D*_A(u) = U*_A(u) - H(A).
#   (G1) glue  A u_v B  (identify u in A with w in B, A-u -> B-w):
#        H = H_A H_B + (H(A-u) - E_A)(H(B-w) - S_B) + D_A (H_B - H(B-w)) + D*_B (H_A - H(A-u)) + D_A D*_B
#   (G2) bridge (join A -> B with the single arc u -> w reversed):  H = H_A H_B - E_A S_B + U_A U*_B
#   (G3) tail vertex q (q -> u only, A - u -> q):                  H = 2 H_A + D_A - E_A
#   (G4) middle vertex x on the bridge (A -> x -> B):               H = H_A H_B + 2 U_A U*_B
#   (G5) D_A(u) >= H(A-u) - E_A(u)
import random, itertools, sys
def ham(adj, vs, start=None, end=None):
    vs = list(vs); n = len(vs)
    if n == 0: return 1
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        if start is None or vs[i] == start: dp[1 << i][i] = 1
    for S in range(1, 1 << n):
        for i in range(n):
            c = dp[S][i]
            if c:
                for j in range(n):
                    if not (S >> j) & 1 and adj[vs[i]][vs[j]]: dp[S | (1 << j)][j] += c
    full = (1 << n) - 1
    return sum(dp[full][i] for i in range(n) if end is None or vs[i] == end)
def rand_t(n, rng):
    adj = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < 0.5: adj[i][j] = 1
            else: adj[j][i] = 1
    return adj
def port(adj, n, u):
    vs = list(range(n)); H = ham(adj, vs); E = ham(adj, vs, end=u); S = ham(adj, vs, start=u)
    others = [v for v in vs if v != u]; U = 0; Us = 0
    for r in range(len(others) + 1):
        for P1 in itertools.combinations(others, r):
            P2 = [v for v in others if v not in P1]
            U += ham(adj, P1) * ham(adj, [u] + P2, start=u)       # (P1, u P2)
            Us += ham(adj, list(P1) + [u], end=u) * ham(adj, P2)   # (P1 u, P2)
    Hm = ham(adj, others)
    return dict(H=H, E=E, S=S, U=U, Us=Us, D=U - H, Ds=Us - H, Hm=Hm)
rng = random.Random(20261115); bad = 0; tests = 0
for trial in range(400):
    na = rng.randint(1, 5); nb = rng.randint(1, 5)
    A = rand_t(na, rng); B = rand_t(nb, rng); u = rng.randrange(na); w = rng.randrange(nb)
    pa = port(A, na, u); pb = port(B, nb, w)
    assert pa['D'] >= pa['Hm'] - pa['E'] >= 0 or na == 1, (pa,)
    # bridge: vertices A: 0..na-1, B: na..na+nb-1
    N = na + nb; T = [[0] * (N + 1) for _ in range(N + 1)]
    for i in range(na):
        for j in range(na): T[i][j] = A[i][j]
    for i in range(nb):
        for j in range(nb): T[na + i][na + j] = B[i][j]
    for i in range(na):
        for j in range(nb): T[i][na + j] = 1
    T[u][na + w] = 0; T[na + w][u] = 1
    hb = ham(T, range(N)); pred = pa['H'] * pb['H'] - pa['E'] * pb['S'] + pa['U'] * pb['Us']
    tests += 1; bad += (hb != pred)
    # middle vertex x = N: A -> x -> B
    for i in range(na): T[i][N] = 1
    for j in range(nb): T[N][na + j] = 1
    hm = ham(T, range(N + 1)); pred = pa['H'] * pb['H'] + 2 * pa['U'] * pb['Us']
    tests += 1; bad += (hm != pred)
    # tail vertex on A: q -> u only
    T2 = [[0] * (na + 1) for _ in range(na + 1)]
    for i in range(na):
        for j in range(na): T2[i][j] = A[i][j]
    for i in range(na):
        if i != u: T2[i][na] = 1
    T2[na][u] = 1
    ht = ham(T2, range(na + 1)); pred = 2 * pa['H'] + pa['D'] - pa['E']
    tests += 1; bad += (ht != pred)
    # glue: identify u (in A) with w (in B); A-u -> B-w
    Aoth = [i for i in range(na) if i != u]; Both = [j for j in range(nb) if j != w]
    idx = {('a', i): k for k, i in enumerate(Aoth)}; v = len(Aoth); idx.update({('b', j): v + 1 + k for k, j in enumerate(Both)})
    M = len(Aoth) + 1 + len(Both); G = [[0] * M for _ in range(M)]
    def ia(i): return v if i == u else idx[('a', i)]
    def ib(j): return v if j == w else idx[('b', j)]
    for i in range(na):
        for j in range(na):
            if A[i][j]: G[ia(i)][ia(j)] = 1
    for i in range(nb):
        for j in range(nb):
            if B[i][j]: G[ib(i)][ib(j)] = 1
    for i in Aoth:
        for j in Both: G[ia(i)][ib(j)] = 1
    hg = ham(G, range(M))
    pred = pa['H'] * pb['H'] + (pa['Hm'] - pa['E']) * (pb['Hm'] - pb['S']) + pa['D'] * (pb['H'] - pb['Hm']) + pb['Ds'] * (pa['H'] - pa['Hm']) + pa['D'] * pb['Ds']
    tests += 1; bad += (hg != pred)
print("identities tested:", tests, "failures:", bad)
