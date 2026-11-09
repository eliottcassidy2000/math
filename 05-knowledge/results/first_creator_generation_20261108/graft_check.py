"""Check H(A => B with the arc a->b reversed) = H(A)H(B) - F_A(a) G_B(b) + U_A(a) U'_B(b) on random small cases."""
import random, itertools
def hp_counts(adj, verts):
    """return dict: (S frozenset) -> {end: count} of Hamiltonian paths of the induced subtournament on S"""
    verts = list(verts); idx = {v: i for i, v in enumerate(verts)}; m = len(verts)
    f = {}
    for i in range(m): f[(1 << i, i)] = 1
    for S in range(1, 1 << m):
        for i in range(m):
            x = f.get((S, i), 0)
            if not x: continue
            for j in range(m):
                if not S >> j & 1 and adj[verts[i]][verts[j]]:
                    f[(S | 1 << j, j)] = f.get((S | 1 << j, j), 0) + x
    return f, verts
def H(adj, verts):
    f, vs = hp_counts(adj, verts); m = len(vs)
    return sum(f.get(((1 << m) - 1, i), 0) for i in range(m)) if m else 1
def paths_end(adj, S, e):   # Ham paths of T[S] ending at e
    S = list(S)
    if not S: return 0
    f, vs = hp_counts(adj, S); m = len(vs)
    return f.get(((1 << m) - 1, vs.index(e)), 0)
def paths_start(adj, S, s):
    S = list(S)
    if not S: return 0
    radj = [[adj[j][i] for j in range(len(adj))] for i in range(len(adj))]
    return paths_end(radj, S, s)
def rand_tour(vs, adj):
    for i, j in itertools.combinations(vs, 2):
        if random.random() < 0.5: adj[i][j] = 1; adj[j][i] = 0
        else: adj[j][i] = 1; adj[i][j] = 0
random.seed(7); bad = 0
for trial in range(300):
    na = random.randint(1, 4); nb = random.randint(1, 4); n = na + nb
    adj = [[0]*n for _ in range(n)]
    A = list(range(na)); B = list(range(na, n))
    rand_tour(A, adj); rand_tour(B, adj)
    for x in A:
        for y in B: adj[x][y] = 1; adj[y][x] = 0
    a = random.choice(A); b = random.choice(B)
    adj[a][b] = 0; adj[b][a] = 1
    direct = H(adj, range(n))
    HA = H(adj, A); HB = H(adj, B)
    FA = paths_end(adj, A, a); GB = paths_start(adj, B, b)
    # U_A(a): ordered (P1, P2) covering A, P2 nonempty starting at a, P1 possibly empty
    UA = 0
    others = [v for v in A if v != a]
    for r in range(len(others) + 1):
        for P1 in itertools.combinations(others, r):
            A2 = [v for v in A if v not in P1]
            UA += (H(adj, P1) if P1 else 1) * paths_start(adj, A2, a)
    UB = 0
    othersb = [v for v in B if v != b]
    for r in range(len(othersb) + 1):
        for Q2 in itertools.combinations(othersb, r):
            B1 = [v for v in B if v not in Q2]
            UB += paths_end(adj, B1, b) * (H(adj, Q2) if Q2 else 1)
    pred = HA * HB - FA * GB + UA * UB
    if pred != direct: bad += 1
print('trials 300, mismatches with H(A)H(B) - F_A(a)G_B(b) + U_A(a)U_B(b):', bad)
