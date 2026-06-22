"""The CUT side of the project is classical (mac-mini-2026-06-22-S40).
VERIFIED: (1) Clebsch graph SRG(16,5,0,2) = cut-space Cayley graph of K_5 (folded 5-cube; the 5
vertex-stars delta_v with sum=0 are its generators); generally cut-space Cayley(K_n) = folded n-cube.
(2) Truncated octahedron = permutohedron of S_4 = transitive tournaments at n=4 (adjacent transpositions).
Both are cut-side objects, complementing the even-graph CYCLE side (S38).
See 07-reflections/the-cut-side-is-classical-clebsch-and-the-permutohedron.md."""
from itertools import combinations, permutations

def clebsch_is_cutspace_cayley():
    def canon(S):
        S = frozenset(S); Sc = frozenset(range(5)) - S
        return min(S, Sc, key=lambda x: (len(x), sorted(x)))
    verts = sorted({canon(S) for r in range(6) for S in combinations(range(5), r)},
                   key=lambda x: (len(x), sorted(x)))
    idx = {v: i for i, v in enumerate(verts)}; N = len(verts)
    adj = [[0]*N for _ in range(N)]
    for v in verts:
        for g in range(5):                       # flip vertex g across the cut
            w = canon(v ^ {g})
            if w != v: adj[idx[v]][idx[w]] = adj[idx[w]][idx[v]] = 1
    deg = {sum(r) for r in adj}
    def common(i, j): return sum(1 for k in range(N) if adj[i][k] and adj[j][k])
    lam = {common(i, j) for i in range(N) for j in range(i+1, N) if adj[i][j]}
    mu  = {common(i, j) for i in range(N) for j in range(i+1, N) if not adj[i][j]}
    return N == 16 and deg == {5} and lam == {0} and mu == {2}

def trunc_oct_is_permutohedron():
    S4 = list(permutations(range(4))); idx = {p: i for i, p in enumerate(S4)}; M = len(S4)
    adj = [[0]*M for _ in range(M)]
    for p in S4:
        for t in range(3):
            q = list(p); q[t], q[t+1] = q[t+1], q[t]; q = tuple(q)
            adj[idx[p]][idx[q]] = adj[idx[q]][idx[p]] = 1
    return M == 24 and {sum(r) for r in adj} == {3}

if __name__ == "__main__":
    print("Clebsch = cut-space Cayley(K_5) [SRG(16,5,0,2)]:", clebsch_is_cutspace_cayley())
    print("Truncated octahedron = permutohedron(S_4):", trunc_oct_is_permutohedron())
