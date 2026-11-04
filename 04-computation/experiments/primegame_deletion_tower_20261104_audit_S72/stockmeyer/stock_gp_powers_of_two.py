# Reconstruct Stockmeyer's M_p, M*_p (arXiv:2110.06158, Def. 2 / Fig. 1) and test G_p vs G*_p.
import itertools, sys
import networkx as nx
from networkx.algorithms.isomorphism import DiGraphMatcher

M4  = [[0,1,2,3],[-1,0,3,-2],[-2,-3,0,1],[-3,2,-1,0]]
M4s = [[0,-2,-3,-1],[2,0,1,-3],[3,-1,0,2],[1,3,-2,0]]

def classify(d):
    # returns ('0',) or ('+1',) / ('-1',) for d = +-1 mod 4 odd, or ('x', x, ysign) for d = y*2^x
    if d == 0: return ('zero',)
    if d % 2 != 0:
        return ('odd', 1 if d % 4 == 1 else -1)
    x = 0; y = d
    while y % 2 == 0:
        y //= 2; x += 1
    return ('even', x, 1 if y % 4 == 1 else -1)

def build(p, star):
    base = M4s if star else M4
    nb = p // 4
    M = [[0]*p for _ in range(p)]
    for bi in range(nb):
        for bj in range(nb):
            c = classify(bj - bi)
            for a in range(4):
                for b in range(4):
                    v = base[a][b]; I = 1 if a == b else 0
                    if c[0] == 'zero': w = v
                    elif c[0] == 'odd':
                        s = c[1]
                        # M: -M4 + 4I if +1, -M4 - 4I if -1 ; M*: -M4* - 4I if +1, -M4* + 4I if -1
                        w = -v + (4*I if s == 1 else -4*I) if not star else -v + (-4*I if s == 1 else 4*I)
                    else:
                        x, s = c[1], c[2]
                        if not star: w = v + ((x+4)*I if s == 1 else -(x+4)*I)
                        else:        w = v + (-(x+4)*I if s == 1 else (x+4)*I)
                    M[4*bi+a][4*bj+b] = w
    return M

def tour(M):
    p = len(M); G = nx.DiGraph(); G.add_nodes_from(range(p))
    for i in range(p):
        for j in range(p):
            if i != j and M[i][j] > 0: G.add_edge(i, j)
            if i != j: assert M[i][j] == -M[j][i] and M[i][j] != 0
    return G

def iso(A, B):
    return DiGraphMatcher(A, B).is_isomorphic()

def deck_equal(A, B):
    ca = [A.subgraph([v for v in A if v != x]).copy() for x in A]
    cb = [B.subgraph([v for v in B if v != x]).copy() for x in B]
    used = [False]*len(cb)
    for c in ca:
        found = False
        for k, d in enumerate(cb):
            if not used[k] and sorted(dict(c.out_degree()).values()) == sorted(dict(d.out_degree()).values()) and iso(c, d):
                used[k] = True; found = True; break
        if not found: return False
    return True

def ham_paths_cycles(G):
    n = G.number_of_nodes(); adj = [0]*n
    for u, v in G.edges(): adj[u] |= 1 << v
    # paths: dp[mask][v] = # paths covering mask ending at v
    full = (1 << n) - 1
    dp = [dict() for _ in range(1 << n)] if n <= 12 else None
    # use list-of-lists for speed
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        row = dp[mask]
        for v in range(n):
            c = row[v]
            if c:
                nb = adj[v] & ~mask
                while nb:
                    lb = nb & -nb; w = lb.bit_length() - 1
                    dp[mask | lb][w] += c
                    nb ^= lb
    H = sum(dp[full])
    # cycles: paths starting at vertex 0
    dp2 = [[0]*n for _ in range(1 << n)]
    dp2[1][0] = 1
    for mask in range(1, 1 << n, 2):
        row = dp2[mask]
        for v in range(n):
            c = row[v]
            if c:
                nb = adj[v] & ~mask
                while nb:
                    lb = nb & -nb; w = lb.bit_length() - 1
                    dp2[mask | lb][w] += c
                    nb ^= lb
    hc = sum(dp2[full][v] for v in range(n) if (adj[v] >> 0) & 1)
    return H, hc

for p in ([int(a) for a in sys.argv[1:]] if __name__ == '__main__' else []):
    M, Ms = build(p, False), build(p, True)
    G, Gs = tour(M), tour(Ms)
    Gr = G.reverse(copy=True)
    print(f"p={p}: scores G={sorted(d for _,d in G.out_degree())} G*={sorted(d for _,d in Gs.out_degree())}")
    print(f"  G~G*: {iso(G,Gs)}  G*~converse(G): {iso(Gs,Gr)}  G self-converse: {iso(G,Gr)}  G* self-converse: {iso(Gs, Gs.reverse(copy=True))}")
    print(f"  same deck: {deck_equal(G, Gs)}")
    if p <= 16:
        print("  (H, hc) G:", ham_paths_cycles(G), " G*:", ham_paths_cycles(Gs))
    if p in (8, 16):
        print("  row1 M_p:", M[0]); print("  row1 M*_p:", Ms[0])

if __name__ == "__main__" and len(sys.argv) == 1:
    pass
