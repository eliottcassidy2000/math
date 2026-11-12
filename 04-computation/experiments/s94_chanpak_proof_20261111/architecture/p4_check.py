# Lane architecture (S94): exact brute-force check of PROOFS.txt P4 (one-exit extension in source coordinates).
# For every labelled tournament B' on n <= 5 vertices and every vertex x:  B = B' + y with y -> x and every other
# vertex of B' beating y.  Claim: H(B) = h - e + u, U_B(y) = h + u, End_B(y) = h - e, where h = H(B'),
# e = End_B'(x), u = U_B'(x) (ordered covers (P1, P2) by disjoint paths, P2 starting at x, P1 possibly empty).
import itertools, sys
def ham_paths(n, adj, verts):
    # all Hamiltonian paths (as tuples) of the subtournament induced on verts
    out = []
    for p in itertools.permutations(verts):
        if all(adj[p[i]][p[i + 1]] for i in range(len(p) - 1)):
            out.append(p)
    return out
def H(n, adj, verts):
    return len(ham_paths(n, adj, verts)) if verts else 1
def End(n, adj, verts, x):
    return sum(1 for p in ham_paths(n, adj, verts) if p[-1] == x)
def U(n, adj, verts, x):
    # ordered covers (P1, P2), P2 a path starting at x, P1 a possibly empty path on the remaining vertices
    tot = 0
    others = [v for v in verts if v != x]
    for k in range(len(others) + 1):
        for S in itertools.combinations(others, k):
            # P2 on {x} + S starting at x ; P1 on the rest (H(empty) = 1)
            p2 = sum(1 for p in ham_paths(n, adj, (x,) + S) if p[0] == x)
            if p2:
                rest = tuple(v for v in others if v not in S)
                tot += p2 * (H(n, adj, rest) if rest else 1)
    return tot
bad = 0; cnt = 0
for n in range(1, 6):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for bits in range(1 << len(pairs)):
        adj = [[0] * (n + 1) for _ in range(n + 1)]
        for k, (i, j) in enumerate(pairs):
            if bits >> k & 1: adj[i][j] = 1
            else: adj[j][i] = 1
        V = tuple(range(n))
        h = H(n, adj, V)
        for x in range(n):
            e = End(n, adj, V, x); u = U(n, adj, V, x)
            y = n
            for v in range(n):
                adj[v][y] = 0; adj[y][v] = 0
            for v in range(n):
                if v == x: adj[y][v] = 1
                else: adj[v][y] = 1
            W = V + (y,)
            HB = H(n + 1, adj, W); UB = U(n + 1, adj, W, y); EB = End(n + 1, adj, W, y)
            cnt += 1
            if not (HB == h - e + u and UB == h + u and EB == h - e and u >= h):
                bad += 1
                if bad < 5: print("FAIL", n, bits, x, (h, e, u), (HB, UB, EB))
print("rooted labelled tournaments tested (orders 1..5):", cnt, " failures:", bad)
print("CHECK", "PASS" if bad == 0 else "FAIL")
