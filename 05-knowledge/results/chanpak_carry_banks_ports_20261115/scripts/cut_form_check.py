# The cut form.  For a tournament B and a new vertex x with out-set X:
#     H(B + x) = H(B) + 2 * sum_{a in X, b not in X} W(a, b),
#     W(a, b) = sum over directed paths P from a to b in B with an EVEN number of vertices of H(B - V(P))
# (Claim A: H(T) - H(T - v) = 2 sum over odd cycles C through v of H(T - V(C)); a cycle through x is x -> a ~> b -> x).
# Consequences checked: second differences are -2 (W(e,f) + W(f,e)) and W(e,f) + W(f,e) >= H(B - e - f) >= 1.
import random, itertools
def ham(adj, vs):
    vs = list(vs); n = len(vs)
    if n == 0: return 1
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n): dp[1 << i][i] = 1
    for S in range(1, 1 << n):
        for i in range(n):
            c = dp[S][i]
            if c:
                for j in range(n):
                    if not (S >> j) & 1 and adj[vs[i]][vs[j]]: dp[S | (1 << j)][j] += c
    return sum(dp[(1 << n) - 1])
def Wmat(adj, n):
    W = [[0] * n for _ in range(n)]
    def ext(path, used):
        a, last = path[0], path[-1]
        if len(path) % 2 == 0: W[a][last] += ham(adj, [v for v in range(n) if v not in used])
        for v in range(n):
            if v not in used and adj[last][v]: ext(path + [v], used | {v})
    for a in range(n): ext([a], {a})
    return W
rng = random.Random(7); tests = 0; bad = 0; minratio = None
for trial in range(120):
    n = rng.randint(2, 6); adj = [[0] * (n + 1) for _ in range(n + 1)]
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < 0.5: adj[i][j] = 1
            else: adj[j][i] = 1
    W = Wmat(adj, n); hB = ham(adj, range(n))
    for X in range(1 << n):
        for v in range(n): adj[n][v] = 1 if X >> v & 1 else 0; adj[v][n] = 0 if X >> v & 1 else 1
        h = ham(adj, range(n + 1))
        pred = hB + 2 * sum(W[a][b] for a in range(n) for b in range(n) if (X >> a & 1) and not (X >> b & 1))
        tests += 1; bad += (h != pred)
    for e in range(n):
        for f in range(e + 1, n):
            low = ham(adj, [v for v in range(n) if v not in (e, f)])
            assert W[e][f] + W[f][e] >= low >= 1
print("cut form verified on %d (base, out-set) pairs; failures: %d; pair interaction >= H(B - e - f) on every pair" % (tests, bad))
