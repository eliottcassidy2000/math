"""Wave 44: random check of the hub-seam identity.
T' any tournament, A any vertex.  T = triangle x->y->z->x placed over T' (triangle beats all of T') except A -> x.
T'_+ = T' plus a vertex b with A -> b and b beating every other vertex of T'.
Claim (proved in the README):  H(T) = 3 H(T'_+) + 2 Psi,  Psi = sum over Q of (-1)^|Q| H(T' - A - Q),
Q = the empty path or a nonempty path of T' - A ending at an in-neighbour of A.
Corollary: if T' - A has a source s, Psi = [A -> s] H(T' - A)."""
import random, itertools, sys

def ham(n, adj):
    if n == 0:
        return 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        row = dp[S]
        for v in range(n):
            c = row[v]
            if not c:
                continue
            for u in range(n):
                if not (S >> u) & 1 and adj[v][u]:
                    dp[S | 1 << u][u] += c
    return sum(dp[(1 << n) - 1])

def sub(adj, keep):
    return len(keep), [[adj[a][b] for b in keep] for a in keep]

def rand_tour(n, rng):
    adj = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < 0.5:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
    return adj

def paths_ending(adj, verts, targets):
    """all nonempty simple paths (as vertex sets with orientation irrelevant for the remainder) in the subtournament
    on verts that end at a vertex of targets; returns list of (length, vertex set)"""
    out = []
    def ext(path, used):
        if path[-1] in targets:
            out.append((len(path), frozenset(path)))
        for w in verts:
            if w not in used and adj[path[-1]][w]:
                ext(path + [w], used | {w})
    for v in verts:
        ext([v], {v})
    return out

def check(n, rng, force_source):
    adj = rand_tour(n, rng)
    A = 0
    rest = list(range(1, n))
    if force_source:          # make vertex 1 a source of T' - A, random arc to A
        s = 1
        for w in rest:
            if w != s:
                adj[s][w], adj[w][s] = 1, 0
    # T: triangle on n, n+1, n+2
    N = n + 3
    T = [[0] * N for _ in range(N)]
    for i in range(n):
        for j in range(n):
            T[i][j] = adj[i][j]
    x, y, z = n, n + 1, n + 2
    T[x][y] = T[y][z] = T[z][x] = 1
    for t in (x, y, z):
        for i in range(n):
            T[t][i] = 1
    T[x][A] = 0; T[A][x] = 1
    H_T = ham(N, T)
    # T'_+
    P = [[0] * (n + 1) for _ in range(n + 1)]
    for i in range(n):
        for j in range(n):
            P[i][j] = adj[i][j]
    b = n
    P[A][b] = 1
    for w in rest:
        P[b][w] = 1
    H_P = ham(n + 1, P)
    # Psi
    inA = {w for w in rest if adj[w][A]}
    psi = ham(*sub(adj, rest))
    for (k, Q) in paths_ending(adj, rest, inA):
        keep = [w for w in rest if w not in Q]
        psi += (-1) ** k * ham(*sub(adj, keep))
    ok = (H_T == 3 * H_P + 2 * psi)
    cor = None
    if force_source:
        cor = (psi == (ham(*sub(adj, rest)) if adj[A][1] else 0))
    return ok, cor, H_T, H_P, psi

rng = random.Random(44)
bad = 0; badc = 0; tests = 0
for trial in range(300):
    n = rng.randint(2, 7)
    fs = trial % 2 == 1
    ok, cor, a, b_, psi = check(n, rng, fs)
    tests += 1
    if not ok:
        bad += 1
        print("IDENTITY FAILS", n, a, b_, psi)
    if fs and not cor:
        badc += 1
        print("COROLLARY FAILS", n, a, b_, psi)
print(f"hub-seam identity: {tests} random tournaments (orders 2..7 plus the triangle), {bad} failures; source corollary failures: {badc}")
