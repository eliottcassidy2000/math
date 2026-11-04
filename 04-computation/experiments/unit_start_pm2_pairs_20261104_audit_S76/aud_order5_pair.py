"""aud_order5_pair.py - literal check of the order-5 Session189-wrapper +-2 pair (no unit-start condition).

Finds labelled representatives (root 0) of the strong rooted profiles (h, s, F) = (13, 3, 15) and (15, 3, 15) at
order 5, builds the literal wrappers D_m (TT_m block Z and root w: w -> Z -> U, r -> w -> U - r) for m = 1..5, and
counts their Hamiltonian paths directly.  Also reports the arc-reversal distance between the two rooted classes
through strong labelled tournaments (root fixed).
"""
from collections import deque


def ham_count(n, adj):
    dp = [[0] * n for _ in range(1 << n)]
    for x in range(n):
        dp[1 << x][x] = 1
    for S in range(1, 1 << n):
        for x in range(n):
            c = dp[S][x]
            if not c:
                continue
            for y in range(n):
                if not (S >> y) & 1 and adj[x][y]:
                    dp[S | (1 << y)][y] += c
    return sum(dp[(1 << n) - 1])


def profile(n, out):
    full = (1 << n) - 1
    end = [[0] * n for _ in range(1 << n)]
    beg0 = [[0] * n for _ in range(1 << n)]
    for x in range(n):
        end[1 << x][x] = 1
    beg0[1][0] = 1
    for S in range(1, full + 1):
        for x in range(n):
            a, b = end[S][x], beg0[S][x]
            if not a and not b:
                continue
            for y in range(n):
                if not (S >> y) & 1 and (out[x] >> y) & 1:
                    end[S | (1 << y)][y] += a
                    beg0[S | (1 << y)][y] += b
    HS = [1] + [sum(end[S]) for S in range(1, full + 1)]
    return HS[full], sum(beg0[full]), sum(end[S][0] * HS[full ^ S] for S in range(1, full + 1, 2))


def strong(n, out):
    full = (1 << n) - 1
    f = 1
    for _ in range(n):
        for x in range(n):
            if (f >> x) & 1:
                f |= out[x]
    b = 1
    for _ in range(n):
        for y in range(n):
            if out[y] & b:
                b |= 1 << y
    return f == full and b == full


n = 5
pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]


def decode(code):
    out = [0] * n
    for k, (i, j) in enumerate(pairs):
        if (code >> k) & 1:
            out[i] |= 1 << j
        else:
            out[j] |= 1 << i
    return out


prof = {}
for code in range(1 << len(pairs)):
    out = decode(code)
    if strong(n, out):
        prof[code] = profile(n, out)
repA = next(c for c, p in prof.items() if p == (13, 3, 15))
repB = next(c for c, p in prof.items() if p == (15, 3, 15))


def wrapperH(out, m):
    N = n + m + 1
    adj = [[0] * N for _ in range(N)]
    for a in range(n):
        for b in range(n):
            adj[a][b] = (out[a] >> b) & 1
    w = n + m
    for i in range(m):
        z = n + i
        for v in range(n):
            adj[z][v] = 1
        for j in range(i + 1, m):
            adj[z][n + j] = 1
        adj[w][z] = 1
    for v in range(1, n):
        adj[w][v] = 1
    adj[0][w] = 1
    return ham_count(N, adj)


for nm, c in (("A5", repA), ("B5", repB)):
    out = decode(c)
    print(nm, "out-rows (root 0):", out, "profile (h,s,F):", prof[c], "H(D_m), m=1..5:", [wrapperH(out, m) for m in range(1, 6)])

# distance between the two rooted classes through strong labelled tournaments, root fixed
targets = {c for c, p in prof.items() if p == (15, 3, 15)}
dist = {repA: 0}
dq = deque([repA])
D = None
while dq:
    t = dq.popleft()
    if D is not None and dist[t] >= D:
        continue
    for k in range(len(pairs)):
        u = t ^ (1 << k)
        if u in dist or u not in prof:
            continue
        dist[u] = dist[t] + 1
        if u in targets and D is None:
            D = dist[u]
        dq.append(u)
print("strong labelled tournaments with profile (15,3,15):", len(targets), "; distance from A5 to that profile set:", D)
