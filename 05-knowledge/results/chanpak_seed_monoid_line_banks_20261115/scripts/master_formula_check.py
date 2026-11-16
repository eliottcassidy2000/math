# One head over  (upper tape) -> C -> (lower tape).  Master formula, checked against brute force.
#   upper tape: heights 1 (next to the block) .. i ; higher beats lower ; tape -> C -> lower tape (depth 1 next to block).
#   head x: contacts U (upper), X_C (block), L (lower); Z = lower non-contacts.
#   alpha = sum_{a in U, h >= 2} 2^(h-2), u1 = [height 1 in U];  beta = sum_{b in Z, depth >= 2} 2^(depth-2), z1 = [depth 1 in Z].
#   H = h(C) (1 + 2 N_up + 2 N_low) + 2 [ cut_C(X_C) + In + Out + Cross ],
#     In    = sum_{s in C - X_C} ( lam_odd[s] (alpha + u1) + lam_even[s] alpha ),
#     Out   = sum_{s in X_C}     ( kap_odd[s] (beta + z1)  + kap_even[s] beta ),
#     Cross = G_even ( (alpha+u1)(beta+z1) + alpha beta ) + G_odd ( (alpha+u1) beta + alpha (beta+z1) ),
#   lam_par[s] = sum over paths Q in C ending at s with |Q| of that parity of H(C - Q); kap likewise for paths starting at s;
#   G_par = sum over all paths Q in C (the empty path is even) of H(C - Q); N_up, N_low = one-apex half-counts of the tapes.
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
def block_stats(C, m):
    full = (1 << m) - 1
    h = [ham(C, [v for v in range(m) if S >> v & 1]) for S in range(1 << m)]
    W = [[[0, 0] for _ in range(m)] for _ in range(m)]          # W[a][b][parity of |Q|], paths a ~> b (a = b: single vertex)
    G = [h[full], 0]
    for s in range(m):
        stack = [(s, 1 << s)]
        while stack:
            v, S = stack.pop(); k = bin(S).count('1'); w = h[full ^ S]
            W[s][v][k % 2] += w; G[k % 2] += w
            for nx in range(m):
                if not S >> nx & 1 and C[v][nx]: stack.append((nx, S | (1 << nx)))
    lam = [[sum(W[t][s][p] for t in range(m)) for p in (0, 1)] for s in range(m)]   # [even, odd]
    kap = [[sum(W[s][t][p] for t in range(m)) for p in (0, 1)] for s in range(m)]
    return h, W, G, lam, kap
def one_apex_half(marks, order):
    """order: list of tape vertices from top to bottom; marks: set of contacts.  sum over a in marks, b not, a above b of c(gap)."""
    n = 0
    for i, a in enumerate(order):
        if a in marks:
            for j in range(i + 1, len(order)):
                if order[j] not in marks: g = j - i; n += 1 if g == 1 else 2 ** (g - 2)
    return n
def formula(C, m, i, j, U, XC, L):
    h, W, G, lam, kap = block_stats(C, m); full = (1 << m) - 1
    alpha = sum(2 ** (hh - 2) for hh in U if hh >= 2); u1 = 1 if 1 in U else 0
    Z = [d for d in range(1, j + 1) if d not in L]
    beta = sum(2 ** (d - 2) for d in Z if d >= 2); z1 = 1 if 1 in Z else 0
    Nup = one_apex_half(set(U), list(range(i, 0, -1))); Nlow = one_apex_half(set(L), list(range(1, j + 1)))
    cutC = sum(W[a][b][0] for a in range(m) for b in range(m) if a != b and (XC >> a & 1) and not (XC >> b & 1))
    In = sum(lam[s][1] * (alpha + u1) + lam[s][0] * alpha for s in range(m) if not XC >> s & 1)
    Out = sum(kap[s][1] * (beta + z1) + kap[s][0] * beta for s in range(m) if XC >> s & 1)
    Cross = G[0] * ((alpha + u1) * (beta + z1) + alpha * beta) + G[1] * ((alpha + u1) * beta + alpha * (beta + z1))
    return h[full] * (1 + 2 * Nup + 2 * Nlow) + 2 * (cutC + In + Out + Cross)
def brute(C, m, i, j, U, XC, L):
    # vertex order: upper heights i..1, block 0..m-1, lower depths 1..j, head
    names = [('u', hh) for hh in range(i, 0, -1)] + [('c', v) for v in range(m)] + [('l', d) for d in range(1, j + 1)] + [('x', 0)]
    n = len(names); adj = [[0] * n for _ in range(n)]
    def rank(nm): return (0, -nm[1]) if nm[0] == 'u' else ((1, 0) if nm[0] == 'c' else (2, nm[1]))
    for p in range(n - 1):
        for q in range(n - 1):
            if p == q: continue
            a, b = names[p], names[q]
            if a[0] == 'c' and b[0] == 'c': adj[p][q] = C[a[1]][b[1]]
            else: adj[p][q] = 1 if rank(a) < rank(b) else 0
    for p in range(n - 1):
        a = names[p]
        con = (a[0] == 'u' and a[1] in U) or (a[0] == 'c' and (XC >> a[1] & 1)) or (a[0] == 'l' and a[1] in L)
        if con: adj[n - 1][p] = 1
        else: adj[p][n - 1] = 1
    return ham(adj, range(n))
if __name__ == '__main__':
    rng = random.Random(11); bad = 0; tests = 0
    for trial in range(60):
        m = rng.randint(1, 4); C = [[0] * m for _ in range(m)]
        for a in range(m):
            for b in range(a + 1, m):
                if rng.random() < 0.5: C[a][b] = 1
                else: C[b][a] = 1
        i = rng.randint(0, 3); j = rng.randint(0, 3)
        U = {hh for hh in range(1, i + 1) if rng.random() < 0.5}; L = {d for d in range(1, j + 1) if rng.random() < 0.5}
        XC = rng.randrange(1 << m)
        f = formula(C, m, i, j, U, XC, L); b = brute(C, m, i, j, U, XC, L); tests += 1
        if f != b: bad += 1; print("MISMATCH", m, i, j, U, XC, L, f, b)
    print("master formula: %d random instances, %d mismatches" % (tests, bad))
