"""Verification of the sole-exit tower law (PROOF.md).

For a rooted tournament (T, t) let T_k = T + x_1 + ... + x_k where x_1 -> t, x_j -> x_{j-1}, and every other vertex
beats x_j.  With h = h(T), e = e_t(T) (Hamiltonian paths ending at t) and
B = B_t(T) = sum_{A subset V - t} h(T[A]) * s_t(T[V - A])  (root response; h(T[empty]) = 1),
the claims are  h_1 = h - e + B,  h_2 = h + 2B,  h_k = h_{k-1} + h_{k-2} + h_{k-3} (k >= 3),
e_k = h_{k-1} - e_{k-1},  B_k = B_{k-1} + h_{k-1},  and the closed form h_k = h*A_k - e*C_k + B*D_k.
Checked exhaustively on every rooted tournament of order <= 4 (towers to height 5) and order 5 (height 3), and on 40 seeded
random rooted tournaments of orders 6-7 (height 3); the B-recursion is checked while the order stays <= 6 (its subset
sum is exponential in pure Python), the h- and e-recursions and the closed form at every height.  Pure Python subset DP; independent of the C census."""
import random, itertools

def ham_paths_ending(rows, verts):
    verts = list(verts); idx = {v: k for k, v in enumerate(verts)}; m = len(verts)
    if m == 0: return {}
    adj = [sum(1 << idx[u] for u in verts if rows[v] >> u & 1) for v in verts]
    full = (1 << m) - 1
    dp = [[0] * m for _ in range(1 << m)]
    for k in range(m): dp[1 << k][k] = 1
    for S in range(1, 1 << m):
        row = dp[S]
        for k in range(m):
            c = row[k]
            if not c: continue
            nxt = adj[k] & ~S
            while nxt:
                b = nxt & -nxt; nxt ^= b; j = b.bit_length() - 1
                dp[S | b][j] += c
    return {verts[k]: dp[full][k] for k in range(m)}

def ham_paths_starting(rows, verts):
    verts = list(verts); n = len(rows); rev = [0] * n
    for i in verts:
        for j in verts:
            if rows[j] >> i & 1: rev[i] |= 1 << j
    return ham_paths_ending(rev, verts)

def h_of(rows, verts):
    verts = list(verts)
    return 1 if not verts else sum(ham_paths_ending(rows, verts).values())

def B_of(rows, verts, t):
    others = [v for v in verts if v != t]; tot = 0
    for mask in range(1 << len(others)):
        A = [others[i] for i in range(len(others)) if mask >> i & 1]
        W = [v for v in verts if v not in A]
        tot += h_of(rows, A) * ham_paths_starting(rows, W)[t]
    return tot

def tower(rows, t0, K, withB=True):
    rows = list(rows); n = len(rows); target = t0
    hs = [h_of(rows, range(n))]; es = [ham_paths_ending(rows, range(n))[t0]]; Bs = [B_of(rows, range(n), t0)]
    for k in range(1, K + 1):
        x = n; rows.append(0)
        for v in range(n):
            if v == target: rows[x] |= 1 << v
            else: rows[v] |= 1 << x
        n += 1
        hs.append(h_of(rows, range(n))); es.append(ham_paths_ending(rows, range(n))[x])
        Bs.append(B_of(rows, range(n), x) if (withB and n <= 6) else None)
        target = x
    return hs, es, Bs

def trib(a0, a1, a2, K):
    s = [a0, a1, a2]
    while len(s) <= K: s.append(s[-1] + s[-2] + s[-3])
    return s

def check(rows, t0, K):
    hs, es, Bs = tower(rows, t0, K)
    h, e, B = hs[0], es[0], Bs[0]
    A, C, D = trib(1, 1, 1, K), trib(0, 1, 0, K), trib(0, 1, 2, K)
    ok = hs[1] == h - e + B and hs[2] == h + 2 * B
    ok &= all(hs[k] == hs[k-1] + hs[k-2] + hs[k-3] for k in range(3, K + 1))
    ok &= all(es[k] == hs[k-1] - es[k-1] for k in range(1, K + 1))
    ok &= all(Bs[k] == Bs[k-1] + hs[k-1] for k in range(1, K + 1) if Bs[k] is not None)
    ok &= all(hs[k] == h * A[k] - e * C[k] + B * D[k] for k in range(K + 1))
    return ok, (h, e, B, hs)

def all_tournaments(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for bits in itertools.product((0, 1), repeat=len(pairs)):
        rows = [0] * n
        for (i, j), b in zip(pairs, bits):
            if b: rows[i] |= 1 << j
            else: rows[j] |= 1 << i
        yield rows

total = bad = 0
for n in range(1, 6):
    for rows in all_tournaments(n):
        for t0 in range(n):
            ok, info = check(rows, t0, 5 if n <= 4 else 3); total += 1
            if not ok: bad += 1; print("FAIL", n, t0, rows, info)
print(f"exhaustive: every rooted tournament of order <= 4 (height 5) and order 5 (height 3): {total} towers, failures {bad}")
random.seed(20261104); rtotal = rbad = 0
for trial in range(40):
    n = random.choice((6, 7)); rows = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if random.random() < 0.5: rows[i] |= 1 << j
            else: rows[j] |= 1 << i
    t0 = random.randrange(n); ok, info = check(rows, t0, 3); rtotal += 1
    if not ok: rbad += 1; print("FAIL", n, t0, rows, info)
print(f"random: {rtotal} rooted tournaments of orders 6-7, towers to height 3: failures {rbad}")
# the tower over a single vertex is the full tiling tournament T*_n of THM-4731 (6): h_k = A000213(k+1)
hs, es, Bs = tower([0], 0, 11, withB=False)
print("tower over one vertex (T*_n, n = k+1):", hs, "= A000213(2..13)")
print("A_k =", trib(1,1,1,10), " C_k =", trib(0,1,0,10), " D_k =", trib(0,1,2,10))
