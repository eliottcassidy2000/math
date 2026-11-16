# Literal check of the near-floor family.  phi_0 = 1 (single vertex), phi_a = 2^a + 1 (fan U(1 0^a), apex port),
# dual fan U(1^b 0).  T(a,b): bridge of the fan and the dual fan (all arcs fan -> dual fan, except apex(dual) -> apex(fan)).
# T'(a,b) = T(a,b) + x with  fan -> x -> dual fan.   Claim: H(T) = 2 phi_a phi_b - 1,  H(T') = 3 phi_a phi_b,
# so 2 H(T') - 3 H(T) = 3: the least value an odd count >= 5 allows.
def ham(adj):
    n = len(adj); dp = [[0] * n for _ in range(1 << n)]
    for i in range(n): dp[1 << i][i] = 1
    for S in range(1, 1 << n):
        for i in range(n):
            c = dp[S][i]
            if c:
                for j in range(n):
                    if not (S >> j) & 1 and adj[i][j]: dp[S | (1 << j)][j] += c
    return sum(dp[(1 << n) - 1])
def one_apex(word):
    """vertices: 0 = apex, 1..m = chain (i -> j for i < j); apex -> chain vertex i iff word[i-1] == '1'."""
    m = len(word); n = m + 1; adj = [[0] * n for _ in range(n)]
    for i in range(1, n):
        for j in range(i + 1, n): adj[i][j] = 1
        if word[i - 1] == '1': adj[0][i] = 1
        else: adj[i][0] = 1
    return adj
def build(a, b, middle):
    A = one_apex('1' + '0' * a) if a else [[0]]
    B = one_apex('1' * b + '0') if b else [[0]]
    na, nb = len(A), len(B); n = na + nb + (1 if middle else 0); T = [[0] * n for _ in range(n)]
    for i in range(na):
        for j in range(na): T[i][j] = A[i][j]
    for i in range(nb):
        for j in range(nb): T[na + i][na + j] = B[i][j]
    for i in range(na):
        for j in range(nb): T[i][na + j] = 1
    T[0][na] = 0; T[na][0] = 1                      # reversed bridge between the two apexes (ports)
    if middle:
        x = n - 1
        for i in range(na): T[i][x] = 1
        for j in range(nb): T[x][na + j] = 1
    return T
ok = True
for a in range(0, 5):
    for b in range(0, 5):
        if a + b == 0: continue
        pa = 1 if a == 0 else 2 ** a + 1; pb = 1 if b == 0 else 2 ** b + 1; m = pa * pb
        h0 = ham(build(a, b, False)); h1 = ham(build(a, b, True))
        good = (h0 == 2 * m - 1 and h1 == 3 * m); ok &= good
        print("a=%d b=%d  m=%4d  H(T)=%5d (2m-1=%5d)  H(T')=%5d (3m=%5d)  n: %d -> %d -> %d  %s" % (a, b, m, h0, 2 * m - 1, h1, 3 * m, (m - 1) // 2, m - 1, (3 * m - 1) // 2, "ok" if good else "FAIL"))
print("ALL OK" if ok else "FAILURE")
