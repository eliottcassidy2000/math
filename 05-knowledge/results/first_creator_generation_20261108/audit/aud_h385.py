"""aud_h385.py - structure of the order-9 REST class with no pendant vertex (H = 385), auditor's code format
(aud_enum.c bit order: pair (0,1) is the most significant of C(n,2) bits; bit set = i -> j)."""
import itertools
def decode(code, n):
    A = [[0] * n for _ in range(n)]; b = n * (n - 1) // 2 - 1
    for i in range(n):
        for j in range(i + 1, n):
            if code >> b & 1: A[i][j] = 1
            else: A[j][i] = 1
            b -= 1
    return A
def H(A, V):
    V = list(V); idx = {v: k for k, v in enumerate(V)}; m = len(V)
    f = [[0] * m for _ in range(1 << m)]
    for k in range(m): f[1 << k][k] = 1
    for S in range(1, 1 << m):
        for k in range(m):
            if f[S][k]:
                for l in range(m):
                    if not S >> l & 1 and A[V[k]][V[l]]: f[S | 1 << l][l] += f[S][k]
    return sum(f[(1 << m) - 1])
def strong(A, V):
    V = list(V)
    for X in range(1, (1 << len(V)) - 1):
        Xs = [V[k] for k in range(len(V)) if X >> k & 1]; Ys = [v for v in V if v not in Xs]
        if all(A[x][y] for x in Xs for y in Ys): return False
    return True
n = 9; A = decode(43151540450, n)
print("H =", H(A, range(n)), " scores:", sorted(sum(r) for r in A))
for v in range(n):
    R = [u for u in range(n) if u != v]
    # ordinal sum P => Q of R with |P| = |Q| = 4
    for P in itertools.combinations(R, 4):
        Q = [u for u in R if u not in P]
        if all(A[p][q] for p in P for q in Q):
            sp, sq = strong(A, P), strong(A, Q)
            print("delete v=%d: T-v = P=>Q, P=%s (strong %s, H %d), Q=%s (strong %s, H %d); v beats %d of P and %d of Q; H(T-v)=%d"
                  % (v, P, sp, H(A, P), Q, sq, H(A, Q), sum(A[v][p] for p in P), sum(A[v][q] for q in Q), H(A, R)))
