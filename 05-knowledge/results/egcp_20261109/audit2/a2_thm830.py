"""AUDIT2: THM-830 (1.1)/(1.3) conventions vs THM-5041 (1d).
Tiling t on S_n = {(a,b): 1 <= b < a <= n, a - b >= 2}; G_t on gaps [N], N = n - 1, with (G_t)_(i,j) = t_(i+1, j)
for j < i (matrix entry: arc i -> j iff 1).  Staircase reflection sigma: (a, b) -> (n+1-b, n+1-a).
Checks for all tilings at n = 4, 5, 6: (1.3) Gamma(kappa t) = G_t^op and Gamma(sigma t) = rho(G_t^op);
then THM-5041 (1d) for the tournament A = G_t: U(A), L(A) as arrays on gap pairs {j < i}; the reversal relabelling
rho(A) has (U, L) = (s L, s U) with s = sigma transported to gap pairs, i.e. (sX)_{j,i} = X_{N+1-i, N+1-j}."""
import itertools
def check(n):
    S = [(a, b) for a in range(1, n + 1) for b in range(1, a) if a - b >= 2]
    N = n - 1; bad = 0
    for bits in range(2 ** len(S)):
        t = {s: bits >> k & 1 for k, s in enumerate(S)}
        def G(t):
            A = [[0] * (N + 1) for _ in range(N + 1)]
            for i in range(1, N + 1):
                for j in range(1, i):
                    if t[(i + 1, j)]: A[i][j] = 1
                    else: A[j][i] = 1
            return A
        A = G(t)
        op = [[A[j][i] for j in range(N + 1)] for i in range(N + 1)]
        rho = lambda M: [[M[N + 1 - i][N + 1 - j] if i and j else 0 for j in range(N + 1)] for i in range(N + 1)]
        kt = {s: 1 - v for s, v in t.items()}
        st = {(a, b): t[(n + 1 - b, n + 1 - a)] for (a, b) in S}
        bad += G(kt) != op
        bad += G(st) != rho(op)
        # THM-5041 (1d) on gap coordinates
        U = {(j, i): A[j][i] for i in range(1, N + 1) for j in range(1, i)}
        L = {(j, i): A[i][j] for i in range(1, N + 1) for j in range(1, i)}
        s = lambda X: {(j, i): X[(N + 1 - i, N + 1 - j)] for (j, i) in X}
        R = rho(A)
        Ur = {(j, i): R[j][i] for i in range(1, N + 1) for j in range(1, i)}
        Lr = {(j, i): R[i][j] for i in range(1, N + 1) for j in range(1, i)}
        bad += not (Ur == s(L) and Lr == s(U))
        # under THM-830 the tournament A = G_t has L(A) = t (tile (i+1, j) <-> gap pair (j, i))
        bad += any(L[(j, i)] != t[(i + 1, j)] for (j, i) in L)
    return 2 ** len(S), bad
for n in (4, 5, 6):
    print("n=%d tilings %d failures %d" % ((n,) + check(n)))
