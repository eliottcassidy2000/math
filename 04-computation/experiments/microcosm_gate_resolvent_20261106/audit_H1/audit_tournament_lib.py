# audit_H1 helper library (independent of the packet's scripts): tournaments, gate/padding/double/join constructions, marked pair via the determinant lemma, exact rational solves, path/cycle/cover counts.
import random
from fractions import Fraction
import subprocess

GENTOURNG = "/opt/homebrew/bin/gentourng"


def rand_tour(n, rng):
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


def from_upper(n, st):
    """gentourng default ascii: upper triangle row by row, '1' at (i,j), i<j, means i->j."""
    A = [[0] * n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if st[k] == '1':
                A[i][j] = 1
            else:
                A[j][i] = 1
            k += 1
    return A


def classes(n, strong=False):
    if n == 1:
        return [[[0]]]
    args = [GENTOURNG, "-q"] + (["-c"] if strong else []) + [str(n)]
    out = subprocess.run(args, capture_output=True, text=True, check=True).stdout.split()
    return [from_upper(n, st) for st in out]


def transitive(n):
    return [[1 if i < j else 0 for j in range(n)] for i in range(n)]


def dominance(A):
    n = len(A)
    return [[A[i][j] - A[j][i] for j in range(n)] for i in range(n)]


def scores(A):
    return [sum(r) for r in A]


def relabel(A, perm):
    """perm[i] = new label of old vertex i."""
    n = len(A)
    B = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            B[perm[i]][perm[j]] = A[i][j]
    return B


def _embed(A, extra):
    n = len(A)
    N = n + extra
    B = [[0] * N for _ in range(N)]
    for i in range(n):
        for j in range(n):
            B[i][j] = A[i][j]
    return B


def gate(A, rng=None):
    """Gamma(T): new s beats every core vertex, every core vertex beats new t, and t -> s.
    s = index n, t = index n+1 before an optional random relabelling."""
    n = len(A)
    B = _embed(A, 2)
    s, t = n, n + 1
    for v in range(n):
        B[s][v] = 1
        B[v][t] = 1
    B[t][s] = 1
    if rng is not None:
        perm = list(range(n + 2))
        rng.shuffle(perm)
        B = relabel(B, perm)
    return B


def pad(A, rng=None):
    """E(T): s -> core -> t, and s -> t."""
    n = len(A)
    B = _embed(A, 2)
    s, t = n, n + 1
    for v in range(n):
        B[s][v] = 1
        B[v][t] = 1
    B[s][t] = 1
    if rng is not None:
        perm = list(range(n + 2))
        rng.shuffle(perm)
        B = relabel(B, perm)
    return B


def double(A):
    """Skew-Sylvester double with dominance [[M, M+I], [M-I, -M]] (THM-447)."""
    n = len(A)
    M = dominance(A)
    D = [[0] * (2 * n) for _ in range(2 * n)]
    for i in range(n):
        for j in range(n):
            I = 1 if i == j else 0
            D[i][j] = M[i][j]
            D[i][n + j] = M[i][j] + I
            D[n + i][j] = M[i][j] - I
            D[n + i][n + j] = -M[i][j]
    B = [[1 if D[i][j] == 1 else 0 for j in range(2 * n)] for i in range(2 * n)]
    # sanity: a tournament
    for i in range(2 * n):
        assert B[i][i] == 0
        for j in range(i + 1, 2 * n):
            assert B[i][j] + B[j][i] == 1
    return B


def join(A, B):
    """A => B: every vertex of A beats every vertex of B."""
    a, b = len(A), len(B)
    C = [[0] * (a + b) for _ in range(a + b)]
    for i in range(a):
        for j in range(a):
            C[i][j] = A[i][j]
        for j in range(b):
            C[i][a + j] = 1
    for i in range(b):
        for j in range(b):
            C[a + i][a + j] = B[i][j]
    return C


# ---------- exact rational linear algebra ----------

def solve_frac(Mat, rhs):
    """Solve Mat x = rhs exactly over Q (Mat: list of lists of ints/Fractions)."""
    n = len(Mat)
    a = [[Fraction(x) for x in row] + [Fraction(rhs[i])] for i, row in enumerate(Mat)]
    for c in range(n):
        p = next(r for r in range(c, n) if a[r][c] != 0)
        a[c], a[p] = a[p], a[c]
        piv = a[c][c]
        rowc = a[c]
        for r in range(n):
            if r != c and a[r][c] != 0:
                f = a[r][c] / piv
                rr = a[r]
                for k in range(c, n + 1):
                    rr[k] -= f * rowc[k]
    return [a[i][n] / a[i][i] for i in range(n)]


def det_frac(Mat):
    n = len(Mat)
    a = [[Fraction(x) for x in row] for row in Mat]
    d = Fraction(1)
    for c in range(n):
        p = next((r for r in range(c, n) if a[r][c] != 0), None)
        if p is None:
            return Fraction(0)
        if p != c:
            a[c], a[p] = a[p], a[c]
            d = -d
        piv = a[c][c]
        d *= piv
        for r in range(c + 1, n):
            if a[r][c] != 0:
                f = a[r][c] / piv
                for k in range(c, n):
                    a[r][k] -= f * a[c][k]
    return d


def F_at(A, zval):
    """F_T(z) = u^T (zI - M)^{-1} u, exact for rational z."""
    M = dominance(A)
    n = len(A)
    Mat = [[(zval if i == j else 0) - M[i][j] for j in range(n)] for i in range(n)]
    x = solve_frac(Mat, [1] * n)
    return sum(x)


def phi(A):
    return F_at(A, Fraction(1))


def phi_gram(A):
    """u^T (I + M^T M)^{-1} u, computed independently."""
    M = dominance(A)
    n = len(A)
    G = [[(1 if i == j else 0) + sum(M[k][i] * M[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    x = solve_frac(G, [1] * n)
    return sum(x)


# ---------- path / cycle counts ----------

def ham_paths_all_subsets(A):
    """h[S] = number of Hamiltonian paths of T[S] (h[0] = 1)."""
    n = len(A)
    out = [[0] * n for _ in range(n)]
    for i in range(n):
        out[i] = sum(1 << j for j in range(n) if A[i][j])
    # dp[S][v] = paths covering S ending at v
    size = 1 << n
    dp = [[0] * n for _ in range(size)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1, size):
        row = dp[S]
        for v in range(n):
            c = row[v]
            if c == 0:
                continue
            nxt = out[v] & ~S
            while nxt:
                w = (nxt & -nxt).bit_length() - 1
                nxt &= nxt - 1
                dp[S | (1 << w)][w] += c
    h = [sum(dp[S]) for S in range(size)]
    h[0] = 1
    return h


def ham_paths(A):
    n = len(A)
    if n == 0:
        return 1
    return ham_paths_all_subsets(A)[(1 << n) - 1]


def ham_cycles(A):
    """Directed Hamiltonian cycles modulo rotation (start fixed at vertex 0)."""
    n = len(A)
    if n < 3:
        return 0
    size = 1 << n
    dp = [[0] * n for _ in range(size)]
    dp[1][0] = 1
    for S in range(1, size, 2):
        for v in range(n):
            c = dp[S][v]
            if c == 0:
                continue
            for w in range(n):
                if not (S >> w) & 1 and A[v][w]:
                    dp[S | (1 << w)][w] += c
    full = size - 1
    return sum(dp[full][v] for v in range(1, n) if A[v][0])


def c3(A):
    n = len(A)
    cnt = 0
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                if (A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[k][j] and A[i][k]):
                    cnt += 1
    return cnt


def cover_poly_values(A, qmax):
    """P_T(q) for q = 0..qmax by the definition: ordered q-tuples of (possibly empty) directed paths
    partitioning V, computed as the q-fold subset convolution of h (h(empty) = 1)."""
    n = len(A)
    h = ham_paths_all_subsets(A)
    full = (1 << n) - 1
    size = 1 << n
    vals = []
    cur = [0] * size
    cur[0] = 1  # 0-fold convolution: delta at empty set
    vals.append(cur[full] if n > 0 else 1)
    for q in range(1, qmax + 1):
        new = [0] * size
        for S in range(size):
            # sum over subsets B of S: cur[S \ B] * h[B]
            B = S
            acc = 0
            while True:
                acc += cur[S ^ B] * h[B]
                if B == 0:
                    break
                B = (B - 1) & S
            new[S] = acc
        cur = new
        vals.append(cur[full])
    return vals


def interpolate_int_poly(vals):
    """Coefficients (low to high) of the unique polynomial of degree < len(vals) with p(q) = vals[q]."""
    import sympy as sp
    q = sp.symbols('q')
    poly = sp.interpolate([(i, v) for i, v in enumerate(vals)], q)
    poly = sp.Poly(sp.expand(poly), q)
    coeffs = poly.all_coeffs()[::-1]
    return [sp.Rational(c) for c in coeffs]


def cover_poly(A):
    n = len(A)
    vals = cover_poly_values(A, n + 1)
    co = interpolate_int_poly(vals)
    co = [int(c) for c in co]
    while len(co) < n + 1:
        co.append(0)
    return co  # co[k] = coefficient of q^k
