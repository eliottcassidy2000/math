"""AUDIT helpers (independent of the author's modules)."""
import itertools
import subprocess
import sys

P32 = 4294967291  # largest prime below 2^32


def gf_field(q):
    """Return (elements, add, mul, neg, squares) for GF(q), q = p or p^3 (p = 3 only for the cube)."""
    # prime field
    def is_prime(x):
        return x > 1 and all(x % d for d in range(2, int(x ** 0.5) + 1))
    if is_prime(q):
        els = list(range(q))
        add = lambda a, b: (a + b) % q
        mul = lambda a, b: (a * b) % q
        neg = lambda a: (-a) % q
        return els, add, mul, neg
    if q == 27:
        # GF(27) = F_3[x]/(x^3 - x - 1)  (x^3 + 2x + 2 is irreducible over F_3: no root 0,1,2)
        def enc(c):
            return c[0] + 3 * c[1] + 9 * c[2]
        def dec(a):
            return [a % 3, (a // 3) % 3, a // 9]
        def add(a, b):
            x, y = dec(a), dec(b)
            return enc([(x[i] + y[i]) % 3 for i in range(3)])
        def neg(a):
            x = dec(a)
            return enc([(-x[i]) % 3 for i in range(3)])
        def mul(a, b):
            x, y = dec(a), dec(b)
            prod = [0] * 5
            for i in range(3):
                for j in range(3):
                    prod[i + j] = (prod[i + j] + x[i] * y[j]) % 3
            # x^3 = x + 1, x^4 = x^2 + x
            c = prod[:3]
            c[0] = (c[0] + prod[3]) % 3
            c[1] = (c[1] + prod[3]) % 3
            c[1] = (c[1] + prod[4]) % 3
            c[2] = (c[2] + prod[4]) % 3
            return enc(c)
        return list(range(27)), add, mul, neg
    raise ValueError("field not implemented")


def paley(q):
    """Adjacency matrix (list of lists of 0/1) of the Paley tournament: a -> b iff b - a is a non-zero square."""
    els, add, mul, neg = gf_field(q)
    squares = set(mul(a, a) for a in els if a != 0)
    assert len(squares) == (q - 1) // 2
    A = [[1 if (a != b and add(b, neg(a)) in squares) else 0 for b in els] for a in els]
    check_tournament(A)
    return A


def check_tournament(A):
    n = len(A)
    for i in range(n):
        assert A[i][i] == 0
        for j in range(i + 1, n):
            assert A[i][j] + A[j][i] == 1, (i, j)


def from_upper_code(code):
    """Upper-triangle code: '1' at (i, j), i < j, iff i -> j (rows i = 0.., j = i+1..)."""
    L = len(code)
    n = 2
    while n * (n - 1) // 2 < L:
        n += 1
    assert n * (n - 1) // 2 == L, L
    A = [[0] * n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if code[k] == '1':
                A[i][j] = 1
            else:
                A[j][i] = 1
            k += 1
    return A


def is_doubly_regular(A):
    n = len(A)
    if n % 4 != 3:
        return False
    t = (n - 3) // 4
    for i in range(n):
        if sum(A[i]) != 2 * t + 1:
            return False
    for i in range(n):
        for j in range(i + 1, n):
            if sum(A[i][x] & A[j][x] for x in range(n)) != t:
                return False
    return True


def ham_exact(A, binary):
    """Exact number of Hamiltonian paths through the audit C program (two moduli + CRT)."""
    n = len(A)
    text = str(n) + "\n" + "\n".join("".join(str(x) for x in row) for row in A) + "\n"
    r1 = int(subprocess.run([binary, "0"], input=text, capture_output=True, text=True, check=True).stdout)
    r2 = int(subprocess.run([binary, str(P32)], input=text, capture_output=True, text=True, check=True).stdout)
    # CRT: h = r1 mod 2^32, h = r2 mod P32
    M1 = 1 << 32
    inv = pow(M1, -1, P32)
    k = ((r2 - r1) * inv) % P32
    return r1 + M1 * k


def ham_python(A):
    """Exact number of Hamiltonian paths, pure Python subset DP (n <= ~15)."""
    n = len(A)
    out = [sum(1 << j for j in range(n) if A[i][j]) for i in range(n)]
    full = (1 << n) - 1
    dp = [dict() for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        d = dp[mask]
        if not d:
            continue
        rem = full & ~mask
        for u, x in d.items():
            t = out[u] & rem
            while t:
                w = (t & -t).bit_length() - 1
                t &= t - 1
                nd = dp[mask | (1 << w)]
                nd[w] = nd.get(w, 0) + x
        if mask != full:
            dp[mask] = None
    return sum(dp[full].values())


def anti_automorphisms(A, limit=None, involutions_only=False):
    """All bijections psi with a -> b iff psi(b) -> psi(a) (backtracking). Returns a list of tuples."""
    n = len(A)
    # anti-automorphism of A = isomorphism from A to its converse B (B[i][j] = A[j][i])
    B = [[A[j][i] for j in range(n)] for i in range(n)]
    # simple invariants to prune: (out-degree, sorted out-degrees of out-neighbours)
    def inv(M, v):
        return (sum(M[v]), tuple(sorted(sum(M[w]) for w in range(n) if M[v][w])))
    ia = [inv(A, v) for v in range(n)]
    ib = [inv(B, v) for v in range(n)]
    res = []
    img = [-1] * n
    used = [False] * n

    def rec(v):
        if limit is not None and len(res) >= limit:
            return
        if v == n:
            p = tuple(img)
            if involutions_only and any(p[p[i]] != i for i in range(n)):
                return
            res.append(p)
            return
        for w in range(n):
            if used[w] or ia[v] != ib[w]:
                continue
            ok = True
            for u in range(v):
                if A[u][v] != B[img[u]][w]:
                    ok = False
                    break
            if ok:
                img[v] = w
                used[w] = True
                rec(v + 1)
                used[w] = False
                img[v] = -1

    rec(0)
    return res


def fold(A, psi):
    """The fold (Z, D1, plus) of the note from its definition. psi: tuple, anti-involution with one fixed vertex."""
    n = len(A)
    fixed = [v for v in range(n) if psi[v] == v]
    assert len(fixed) == 1
    v0 = fixed[0]
    for a in range(n):
        assert psi[psi[a]] == a
        for b in range(n):
            if a != b:
                assert A[a][b] == A[psi[b]][psi[a]]
    plus = [a for a in range(n) if a != v0 and A[v0][a]]
    m = (n - 1) // 2
    assert len(plus) == m
    minus = [psi[a] for a in plus]
    assert all(A[b][v0] for b in minus)
    D1 = [A[plus[i]][minus[i]] for i in range(m)]
    Z = [[0] * m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            if i != j:
                Z[i][j] = 1 if (A[plus[i]][plus[j]] + A[plus[i]][minus[j]]) == 1 else 0
    check_tournament(Z)
    return Z, D1, plus, minus, v0


def paths_ending_in(Z, D1):
    """Number of k-vertex directed paths of Z whose last vertex is in D1, for k = 1..m (subset DP)."""
    m = len(Z)
    out = [sum(1 << j for j in range(m) if Z[i][j]) for i in range(m)]
    full = (1 << m) - 1
    dp = {}
    for v in range(m):
        dp[(1 << v, v)] = 1
    res = [0] * (m + 1)
    layer = dp
    for k in range(1, m + 1):
        nxt = {}
        for (mask, u), x in layer.items():
            if D1[u]:
                res[k] += x
            t = out[u] & ~mask & full
            while t:
                w = (t & -t).bit_length() - 1
                t &= t - 1
                key = (mask | (1 << w), w)
                nxt[key] = nxt.get(key, 0) + x
        layer = nxt
    return res


def cd_sign(k):
    """Sign table s[i][j] of the Cayley-Dickson algebra of dimension 2^k: e_i e_j = s e_(i xor j).
    Doubling rule (a,b)(c,d) = (ac - conj(d) b, d a + b conj(c))."""
    s = [[1]]
    N = 1
    for _ in range(k):
        kap = [1] + [-1] * (N - 1)
        t = [[0] * (2 * N) for _ in range(2 * N)]
        for i in range(N):
            for j in range(N):
                t[i][j] = s[i][j]
                t[i][N + j] = s[j][i]
                t[N + i][j] = s[i][j] * kap[j]
                t[N + i][N + j] = -kap[j] * s[j][i]
        s = t
        N *= 2
    return s


def cd_tournament(k):
    s = cd_sign(k)
    N = 1 << k
    A = [[0] * (N - 1) for _ in range(N - 1)]
    for a in range(1, N):
        for b in range(1, N):
            if a != b:
                A[a - 1][b - 1] = 1 if s[a][b] == 1 else 0
    check_tournament(A)
    return A


def double_tournament(T):
    """The double of the note (section 5.3): V_1 + {oo} + V_2."""
    n = len(T)
    N = 2 * n + 1
    A = [[0] * N for _ in range(N)]
    oo = n
    for a in range(n):
        for b in range(n):
            if a != b:
                A[a][b] = T[a][b]                 # T on V_1
                A[n + 1 + a][n + 1 + b] = T[b][a]  # converse on V_2
                A[a][n + 1 + b] = T[a][b]          # a_1 -> b_2 iff a -> b
                A[n + 1 + b][a] = 1 - T[a][b]
        A[oo][a] = 1                               # oo -> V_1
        A[n + 1 + a][oo] = 1                       # V_2 -> oo
        A[a][n + 1 + a] = 1                        # a_1 -> a_2
    check_tournament(A)
    return A
