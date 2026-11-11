"""a8_thm448.py (AUDIT) - is the Cayley-Dickson sign tournament CD_n the doubly regular 'Mersenne doubling tower'
T_(2^n - 1) of THM-448 (kind-pasteur-2026-06-09-S1; S -> [[S, S], [S - 2I, 2I - S]], closed form
S[i, j] = (-1)^(popcount((i & j) >> (v + 1)) + bit_v(i)), v = lsb(i xor j);  H(T_15) = 198,335,025, |Aut| = 21)?

Test: explicit isomorphism search (backtracking with invariants) between CD_n (own doubling) and T_(2^n-1)
(both the block recursion and the closed form of THM-448), n = 3, 4, 5; also against the converse."""
import sys, os, itertools
from aud_common import ham_c, HERE

out = []


def P(s):
    print(s)
    out.append(s)
    sys.stdout.flush()


def conj(x):
    if len(x) == 1:
        return x[:]
    h = len(x) // 2
    return conj(x[:h]) + [-t for t in x[h:]]


def mul(x, y):
    if len(x) == 1:
        return [x[0] * y[0]]
    h = len(x) // 2
    a, b, c, d = x[:h], x[h:], y[:h], y[h:]
    return [p - q for p, q in zip(mul(a, c), mul(conj(d), b))] + [p + q for p, q in zip(mul(d, a), mul(b, conj(c)))]


def cd_tournament(n):
    N = 1 << n
    out_ = [set() for _ in range(N - 1)]
    for i in range(1, N):
        for j in range(1, N):
            if i == j:
                continue
            x = [0] * N
            y = [0] * N
            x[i] = 1
            y[j] = 1
            if mul(x, y)[i ^ j] == 1:
                out_[i - 1].add(j - 1)
    return out_


def tower_matrix(k):
    """THM-447/448 block recursion from S_1 = [1]."""
    S = [[1]]
    for _ in range(k):
        m = len(S)
        I2 = [[2 if i == j else 0 for j in range(m)] for i in range(m)]
        top = [S[i] + S[i] for i in range(m)]
        bot = [[S[i][j] - I2[i][j] for j in range(m)] + [I2[i][j] - S[i][j] for j in range(m)] for i in range(m)]
        S = top + bot
    return S


def tower_closed(k):
    N = 1 << k
    S = [[1] * N for _ in range(N)]
    for i in range(N):
        for j in range(N):
            if i == j:
                continue
            x = i ^ j
            v = (x & -x).bit_length() - 1
            e = bin((i & j) >> (v + 1)).count("1") + ((i >> v) & 1)
            S[i][j] = -1 if e & 1 else 1
    return S


def core(S):
    N = len(S)
    assert all(S[0][j] == 1 for j in range(N))
    return [set(j - 1 for j in range(1, N) if j != i and S[i][j] == 1) for i in range(1, N)]


def invariants(T):
    n = len(T)
    inv = []
    for v in range(n):
        # number of cyclic triangles through v, and the multiset of triple intersections with v
        outs = T[v]
        c3 = sum(1 for a in outs for b in range(n) if b != v and b != a and b in T[a] and v in T[b])
        tri = sorted(len(T[v] & T[a] & T[b]) for a in range(n) for b in range(a + 1, n) if a != v and b != v)
        inv.append((len(outs), c3, tuple(tri)))
    return inv


def isomorphic(A, B):
    n = len(A)
    ia, ib = invariants(A), invariants(B)
    if sorted(ia) != sorted(ib):
        return False
    img = [-1] * n
    used = [False] * n
    sys.setrecursionlimit(10000)

    def rec(k):
        if k == n:
            return True
        for v in range(n):
            if used[v] or ib[v] != ia[k]:
                continue
            if all((k in A[i]) == (v in B[img[i]]) for i in range(k)):
                img[k] = v
                used[v] = True
                if rec(k + 1):
                    return True
                used[v] = False
        return False
    return rec(0)


def converse(T):
    n = len(T)
    return [set(j for j in range(n) if j != i and i in T[j]) for i in range(n)]


for k in (2, 3, 4, 5):
    Sb = tower_matrix(k)
    Sc = tower_closed(k)
    same = (Sb == Sc)
    Tt = core(Sb)
    Tc = cd_tournament(k)
    N = 1 << k
    had = all(sum(Sb[i][t] * Sb[j][t] for t in range(N)) == (N if i == j else 0) for i in range(N) for j in range(N))
    iso = isomorphic(Tc, Tt)
    isoc = isomorphic(Tc, converse(Tt))
    hs = ham_c([Tc, Tt]) if k <= 4 else ("-", "-")
    P("n = %d (%d vertices): THM-448 block recursion == its closed form: %s; skew-Hadamard: %s; H(CD_n) = %s, H(T) = %s; CD_n isomorphic to T: %s; to the converse of T: %s"
      % (k, N - 1, same, had, hs[0], hs[1], iso, isoc))
open(os.path.join(HERE, "a8_thm448.out"), "w").write("\n".join(out) + "\n")
