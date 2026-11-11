"""aud_core.py - independent counters written by the auditor of kind-pasteur-2026-11-10-S22.

Nothing here imports the session's towers.py / tower_build.py.

Two unrelated counting methods:
  (1) brute_*  : plain enumeration of all vertex permutations (ground truth, small orders);
  (2) dp_*     : a "pull" subset recursion over (vertex set, last vertex) written from scratch.

Conventions (as in THM-5130):
  A[i][j] = 1 iff i -> j.
  Q_B(e,f): ordered pairs (P1,P2) of disjoint NONEMPTY paths covering B, P1 ending at e, P2 starting at f.
  W_B(O) = sum over ordered pairs e != f in O of Q_B(e,f).
  U_B(o): ordered pairs (P1,P2) covering B, P1 possibly empty, P2 starting at o.
  state(B,O) = (r, M, h, E) = (H(B) - W_B(O), U_B(O), H(B), #Hamiltonian paths of B ending in O).
"""
import itertools


def adj_from_code(m, code):
    """upper-triangle code: bit b (b-th pair i<j in lexicographic order) is 1 iff i -> j."""
    A = [[0] * m for _ in range(m)]
    b = 0
    for i in range(m):
        for j in range(i + 1, m):
            if (code >> b) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            b += 1
    return A


def is_tournament(A):
    n = len(A)
    for i in range(n):
        if A[i][i]:
            return False
        for j in range(i + 1, n):
            if A[i][j] + A[j][i] != 1:
                return False
    return True


def transitive(a):
    """TT_a on 0..a-1 with i -> j for i < j (vertex k of the theorem = index k-1)."""
    return [[1 if i < j else 0 for j in range(a)] for i in range(a)]


def close(A, O, X):
    """C(B, O; X): B => X => v, v -> O, B - O -> v.  Vertices: B (0..n-1), X (n..n+k-1), v (n+k)."""
    n, k = len(A), len(X)
    N = n + k + 1
    T = [[0] * N for _ in range(N)]
    Oset = set(O)
    for i in range(n):
        for j in range(n):
            T[i][j] = A[i][j]
        for j in range(k):
            T[i][n + j] = 1
        if i in Oset:
            T[N - 1][i] = 1
        else:
            T[i][N - 1] = 1
    for i in range(k):
        for j in range(k):
            T[n + i][n + j] = X[i][j]
        T[n + i][N - 1] = 1
    assert is_tournament(T)
    return T


# ---------------------------------------------------------------- brute force (permutations)

def _breaks(A, perm):
    """positions k (1..n-1) where perm[k-1] does not beat perm[k]."""
    return [k for k in range(1, len(perm)) if not A[perm[k - 1]][perm[k]]]


def brute_H(A):
    n = len(A)
    if n == 0:
        return 1
    return sum(1 for p in itertools.permutations(range(n)) if not _breaks(A, p))


def brute_state(A, O):
    """(r, M, h, E) by enumerating every permutation and every cut position."""
    n = len(A)
    Oset = set(O)
    h = W = M = E = 0
    for p in itertools.permutations(range(n)):
        br = _breaks(A, p)
        if len(br) >= 2:
            continue
        if not br:
            h += 1
            if p[-1] in Oset:
                E += 1
            # every cut is valid
            for k in range(0, n):           # P1 = p[:k] (possibly empty), P2 = p[k:] nonempty
                if p[k] in Oset:
                    M += 1
                if k >= 1 and p[k - 1] in Oset and p[k] in Oset:
                    W += 1
        else:
            k = br[0]                       # the only admissible cut
            if p[k] in Oset:
                M += 1
                if p[k - 1] in Oset:
                    W += 1
    return (h - W, M, h, E)


def brute_Q(A, e, f):
    n = len(A)
    c = 0
    for p in itertools.permutations(range(n)):
        br = _breaks(A, p)
        if len(br) >= 2:
            continue
        ks = range(1, n) if not br else br
        for k in ks:
            if p[k - 1] == e and p[k] == f:
                c += 1
    return c


def brute_gadget(X, P):
    """constants of a gadget (X, P) by permutation enumeration with one or two cuts.
    returns dict with H, Pi, w, w3, sigma, u, eps, nu, u3."""
    n = len(X)
    Ps = set(P)
    HX = Pi = w = w3 = sigma = u = eps = nu = u3 = 0
    for p in itertools.permutations(range(n)):
        br = _breaks(X, p)
        nb = len(br)
        if nb == 0:
            HX += 1
            Pi += n + 1
            if p[0] in Ps:
                sigma += 1
            if p[-1] in Ps:
                eps += 1
        elif nb == 1:
            Pi += 1
        # ordered pairs: cut k in 0..n ; (D1, D2) = (p[:k], p[k:])
        if nb <= 1:
            ks = range(0, n + 1) if nb == 0 else br
            for k in ks:
                # u: D2 nonempty starting in P, D1 possibly empty
                if k <= n - 1 and p[k] in Ps:
                    u += 1
                # nu: D1 nonempty ending in P, D2 possibly empty
                if k >= 1 and p[k - 1] in Ps:
                    nu += 1
                # w: both nonempty, D1 ends at e in P, D2 starts at f in P
                if 1 <= k <= n - 1 and p[k - 1] in Ps and p[k] in Ps:
                    w += 1
        # ordered triples of nonempty paths: cuts 1 <= k1 < k2 <= n-1, breaks inside {k1, k2}
        if nb <= 2:
            for k1 in range(1, n - 1):
                for k2 in range(k1 + 1, n):
                    if all(b in (k1, k2) for b in br):
                        # w3: D1 ends at e, D2 starts at f, e != f in P
                        if p[k1 - 1] in Ps and p[k1] in Ps:
                            w3 += 1
                        # u3: D1 starts at o in P
                        if p[0] in Ps:
                            u3 += 1
    return dict(H=HX, Pi=Pi, w=w, w3=w3, sigma=sigma, u=u, eps=eps, nu=nu, u3=u3)


def gadget_matrix(g):
    """the 4x4 matrix of Theorem A acting on (r, M, h, E), built from gadget constants."""
    rX = g["H"] - g["w"]
    kap = g["Pi"] - 2 * g["w"] - 2 * g["w3"]
    return [[rX, kap, 0, 0],
            [g["u"], 4 * g["u"] - 2 * g["sigma"] + 2 * g["u3"], g["u"], -g["sigma"]],
            [g["H"], g["Pi"], 0, 0],
            [g["eps"], g["nu"], -g["eps"], 0]]


def mat_apply(G, s):
    return tuple(sum(G[i][j] * s[j] for j in range(4)) for i in range(4))


# ---------------------------------------------------------------- subset recursion ("pull" form)

def dp_tables(A):
    """endat[S][v] = #Hamiltonian paths of A[S] ending at v; startat[S][v] = ... starting at v; tot[S] = H(A[S])."""
    n = len(A)
    N = 1 << n
    inmask = [0] * n
    outmask = [0] * n
    for i in range(n):
        for j in range(n):
            if A[i][j]:
                outmask[i] |= 1 << j
                inmask[j] |= 1 << i
    endat = [None] * N
    startat = [None] * N
    tot = [0] * N
    tot[0] = 1
    endat[0] = [0] * n
    startat[0] = [0] * n
    for S in range(1, N):
        e = [0] * n
        s = [0] * n
        if S & (S - 1) == 0:
            v = S.bit_length() - 1
            e[v] = 1
            s[v] = 1
        else:
            T = S
            while T:
                low = T & -T
                v = low.bit_length() - 1
                T ^= low
                R = S ^ low
                pe = endat[R]
                ps = startat[R]
                # last vertex v, previous vertex x with x -> v
                c = 0
                Xs = R & inmask[v]
                while Xs:
                    l2 = Xs & -Xs
                    c += pe[l2.bit_length() - 1]
                    Xs ^= l2
                e[v] = c
                # first vertex v, next vertex y with v -> y
                c = 0
                Ys = R & outmask[v]
                while Ys:
                    l2 = Ys & -Ys
                    c += ps[l2.bit_length() - 1]
                    Ys ^= l2
                s[v] = c
        endat[S] = e
        startat[S] = s
        tot[S] = sum(e)
    return endat, startat, tot


def dp_H(A):
    if not A:
        return 1
    return dp_tables(A)[2][(1 << len(A)) - 1]


def dp_state(A, O):
    n = len(A)
    full = (1 << n) - 1
    endat, startat, tot = dp_tables(A)
    h = tot[full]
    E = sum(endat[full][o] for o in O)
    M = 0
    W = 0
    for S in range(0, full + 1):       # S = vertex set of P1
        R = full ^ S
        if R == 0:
            continue
        sR = sum(startat[R][o] for o in O if (R >> o) & 1)
        if sR == 0:
            continue
        M += tot[S] * sR
        if S:
            eS = sum(endat[S][o] for o in O if (S >> o) & 1)
            W += eS * sR
    return (h - W, M, h, E)


def dp_pi2(A):
    n = len(A)
    full = (1 << n) - 1
    tot = dp_tables(A)[2]
    return sum(tot[S] * tot[full ^ S] for S in range(full + 1))
