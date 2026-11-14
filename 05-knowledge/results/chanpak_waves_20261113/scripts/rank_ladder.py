#!/usr/bin/env python3
"""Wave 23: the rank ladder of Chan-Pak transfer systems.

For each tournament family that is computed by a finite linear transfer (one matrix per letter), compute the
associative algebra A spanned by the letter words, its semisimple quotient A/rad(A) (rank of the trace form
tr(ab), exact over Q via two large primes), and whether the letter commutators lie in rad(A), i.e. whether the
whole family is simultaneously triangularizable (McCoy: triangularizable iff A/rad(A) is commutative).

Families:
  one-apex  O/Z     4-dim state (H,E,L,h)          letters O, Z
  root-change O/Z/J 4-dim state (H,E,L,h)          letters O, Z, J
  one-hub           6-state segment automaton       letters x_i in {0,1}
  two-hub           36-state segment automaton      letters (x_i, y_i) in {0,1}^2 (scripts/twohub.c)

Then the Zariski-density certificate for the invertible root-change blocks acting on (H,E,lam):
the characteristic polynomial of J is x^3 - x^2 - x - 1 (irreducible, discriminant -44, Galois group S3), and an
invertible block B with B^6 J != J B^6 shows that the Zariski closure is all of GL3 (proof in the README).
"""
import itertools
import sys

PRIMES = (2305843009213693951, 4611686018427387847)  # 2^61-1 and a second large prime


def matmul(A, B, p=None):
    n, m, k = len(A), len(B), len(B[0])
    out = [[0] * k for _ in range(n)]
    for i in range(n):
        Ai = A[i]
        oi = out[i]
        for t in range(m):
            a = Ai[t]
            if a:
                Bt = B[t]
                for j in range(k):
                    oi[j] += a * Bt[j]
        if p:
            for j in range(k):
                oi[j] %= p
    return out


def identity(n):
    return [[int(i == j) for j in range(n)] for i in range(n)]


def flat(M):
    return [x for row in M for x in row]


def unflat(v, n):
    return [v[i * n:(i + 1) * n] for i in range(n)]


class ModBasis:
    """Incremental row-echelon basis mod p."""

    def __init__(self, p):
        self.p = p
        self.rows = {}  # pivot -> row (normalised, pivot entry 1)

    def reduce(self, v):
        p = self.p
        v = [x % p for x in v]
        for piv, r in self.rows.items():
            c = v[piv]
            if c:
                v = [(x - c * y) % p for x, y in zip(v, r)]
        return v

    def add(self, v):
        v = self.reduce(v)
        for i, x in enumerate(v):
            if x:
                inv = pow(x, self.p - 2, self.p)
                v = [(y * inv) % self.p for y in v]
                # keep fully reduced
                for piv in list(self.rows):
                    r = self.rows[piv]
                    c = r[i]
                    if c:
                        self.rows[piv] = [(a - c * b) % self.p for a, b in zip(r, v)]
                self.rows[i] = v
                return True
        return False


def rank_mod(rows, p):
    rows = [[x % p for x in r] for r in rows]
    rank = 0
    ncols = len(rows[0]) if rows else 0
    for c in range(ncols):
        piv = None
        for i in range(rank, len(rows)):
            if rows[i][c]:
                piv = i
                break
        if piv is None:
            continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        inv = pow(rows[rank][c], p - 2, p)
        rows[rank] = [(x * inv) % p for x in rows[rank]]
        for i in range(len(rows)):
            if i != rank and rows[i][c]:
                f = rows[i][c]
                rows[i] = [(x - f * y) % p for x, y in zip(rows[i], rows[rank])]
        rank += 1
    return rank


def algebra(letters, p):
    """Span of all words in the letter matrices (unital), mod p. Returns list of basis matrices (mod p)."""
    n = len(letters[0])
    basis = ModBasis(p)
    mats = []
    frontier = [identity(n)]
    basis.add(flat(identity(n)))
    mats.append(identity(n))
    while frontier:
        new = []
        for W in frontier:
            for L in letters:
                P = matmul(W, L, p)
                if basis.add(flat(P)):
                    mats.append(P)
                    new.append(P)
        frontier = new
    return mats


def trace_prod(A, B, p):
    n = len(A)
    s = 0
    for i in range(n):
        Ai = A[i]
        for t in range(n):
            s += Ai[t] * B[t][i]
    return s % p


def analyse(name, letters):
    out = []
    for p in PRIMES:
        mats = algebra(letters, p)
        d = len(mats)
        G = [[trace_prod(a, b, p) for b in mats] for a in mats]
        rs = rank_mod(G, p)
        # commutators of letters in rad(A)?
        tri = True
        for x, y in itertools.combinations(range(len(letters)), 2):
            C = [[(u - v) % p for u, v in zip(r1, r2)]
                 for r1, r2 in zip(matmul(letters[x], letters[y], p), matmul(letters[y], letters[x], p))]
            if any(trace_prod(C, b, p) for b in mats):
                tri = False
        # centre dimension of A/rad: elements z with tr([z,a] b)=0 ... (skipped; report rs only)
        out.append((d, rs, tri))
    (d1, r1, t1), (d2, r2, t2) = out
    assert (d1, r1, t1) == (d2, r2, t2), (name, out)
    n = len(letters[0])
    print(f"{name:28s} state dim {n:3d}  dim A = {d1:5d}  dim A/rad A = {r1:4d}  triangularizable: {t1}")
    return d1, r1, t1


# ---------------- families ----------------
O4 = [[1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 2, 0], [0, 0, 0, 1]]
Z4 = [[1, -1, 1, 0], [0, 0, 0, 1], [0, 0, 2, -1], [0, 0, 0, 1]]
J4 = [[1, -1, 1, 0], [1, -1, 0, 0], [1, 0, 1, 0], [1, 0, 0, 0]]


def check_ozj_rules():
    # the matrices must reproduce caterpillar.c's rules on random states
    import random
    for _ in range(200):
        H, E, L, h = (random.randint(-50, 50) for _ in range(4))
        v = [[H], [E], [L], [h]]
        assert [r[0] for r in matmul(O4, v)] == [H, 0, 2 * L, h]
        assert [r[0] for r in matmul(Z4, v)] == [H + L - E, h, 2 * L - h, h]
        assert [r[0] for r in matmul(J4, v)] == [H + L - E, H - E, H + L, H]


def onehub_letters():
    # states (s0 in 0..2, s1 in 0..1); index = s0*2 + s1 ; row vector convention: new = old * M
    idx = lambda s0, s1: s0 * 2 + s1
    Ms = []
    for x in (0, 1):
        M = [[0] * 6 for _ in range(6)]
        for s0 in range(3):
            for s1 in range(2):
                a = idx(s0, s1)
                # i into S0
                M[a][idx(2 if x else 1, s1)] += 1
                # i into S1
                if s1 == 0:
                    if x:
                        M[a][idx(s0, 1)] += 1
                else:
                    M[a][idx(s0, 1)] += 1
        Ms.append(M)
    return Ms


def onehub_count(word):
    Ms = onehub_letters()
    v = [0] * 6
    v[0] = 1
    for x in word:
        v = [sum(v[a] * Ms[x][a][b] for a in range(6)) for b in range(6)]
    # final: s0 != 2
    return sum(v[s0 * 2 + s1] for s0 in (0, 1) for s1 in (0, 1))


def onehub_brute(word):
    # spine 1..m, i->j for i<j, hub u with u->i iff x_i = 1 ; count Hamiltonian paths
    m = len(word)
    n = m + 1
    U = m
    adj = [[False] * n for _ in range(n)]
    for i in range(m):
        for j in range(i + 1, m):
            adj[i][j] = True
        if word[i]:
            adj[U][i] = True
        else:
            adj[i][U] = True
    cnt = 0
    for perm in itertools.permutations(range(n)):
        if all(adj[perm[k]][perm[k + 1]] for k in range(n - 1)):
            cnt += 1
    return cnt


def twohub_letters():
    # port of scripts/twohub.c step(): state c[o][s0][s1][s2], o in 0..1, s0 in 0..2, s1 in 0..2, s2 in 0..1
    def idx(o, s0, s1, s2):
        return ((o * 3 + s0) * 3 + s1) * 2 + s2
    Ms = []
    for xi in (0, 1):
        for yi in (0, 1):
            M = [[0] * 36 for _ in range(36)]
            for o in range(2):
                h1beats = xi if o == 0 else yi
                h2beats = yi if o == 0 else xi
                for s0 in range(3):
                    for s1 in range(3):
                        for s2 in range(2):
                            a = idx(o, s0, s1, s2)
                            M[a][idx(o, 2 if h1beats else 1, s1, s2)] += 1
                            if s1 == 0:
                                if h1beats:
                                    M[a][idx(o, s0, 2 if h2beats else 1, s2)] += 1
                            else:
                                M[a][idx(o, s0, 2 if h2beats else 1, s2)] += 1
                            if s2 == 0:
                                if h2beats:
                                    M[a][idx(o, s0, s1, 1)] += 1
                            else:
                                M[a][idx(o, s0, s1, 1)] += 1
            Ms.append(M)
    return Ms


def twohub_count(xw, yw, e):
    Ms = twohub_letters()
    def idx(o, s0, s1, s2):
        return ((o * 3 + s0) * 3 + s1) * 2 + s2
    v = [0] * 36
    v[idx(0, 0, 0, 0)] = 1
    v[idx(1, 0, 0, 0)] = 1
    for xi, yi in zip(xw, yw):
        M = Ms[xi * 2 + yi]
        v = [sum(v[a] * M[a][b] for a in range(36)) for b in range(36)]
    H = 0
    for o in range(2):
        h1toh2 = e if o == 0 else 1 - e
        for s0 in range(3):
            for s1 in range(3):
                for s2 in range(2):
                    c = v[idx(o, s0, s1, s2)]
                    if not c or s0 == 2 or s1 == 2:
                        continue
                    if s1 == 0 and not h1toh2:
                        continue
                    H += c
    return H


def twohub_brute(xw, yw, e):
    m = len(xw)
    n = m + 2
    u1, u2 = m, m + 1
    adj = [[False] * n for _ in range(n)]
    for i in range(m):
        for j in range(i + 1, m):
            adj[i][j] = True
        if xw[i]:
            adj[u1][i] = True
        else:
            adj[i][u1] = True
        if yw[i]:
            adj[u2][i] = True
        else:
            adj[i][u2] = True
    if e:
        adj[u1][u2] = True
    else:
        adj[u2][u1] = True
    # Hamiltonian paths by DP over subsets
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1 << n):
        for v in range(n):
            c = dp[S][v]
            if not c:
                continue
            for w in range(n):
                if not (S >> w) & 1 and adj[v][w]:
                    dp[S | (1 << w)][w] += c
    return sum(dp[(1 << n) - 1])


def main():
    check_ozj_rules()
    import random
    random.seed(23)
    for _ in range(300):
        m = random.randint(1, 7)
        xw = [random.randint(0, 1) for _ in range(m)]
        yw = [random.randint(0, 1) for _ in range(m)]
        e = random.randint(0, 1)
        assert twohub_count(xw, yw, e) == twohub_brute(xw, yw, e), (xw, yw, e)
    print("two-hub transfer (port of scripts/twohub.c) matches subset-DP brute force on 300 random words")
    for w in itertools.product((0, 1), repeat=6):
        assert onehub_count(list(w)) == onehub_brute(list(w)), w
    print("one-hub transfer matches brute force on all 64 words of length 6")
    print("P46 rank ladder (exact algebra, two primes):")
    res = {}
    res['oneapex'] = analyse("one-apex O/Z (H,E,L,h)", [O4, Z4])
    res['ozj'] = analyse("root-change O/Z/J (H,E,L,h)", [O4, Z4, J4])
    res['onehub'] = analyse("one-hub (6 states)", onehub_letters())
    res['twohub'] = analyse("two-hub (36 states)", twohub_letters())
    ok = res['oneapex'][2] and res['onehub'][2] and not res['ozj'][2] and not res['twohub'][2]
    print("P46 prediction (a) and (b) both hold:", ok)

    # ---------- Zariski density certificate for the invertible root-change blocks (3x3 on (H,E,lam)) -------------
    import sympy as sp
    x = sp.symbols('x')

    def step(s, m):
        H, E, L, h = s
        if m == 'O':
            return (H, 0, 2 * L, h)
        if m == 'Z':
            return (H + L - E, h, 2 * L - h, h)
        return (H + L - E, H - E, H + L, H)

    def blockmat(w):
        cols = []
        for e in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
            s = step((e[0], e[1], e[2], 0), 'J')
            for ch in w:
                s = step(s, ch)
            cols.append(s[:3])
        return sp.Matrix([[cols[c][r] for c in range(3)] for r in range(3)])

    J = blockmat('')
    f = J.charpoly(x).as_expr()
    print("char poly of J on (H,E,lam):", sp.expand(f), " discriminant:", sp.discriminant(f, x),
          " irreducible over Q:", sp.Poly(f, x).is_irreducible)
    for w in ['OZ', 'OOZ', 'OZZ', 'OOOZ', 'OZOZ', 'OZZZ']:
        B = blockmat(w)
        B6 = B ** 6
        print(f"  block J{w:5s} det {int(B.det()):3d}   B^6 J == J B^6 ? {B6 * J == J * B6}")


if __name__ == '__main__':
    main()
