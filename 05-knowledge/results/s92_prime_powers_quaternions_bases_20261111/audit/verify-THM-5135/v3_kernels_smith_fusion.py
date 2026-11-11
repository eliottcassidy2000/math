"""Verifier of THM-5135 (S92 audit follow-up): statement (4): (a) unit determinant, (b) the Smith-form law (2 x 2 and
3 x 3, singular matrices included, the prime 2 included), the quoted numbers, (c) the two Gram matrices of rank 4, the
fusion law of Session264 and the parities of its states.  Own code, exact.
Run: python v3_kernels_smith_fusion.py > v3_kernels_smith_fusion.out
"""
import itertools
import math
import random
from fractions import Fraction as F

import numpy as np

FAILS = []
INF = 10 ** 9


def check(cond, msg):
    print(("  ok   " if cond else "  FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def vp(n, p):
    if n == 0:
        return INF
    n = abs(n)
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def scaled_table(q, scale):
    d = np.zeros(q, dtype=object)
    for a in range(q):
        for b in range(q):
            d[(scale * a * b) % q] += 1
    return d


def conv(a, b):
    q = len(a)
    out = np.zeros(q, dtype=object)
    for i in range(q):
        if a[i]:
            for j in range(q):
                if b[j]:
                    out[(i + j) % q] += a[i] * b[j]
    return out


def law_brute(G, q):
    """class counts of u^T G w over all (u, w) in (Z/q)^r x (Z/q)^r, direct enumeration."""
    r = len(G)
    g = np.arange(q, dtype=np.int64)
    grids = np.meshgrid(*[g] * (2 * r), indexing="ij")
    val = np.zeros_like(grids[0])
    for i in range(r):
        for j in range(r):
            c = G[i][j] % q
            if c:
                val = (val + c * grids[i] * grids[r + j]) % q
    return np.bincount(val.ravel(), minlength=q).astype(object)


def det_exact(M):
    n = len(M)
    if n == 1:
        return M[0][0]
    if n == 2:
        return M[0][0] * M[1][1] - M[0][1] * M[1][0]
    return sum((-1) ** j * M[0][j] * det_exact([row[:j] + row[j + 1:] for row in M[1:]]) for j in range(n))


def minors_gcd(G, size):
    r = len(G)
    g = 0
    for rows in itertools.combinations(range(r), size):
        for cols in itertools.combinations(range(r), size):
            g = math.gcd(g, det_exact([[G[i][j] for j in cols] for i in rows]))
    return g


def smith_exponents(G, p, k):
    """min(e_i, k) for the elementary divisors p^(e_i) of the integer matrix G over Z_p (determinantal divisors)."""
    r = len(G)
    prev = 0
    out = []
    for size in range(1, r + 1):
        v = vp(minors_gcd(G, size), p)
        if v >= INF or prev >= INF:
            e = INF
            prev = INF
        else:
            e = v - prev
            prev = v
        out.append(min(e, k))
    return out


def law_smith(G, p, k):
    q = p ** k
    es = smith_exponents(G, p, k)
    out = None
    for e in es:
        t = scaled_table(q, p ** e % q if e < k else 0)
        out = t if out is None else conv(out, t)
    return out, es


def by_valuation(cnt, p, k):
    q = p ** k
    out = []
    for v in range(k + 1):
        c = (p ** v) % q
        out.append(int(cnt[c]))
    return out


# ---------------------------------------------------------------- [1] Smith law, 2 x 2
print("[1] (4)(b): Smith-form law, all 2 x 2 matrices modulo 9 and modulo 8 (direct enumeration of (u, w))")
for (p, k) in ((3, 2), (2, 3)):
    q = p ** k
    bad = 0
    classes = {}
    pred_cache = {}
    g = np.arange(q, dtype=np.int64)
    U1, U2, W1, W2 = np.meshgrid(g, g, g, g, indexing="ij")
    m11, m12, m21, m22 = (U1 * W1).ravel(), (U1 * W2).ravel(), (U2 * W1).ravel(), (U2 * W2).ravel()
    for a, b, c, d in itertools.product(range(q), repeat=4):
        cnt = np.bincount((a * m11 + b * m12 + c * m21 + d * m22) % q, minlength=q)
        es = tuple(smith_exponents([[a, b], [c, d]], p, k))
        if es not in pred_cache:
            pred_cache[es] = law_smith([[a, b], [c, d]], p, k)[0]
        if [int(x) for x in cnt] != [int(x) for x in pred_cache[es]]:
            bad += 1
        classes[es] = classes.get(es, 0) + 1
    check(bad == 0, f"modulo {q}: {q ** 4} matrices, {len(classes)} Smith classes {dict(sorted(classes.items()))}, mismatches {bad}")

print("[1b] (4)(b): random 3 x 3 matrices (singular ones included) and non-symmetric ones, moduli 4, 8, 9, 27 (r = 2), 25 (r = 2)")
random.seed(5135)
bad = 0
n = 0
seen = set()
for (p, k, r, trials, lo, hi) in ((2, 2, 3, 60, -4, 4), (2, 3, 3, 25, -8, 8), (3, 2, 3, 30, -9, 9), (3, 3, 2, 60, -13, 13), (5, 2, 2, 60, -12, 12)):
    q = p ** k
    for t in range(trials):
        G = [[random.randint(lo, hi) for _ in range(r)] for _ in range(r)]
        if t % 4 == 0:                                   # force a singular matrix: last row = combination of the others
            G[r - 1] = [sum(random.randint(-2, 2) * G[i][j] for i in range(r - 1)) for j in range(r)]
        if t % 7 == 1:                                   # force divisibility by p
            G = [[p * x for x in row] for row in G]
        pred, es = law_smith(G, p, k)
        seen.add((p, k, r, tuple(es)))
        n += 1
        if [int(x) for x in law_brute(G, q)] != [int(x) for x in pred]:
            bad += 1
check(bad == 0, f"{n} matrices, {len(seen)} distinct (modulus, rank, Smith class) cases, mismatches {bad}")
print(f"       Smith classes met: {sorted(seen)}")

# ---------------------------------------------------------------- [2] the quoted numbers
print("[2] (4): the quoted kernels")


def rank_law_counts(p, k, r):
    q = p ** k
    t = scaled_table(q, 1)
    out = t
    for _ in range(r - 1):
        out = conv(out, t)
    return out


G1 = [[5, 2], [2, 1]]
G2 = [[11, 6], [6, 3]]
G3 = [[1, 0], [0, 2]]
check(5 * 1 - 2 * 2 == 1 and 11 * 3 - 36 == -3 and 149 * 37 - 72 * 72 == 329 == 7 * 47 and 1 * 5 - 2 * 3 == -1, "determinants 1, -3, 329 = 7 * 47, -1")
for name, G in (("[[5,2],[2,1]]", G1), ("[[11,6],[6,3]]", G2), ("diag(1,2)", G3)):
    res = []
    for (p, k) in ((3, 2), (5, 2), (7, 2)):
        res.append([int(x) for x in law_brute(G, p ** k)] == [int(x) for x in rank_law_counts(p, k, 2)])
    want = [name != "[[11,6],[6,3]]", True, True]
    check(res == want, f"{name}: rank-2 law modulo 9, 25, 49: {res}")
c9 = law_brute(G2, 9)
check(by_valuation(c9, 3, 2) == [486, 1134, 1377] and by_valuation(rank_law_counts(3, 2, 2), 3, 2) == [648, 864, 945]
      and smith_exponents(G2, 3, 2) == [0, 1] and [int(x) for x in c9] == [int(x) for x in conv(scaled_table(9, 1), scaled_table(9, 3))],
      f"[[11,6],[6,3]] modulo 9: counts by valuation {by_valuation(c9, 3, 2)} (law of a1 b1 + 3 a2 b2), rank-2 law {by_valuation(rank_law_counts(3, 2, 2), 3, 2)}")
K = [[149, 72], [72, 37]]
c49 = law_brute(K, 49)
fac = [F(x * 49, 49 ** 4) for x in by_valuation(c49, 7, 2)]
check(fac == [F(6, 7), F(90, 49), F(97, 49)] and [int(x) for x in c49] == [int(x) for x in conv(scaled_table(49, 1), scaled_table(49, 7))],
      f"Session261 kernel modulo 49: factors {[str(f) for f in fac]} at the valuations 0, 1, 2; law of a1 b1 + 7 a2 b2")
for p in (7, 47):
    cp = law_brute(K, p)
    check([int(x) for x in cp] == [p * p * int(x) for x in scaled_table(p, 1)] and smith_exponents(K, p, 1) == [0, 1],
          f"Session261 kernel modulo {p}: p^2 times the multiplication table (elementary divisors 1, p); class 0 factor {F(int(cp[0]) * p, p ** 4)}")
check([x % 7 for row in K for x in row] == [2, 2, 2, 2], "Session261 kernel modulo 7 is 2 [[1,1],[1,1]]")
check(all([int(x) for x in law_brute(K, p)] == [int(x) for x in rank_law_counts(p, 1, 2)] for p in (3, 5, 11, 13)), "Session261 kernel: rank-2 law modulo 3, 5, 11, 13")
check([int(x) for x in law_brute([[1, 2], [3, 5]], 27)] == [int(x) for x in rank_law_counts(3, 3, 2)], "[[1,2],[3,5]] (det -1): rank-2 law modulo 27")

# ---------------------------------------------------------------- [3] the two rank-4 forms
print("[3] (4)(c): Gram matrices of rank 4")


def det_int(M):
    n = len(M)
    A = [[F(x) for x in row] for row in M]
    det = F(1)
    for i in range(n):
        piv = next((r for r in range(i, n) if A[r][i] != 0), None)
        if piv is None:
            return 0
        if piv != i:
            A[i], A[piv] = A[piv], A[i]
            det = -det
        det *= A[i][i]
        for r in range(i + 1, n):
            f = A[r][i] / A[i][i]
            for c in range(i, n):
                A[r][c] -= f * A[i][c]
    return int(det)


GK = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, -1, 1], [0, 0, 1, 0]]
GF = [[1, 1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
GM = [[1, 1, 0, 0], [1, 0, 0, 0], [0, 0, 2, 0], [0, 0, 0, 2]]
check(det_int(GK) == -1 and det_int(GF) == -1 and det_int(GM) == -4, f"determinants: K-form {det_int(GK)}, fusion corner {det_int(GF)}, path-mass coordinates {det_int(GM)}")
random.seed(264)
okK = okF = okM = True
for _ in range(2000):
    U, E, d, a, U2, E2, d2, a2 = [random.randint(-50, 50) for _ in range(8)]
    lhs = U * U2 + E * E2 - d * d2 + a * d2 + a2 * d
    rhs = sum(GK[i][j] * x * y for i, x in enumerate((U, E, d, a)) for j, y in enumerate((U2, E2, d2, a2)))
    okK &= lhs == rhs
    h, H, Uu, dd, h2, H2, V2, e2 = [random.randint(-50, 50) for _ in range(8)]
    s264 = H * h2 + h * H2 - h * h2 + (Uu - h) * (V2 - h2) + (dd - h) * (e2 - h2)
    disp = h * h2 + (H - h) * h2 + h * (H2 - h2) + (Uu - h) * (V2 - h2) + (dd - h) * (e2 - h2)
    left = (h, H - h, Uu - h, dd - h)
    right = (h2, H2 - h2, V2 - h2, e2 - h2)
    okF &= s264 == disp == sum(GF[i][j] * left[i] * right[j] for i in range(4) for j in range(4))
    ao, ae, bo, be = [random.randint(-50, 50) for _ in range(4)]
    okM &= (ao + ae) * (bo + be) + (ao - ae) * (bo - be) == 2 * (ao * bo + ae * be)
check(okK, "K-form U U' + E E' - d d' + a d' + a' d = u^T G w with the displayed matrix (2000 random integer points)")
check(okF, "Session264 eq. (3) = the form displayed in (4)(c) = u^T G w with [[1,1,0,0],[1,0,0,0],[0,0,1,0],[0,0,0,1]] (2000 points)")
check(okM, "(a_o + a_e)(b_o + b_e) + (a_o - a_e)(b_o - b_e) = 2 (a_o b_o + a_e b_e)")
for name, G in (("K-form", GK), ("fusion corner", GF)):
    res = [[int(x) for x in law_brute(G, q)] == [int(x) for x in rank_law_counts(p, k, 4)] for (p, k, q) in ((2, 1, 2), (3, 1, 3), (2, 2, 4), (5, 1, 5))]
    check(all(res), f"{name}, uniform states: law (1) with r = 4 modulo 2, 3, 4, 5: {res}")
res = {}
for (p, k) in ((2, 1), (2, 2), (3, 1), (5, 1)):
    q = p ** k
    res[q] = ([int(x) for x in law_brute(GM, q)] == [int(x) for x in rank_law_counts(p, k, 4)], smith_exponents(GM, p, k))
check(res[2][0] is False and res[4][0] is False and res[3][0] is True and res[5][0] is True and res[4][1] == [0, 0, 1, 1],
      f"path-mass Gram matrix (determinant -4): law (1) with r = 4 modulo 2, 4, 3, 5: {[res[q][0] for q in (2, 4, 3, 5)]} (unimodular at the odd primes only; "
      f"Smith exponents at 2: {res[4][1]})")
check(all([int(x) for x in law_brute(GM, 2 ** k)] == [int(x) for x in law_smith(GM, 2, k)[0]] for k in (1, 2)),
      "path-mass Gram matrix at 2: its law modulo 2 and 4 is the convolution of the tables a b, a b, 2 a b, 2 a b (Smith-form law, rank 4)")
check(1 + F(7 * 2, 2 * 27) == F(34, 27) and 1 + F(1, 27) - F(1, 81) == F(83, 81) and F(15, 16) + F(1, 8) == F(17, 16),
      "34/27 = 1 + (3p-2)(p-1)/(2p^3) at p = 3; 83/81 = 1 + p^-3 - p^-4; 17/16 = 15/16 + 1/8")

# ---------------------------------------------------------------- [4] fusion law and parities on tournaments
print("[4] (4)(c): Session264's fusion law (3) and the parities of the fusion states (pointed tournaments, root 0)")


def ham_table(adj, n):
    """hp[mask] = number of Hamiltonian paths of the sub-tournament induced on mask (hp[0] = 1)."""
    full = 1 << n
    dp = [[0] * n for _ in range(full)]
    hp = [0] * full
    hp[0] = 1
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, full):
        tot = 0
        for v in range(n):
            c = dp[mask][v]
            if c:
                tot += c
                for w in range(n):
                    if not (mask >> w) & 1 and adj[v][w]:
                        dp[mask | (1 << w)][w] += c
        hp[mask] = tot
    return hp


def profile(adj, n):
    """(H, h, U, V, d, e) for the tournament pointed at 0, as in Session264 PROOF.md section 2."""
    hp = ham_table(adj, n)
    full = (1 << n) - 1
    H = hp[full]
    h = hp[full & ~1]
    U = V = d = e = 0
    a_o = a_e = 0
    # paths starting at 0 (for U, d) and ending at 0 (for V, e): depth-first enumeration
    def walk(v, mask, length, forward):
        nonlocal U, V, d, e, a_o, a_e
        comp = hp[full & ~mask]
        sgn = 1 if length % 2 else -1
        if forward:
            U += comp
            d += sgn * comp
            if length >= 2:
                if length % 2:
                    a_o += comp
                else:
                    a_e += comp
        else:
            V += comp
            e += sgn * comp
        for w in range(n):
            if not (mask >> w) & 1:
                if (forward and adj[v][w]) or ((not forward) and adj[w][v]):
                    walk(w, mask | (1 << w), length + 1, forward)
    walk(0, 1, 1, True)
    walk(0, 1, 1, False)
    return H, h, U, V, d, e, a_o, a_e


def all_tournaments(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for bits in range(1 << len(pairs)):
        adj = [[0] * n for _ in range(n)]
        for t, (i, j) in enumerate(pairs):
            if (bits >> t) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
        yield adj


def fuse(A, na, Bm, nb):
    n = na + nb - 1
    adj = [[0] * n for _ in range(n)]
    for i in range(na):
        for j in range(na):
            adj[i][j] = A[i][j]
    def idx(j):
        return 0 if j == 0 else na + j - 1
    for i in range(nb):
        for j in range(nb):
            if i != j and Bm[i][j]:
                adj[idx(i)][idx(j)] = 1
    for i in range(1, na):
        for j in range(1, nb):
            adj[i][na + j - 1] = 1
    return adj, n


profs = {}
for n in range(1, 7):
    lst = []
    par_ok = True
    mass_ok = True
    pairs_mod4 = set()
    for adj in all_tournaments(n):
        H, h, U, V, d, e, a_o, a_e = profile(adj, n)
        par_ok &= (H % 2 == 1 and h % 2 == 1 and (U - d) % 2 == 0 and (V - e) % 2 == 0)
        mass_ok &= (U - h == a_o + a_e and d - h == a_o - a_e and a_o == (U + d) // 2 - h and a_e == (U - d) // 2)
        pairs_mod4.add(((U - h) % 4, (d - h) % 4))
        if n <= 4:
            lst.append((adj, (H, h, U, V, d, e)))
    profs[n] = lst
    check(par_ok and mass_ok, f"order {n}: every labelled pointed tournament has H, h odd, U = d and V = e (mod 2); U - h = a_o + a_e, d - h = a_o - a_e; "
                              f"pairs ((U - h) mod 4, (d - h) mod 4): {sorted(pairs_mod4)}")
viol = 0
npairs = 0
odd_all = True
for na in range(1, 5):
    for nb in range(1, 5):
        for (A, pa) in profs[na]:
            for (Bm, pb) in profs[nb]:
                adj, n = fuse(A, na, Bm, nb)
                Hf = ham_table(adj, n)[(1 << n) - 1]
                H, h, U, V, d, e = pa
                H2, h2, U2, V2, d2, e2 = pb
                pred = h * h2 + (H - h) * h2 + h * (H2 - h2) + (U - h) * (V2 - h2) + (d - h) * (e2 - h2)
                npairs += 1
                odd_all &= (Hf % 2 == 1)
                if Hf != pred:
                    viol += 1
check(viol == 0 and npairs == 5625 and odd_all, f"fusion law on {npairs} ordered pairs of labelled pointed tournaments of orders 1..4: violations {viol}; H(A * B) odd in every case")

print("ALL CHECKS PASSED" if not FAILS else f"{len(FAILS)} FAILURES: {FAILS}")
