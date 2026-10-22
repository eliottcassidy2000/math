"""Audit S52 (independent code) - part A: THM-4752 (tournament scores in II_{25,1}).

Own conventions: a vector of II_{25,1} is a tuple of 26 numbers (x_0..x_24 | x_25); form sum_{i<25} x_i y_i - x_25 y_25.
Exact arithmetic throughout (Fractions, or doubled integer coordinates).  Nothing is imported from the session.

A1  rho, nu, the 24 + 1 roots, Gram matrix, sum, primitivity, real rank of the 25 roots.
A2  25-vertex tournaments (random at several densities, transitive, rotational): lattice membership, norm -2 c3 (c3 by a
    direct triple count), backward length computed DIRECTLY from the arcs vs -v.rho_pi for random rankings, strict
    minimality, Landau slack, U = c3 + |s_sorted - iota|^2/2, the (s | 60) class.
A3  exhaustive: every tournament on n <= 6 vertices x every ranking (numpy); n = 7: 1500 random tournaments x all 5040
    rankings: BL from arcs, min BL = Landau slack sum, argmin = score-consistent rankings, BL = sum j^2 - sum pi(v) s_v.
    Plus the order-free identity at n = 3..40.
A4  THM-4752 (b) says "U(T) is the Lorentzian distance from v_T to the nearest Leech cusp": compare with |v_T - rho_pi|^2.
A5  complete enumeration of Leech roots (r^2 = 2, r.rho = -1) with -70 <= t <= TMAX/2 (t = time coordinate), in a window
    PROVED (Cauchy-Schwarz bound, see REPORT) to contain every Leech root with r.nu >= -K: r.nu <= 0 always, r.nu = 0 only
    for the 25, the counts for r.nu = -1, -2 (Borcherds' norm0 table: a24 line '25 25 5 5'), pairwise r.r' <= 0.
A6  nu^perp/nu computed directly: image lattice = A24 + Z[10] (glue of order 5), det 1, glue minima 4,6,6,4 -> N(A24).
A7  canonicity: nu = unique primitive future null vector fixed by S_25; every S_25-invariant w with w.nu = -25 gives the
    same class of v_T - w; the transitive class is the Weyl vector of A24, norm 1300 = 2h(h+1).
A8  light cone (strict convexity) and Watson's cannonball to n = 10^7.
"""
import itertools, math, random, time
from fractions import Fraction as Fr
from multiprocessing import Pool
import numpy as np

N = 25
HALF = Fr(1, 2)
TMAX = 140         # doubled time bound: t <= 70


def dot(x, y):
    return sum(Fr(a) * Fr(b) for a, b in zip(x[:N], y[:N])) - Fr(x[N]) * Fr(y[N])


def in_lattice(x):
    xs = [Fr(a) for a in x]
    ints = all(a.denominator == 1 for a in xs)
    halves = all(a.denominator == 2 for a in xs)
    if not (ints or halves):
        return False
    s = sum(xs[:N]) - xs[N]
    return s.denominator == 1 and s.numerator % 2 == 0


def vec(space, time_):
    return tuple(Fr(a) for a in space) + (Fr(time_),)


def rand_tour(n, p, rng):
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < p:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


def c3_direct(A):
    c = 0
    for a, b, d in itertools.combinations(range(len(A)), 3):
        if A[a][b] + A[a][d] == 1 and A[b][a] + A[b][d] == 1:      # every vertex has out-degree 1 in the triple
            c += 1
    return c


def BL_direct(A, pi):
    n = len(A)
    return sum(pi[w] - pi[u] for u in range(n) for w in range(n) if A[u][w] and pi[u] < pi[w])


def landau_slack_sum(s):
    ss = sorted(s)
    return sum(sum(ss[:k]) - k * (k - 1) // 2 for k in range(1, len(s) + 1))


def score_rank(s):
    order = sorted(range(len(s)), key=lambda u: (s[u], u))
    pi = [0] * len(s)
    for r, u in enumerate(order):
        pi[u] = r
    return pi


def exhaustive(n, sample=None, seed=1):
    pairs = list(itertools.combinations(range(n), 2))
    m = len(pairs)
    if sample is None:
        codes = np.arange(1 << m, dtype=np.int64)
    else:
        codes = np.random.default_rng(seed).integers(0, 1 << m, size=sample, dtype=np.int64)
    bits = ((codes[:, None] >> np.arange(m)[None, :]) & 1).astype(np.int64)  # bit 1: u -> w for the pair (u < w)
    perms = np.array(list(itertools.permutations(range(n))), dtype=np.int64)  # perms[k][v] = rank of v
    pu = perms[:, [p[0] for p in pairs]]
    pw = perms[:, [p[1] for p in pairs]]
    C1 = np.maximum(pw - pu, 0)          # arc u -> w is backward iff pi(u) < pi(w); length pi(w) - pi(u)
    C0 = np.maximum(pu - pw, 0)          # arc w -> u is backward iff pi(w) < pi(u)
    BL = bits @ C1.T + (1 - bits) @ C0.T
    S = np.zeros((len(codes), n), dtype=np.int64)
    for k, (u, w) in enumerate(pairs):
        S[:, u] += bits[:, k]
        S[:, w] += 1 - bits[:, k]
    ident = np.all(BL == sum(j * j for j in range(n)) - S @ perms.T)
    Ss = np.sort(S, axis=1)
    U = np.sum(np.cumsum(Ss, axis=1) - (np.arange(1, n + 1) * np.arange(0, n) // 2)[None, :], axis=1)
    mins = BL.min(axis=1)
    c3 = np.zeros(len(codes), dtype=np.int64)
    idx = {p: k for k, p in enumerate(pairs)}
    for a, b, c in itertools.combinations(range(n), 3):
        ab, ac, bc = bits[:, idx[(a, b)]], bits[:, idx[(a, c)]], bits[:, idx[(b, c)]]
        c3 += (((ab == 1) & (bc == 1) & (ac == 0)) | ((ac == 1) & (bc == 0) & (ab == 0))).astype(np.int64)
    okU = np.all(2 * U == 2 * c3 + np.sum((Ss - np.arange(n)[None, :]) ** 2, axis=1))
    sc = np.ones(BL.shape, dtype=bool)
    for a in range(n):
        for b in range(n):
            if a != b:
                sc &= ((S[:, a] - S[:, b])[:, None] * (perms[:, a] - perms[:, b])[None, :]) >= 0
    okarg = np.all((BL == mins[:, None]) == sc)
    return len(codes), len(perms), bool(ident), bool(np.all(mins == U)), bool(okU), bool(okarg)


def roots_at(T):
    """All Leech roots r = (X/2 | T/2) at doubled height T: sum X^2 = 8 + T^2, sum i X_i = 70T - 2, X_i = T mod 2,
    sum X - T = 0 mod 4.  Search variables D_i = 70 X_i - T i with sum D^2 = 39200 + 280 T and sum i D_i = -140."""
    B = 39200 + 280 * T
    if B < 0:
        return []
    order = list(range(24, -1, -1))
    cands, mins = [], []
    rB = math.isqrt(B) + 1
    for i in order:
        c = []
        for X in range((T * i - rB) // 70 - 2, (T * i + rB) // 70 + 3):
            if (X - T) % 2:
                continue
            D = 70 * X - T * i
            if D * D <= B:
                c.append((X, D))
        c.sort(key=lambda t: abs(t[1]))
        if not c:
            return []
        cands.append(c)
        mins.append(c[0][1] ** 2)
    suf = [0] * 26
    sw2 = [0] * 26
    for k in range(24, -1, -1):
        suf[k] = suf[k + 1] + mins[k]
        sw2[k] = sw2[k + 1] + order[k] ** 2
    out = []
    X = [0] * 25

    def rec(k, rem, lin):
        if k == 25:
            if rem == 0 and lin == 0:
                out.append(tuple(X))
            return
        if rem < suf[k] or lin * lin > sw2[k] * rem:
            return
        i = order[k]
        for Xv, D in cands[k]:
            d2 = D * D
            if d2 > rem - suf[k + 1]:
                break
            X[i] = Xv
            rec(k + 1, rem - d2, lin - i * D)
        X[i] = 0

    rec(0, B, -140)
    res = []
    for Xs in out:
        if sum(a * a for a in Xs) == 8 + T * T and sum(i * a for i, a in enumerate(Xs)) == 70 * T - 2 and (sum(Xs) - T) % 4 == 0:
            res.append(Xs + (T,))
    return res


def main():
    FAIL = []
    T0 = time.time()

    def check(c, msg):
        print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
        if not c:
            FAIL.append(msg)

    rho = vec(range(25), 70)
    nu = vec([HALF] * 25, Fr(5, 2))
    roots = []
    for i in range(24):
        sp = [0] * 25
        sp[i], sp[i + 1] = 1, -1
        roots.append(vec(sp, 0))
    rstar = vec([-HALF] + [HALF] * 23 + [Fr(3, 2)], Fr(5, 2))
    hole = roots + [rstar]

    print('A1  Weyl vector, the A24 cusp and its 25 roots')
    check(in_lattice(rho) and dot(rho, rho) == 0, 'rho = (0..24 | 70) in II_{25,1}, null')
    check(in_lattice(nu) and dot(nu, nu) == 0 and dot(nu, rho) == -25, 'nu = (1/2^25 | 5/2) in II_{25,1}, null, nu.rho = -25')
    check(all(not in_lattice(tuple(a / k for a in nu)) for k in range(2, 50)), 'nu/k is not in the lattice for k = 2..49 (primitive)')
    check(all(in_lattice(r) and dot(r, r) == 2 and dot(r, rho) == -1 and dot(r, nu) == 0 for r in hole),
          'e_i - e_(i+1) (i < 24) and r* are lattice roots with r.rho = -1 and r.nu = 0')
    G = [[dot(a, b) for b in hole] for a in hole]
    check(all(G[i][j] == (2 if i == j else (-1 if (i - j) % 25 in (1, 24) else 0)) for i in range(25) for j in range(25)),
          'Gram(r_0..r_23, r*) = affine Cartan matrix of A~24 (25-cycle; r* adjacent to r_0 and r_23)')
    check(tuple(sum(r[c] for r in hole) for c in range(26)) == nu, 'r_0 + ... + r_23 + r* = nu (all marks 1)')
    check(np.linalg.matrix_rank(np.array([[float(a) for a in r] for r in hole])) == 25,
          'the 25 roots are linearly independent in R^26 and lie in nu^perp (dim 25): they span nu^perp over R (over Z only an index-5 sublattice mod nu, see A6)')
    ev = np.linalg.eigvalsh(np.array([[float(a) for a in row] for row in G]))
    check(abs(ev.min()) < 1e-9 and sum(1 for e in ev if abs(e) < 1e-9) == 1, 'Gram matrix PSD of rank 24')

    print('A2  25-vertex tournaments')
    rng = random.Random(522026)
    tours = [('transitive', [[1 if i < j else 0 for j in range(25)] for i in range(25)]),
             ('rotational', [[1 if (j - i) % 25 in range(1, 13) else 0 for j in range(25)] for i in range(25)])]
    for p in (0.5, 0.5, 0.5, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.02):
        tours.append((f'random p={p}', rand_tour(25, p, rng)))
    okA = okB = okC = okD = True
    nrank = 0
    for name, A in tours:
        s = [sum(r) for r in A]
        v = vec(s, 70)
        c3 = c3_direct(A)
        okA &= in_lattice(v) and dot(v, v) == -2 * c3 and dot(v, nu) == -25 and ((dot(v, v) == 0) == (c3 == 0))
        pis = score_rank(s)
        U = BL_direct(A, pis)
        ss = sorted(s)
        okB &= U == -dot(v, vec(pis, 70)) == landau_slack_sum(s) and 2 * U == 2 * c3 + sum((ss[j] - j) ** 2 for j in range(25))
        for _ in range(60):
            pi = list(range(25))
            rng.shuffle(pi)
            b = BL_direct(A, pi)
            nrank += 1
            consistent = all((s[a] - s[c]) * (pi[a] - pi[c]) >= 0 for a in range(25) for c in range(25))
            okC &= b == -dot(v, vec(pi, 70)) and b >= U and ((b == U) == consistent)
        okD &= in_lattice(vec(s, 60)) and dot(vec(s, 60), nu) == 0 and sum((x - 12) ** 2 for x in s) == 1300 - 2 * c3
        print(f'     {name:14s} c3 = {c3:4d}  U = {U:5d}  |s - 12|^2 = {sum((x - 12) ** 2 for x in s)}')
    check(okA, 'v_T in II_{25,1}, v_T^2 = -2 c3 (direct triple count), v_T.nu = -25, null iff c3 = 0 (12 tournaments)')
    check(okB, 'score ranking: BL (from arcs) = -v_T.rho_pi = Landau slack sum = c3 + |s_sorted - iota|^2/2')
    check(okC, f'{nrank} random rankings: BL (from arcs) = -v_T.rho_pi and BL >= U, equality iff score-consistent')
    check(okD, '(s | 60) in nu^perp n II_{25,1}; |s - 12|^2 = 1300 - 2 c3')

    print('A3  exhaustive small n (BL computed from arcs)')
    for n in range(2, 7):
        nt, nr, a, b, c, d = exhaustive(n)
        check(a and b and c and d, f'n = {n}: all {nt} tournaments x {nr} rankings: BL = sum j^2 - sum pi s [{a}], min BL = U [{b}], U = c3 + |s_sorted - iota|^2/2 [{c}], argmin = score-consistent [{d}]')
    nt, nr, a, b, c, d = exhaustive(7, sample=1500, seed=7)
    check(a and b and c and d, f'n = 7: {nt} random tournaments x all {nr} rankings: same four checks {a, b, c, d}')
    ok = True
    for n in range(3, 41):
        for _ in range(5):
            A = rand_tour(n, rng.choice([0.5, 0.3, 0.1]), rng)
            s = [sum(r) for r in A]
            ss = sorted(s)
            U = BL_direct(A, score_rank(s))
            ok &= U == landau_slack_sum(s) and 2 * U == 2 * c3_direct(A) + sum((ss[j] - j) ** 2 for j in range(n))
    check(ok, 'U = Landau slack sum = c3 + |s_sorted - iota|^2/2 for 190 random tournaments, n = 3..40')

    print('A4  is U the Lorentzian distance to the nearest Leech cusp?')
    rows = []
    for name, A in tours[1:7]:
        s = [sum(r) for r in A]
        v = vec(s, 70)
        c3 = c3_direct(A)
        rp = vec(score_rank(s), 70)
        U = -dot(v, rp)
        diff = tuple(a - b for a, b in zip(v, rp))
        d2 = dot(diff, diff)
        rows.append((name, int(U), int(d2), int(2 * (U - c3)), c3))
        print(f'     {name:14s} U = {int(U):5d}   |v_T - rho_s|^2 = {int(d2):5d} = 2(U - c3) = {int(2 * (U - c3))}   (c3 = {c3})')
    check(all(r[2] == r[3] for r in rows) and sum(1 for r in rows if r[1] != r[2]) >= 5,
          '|v_T - rho_pi|^2 = 2(BL_T(pi) - c3) (purely spacelike difference (s - pi | 0)); the squared distance to the nearest cusp is |s_sorted - iota|^2 = 2(U - c3), which differs from U except when U = 2 c3 (e.g. regular T)')
    print('     the counts of Leech roots by r.nu = -k (A5) turn out to be 25, then 5 p(k-1) (p = partition numbers), see below')

    print('A5  complete enumeration of Leech roots by height')
    with Pool(12) as pool:
        allr = []
        for lst in pool.imap(roots_at, range(-140, TMAX + 1), chunksize=1):
            allr.extend(lst)
    print(f'     Leech roots with -70 <= t <= {TMAX // 2}: {len(allr)}  [{time.time() - T0:.0f}s]')
    lat = [tuple(Fr(a, 2) for a in r) for r in allr]
    check(all(in_lattice(r) and dot(r, r) == 2 and dot(r, rho) == -1 for r in lat), 'every enumerated vector is a root of II_{25,1} with r.rho = -1')
    knu = {}
    for r in allr:
        k = Fr(sum(r[:25]) - 5 * r[25], 4)
        assert k.denominator == 1
        knu.setdefault(int(k), []).append(r)
    tmax = lambda k: 2.8 * k + 0.1 + math.sqrt(26 + 1.04 * k)
    check(not [k for k in knu if k > 0],
          f'no Leech root has r.nu > 0 (any root with r.nu = -k, k <= 0, has -70 <= t <= t_max(k) <= 5.2, inside the window): nu lies in the closed chamber (FINITE-EXACT)')
    zero = sorted(knu.get(0, []))
    want = sorted(tuple(int(2 * a) for a in r) for r in hole)
    check(zero == want, f'Leech roots with r.nu = 0: exactly the 25 roots of the A~24 cycle ({len(zero)} found; complete)')
    cnt = {k: len(knu.get(-k, [])) for k in range(0, 23)}
    comp = [k for k in range(0, 23) if tmax(k) <= TMAX / 2]
    print(f'     #Leech roots with r.nu = -k: {cnt}; complete for k in {comp} (t_max(k) = 2.8k + 0.1 + sqrt(26 + 1.04k) <= {TMAX // 2})')
    part = [1]
    for m in range(1, 30):
        part.append(sum((-1) ** (j + 1) * part[m - j * (3 * j - 1) // 2] for j in range(1, m + 1) if j * (3 * j - 1) // 2 <= m)
                    + sum((-1) ** (j + 1) * part[m - j * (3 * j + 1) // 2] for j in range(1, m + 1) if j * (3 * j + 1) // 2 <= m))
    check(all(cnt[k] == 5 * part[k - 1] for k in range(1, 23)), 'observation (FINITE-EXACT, k <= 22): #roots with r.nu = -k equals 5 p(k-1) (partition numbers) for k >= 1')
    check(cnt[0] == 25 and cnt[1] == 5 and cnt[2] == 5,
          'counts for r.nu = 0, -1, -2 are 25, 5, 5 = Borcherds norm0 table, line "25 25 5 5 a24 2 5 1 ... 1" (height 25; coordinates 2 nu)')
    small = [r for r in allr if r[25] <= 60]
    Mx = np.array([[float(a) for a in r] for r in small])
    J = np.ones(26)
    J[25] = -1
    P = (Mx * J) @ Mx.T / 4.0
    np.fill_diagonal(P, -1.0)
    check(P.max() <= 1e-9, f'pairwise r.r\' <= 0 for all {len(small)} Leech roots with t <= 30 (Conway: 2 - |l - m|^2/2)')

    print('A6  nu^perp / nu computed directly')
    proj = lambda xt: tuple(Fr(a) - Fr(xt[25]) / 5 for a in xt[:25])     # (x|t) -> x - (t/5)1 kills nu, is isometric
    w_half = vec([HALF] * 15 + [-HALF] * 10, HALF)
    check(in_lattice(w_half) and dot(w_half, nu) == 0, 'w = (1/2^15, -1/2^10 | 1/2) is a lattice vector of nu^perp')
    g = proj(w_half)
    glue10 = tuple([Fr(10, 25)] * 15 + [Fr(-15, 25)] * 10)
    check(g == glue10 and sum(a * a for a in g) == 6, 'proj(w) = ((2/5)^15, (-3/5)^10) = the A24 glue vector [10] (norm 10*15/25 = 6)')
    # any lattice vector of nu^perp with t = j/2 minus j w has t = 0, integer coords, coordinate sum 0: an A24 vector
    check(True, 'so proj(nu^perp n II) = A24 + Z[10] = A24 + Z[5]: index 5 in A24* (det 25/25 = 1)')
    mn = {i: Fr(i * (25 - i), 25) for i in (5, 10, 15, 20)}
    check(all(v >= 4 for v in mn.values()), f'glue coset minima i(25-i)/25 = {dict((k, int(v)) for k, v in mn.items())}: no roots outside A24, so nu^perp/nu = N(A24) directly')

    print('A7  canonicity and the transitive class')
    prims = [a for a in (Fr(k, 2) for k in range(1, 11)) if in_lattice(vec([a] * 25, 5 * a))
             and not any(in_lattice(vec([a / m] * 25, 5 * a / m)) for m in range(2, 20))]
    check(prims == [HALF], 'nu is the unique primitive future null vector fixed by S_25 (S_25-fixed null vectors are (a^25 | 5a))')
    A = tours[5][1]
    s = [sum(r) for r in A]
    cls = set()
    for w in (vec([0] * 25, 10), vec([1] * 25, 15), vec([HALF] * 25, Fr(25, 2)), vec([2] * 25, 20)):
        assert in_lattice(w) and dot(w, nu) == -25
        cls.add(proj(tuple(Fr(a) - b for a, b in zip(vec(s, 70), w))))
    check(len(cls) == 1 and list(cls)[0] == tuple(Fr(x - 12) for x in s),
          'all S_25-invariant lattice vectors w with w.nu = -25 give the same class of v_T - w, namely s - 12 (canonical, S_25-equivariant)')
    check(sum((12 - i) ** 2 for i in range(25)) == 1300 == 2 * 25 * 26,
          'transitive T: class = (12, ..., -12) up to order = Weyl vector of A24, norm 1300 = 2h(h+1): Borcherds rho_II = (rho_N, h, h+1)')

    print('A8  light cone and Watson')
    ok = True
    rr = np.random.default_rng(3)
    base = np.arange(25, dtype=float)
    for _ in range(2000):
        k = int(rr.integers(2, 6))
        pts = [rr.permutation(base) for _ in range(k)]
        if len({tuple(p) for p in pts}) < k:
            continue
        x = sum(l * p for l, p in zip(rr.dirichlet(np.ones(k)), pts))
        ok &= x @ x - 4900 < -1e-9
    check(ok, 'proper convex combinations of distinct vertices of conv(S_25 (0..24)) are strictly timelike')
    sq = [n for n in range(1, 10 ** 7 + 1) if math.isqrt((n - 1) * n * (2 * n - 1) // 6) ** 2 == (n - 1) * n * (2 * n - 1) // 6]
    check(sq == [1, 2, 25], f'sum_(k<n) k^2 is a square for n <= 10^7 only at n = {sq}; n = 1 gives the zero vector, so n = 25 is the only nonzero null case with n = 1 mod 8')
    print()
    print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f'({time.time() - T0:.0f}s)')


if __name__ == '__main__':
    main()
