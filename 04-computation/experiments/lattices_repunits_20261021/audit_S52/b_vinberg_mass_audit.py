"""Audit S52 (independent code) - part B: lattice_cusps.py claims.

B1  Vinberg's algorithm for II_{9,1} and II_{17,1}, written from scratch with the OPPOSITE sign convention to the
    session (polytope P = {x : x.r >= 0}, new roots have negative time coordinate, D_n chamber |x_1| <= x_2 <= ... <= x_n),
    exact integer arithmetic (doubled coordinates), exact positive-definiteness by Bareiss determinants.  Termination is
    certified by Vinberg's finite-volume criterion (1985 survey, Prop. 4.2) implemented independently.  Cusps =
    parabolic subdiagrams of rank n-1; null vectors from exact kernels.  Weyl vector of E10.
B2  Siegel mass formula in ranks 8, 16, 24 with my own Bernoulli numbers, Weyl group orders and Niemeier table
    (|G1|, |G2| from Conway-Sloane SPLAG Table 16.1 as I recall it), plus the M11-for-A2^12 variant.
B3  E8 holes by brute-force enumeration of all E8 vectors of norm <= 6; Milnor's pair E8+E8 / D16+ through norm 8 from
    Jacobi theta series (independent of the session's closed-form counts).
"""
import itertools, math, time
from fractions import Fraction as Fr
import numpy as np

FAIL = []
T0 = time.time()


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


# ---------------------------------------------------------------- exact linear algebra on small integer matrices
def bareiss_det(M):
    A = [list(r) for r in M]
    n = len(A)
    if n == 0:
        return 1
    sign, prev = 1, 1
    for k in range(n - 1):
        if A[k][k] == 0:
            sw = next((i for i in range(k + 1, n) if A[i][k] != 0), None)
            if sw is None:
                return 0
            A[k], A[sw] = A[sw], A[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) // prev
        prev = A[k][k]
    return sign * A[n - 1][n - 1]


def is_pd(G, idx):
    idx = list(idx)
    return all(bareiss_det([[G[a][b] for b in idx[:m]] for a in idx[:m]]) > 0 for m in range(1, len(idx) + 1))


def comps(G, idx):
    idx = list(idx)
    left, out = set(idx), []
    while left:
        a = left.pop()
        comp, st = [a], [a]
        while st:
            u = st.pop()
            for b in list(left):
                if G[u][b] != 0:
                    left.discard(b)
                    comp.append(b)
                    st.append(b)
        out.append(sorted(comp))
    return out


def is_affine_component(G, comp):
    if len(comp) < 2:
        return False
    if bareiss_det([[G[a][b] for b in comp] for a in comp]) != 0:
        return False
    return is_pd(G, comp[1:])          # proper subdiagrams of an affine diagram are elliptic; with det 0 -> PSD nullity 1


def kernel_vector(G, comp):
    """Exact kernel of the (affine) Gram submatrix, scaled to a primitive positive integer vector."""
    m = len(comp)
    A = [[Fr(G[a][b]) for b in comp] for a in comp]
    # solve with last coordinate = 1 using the PD leading block
    B = [row[:-1] for row in A[:-1]]
    rhs = [-row[-1] for row in A[:-1]]
    n = m - 1
    M = [B[i][:] + [rhs[i]] for i in range(n)]
    for c in range(n):
        p = next(i for i in range(c, n) if M[i][c] != 0)
        M[c], M[p] = M[p], M[c]
        for i in range(n):
            if i != c and M[i][c] != 0:
                f = M[i][c] / M[c][c]
                M[i] = [x - f * y for x, y in zip(M[i], M[c])]
    sol = [M[i][n] / M[i][i] for i in range(n)] + [Fr(1)]
    den = 1
    for x in sol:
        den = den * x.denominator // math.gcd(den, x.denominator)
    v = [int(x * den) for x in sol]
    g = 0
    for x in v:
        g = math.gcd(g, abs(x))
    v = [x // g for x in v]
    if v[0] < 0:
        v = [-x for x in v]
    return v


# ---------------------------------------------------------------- II_{n,1} in doubled coordinates
# a vector is (X_1..X_n | X_0) with X = 2x; form x.y = (sum X_i Y_i - X_0 Y_0)/4; lattice: all X even or all odd,
# (sum X_i - X_0)/2 even.
def ip(X, Y):
    v = sum(a * b for a, b in zip(X[:-1], Y[:-1])) - X[-1] * Y[-1]
    assert v % 4 == 0
    return v // 4


def in_II(X):
    par = {a % 2 for a in X}
    if len(par) != 1:
        return False
    s = sum(X[:-1]) - X[-1]
    return s % 4 == 0


def anti_dominant(n, K):
    """Doubled vectors Y (all = K mod 2) with sum Y^2 = 8 + K^2 and Y_n <= ... <= Y_2 <= -|Y_1| (my convention)."""
    target = 8 + K * K
    par = K % 2
    out = []
    # write Y_i = -Z_i with Z_n >= ... >= Z_2 >= |Z_1| >= 0, Z_i = par mod 2
    vals = [z for z in range(0, math.isqrt(target) + 2) if z % 2 == par]

    def rec(pos, prefix, rem, cap):
        # fill Z_n, Z_(n-1), ..., Z_2 (non-increasing), then Z_1 with |Z_1| <= Z_2
        if pos == 1:
            for z in vals:
                if z <= cap and z * z == rem:
                    for sgn in ((1, -1) if z else (1,)):
                        Z = [sgn * z] + prefix[::-1]
                        out.append(tuple(-t for t in Z))
            return
        for z in vals:
            if z > cap or z * z > rem:
                continue
            if rem - z * z < (pos - 1) * vals[0] ** 2:
                continue
            rec(pos - 1, prefix + [z], rem - z * z, z)

    rec(n, [], target, vals[-1])
    return out


def vinberg(n, Kmax=12):
    roots = []
    # D_n simple roots for the chamber |x_1| <= x_2 <= ... <= x_n, convention P = {x : x.r >= 0}
    r = [0] * (n + 1); r[0] = 2; r[1] = 2; roots.append(tuple(r))            # e_1 + e_2
    for i in range(1, n):
        r = [0] * (n + 1); r[i] = 2; r[i - 1] = -2; roots.append(tuple(r))   # e_(i+1) - e_i
    assert all(in_II(r) and ip(r, r) == 2 for r in roots)
    levels = []
    for K in range(1, Kmax + 1):            # time coordinate -K/2 (doubled -K)
        cand = []
        for Y in anti_dominant(n, K):
            X = Y + (-K,)
            if not in_II(X):
                continue
            if ip(X, X) != 2:
                continue
            if all(ip(X, s) <= 0 for s in roots):
                cand.append(X)
        # order-independence: all accepted-at-this-level candidates must be pairwise <= 0
        pair_ok = all(ip(a, b) <= 0 for a, b in itertools.combinations(cand, 2))
        roots += cand
        levels.append((K, len(cand), pair_ok))
        G = [[ip(a, b) for b in roots] for a in roots]
        fv, par = finite_volume(G, n)
        if fv:
            return roots, G, par, levels
    return roots, None, None, levels


def finite_volume(G, n):
    m = len(G)
    nodes = range(m)
    parab = []
    for size in range(n, m + 1):
        for Q in itertools.combinations(nodes, size):
            cs = comps(G, Q)
            if all(is_affine_component(G, c) for c in cs) and size - len(cs) == n - 1:
                parab.append((Q, cs))
    ell_n = [S for S in itertools.combinations(nodes, n) if is_pd(G, S)]
    if not parab and not ell_n:
        return False, parab
    for S in itertools.combinations(nodes, n - 1):
        if not is_pd(G, S):
            continue
        ext = sum(1 for v in nodes if v not in S and is_pd(G, S + (v,)))
        ext += sum(1 for Q, _ in parab if set(S) <= set(Q))
        if ext != 2:
            return False, parab
    return True, parab


def affine_name(G, comp):
    m = len(comp)
    deg = {a: sum(1 for b in comp if b != a and G[a][b] != 0) for a in comp}
    if any(abs(G[a][b]) >= 2 for a in comp for b in comp if a != b):
        return 'A~1' if m == 2 else '?'
    if all(d == 2 for d in deg.values()):
        return f'A~{m - 1}'
    b3 = [a for a in comp if deg[a] == 3]
    if len(b3) == 2:
        return f'D~{m - 1}'
    if len(b3) == 1:
        c = b3[0]
        arms = []
        for b in comp:
            if b != c and G[c][b] != 0:
                L, prev, cur = 1, c, b
                while True:
                    nx = [x for x in comp if x not in (prev, cur) and G[cur][x] != 0]
                    if not nx:
                        break
                    prev, cur = cur, nx[0]
                    L += 1
                arms.append(L)
        return {(1, 1, m - 3): f'D~{m - 1}', (2, 2, 2): 'E~6', (1, 3, 3): 'E~7', (1, 2, 5): 'E~8'}.get(tuple(sorted(arms)), str(arms))
    return '?'


print('B1  Vinberg diagrams of II_{9,1}, II_{17,1} (independent implementation)')
for n, exp_roots, exp_cusps in ((9, 10, [['E~8']]), (17, 19, [['D~16'], ['E~8', 'E~8']])):
    roots, G, par, levels = vinberg(n)
    check(G is not None and len(roots) == exp_roots and all(lv[2] for lv in levels),
          f'II_{{{n},1}}: {len(roots)} simple roots; levels (2|time|, #new, same-level pairwise <= 0) = {levels}; finite volume certified (Vinberg criterion)')
    names = sorted(sorted(affine_name(G, c) for c in cs) for _, cs in par)
    check(names == sorted(exp_cusps), f'II_{{{n},1}}: cusps (parabolic subdiagrams of rank {n - 1}) = {names}')
    for Q, cs in par:
        zs = []
        for c in cs:
            marks = kernel_vector(G, c)
            z = tuple(sum(mk * roots[a][i] for mk, a in zip(marks, c)) for i in range(n + 1))
            zs.append((affine_name(G, c), sum(marks), z))
        same = len({z for _, _, z in zs}) == 1
        z = zs[0][2]
        zf = tuple(-a for a in z) if z[-1] < 0 else z
        # primitive in the lattice?  z (doubled coords) / 2 = actual vector
        actual = tuple(Fr(a, 2) for a in zf)
        check(same and ip(z, z) == 0, f'II_{{{n},1}} cusp {[t for t, _, _ in zs]}: marks sum = {[h for _, h, _ in zs]}, common null vector (future) = ({",".join(str(a) for a in actual[:-1])} | {actual[-1]})')
    if n == 9:
        A = np.array(G, dtype=float)
        c = np.linalg.solve(A, np.ones(len(roots)))          # rho . r = +1 in my convention (rho in P)
        rho2 = float(np.ones(len(roots)) @ c)
        Gi = [[Fr(x) for x in row] for row in G]
        check(abs(rho2 + 1240) < 1e-7 and round(np.linalg.det(A)) == -1, f'E10 Gram det = {round(np.linalg.det(A))}, Weyl vector norm rho^2 = {rho2:.6f} (expected -1240)')
        z = par[0][1][0]
        marks = kernel_vector(G, z)
        check(sum(marks) == 30, f'rho . z = sum of marks of E~8 = {sum(marks)} = h(E8) (up to the sign convention)')

print('B2  Siegel mass formula')


def bernoulli(m):
    B = [Fr(1)]
    for k in range(1, m + 1):
        B.append(-sum(math.comb(k + 1, j) * B[j] for j in range(k)) / (k + 1))
    return B


Bn = bernoulli(30)


def mass(dim):
    k = dim // 2
    M = abs(Bn[k]) / dim
    for j in range(1, k):
        M *= abs(Bn[2 * j]) / (4 * j)
    return M


fac = math.factorial
WE = {6: 51840, 7: 2903040, 8: 696729600}


def Wo(t, r):
    return fac(r + 1) if t == 'A' else (2 ** (r - 1) * fac(r) if t == 'D' else WE[r])


def hcox(t, r):
    return r + 1 if t == 'A' else (2 * r - 2 if t == 'D' else {6: 12, 7: 18, 8: 30}[r])


def nroot(t, r):
    return r * (r + 1) if t == 'A' else (2 * r * (r - 1) if t == 'D' else {6: 72, 7: 126, 8: 240}[r])


# (root system as string, |G1|, |G2|), my recollection of SPLAG Table 16.1
TABLE = [('D24', 1, 1), ('D16 E8', 1, 1), ('E8^3', 1, 6), ('A24', 2, 1), ('D12^2', 1, 2), ('A17 E7', 2, 1),
         ('D10 E7^2', 1, 2), ('A15 D9', 2, 1), ('D8^3', 1, 6), ('A12^2', 2, 2), ('A11 D7 E6', 2, 1), ('E6^4', 2, 24),
         ('A9^2 D6', 2, 2), ('D6^4', 1, 24), ('A8^3', 2, 6), ('A7^2 D5^2', 2, 4), ('A6^4', 2, 12), ('A5^4 D4', 2, 24),
         ('D4^6', 3, 720), ('A4^6', 2, 120), ('A3^8', 2, 1344), ('A2^12', 2, 95040), ('A1^24', 1, 244823040)]
CO0 = 2 ** 22 * 3 ** 9 * 5 ** 4 * 7 ** 2 * 11 * 13 * 23


def parse(rs):
    out = []
    for tok in rs.split():
        base, _, mult = tok.partition('^')
        out.append((base[0], int(base[1:]), int(mult or 1)))
    return out


def niemeier_mass(table):
    tot = Fr(1, CO0)
    good = True
    for rs, g1, g2 in table:
        R = parse(rs)
        w = 1
        for t, r, m in R:
            w *= Wo(t, r) ** m
        tot += Fr(1, w * g1 * g2)
        hs = {hcox(t, r) for t, r, m in R}
        good &= len(hs) == 1 and sum(r * m for t, r, m in R) == 24 and sum(nroot(t, r) * m for t, r, m in R) == 24 * min(hs)
    return tot, good


check(mass(8) == Fr(1, WE[8]), f'rank 8: Siegel mass = 1/{1 / mass(8)} = 1/|W(E8)|')
m16 = Fr(1, 2 * WE[8] ** 2) + Fr(1, 2 ** 15 * fac(16))
check(mass(16) == m16, f'rank 16: Siegel mass {mass(16)} = 1/|Aut(E8^2)| + 1/|Aut(D16+)| (|Aut(D16+)| = 2^15 16!)')
tot, good = niemeier_mass(TABLE)
check(good and CO0 == 8315553613086720000, '23 root systems: rank 24, single Coxeter number h, 24h roots; |Co_0| = 8315553613086720000')
check(tot == mass(24), f'rank 24: sum of 1/|Aut| over the 24 Niemeier lattices = Siegel mass exactly ({mass(24)})')
bad = [(rs, g1, (7920 if rs == 'A2^12' else g2)) for rs, g1, g2 in TABLE]
tb, _ = niemeier_mass(bad)
check(tb != mass(24), 'with |G2(A2^12)| = |M11| = 7920 instead of |M12| = 95040 the mass identity fails (M12 is forced)')
# sensitivity: perturbing any single |G1| or |G2| by a factor 2 breaks the identity
sens = True
for i in range(len(TABLE)):
    for which in (1, 2):
        T2 = [list(r) for r in TABLE]
        T2[i][which] *= 2
        sens &= niemeier_mass([tuple(r) for r in T2])[0] != mass(24)
check(sens, 'doubling any single |G1| or |G2| breaks the identity (the exact mass check pins every table entry it uses)')

print('B3  E8 holes and Milnor\'s pair')
vecs = []
for x in itertools.product(range(-2, 3), repeat=8):
    if sum(x) % 2 == 0 and sum(t * t for t in x) <= 6:
        vecs.append(tuple(Fr(t) for t in x))
for x in itertools.product((-5, -3, -1, 1, 3, 5), repeat=8):
    if (sum(x) // 2) % 2 == 0 and sum(t * t for t in x) <= 24:
        vecs.append(tuple(Fr(t, 2) for t in x))
cnt = {}
for v in vecs:
    k = sum(t * t for t in v)
    cnt[k] = cnt.get(k, 0) + 1
check(cnt.get(0) == 1 and cnt.get(2) == 240 and cnt.get(4) == 2160 and cnt.get(6) == 6720,
      f'E8 shells by brute force: {dict(sorted((int(k), v) for k, v in cnt.items()))} (1, 240, 2160, 6720 = 240 sigma_3(n))')
for name, c, R2, npts in (('shallow', tuple([Fr(5, 6)] + [Fr(1, 6)] * 7), Fr(8, 9), 9),
                         ('deep', tuple([Fr(1)] + [Fr(0)] * 7), Fr(1), 16)):
    d = [sum((a - b) ** 2 for a, b in zip(v, c)) for v in vecs]
    mn = min(d)
    near = [v for v, x in zip(vecs, d) if x == mn]
    cen = tuple(sum(v[i] for v in near) / len(near) for i in range(8))
    check(mn == R2 and len(near) == npts and cen == c,
          f'{name} hole {tuple(str(t) for t in c)}: squared distance {mn}, {len(near)} nearest points, centroid = hole: {cen == c} (all lattice points of norm <= 6 searched; enough since |c| < 1)')
    if name == 'deep':
        anti = all(tuple(2 * ci - vi for ci, vi in zip(c, v)) in set(near) for v in near)
        check(anti, 'deep-hole vertices come in 8 antipodal pairs about the hole: a cross-polytope')
    else:
        dd = {sum((a - b) ** 2 for a, b in zip(u, v)) for u, v in itertools.combinations(near, 2)}
        check(dd == {2}, 'shallow-hole vertices pairwise at squared distance 2: regular simplex')
# Milnor via Jacobi thetas in powers of q^(1/4) (exponent = 4 * norm)
L = 4 * 8 + 1


def series_mul(a, b):
    c = [0] * L
    for i, x in enumerate(a):
        if x:
            for j in range(L - i):
                c[i + j] += x * b[j]
    return c


def power(a, k):
    r = [1] + [0] * (L - 1)
    for _ in range(k):
        r = series_mul(r, a)
    return r


th3 = [0] * L; th4 = [0] * L; th2 = [0] * L
for m in range(-6, 7):
    if 4 * m * m < L:
        th3[4 * m * m] += 1
        th4[4 * m * m] += (-1) ** abs(m)
for m in range(-7, 7):
    e = (2 * m + 1) ** 2      # 4 (m + 1/2)^2
    if e < L:
        th2[e] += 1
E8th = [(a + b + c) // 2 for a, b, c in zip(power(th3, 8), power(th4, 8), power(th2, 8))]
D16p = [(a + b + c) // 2 for a, b, c in zip(power(th3, 16), power(th4, 16), power(th2, 16))]
E8sq = series_mul(E8th, E8th)
by_norm = lambda s: [s[4 * k] for k in range(0, 9)]
check(by_norm(E8sq) == by_norm(D16p) and by_norm(E8sq)[2] == 480 and by_norm(E8sq)[4] == 61920,
      f'theta(E8)^2 = theta(D16+) through norm 8: {by_norm(E8sq)[:9:2]} (480 roots, 61920 norm-4 vectors) - Milnor\'s pair')
check(all(x == 0 for i, x in enumerate(E8sq) if i % 8) and all(x == 0 for i, x in enumerate(D16p) if i % 8), 'both theta series are supported on even norms (even lattices)')
print()
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f'({time.time() - T0:.0f}s)')
