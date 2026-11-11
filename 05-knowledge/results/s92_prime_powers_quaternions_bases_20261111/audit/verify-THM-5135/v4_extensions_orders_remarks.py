"""Verifier of THM-5135 (S92 audit follow-up): statement (5) (census, groups, algebras, norms, the three orders and
their counts), statement (6) (identity and the sign of the coefficient), Remarks 2 and 3.  Own code, exact.
Run: python v4_extensions_orders_remarks.py > v4_extensions_orders_remarks.out
"""
import itertools
import math
import random
from fractions import Fraction as F

import numpy as np

FAILS = []


def check(cond, msg):
    print(("  ok   " if cond else "  FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


# ---------------------------------------------------------------- [1] the census
print("[1] (5): normalized 2-cocycles on F_2^2 with values in F_2")
V = range(4)
nz = (1, 2, 3)
cochains = 0
cocycles = []
for bits in range(1 << 9):
    cochains += 1
    f = [[0] * 4 for _ in V]
    k = 0
    for a in nz:
        for b in nz:
            f[a][b] = (bits >> k) & 1
            k += 1
    if all((f[a][b] + f[a ^ b][c] + f[b][c] + f[a][b ^ c]) % 2 == 0 for a in V for b in V for c in V):
        cocycles.append(f)
# coboundaries of normalized 1-cochains g (g(0) = 0): (dg)(a, b) = g(a) + g(b) + g(a + b)
cobs = set()
for g1, g2, g3 in itertools.product(range(2), repeat=3):
    g = (0, g1, g2, g3)
    cobs.add(tuple(tuple((g[a] + g[b] + g[a ^ b]) % 2 for b in V) for a in V))
cobs = sorted(cobs)
nontriv = [c for c in cobs if any(any(r) for r in c)]
pairs_changed = [(a, b) for a in nz for b in nz if nontriv[0][a][b]]
check(cochains == 512 and len(cocycles) == 16 and len(cobs) == 2 and len(nontriv) == 1,
      f"{cochains} normalized cochains, {len(cocycles)} cocycles, {len(cobs)} coboundaries (one non-trivial)")
check(sorted(pairs_changed) == [(a, b) for a in nz for b in nz if a != b] and len(pairs_changed) == 6,
      f"the non-trivial coboundary is 1 exactly on the six ordered pairs a != b non-zero: {sorted(pairs_changed)}")
classes = {}
for f in cocycles:
    key = tuple(f[a][b] for a in nz for b in nz)
    g2 = tuple((f[a][b] + nontriv[0][a][b]) % 2 for a in nz for b in nz)
    classes.setdefault(min(key, g2), []).append(f)
check(len(classes) == 8 and all(len(v) == 2 for v in classes.values()), "8 classes, two cocycles each")
supports = sorted(tuple(a for a in nz if fs[0][a][a]) for fs in classes.values())
check(len(set(supports)) == 8, f"the supports S (units of square -1) are the 8 subsets of the three points: {supports}")


def algebra_mult(f, x, y):
    z = [0] * 4
    for a in V:
        if x[a]:
            for b in V:
                if y[b]:
                    z[a ^ b] += x[a] * y[b] * (-1) ** f[a][b]
    return z


def left_matrix(f, x):
    M = np.zeros((4, 4), dtype=object)
    for b in V:
        e = [0] * 4
        e[b] = 1
        col = algebra_mult(f, x, e)
        for a in V:
            M[a, b] = col[a]
    return M


def det4(M):
    M = [[F(int(v)) for v in row] for row in M]
    n = 4
    det = F(1)
    for i in range(n):
        piv = next((r for r in range(i, n) if M[r][i] != 0), None)
        if piv is None:
            return F(0)
        if piv != i:
            M[i], M[piv] = M[piv], M[i]
            det = -det
        det *= M[i][i]
        for r in range(i + 1, n):
            fac = M[r][i] / M[i][i]
            for c in range(i, n):
                M[r][c] -= fac * M[i][c]
    return det


rows = {}
rule_ok = True
polar_ok = True
random.seed(5)
for rep, fs in classes.items():
    f = fs[0]
    S = tuple(a for a in nz if f[a][a])
    elems = [(a, z) for a in V for z in (0, 1)]

    def mul(x, y, f=f):
        return (x[0] ^ y[0], (x[1] + y[1] + f[x[0]][y[0]]) % 2)

    def order(x, mul=mul):
        y, kk = x, 1
        while y != (0, 0):
            y = mul(y, x)
            kk += 1
        return kk
    census = sorted(order(x) for x in elems)
    n4 = census.count(4)
    n8 = census.count(8)
    abel = all(mul(x, y) == mul(y, x) for x in elems for y in elems)
    sq_central = all(mul(x, x)[0] == 0 for x in elems)                    # the square of (a, z) is (0, f(a, a))
    sq_value = all(mul((a, z), (a, z)) == (0, f[a][a]) for a in V for z in (0, 1))
    name = {(True, 0): "C_2^3", (True, 4): "C_4 x C_2", (False, 2): "D_4", (False, 6): "Q_8"}.get((abel, n4), "?")
    # rule of squares and polarization, for both cocycles of the class
    for ff in fs:
        for a, b in itertools.permutations(nz, 2):
            c = a ^ b
            s = lambda t: -1 if ff[t][t] else 1
            eps = 1 if ff[a][b] == ff[b][a] else -1
            rule_ok &= (s(a) * s(b) * s(c) == eps)
            polar_ok &= ((ff[c][c] + ff[a][a] + ff[b][b]) % 2 == (ff[a][b] + ff[b][a]) % 2)
    # sign tables
    kinds = []
    for ff in fs:
        anti = all(ff[a][b] != ff[b][a] for a in nz for b in nz if a < b)
        symm = all(ff[a][b] == ff[b][a] for a in nz for b in nz if a < b)
        if anti:
            out = {a: sum(1 for b in nz if b != a and ff[a][b] == 0) for a in nz}
            if sorted(out.values()) == [1, 1, 1]:
                kinds.append("cyclic")
            else:
                mid = [a for a in nz if out[a] == 1][0]
                kinds.append("transitive, middle vertex in S" if mid in S else "transitive, middle vertex NOT in S")
        elif symm:
            kinds.append("signed triangle")
        else:
            kinds.append("mixed")
    conv_pair = None
    if kinds[0] in ("cyclic",) or kinds[0].startswith("transitive"):
        conv_pair = all(fs[0][a][b] != fs[1][a][b] for a in nz for b in nz if a != b)
    # zero divisors with coordinates in {-1, 0, 1}
    zd = False
    vecs = [v for v in itertools.product((-1, 0, 1), repeat=4) if any(v)]
    for x in vecs:
        for y in vecs:
            if not any(algebra_mult(f, x, y)):
                zd = True
                break
        if zd:
            break
    # algebra norm: det of the left multiplication against the candidate norm forms
    normform_ok = True
    for _ in range(40):
        x = [random.randint(-6, 6) for _ in range(4)]
        dL = det4(left_matrix(f, x))
        sgn = [1] + [(-1 if a in S else 1) for a in nz]                 # e_a^2
        if len(S) % 2 == 1:
            # quaternion algebra (a, b): reduced norm x0^2 - s1 x1^2 - s2 x2^2 - s3 x3^2 with the sign of e_3^2 forced
            red = x[0] ** 2 - sgn[1] * x[1] ** 2 - sgn[2] * x[2] ** 2 - sgn[3] * x[3] ** 2
            normform_ok &= (dL == red ** 2)
    rows[S] = (name, n4, n8, abel, sq_central and sq_value, kinds, conv_pair, zd, normform_ok)
check(rule_ok, "rule of squares s_a s_b s_c = eps (eps = +1 iff e_a, e_b commute), all 16 cocycles and all pairs")
check(polar_ok, "additive form: q(a + b) + q(a) + q(b) = commutator (polar form of the squaring map), all 16 cocycles")
for S in sorted(rows, key=lambda t: (len(t), t)):
    name, n4, n8, abel, sqc, kinds, conv_pair, zd, nf = rows[S]
    print(f"       S = {S}: group {name}, elements of order 4: {n4}, of order 8: {n8}, commutative {abel}, sign tables {kinds}, "
          f"converse pair {conv_pair}, zero divisors in the box {zd}")
by = {}
for S, row in rows.items():
    by.setdefault(len(S), set()).add((row[0], row[1], row[3], row[7]))
check(by == {0: {("C_2^3", 0, True, True)}, 1: {("D_4", 2, False, True)}, 2: {("C_4 x C_2", 4, True, True)}, 3: {("Q_8", 6, False, False)}},
      "table of (5): #S = 0, 1, 2, 3 -> C_2^3, D_4, C_4 x C_2, Q_8 with 0, 2, 4, 6 elements of order 4; commutative iff #S even; zero divisors except for #S = 3")
check(all(r[2] == 0 and r[4] for r in rows.values()), "no element of order 8: the square of (a, z) is (0, f(a, a)) in every class (exponent divides 4)")
check(all(("cyclic" in r[5][0]) == (len(S) == 3) and (r[5][0].startswith("transitive, middle vertex in S")) == (len(S) == 1)
          and (r[5][0] == "signed triangle") == (len(S) % 2 == 0) for S, r in rows.items())
      and all(r[6] is True for S, r in rows.items() if len(S) % 2 == 1),
      "sign tables: cyclic for #S = 3, transitive with middle vertex S for #S = 1 (the two cocycles are converse), signed triangles for #S even")
check(all(r[8] for r in rows.values()), "non-commutative classes: det L_x = (reduced norm)^2, reduced norm x0^2 + x1^2 + x2^2 + x3^2 for Q_8 and "
                                        "x0^2 + (one square) - (two squares) for D_4")

# ---------------------------------------------------------------- [2] C + C: idempotent coordinates and the three orders
print("[2] (5): norms and orders")
# class S = {1, 2}: i = e_1, j = e_2, ij = e_3 with (ij)^2 = +1
fcc = [fs for rep, fs in classes.items() if tuple(a for a in nz if fs[0][a][a]) == (1, 2)][0][0]
sign3 = 1 if fcc[1][2] == 0 else -1          # e_1 e_2 = sign3 * e_3
ok_norm = True
ok_coord = True
for _ in range(300):
    x = [random.randint(-7, 7) for _ in range(4)]
    # g + g' j with g = x0 + x1 i, g' = x2 + sign3 * x3 i   (g' j = x2 j + sign3 x3 (i j) and i j = sign3 e_3)
    g = complex(x[0], x[1])
    gp = complex(x[2], sign3 * x[3])
    z1 = g + 1j * gp
    z2 = g - 1j * gp
    n = round(abs(z1) ** 2) * round(abs(z2) ** 2)
    ok_norm &= (det4(left_matrix(fcc, x)) == n)
    ok_coord &= (round((z1 - z2).real) % 2 == 0 and round((z1 - z2).imag) % 2 == 0)
check(ok_norm, "C + C (S = {1, 2}): det L_x = |g + i g'|^2 |g - i g'|^2 for x = g + g' j (the algebra norm, degree 4)")
check(ok_coord, "integer coordinates in 1, i, j, ij give z_1 = z_2 (mod 2) in Z[i] x Z[i]")
idx = round(abs(np.linalg.det(np.array([[1, 0, 0, -1], [0, 1, 1, 0], [1, 0, 0, 1], [0, 1, -1, 0]], dtype=float))))
check(idx == 4, f"index of the Z-span of 1, i, j, ij in Z[i] x Z[i]: {idx}")

LIM = 220
r2 = [0] * (LIM + 1)
Bd = int(LIM ** 0.5) + 1
gauss = {}
for a in range(-Bd, Bd + 1):
    for b in range(-Bd, Bd + 1):
        n = a * a + b * b
        if n <= LIM:
            r2[n] += 1
            gauss.setdefault(n, []).append((a, b))
sig1 = [0] * (LIM + 1)
for d in range(1, LIM + 1):
    for kk in range(d, LIM + 1, d):
        sig1[kk] += d
# Lipschitz order: direct count on Z^4 for n <= 60
lip = [0] * 61
for a, b, c, d in itertools.product(range(-8, 9), repeat=4):
    n = a * a + b * b + c * c + d * d
    if n <= 60:
        lip[n] += 1
check(all(lip[n] == sum(r2[a] * r2[n - a] for a in range(n + 1)) for n in range(1, 61)) and all(lip[n] == 8 * sig1[n] for n in range(1, 61, 2)) and lip[1] == 8,
      "Lipschitz order: direct count on Z^4 = sum_(a + b = n) r_2(a) r_2(b) for n <= 60, = 8 sigma_1(n) at odd n; 8 units")
# maximal order Z[i] x Z[i] and the suborder z_1 = z_2 (mod 2): direct enumeration of pairs
Mcount = [0] * (LIM + 1)
Bcount = [0] * (LIM + 1)
for n in range(1, LIM + 1):
    for a in range(1, n + 1):
        if n % a == 0 and a in gauss and (n // a) in gauss:
            for (u1, v1) in gauss[a]:
                for (u2, v2) in gauss[n // a]:
                    Mcount[n] += 1
                    if (u1 - u2) % 2 == 0 and (v1 - v2) % 2 == 0:
                        Bcount[n] += 1
check(all(Mcount[n] == sum(r2[a] * r2[n // a] for a in range(1, n + 1) if n % a == 0) for n in range(1, LIM + 1)) and Mcount[1] == 16,
      f"maximal order Z[i] x Z[i]: count = sum_(a b = n) r_2(a) r_2(b), n <= {LIM}; 16 units")
check(all(2 * Bcount[n] == Mcount[n] for n in range(1, LIM + 1, 2)) and Bcount[1] == 8,
      f"suborder z_1 = z_2 (mod 2): half of the count at every odd n <= {LIM}; 8 units")
check(all(Bcount[n] == 0 for n in range(2, LIM + 1, 4)), f"suborder: 0 at every n = 2 (mod 4), n <= {LIM}")
check(Bcount[4] == 16 and Mcount[4] == 48, f"n = 4: {Bcount[4]} against {Mcount[4]}")
check(all(Mcount[n] == 0 and Bcount[n] == 0 for n in range(3, LIM + 1, 4)), "both counts vanish at n = 3 (mod 4)")
# the suborder count directly in the coordinates (x0, x1, x2, x3) of the table
direct = [0] * 41
for x0, x1, x2, x3 in itertools.product(range(-9, 10), repeat=4):
    n = ((x0 - x3) ** 2 + (x1 + x2) ** 2) * ((x0 + x3) ** 2 + (x1 - x2) ** 2)
    if 1 <= n <= 40:
        direct[n] += 1
check(direct[1:21] == Bcount[1:21], f"the same count in the coordinates of the table, n <= 20: {direct[1:21]}")
print(f"       Z[i] x Z[i], n = 1..20: {Mcount[1:21]}")

# norm form of the unit-spanned order of M_2 against the seam
def table(q):
    d = [0] * q
    for a in range(q):
        for b in range(q):
            d[(a * b) % q] += 1
    return d


def conv(a, b):
    q = len(a)
    out = [0] * q
    for i in range(q):
        for j in range(q):
            out[(i + j) % q] += a[i] * b[j]
    return out


res = {}
for q in (3, 9, 5, 7, 2, 4, 8, 16):
    sq = [0] * q
    for a in range(q):
        sq[(a * a) % q] += 1
    neg = [0] * q
    for a in range(q):
        neg[(-a * a) % q] += 1
    form = conv(conv(sq, sq), conv(neg, neg))
    seam = conv(table(q), table(q))
    res[q] = (form == seam)
check(all(res[q] for q in (3, 9, 5, 7)) and not any(res[q] for q in (2, 4, 8, 16)),
      f"x0^2 + x1^2 - x2^2 - x3^2 against a1 a2 + b1 b2: equal modulo 3, 9, 5, 7, different modulo 2, 4, 8, 16: {res}")

# ---------------------------------------------------------------- [3] statement (6)
print("[3] (6): Brahmagupta identity and the coefficient of sqrt(-2^s)")
import sympy as sp
xA, yA, xB, yB, t = sp.symbols("xA yA xB yB t")
ident = sp.expand((xA ** 2 + t * yA ** 2) * (xB ** 2 + t * yB ** 2) - (xA * xB + t * yA * yB) ** 2 - t * (xA * yB - xB * yA) ** 2)
check(ident == 0, "(x^2 + t y^2)(x'^2 + t y'^2) = (x x' + t y y')^2 + t (x y' - x' y)^2 (polynomial identity)")
w = sp.symbols("w")
prod = sp.expand((xA + w * yA) * (xB - w * yB)).subs(w ** 2, -t)
prod = sp.expand(prod)
real_part = prod.coeff(w, 0)
coef_w = prod.coeff(w, 1)
check(sp.expand(real_part - (xA * xB + t * yA * yB)) == 0 and sp.expand(coef_w + (xA * yB - xB * yA)) == 0,
      "z_A conj(z_B) = M + w (x_B y_A - x_A y_B), w^2 = -t: real part M, coefficient of w = minus the cross-determinant")


def psi_table(n):
    # q(h) = h - psi(h): q(0) = 0, q(2h) = q(h), q(2h + 1) = q(h) + ceil(h/2)   (the recursion used by the audit, a6)
    q = [0] * n
    for h in range(1, n):
        q[h] = q[h // 2] if h % 2 == 0 else q[h // 2] + ((h // 2) + 1) // 2
    return [h - q[h] for h in range(n)]


psi = psi_table(300)
states = [(2 * psi[h] + 1, 2 * h + 1) for h in range(300)]
check(states[:8] == [(1, 1), (3, 3), (5, 5), (5, 7), (9, 9), (9, 11), (11, 13), (9, 15)], f"one-apex states (x, y), h < 8: {states[:8]}")
bad = 0
cases = 0
for s in range(0, 7):
    for (x1, y1) in states[::3]:
        for (x2, y2) in states[::5]:
            cases += 1
            if (x1 * x1 + 2 ** s * y1 * y1) * (x2 * x2 + 2 ** s * y2 * y2) != (x1 * x2 + 2 ** s * y1 * y2) ** 2 + 2 ** s * (x1 * y2 - x2 * y1) ** 2:
                bad += 1
check(bad == 0, f"identity on one-apex states: {cases} cases (s <= 6), failures {bad}")

# ---------------------------------------------------------------- [4] Remarks 2 and 3
print("[4] Remarks 2 and 3")
N4 = [[5, 3, 0], [3, 7, 2], [0, 2, 2]]
N5 = [[15, 13, 2], [13, 14, 4], [2, 4, 2]]
d4 = int(sp.Matrix(N4).det())
d5 = int(sp.Matrix(N5).det())
check(d4 == 32 and d5 == -6, f"determinants of the two 3 x 3 visit matrices: {d4}, {d5}")


def covers(adj, n):
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
                for wv in range(n):
                    if not (mask >> wv) & 1 and adj[v][wv]:
                        dp[mask | (1 << wv)][wv] += c
        hp[mask] = tot
    fullm = full - 1
    H = hp[fullm]
    Pi2 = sum(hp[S] * hp[fullm ^ S] for S in range(full))
    P2 = sum(hp[S] * hp[fullm ^ S] for S in range(1, fullm))
    P3 = 0
    for S in range(1, fullm):
        rest = fullm ^ S
        T = rest
        sub = (T - 1) & T
        while sub:
            P3 += hp[S] * hp[sub] * hp[rest ^ sub]
            sub = (sub - 1) & T
    return H, Pi2, P2, P3


def all_tournaments(n):
    prs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for bits in range(1 << len(prs)):
        adj = [[0] * n for _ in range(n)]
        for tt, (i, j) in enumerate(prs):
            if (bits >> tt) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
        yield adj


random.seed(4756)
tot = 0
ok3 = True
for n in range(1, 7):
    gen = all_tournaments(n) if n <= 5 else None
    if gen is None:
        def sample(n=n):
            prs = [(i, j) for i in range(n) for j in range(i + 1, n)]
            for _ in range(1500):
                adj = [[0] * n for _ in range(n)]
                for (i, j) in prs:
                    if random.random() < 0.5:
                        adj[i][j] = 1
                    else:
                        adj[j][i] = 1
                yield adj
        gen = sample()
    for adj in gen:
        H, Pi2, P2, P3 = covers(adj, n)
        tot += 1
        ok3 &= (P2 == Pi2 - 2 * H)
        ok3 &= ((Pi2 - (-1) ** n * H) % 3 == 0)
        ok3 &= ((P2 - ((-1) ** n + 1) * H) % 3 == 0)
        ok3 &= (P2 % 2 == 0 and P3 % 6 == 0)
check(ok3, f"{tot} tournaments (all of order <= 5, 1500 random of order 6): P_2 = Pi_2 - 2 H; Pi_2 = (-1)^n H and P_2 = ((-1)^n + 1) H (mod 3); "
           f"P_2 even, P_3 divisible by 6")

print("ALL CHECKS PASSED" if not FAILS else f"{len(FAILS)} FAILURES: {FAILS}")
