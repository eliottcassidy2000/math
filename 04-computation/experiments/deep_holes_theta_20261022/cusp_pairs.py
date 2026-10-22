r"""HYP-9368 test: is the swap theorem special to the Leech cusp? (session opus-2026-10-22-S53)

Two cusps of II_{25,1} are taken from Borcherds' norm-0 table: z1 and z2, primitive null vectors adjacent to rho.

Setting.
- m = -(z1 . z2) > 0. M = {z1, z2}^perp has discriminant group (Z/m)^2 when span(z1, z2) is primitive; this is
  checked.
- The two isotropic lines come from the two cusps: H_1 = proj(z1^perp), H_2 = proj(z2^perp).

Tested statement (same certificate as swap_theorem_proof.py).
- Theta_M - sigma^* Theta_M = Delta v for the swap sigma of the two lines (normalised so that it maps the
  z1-line {x . z1 = -1} onto the z2-line).
- Consequently theta(L1 + p1) - theta(L2 + p2) = (sum over the z1-line of v) Delta, where the coset L_i + p_i is the
  projection of {x . z_i = -1}.

Certificate: M^* vectors with |x|^2 <= 2 by class; E = 0 below norm 2; v Weil-invariant (exact).
Needs cypari, numba, sympy. Usage: python cusp_pairs.py name1 name2 [name1 name2 ...]"""
import sys, time, math
import numpy as np
from fractions import Fraction
from cypari import pari
from numba import njit
from sympy import Poly, symbols, cyclotomic_poly, rem, ZZ

PAIRS = list(zip(sys.argv[1::2], sys.argv[2::2])) or [('a1^24', 'a24'), ('a2^12', 'd24'), ('a24', 'e8^3'), ('d4^6', 'a6^4')]
pari.allocatemem(2 * 10 ** 9)
T0 = time.time()
FAIL = []


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


src = open(__file__.replace('cusp_pairs.py', 'all_holes.py'), encoding='utf-8').read()
TABLE = src.split('TABLE = """\\\n')[1].split('"""')[0]
n = 25


def lor2(x, y):
    return int(sum(int(a) * int(b) for a, b in zip(x[:-1], y[:-1])) - int(x[-1]) * int(y[-1]))


rows = {}
for line in TABLE.strip('\n').split('\n'):
    t = line.split()
    co = [int(x) for x in t[6:]]
    rows[t[4]] = (int(t[0]), int(t[1]), [co[1 + 24 - i] for i in range(n)] + [co[0]])
rows['leech'] = (0, 0, [2 * i for i in range(n)] + [140])

gens = []
for i in range(n - 1):
    v = [0] * (n + 1); v[i] = 2; v[i + 1] = -2; gens.append(v)
v = [0] * (n + 1); v[n - 2] = 2; v[n - 1] = 2; gens.append(v)
gens.append([1] * (n + 1))
v = [0] * (n + 1); v[n - 1] = 2; v[n] = 2; gens.append(v)
v = [0] * (n + 1); v[0] = 4; gens.append(v)
Am = pari.matrix(n + 1, len(gens), [gens[j][i] for i in range(n + 1) for j in range(len(gens))])
Hb = pari.mathnf(Am)
basis = [[int(Hb[i, j]) for i in range(n + 1)] for j in range(n + 1)]


def kernel(forms, dim):
    F = pari.matrix(len(forms), n + 1, [lor2(b, f) // 4 for f in forms for b in basis])
    K = pari.matkerint(F)
    return [[sum(int(K[i, j]) * basis[i][t] for i in range(n + 1)) for t in range(n + 1)] for j in range(dim)]


def bez(vals):
    g, coef = 0, [0] * len(vals)
    for j, val in enumerate(vals):
        if val == 0:
            continue
        if g == 0:
            g, coef = val, [0] * len(vals); coef[j] = 1; continue
        a0, b0, a1, b1, r0_, r1_ = 1, 0, 0, 1, g, val
        while r1_:
            qq = r0_ // r1_
            r0_, r1_ = r1_, r0_ - qq * r1_
            a0, a1 = a1, a0 - qq * a1
            b0, b1 = b1, b0 - qq * b1
        coef = [a0 * c for c in coef]; coef[j] += b0; g = r0_
    return [-c * g for c in coef], abs(g)


@njit(cache=True)
def enum_dual(Q, Bnd):
    d = Q.shape[0]
    q = Q.astype(np.float64).copy()
    for i in range(d):
        for j in range(i + 1, d):
            q[j, i] = q[i, j]
            q[i, j] = q[i, j] / q[i, i]
        for k in range(i + 1, d):
            for l in range(k, d):
                q[k, l] -= q[k, i] * q[i, l]
    out = []
    a = np.zeros(d, dtype=np.int64)
    T = np.zeros(d + 1)
    cen = np.zeros(d)
    hi = np.zeros(d, dtype=np.int64)
    i = d - 1
    rad = np.sqrt(Bnd / q[i, i])
    a[i] = np.int64(np.ceil(-rad - 1e-9))
    hi[i] = np.int64(np.floor(rad + 1e-9))
    while True:
        if a[i] > hi[i]:
            i += 1
            if i == d:
                break
            a[i] += 1
            continue
        t = a[i] - cen[i]
        T[i] = T[i + 1] + q[i, i] * t * t
        if T[i] > Bnd * (1 + 1e-9) + 1e-9:
            a[i] += 1
            continue
        if i == 0:
            out.append(a.copy())
            a[i] += 1
            continue
        j = i - 1
        s = 0.0
        for l in range(j + 1, d):
            s += q[j, l] * a[l]
        cen[j] = -s
        rad = np.sqrt(max(Bnd * (1 + 1e-9) + 1e-9 - T[i], 0.0) / q[j, j])
        a[j] = np.int64(np.ceil(cen[j] - rad))
        hi[j] = np.int64(np.floor(cen[j] + rad))
        T[j + 1] = T[i]
        i = j
    return out


X = symbols('X')
for n1, n2 in PAIRS:
    t1 = time.time()
    h1, k1, z1 = rows[n1]
    h2, k2, z2 = rows[n2]
    m = -(lor2(z1, z2) // 4)
    if m <= 0:
        print(f'== {n1}, {n2}: z1 . z2 = {-m} (not a timelike pair); skipped'); continue
    Mv = kernel([z1, z2], 24)
    GM = [[Fraction(lor2(a, b), 4) for b in Mv] for a in Mv]
    GMp = pari.matrix(24, 24, [int(x) for row in GM for x in row])
    detM = int(pari.matdet(GMp))
    Gi = pari.matsolve(GMp, pari.matid(24))
    Ginv = [[Fraction(str(Gi[i, j])) for j in range(24)] for i in range(24)]
    den = 1
    for row in Ginv:
        for x in row:
            den = den * x.denominator // math.gcd(den, x.denominator)
    QD = np.array([[int(x * den) for x in row] for row in Ginv], dtype=np.int64)

    def projM(x2):
        g = [Fraction(lor2(a, x2), 4) for a in Mv]
        return [sum(Ginv[i][j] * g[j] for j in range(24)) for i in range(24)]

    def order(c):
        o = 1
        for t in c:
            o = o * t.denominator // math.gcd(o, t.denominator)
        return o
    g1 = max((projM(x) for x in kernel([z1], 25)), key=order)
    g2 = max((projM(x) for x in kernel([z2], 25)), key=order)
    e12 = sum(g1[i] * sum(GM[i][j] * g2[j] for j in range(24)) for i in range(24))
    e = int((e12 * m) % m)
    if detM != m * m or order(g1) != m or order(g2) != m or math.gcd(e, m) != 1:
        print(f'== {n1}, {n2}: m = {m}, det M = {detM} (D not (Z/m)^2 in this normalisation); skipped'); continue

    def cls_dual(y):
        return (int((sum(Fraction(int(y[i])) * g2[i] for i in range(24)) * m) % m),
                int((sum(Fraction(int(y[i])) * g1[i] for i in range(24)) * m) % m))

    def cls_M(c):
        x2 = sum(c[i] * sum(GM[i][j] * g2[j] for j in range(24)) for i in range(24))
        x1 = sum(c[i] * sum(GM[i][j] * g1[j] for j in range(24)) for i in range(24))
        return (int((x2 * m) % m), int((x1 * m) % m))
    einv = pow(e, -1, m)
    qD = lambda g: Fraction((einv * g[0] * g[1]) % m, m)
    bD = lambda g, d: ((einv * (g[0] * d[1] + g[1] * d[0])) % m)
    # lines: projections of {x . z1 = -1} and {x . z2 = -1}
    def line(zz, other_gen_lines):
        cf, gg = bez([lor2(b, zz) // 4 for b in basis])
        y0 = [sum(cf[j] * basis[j][t] for j in range(n + 1)) for t in range(n + 1)]
        assert gg == 1 and lor2(y0, zz) // 4 == -1
        H = {cls_M(projM(x)) for x in kernel([zz], 25)}
        sub, frontier = set(), {(0, 0)}
        while frontier:
            nxt = set()
            for a_ in frontier:
                for b_ in H:
                    c_ = ((a_[0] + b_[0]) % m, (a_[1] + b_[1]) % m)
                    if c_ not in sub:
                        sub.add(c_); nxt.add(c_)
            frontier = nxt
        b0 = cls_M(projM(y0))
        return {((b0[0] + s_[0]) % m, (b0[1] + s_[1]) % m) for s_ in sub}, sub, b0
    A1, H1, a0 = line(z1, None)
    A2, H2, b0 = line(z2, None)
    # H1 = {(a, 0)} (classes of z1^perp: pair integrally with g1 ... ) -- determine orientation
    horiz1 = all(g[1] == 0 for g in H1)
    if horiz1:
        lam = (b0[0] * pow(a0[1], -1, m)) % m if math.gcd(a0[1], m) == 1 else None
    else:
        lam = None
    if lam is None:
        print(f'== {n1}, {n2}: m = {m}: line positions not invertible ({a0}, {b0}); skipped'); continue
    laminv = pow(lam, -1, m)
    sigma = lambda g: ((lam * g[1]) % m, (laminv * g[0]) % m)
    lines_ok = A2 == {sigma(g) for g in A1} and {sigma(g) for g in H1} == H2
    vecs = enum_dual(QD, 2.0 * den)
    table = {}
    for y in vecs:
        nrm = Fraction(int(y @ QD @ y), den)
        if nrm > 2:
            continue
        g = cls_dual(y)
        table.setdefault(g, {}).setdefault(nrm, 0)
        table[g][nrm] += 1
    classes = [(a, b) for a in range(m) for b in range(m)]
    below_ok, v = True, {}
    for g in classes:
        tg, ts = table.get(g, {}), table.get(sigma(g), {})
        for nr in set(tg) | set(ts):
            if nr < 2 and tg.get(nr, 0) != ts.get(nr, 0):
                below_ok = False
        v[g] = tg.get(Fraction(2), 0) - ts.get(Fraction(2), 0)
    T_ok = all(v[g] == 0 or qD(g) == 0 for g in classes)
    Phi = Poly(cyclotomic_poly(m, X), X, domain=ZZ) if m > 1 else None
    S_ok = True
    supp = [d for d in classes if v[d] != 0]
    for g in classes:
        w = [0] * m
        for d in supp:
            w[(-bD(g, d)) % m] += v[d]
        w[0] -= m * v[g]
        if any(w):
            P = Poly(sum(c * X ** j for j, c in enumerate(w)), X, domain=ZZ)
            if not rem(P, Phi).is_zero:
                S_ok = False
                break
    G = sum(v[g] for g in A1)
    ok = below_ok and T_ok and S_ok and lines_ok
    check(ok, f'{n1} (h={h1}) x {n2} (h={h2}): m = -z1.z2 = {m}; |M^* norm<=2| = {sum(sum(t.values()) for t in table.values())}; '
              f'E = 0 below norm 2: {below_ok}; v invariant: {T_ok and S_ok}; => theta(L1 + p1) - theta(L2 + p2) = {G} Delta [{time.time() - t1:.0f}s]')
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f' ({time.time() - T0:.0f}s)')
