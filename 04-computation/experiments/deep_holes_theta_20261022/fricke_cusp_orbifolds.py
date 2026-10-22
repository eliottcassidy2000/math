"""Orbifolds of V_Lambda along the shell-s cusp points and their multiples (session opus-2026-10-22-S53; HYP-9368).

For a Leech point u of order k with k u^2/2 in Z, the cyclic orbifold of V_Lambda by e^{2 pi i u(0)} is the lattice
VOA of the Kneser neighbour Lambda' = Lambda_u + Z w (Lambda_u = {x : x.u in Z}, w in u + Lambda of even norm). Its roots
are the norm-2 vectors x of Lambda + j u (j = 1..k-1) with (x - j w).u in Z. For each class of shell-s cusp points v
of the given types (one representative per class of the cusp certificate), this prints the neighbour type at v and at
every multiple cv, c | m, 1 < c < m, together with the minimal norm and multiplicity of Lambda + cv.
Usage: python fricke_cusp_orbifolds.py SHELL CAP types..."""
import sys, time, math, re
import numpy as np
from cypari import pari
from numba import njit

SMAX = int(sys.argv[1])
CAP = int(sys.argv[2])
QSER = 0
ONLY = sys.argv[3:]
NORM = 6
pari.allocatemem(2 * 10 ** 9)
FAIL = []


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


src = open(__file__.replace('fricke_cusp_orbifolds.py', 'all_holes.py'), encoding='utf-8').read()
TABLE = src.split('TABLE = """\\\n')[1].split('"""')[0]
n = 25


def lor2(x, y):
    return int(sum(int(a) * int(b) for a, b in zip(x[:-1], y[:-1])) - int(x[-1]) * int(y[-1]))


rows = []
for line in TABLE.strip('\n').split('\n'):
    t = line.split()
    co = [int(x) for x in t[6:]]
    rows.append((t[4], int(t[0]), int(t[1]), [co[1 + 24 - i] for i in range(n)] + [co[0]]))
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
rho2 = [2 * i for i in range(n)] + [140]


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
    return [-c * g for c in coef]


def perp_basis(u2, w2):
    F = pari.matrix(2, n + 1, [lor2(b, u2) // 4 for b in basis] + [lor2(b, w2) // 4 for b in basis])
    K = pari.matkerint(F)
    vecs = [[sum(int(K[i, j]) * basis[i][t] for i in range(n + 1)) for t in range(n + 1)] for j in range(24)]
    G = pari.matrix(24, 24, [lor2(a, b) // 4 for a in vecs for b in vecs])
    U = pari.qflllgram(G)
    vecs = [[sum(int(U[i, j]) * vecs[i][t] for i in range(24)) for t in range(n + 1)] for j in range(24)]
    return vecs, [[lor2(a, b) // 4 for b in vecs] for a in vecs]


@njit(cache=True)
def shell_vecs(Gf, Gi, Yd, den, R2):
    """all a with |a + Yd/den|^2 <= R2: returns den^2 |a + Yd/den|^2 and a."""
    d = Gf.shape[0]
    q = Gf.copy()
    for i in range(d):
        for j in range(i + 1, d):
            q[j, i] = q[i, j]
            q[i, j] = q[i, j] / q[i, i]
        for k in range(i + 1, d):
            for l in range(k, d):
                q[k, l] -= q[k, i] * q[i, l]
    y = -Yd.astype(np.float64) / den
    vals = []
    vecs = []
    a = np.zeros(d, dtype=np.int64)
    T = np.zeros(d + 1)
    cen = np.zeros(d)
    hi = np.zeros(d, dtype=np.int64)
    i = d - 1
    cen[i] = y[i]
    rad = np.sqrt(R2 + 1e-9) / np.sqrt(q[i, i])
    a[i] = np.int64(np.ceil(cen[i] - rad))
    hi[i] = np.int64(np.floor(cen[i] + rad))
    while True:
        if a[i] > hi[i]:
            i += 1
            if i == d:
                break
            a[i] += 1
            continue
        t = a[i] - cen[i]
        T[i] = T[i + 1] + q[i, i] * t * t
        if T[i] > R2 + 1e-9:
            a[i] += 1
            continue
        if i == 0:
            s = 0
            for u in range(d):
                zu = den * a[u] + Yd[u]
                if zu == 0:
                    continue
                acc = 0
                for w in range(d):
                    acc += Gi[u, w] * (den * a[w] + Yd[w])
                s += zu * acc
            vals.append(s)
            vecs.append(a.copy())
            a[i] += 1
            continue
        j = i - 1
        sacc = 0.0
        for l in range(j + 1, d):
            sacc += q[j, l] * (a[l] - y[l])
        cen[j] = y[j] - sacc
        rad = np.sqrt(max(R2 + 1e-9 - T[i], 0.0) / q[j, j])
        a[j] = np.int64(np.ceil(cen[j] - rad))
        hi[j] = np.int64(np.floor(cen[j] + rad))
        T[j + 1] = T[i]
        i = j
    return vals, vecs



rhop2 = list(rho2); rhop2[0] += 2; rhop2[1] -= 2
F = pari.matrix(2, n + 1, [lor2(b, rho2) // 4 for b in basis] + [lor2(b, rhop2) // 4 for b in basis])
K = pari.matkerint(F)
lam = [[sum(int(K[i, j]) * basis[i][t] for i in range(n + 1)) for t in range(n + 1)] for j in range(24)]
GL = pari.matrix(24, 24, [lor2(a, b) // 4 for a in lam for b in lam])
U = pari.qflllgram(GL)
lam = [[sum(int(U[i, j]) * lam[i][t] for i in range(24)) for t in range(n + 1)] for j in range(24)]
GmL = [[lor2(a, b) // 4 for b in lam] for a in lam]
GpL = pari.matrix(24, 24, [x for row in GmL for x in row])
GL_i = np.array(GmL, dtype=np.int64)
GL_f = GL_i.astype(np.float64)
import cmath


@njit(cache=True)
def twisted_counts(Gf, Gi, GY, h, Bnd):
    d = Gf.shape[0]
    q = Gf.copy()
    for i in range(d):
        for j in range(i + 1, d):
            q[j, i] = q[i, j]
            q[i, j] = q[i, j] / q[i, i]
        for k in range(i + 1, d):
            for l in range(k, d):
                q[k, l] -= q[k, i] * q[i, l]
    M = Bnd // 2
    out = np.zeros((M + 1, h), dtype=np.int64)
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
        if T[i] > Bnd + 1e-7:
            a[i] += 1
            continue
        if i == 0:
            s = 0
            for u in range(d):
                if a[u] == 0:
                    continue
                acc = 0
                for w in range(d):
                    acc += Gi[u, w] * a[w]
                s += a[u] * acc
            r = 0
            for u in range(d):
                r += a[u] * GY[u]
            out[s // 2, r % h] += 1
            a[i] += 1
            continue
        j = i - 1
        sacc = 0.0
        for l in range(j + 1, d):
            sacc += q[j, l] * a[l]
        cen[j] = -sacc
        rad = np.sqrt(max(Bnd + 1e-7 - T[i], 0.0) / q[j, j])
        a[j] = np.int64(np.ceil(cen[j] - rad))
        hi[j] = np.int64(np.floor(cen[j] + rad))
        T[j + 1] = T[i]
        i = j
    return out


def inv_delta(L, step=1):
    s = [0] * L; s[0] = 1
    for mm in range(1, L // step + 2):
        for _ in range(24):
            for t in range(step * mm, L):
                s[t] += s[t - step * mm]
    return s



def _merge(cnt_b):
    out = {}
    for (t_, a_), c_ in cnt_b.items():
        out[t_] = out.get(t_, 0) + c_
    return out


def cusp_data(Ym, m):
    """Conditions (i)-(iii) for v = Ym/m (Leech coordinates); returns the invariant key and the verdict."""
    vv = int(Ym @ GL_i @ Ym)                                     # m^2 v^2
    type0 = vv % (2 * m) == 0                                    # m v^2/2 in Z
    svL, saL = shell_vecs(GL_f, GL_i, Ym, m, 2.0 + 1e-6)
    vals = [int(s) for s in svL]
    mn = min(vals)
    nmin = vals.count(mn)
    n2 = vals.count(2 * m * m)
    cond1 = mn == 2 * m * m - 2 * m and nmin == 1
    cond2 = n2 == 24
    a0 = (vv // (2 * m)) % m if type0 else 0
    lab = {}
    for b in range(1, m):
        svb, sab = shell_vecs(GL_f, GL_i, b * Ym, m, 2.0 + 1e-6)
        cnt_b = {}
        for s_, a_ in zip(svb, sab):
            R = (int(Ym @ GL_i @ a_) + b * vv // m) % m
            k_ = (int(s_), (R - a0 * b) % m)
            cnt_b[k_] = cnt_b.get(k_, 0) + 1
        lab[b] = cnt_b
    prof = {b: sorted((t_, c_) for t_, c_ in _merge(lab[b]).items()) for b in lab}
    irr = any(prof[d] != prof[1] for d in range(2, m) if math.gcd(d, m) == 1)
    cond0 = type0
    for d in [d for d in range(2, m) if math.gcd(d, m) == 1]:
        di = pow(d, -1, m)
        for b in range(1, m):
            tgt = lab[(di * b) % m]
            mapped = {(t_, (d * a_) % m): c_ for (t_, a_), c_ in lab[b].items()}
            if mapped != tgt:
                cond0 = False
                break
        if not cond0:
            break
    cusps = []
    cond3 = True
    for c in [c for c in range(2, m) if m % c == 0]:
        svc, sac = shell_vecs(GL_f, GL_i, c * Ym, m, 2.0 - 1e-6)
        shells = {}
        for s_, a_ in zip(svc, sac):
            R = (int(Ym @ GL_i @ a_) + c * vv // m) % m              # e (v.x), x = a + c Ym/m
            shells.setdefault(int(s_), []).append(R)
        ok_c = True
        desc = []
        for t_ in sorted(shells):
            Rs = shells[t_]
            step = m // c
            r0 = Rs[0] % step
            coset = all(R % step == r0 for R in Rs)
            P = [0] * c
            for R in Rs:
                P[((R - r0) // step) % c] += 1
            van = coset and int(pari(f'(Pol({P[::-1]}, T) % polcyclo({c}, T)) == 0')) == 1
            ok_c = ok_c and van
            g_ = math.gcd(t_, m * m)
            desc.append(f'{t_ // g_}/{m * m // g_}:{len(Rs)}{"" if van else "!"}')
        cond3 = cond3 and ok_c
        cusps.append(f'c={c}[' + (' '.join(desc) if desc else 'deep') + ']')
    key = (m, cond0, cond1, cond2, cond3, ' '.join(cusps), irr)
    return key, cond0 and cond1 and cond2 and cond3


def qseries(Ym, e):
    cnt = twisted_counts(GL_f, GL_i, (GL_i @ Ym).astype(np.int64), e, NORM)
    M = NORM // 2
    th = [sum(int(cnt[k][r]) * cmath.exp(2j * math.pi * r / e) for r in range(e)) for k in range(M + 1)]
    th = [round(zz.real) for zz in th]
    idl = inv_delta(M + 1)
    return [sum(th[i] * idl[k - i] for i in range(0, k + 1)) for k in range(M + 1)]




def neighbour(Yu, k):
    """Root system of the Kneser neighbour of Lambda at u = Yu/k (Leech coordinates)."""
    from fractions import Fraction as Fr
    uu = int(Yu @ GL_i @ Yu)                           # k^2 u^2
    if uu % (2 * k):
        return 'not type 0'
    kk = uu // (2 * k)                                 # k u^2 / 2
    # lambda with k(u.lambda) = -kk mod k: Bezout on R_i = (Yu G)_i
    Ri = [int(t) % k for t in (Yu @ GL_i)]
    g_, coef = 0, [0] * 24
    for i, r in enumerate(Ri):
        if r == 0:
            continue
        if g_ == 0:
            g_, coef = r, [0] * 24; coef[i] = 1; continue
        a0, b0, a1, b1, r0_, r1_ = 1, 0, 0, 1, g_, r
        while r1_:
            qq = r0_ // r1_
            r0_, r1_ = r1_, r0_ - qq * r1_
            a0, a1 = a1, a0 - qq * a1
            b0, b1 = b1, b0 - qq * b1
        coef = [a0 * c for c in coef]; coef[i] += b0; g_ = r0_
    assert math.gcd(g_, k) == 1, 'u not of exact order k'
    t_ = (-kk * pow(g_, -1, k)) % k
    lam_c = np.array([(c * t_) for c in coef], dtype=np.int64)
    W = Yu + k * lam_c
    assert int(W @ GL_i @ W) % (2 * k * k) == 0
    roots = []
    for j in range(1, k):
        sv_, sa_ = shell_vecs(GL_f, GL_i, j * Yu, k, 2.0 + 1e-6)
        for s_, a_ in zip(sv_, sa_):
            if int(s_) != 2 * k * k:
                continue
            X = k * a_ + j * Yu                        # k x
            if int((X - j * W) @ GL_i @ Yu) % (k * k) == 0:
                roots.append(X)
    nr = len(roots)
    if nr == 0:
        return 'Leech'
    if nr % 24:
        return f'?{nr} roots'
    hh = nr // 24
    Rm = np.array(roots, dtype=np.int64)
    IP = Rm @ GL_i @ Rm.T
    adj = IP != 0
    seen = [False] * nr
    comps = []
    for s0 in range(nr):
        if seen[s0]:
            continue
        stack = [s0]; seen[s0] = True; c = 0
        while stack:
            u_ = stack.pop(); c += 1
            for w_ in np.nonzero(adj[u_])[0]:
                if not seen[w_]:
                    seen[w_] = True; stack.append(int(w_))
        comps.append(c)
    out = []
    for s in sorted(comps, reverse=True):
        if s == hh * (hh - 1):
            out.append(f'A{hh - 1}')
        elif hh % 2 == 0 and s == 2 * (hh // 2 + 1) * (hh // 2):
            out.append(f'D{hh // 2 + 1}')
        elif (s, hh) in ((72, 12), (126, 18), (240, 30)):
            out.append({12: 'E6', 18: 'E7', 30: 'E8'}[hh])
        else:
            out.append(f'?({s})')
    from collections import Counter
    cc = Counter(out)
    return ''.join(f'{a}^{b}' if b > 1 else a for a, b in sorted(cc.items()))


for name, h, n0, z2 in rows:
    if ONLY and name not in ONLY:
        continue
    zv = [lor2(b, z2) // 4 for b in basis]
    cf = bez(zv)
    y2 = [sum(cf[j] * basis[j][t] for j in range(n + 1)) for t in range(n + 1)]
    yy = lor2(y2, y2) // 4
    zp2 = [a + (yy // 2) * b for a, b in zip(y2, z2)]
    Nv, GN = perp_basis(z2, zp2)
    zrho = lor2(z2, rho2) // 4
    rzp = lor2(rho2, zp2) // 4
    rhoN2 = [r + rzp * a + zrho * b for r, a, b in zip(rho2, z2, zp2)]
    GNp = pari.matrix(24, 24, [v for row in GN for v in row])
    sol = pari.matsolve(GNp, pari.matrix(24, 1, [lor2(a, rhoN2) // 4 for a in Nv]))
    Yr = np.array([int(sol[i, 0]) for i in range(24)], dtype=np.int64)
    Gi = np.array(GN, dtype=np.int64)
    Gf = Gi.astype(np.float64)
    sv, sa = shell_vecs(Gf, Gi, -Yr, h, 2.0 * (h + 1) / h + 2.0 * (SMAX + 0.5) / h + 1e-6)
    by_shell = {}
    for s_, a_ in zip(sv, sa):
        by_shell.setdefault(int(s_), []).append(a_)
    lev = sorted(by_shell)[SMAX]
    shell_list = by_shell[lev]
    if len(shell_list) > CAP:
        idx_ = np.random.default_rng(11).choice(len(shell_list), CAP, replace=False)
        shell_list = [shell_list[i_] for i_ in idx_]
    seen_keys = {}
    for a0 in shell_list:
        x2 = [sum(int(a0[i]) * Nv[i][t] for i in range(24)) for t in range(n + 1)]
        xx = lor2(x2, x2) // 4
        zpp2 = [a + b + (xx // 2) * c for a, b, c in zip(zp2, x2, z2)]
        m = -lor2(zpp2, rho2) // 4
        zr = lor2(zpp2, rho2) // 4
        zrp_ = lor2(zpp2, rhop2) // 4
        lz2 = [a + zrp_ * b + zr * c for a, b, c in zip(zpp2, rho2, rhop2)]
        sol = pari.matsolve(GpL, pari.matrix(24, 1, [lor2(a, lz2) // 4 for a in lam]))
        Ym = np.array([int(sol[i, 0]) for i in range(24)], dtype=np.int64)
        key, verdict = cusp_data(Ym, m)
        if key in seen_keys:
            continue
        seen_keys[key] = True
        line = [f'{name:9s} shell {SMAX} m={m} {"HAUPTMODUL" if verdict else "-"}: v -> {neighbour(Ym, m)}']
        for c in [c for c in range(2, m) if m % c == 0]:
            k = m // c                                 # c v = Ym / k, a point of order k
            svc, sac = shell_vecs(GL_f, GL_i, Ym, k, 2.0 + 1e-6)
            vals = [int(s) for s in svc]
            mn = min(vals) if vals else None
            gq = math.gcd(mn, k * k) if mn else 1
            line.append(f'{c}v (order {k}): min {mn // gq}/{k * k // gq} x{vals.count(mn)} -> {neighbour(Ym, k)}')
        print('  ' + '\n      '.join(line), flush=True)
