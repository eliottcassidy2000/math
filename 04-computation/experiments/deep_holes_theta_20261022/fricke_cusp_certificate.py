"""Cusp certificate for Fricke Hauptmoduln of Leech inner automorphisms (session opus-2026-10-22-S53; HYP-9368).

CRITERION (classical theta transformations, no VOA input). Let v in Lambda (x) Q have order e, with e v^2/2 in Z
(type 0), and Z_v = Theta^[v]/eta^24, Theta^[v] = sum_x e^(2 pi i x.v) q^(x^2/2). For gamma = (a b; c d) in SL2(Z),
Z_v o gamma is a root of unity times theta_{cv,dv}/eta^24, where theta_{u,w} = sum_{x in Lambda+u} e^(2 pi i w.x) q^(x^2/2).
Then Z_v - 24 is the normalized Hauptmodul of Gamma_0(e)+e (so X_0(e)+e has genus zero) iff
  (0)   Z_v has rational coefficients, so that it is Gamma_0(e)-invariant (Theta^[v] o gamma = Theta^[dv] there).
        M = {lambda : lambda.v in Z} has M*/M = (Z/e)^2 with q(a,b) = ab/e, Lambda = classes (a,0) and Lambda + bv =
        classes (a,b), a = e(v.x) - a_0 b, a_0 = e v^2/2. Each unit d gives the isometry (a,b) -> (da, b/d), which
        commutes with the Weil representation. By the vector-valued Sturm lemma (weight 12: a holomorphic form
        vanishing to order > 1 at infinity is 0) the theta vector is invariant iff it is through q^1, i.e. iff the
        vectors of norm <= 2 of Lambda + bv, labelled by (norm, a), match those of Lambda + bv/d labelled (norm, a/d)...
        precisely count_b(t, a) = count_{b/d}(t, da) for all b, t <= 2, a. This gives theta_{(a,0)} = theta_{(da,0)}.
  (i)   Lambda + v has a unique vector of norm 2 - 2/e (simple pole at the cusp 0),
  (ii)  Lambda + v has exactly 24 vectors of norm 2 (then Z_v o W_e - Z_v, a holomorphic function on X_0(e), is the
        constant n_2 - 24 = 0),
  (iii) no pole at the other cusps: for every c | e with 1 < c < e and every norm t < 2 of Lambda + cv, the twined
        sum over the norm-t vectors vanishes for every d prime to c. On a fixed norm, e(v.x) runs over one coset
        r_0 + (e/c)Z of Z/e, so the condition is that Phi_c divides P(T) = sum_x T^s(x), s(x) = (e(v.x) - r_0)/(e/c).
        (Integer coefficients: one primitive c-th root suffices, by Galois.)
Tests the cusps w adjacent to the deep-hole cusp z (w.z = -1) in the shells 0..SMAX of N + rho_N/h, i.e. m = -w.rho =
h + 1 + s, for all 23 types; per shell every vector, or CAP random ones. Prints one line per invariant class, with
the verdict HAUPTMODUL ((0)-(iii) hold: PROVED), IRRATIONAL (some unit d gives theta(Lambda + dv) != theta(Lambda + v)
below norm 2, so Theta^[dv] != Theta^[v] and Z_v is not rational: PROVED not a Fricke Hauptmodul), "no" ((0) holds
and (i), (ii) or (iii) fails: PROVED not), or UNDECIDED ((0) fails without a visible irrationality), and, with
QSER = 1, the twining through q^2 for the HAUPTMODUL classes.
Usage: python fricke_cusp_certificate.py SMAX CAP QSER [types...]"""
import sys, time, math, re
import numpy as np
from cypari import pari
from numba import njit

SMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 2
CAP = int(sys.argv[2]) if len(sys.argv) > 2 else 400
QSER = int(sys.argv[3]) if len(sys.argv) > 3 else 1
ONLY = sys.argv[4:]
NORM = 6
pari.allocatemem(2 * 10 ** 9)
FAIL = []


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


src = open(__file__.replace('fricke_cusp_certificate.py', 'all_holes.py'), encoding='utf-8').read()
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
    units = [d for d in range(1, m) if math.gcd(d, m) == 1]
    stab = [d for d in units if prof[d] == prof[1]]
    irr = len(stab) < len(units)
    if irr:
        irr = f'units fixing theta(Lambda+v) below norm 2: {stab} (index {len(units) // len(stab)})'
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


HITS = []
for name, h, n0, z2 in rows:
    if ONLY and name not in ONLY:
        continue
    t1 = time.time()
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
    levels = sorted(by_shell)
    for s_idx, lev in enumerate(levels[:SMAX + 1]):
        shell_list = by_shell[lev]
        tot = len(shell_list)
        if tot > CAP:
            idx_ = np.random.default_rng(11).choice(tot, CAP, replace=False)
            shell_list = [shell_list[i_] for i_ in idx_]
        classes = {}
        for a0 in shell_list:
            x2 = [sum(int(a0[i]) * Nv[i][t] for i in range(24)) for t in range(n + 1)]
            xx = lor2(x2, x2) // 4
            zpp2 = [a + b + (xx // 2) * c for a, b, c in zip(zp2, x2, z2)]
            assert lor2(zpp2, zpp2) == 0 and lor2(zpp2, z2) // 4 == -1
            m = -lor2(zpp2, rho2) // 4
            zr = lor2(zpp2, rho2) // 4
            zrp_ = lor2(zpp2, rhop2) // 4
            lz2 = [a + zrp_ * b + zr * c for a, b, c in zip(zpp2, rho2, rhop2)]
            sol = pari.matsolve(GpL, pari.matrix(24, 1, [lor2(a, lz2) // 4 for a in lam]))
            Ym = np.array([int(sol[i, 0]) for i in range(24)], dtype=np.int64)
            key, verdict = cusp_data(Ym, m)
            if key not in classes:
                classes[key] = [0, verdict, Ym]
            classes[key][0] += 1
        for key, (cnt_, verdict, Ym) in classes.items():
            m, c0, c1, c2, c3, cusps, irr = key
            ser = ''
            if QSER and verdict:
                ser = f' Z-24 = q^-1 + {qseries(Ym, m)[2:]}'
            tag = 'HAUPTMODUL' if verdict else ('IRRATIONAL' if irr else ('no        ' if c0 else 'UNDECIDED '))
            print(f'  {tag} {name:9s} shell {s_idx} m={m:2d} x{cnt_}/{tot}: '
                  f'(0) {c0}, (i) {c1}, (ii) {c2}, (iii) {c3}; {cusps}{ser}'
                  + (f'; {irr}' if irr and c1 and c2 and c3 else ''), flush=True)
            if verdict:
                HITS.append((m, name, s_idx))
    print(f'       {name} done [{time.time() - t1:.0f}s]', flush=True)
print('HAUPTMODUL classes (order, type, shell):', sorted(HITS))
print('orders realised:', sorted(set(h_[0] for h_ in HITS)))
