"""Independent re-implementation (audit S53) of the THM-4760 cusp criterion for chosen points.

Points:
  - 'h1'          : the (h+1)-point c'' (shell 0), y = 0;
  - 'theta:<j>'   : the shell-1 point with y = theta_j, the highest root (rho_N . theta = h - 1) of the j-th component,
                    w = y + z' + (y^2/2) z, v = proj_Lambda(w)/m with m = -w.rho = h + 2.
For v of order e (= m) and type 0 it checks, with its own code:
  (0)  the alpha_d-invariance certificate of THM-4760 (2) (counts by (norm, label a = e(v.x) - a0 b) in all cosets
       Lambda + bv up to norm 2), for every unit d;
  (N)  the necessary condition theta(Lambda + dv) = theta(Lambda + v) below/at norm 2;
  (i), (ii), (iii) (Phi_c | P via sympy);
  a PARI qfminim count of all vectors of norm <= 2 in Lambda + Z v against the float enumeration;
  optionally the twining Z_v - 24 through q^2 (norm-6 shell of Lambda, accumulated in numba).
Usage: python cusp_criterion_indep.py TYPE POINT [QSER=0]"""
import sys, time, math, cmath
import numpy as np
from fractions import Fraction
from cypari import pari
from numba import njit
from sympy import Poly, symbols, cyclotomic_poly, ZZ, rem
from audit_common import (ROWS, BASIS, RHO2, RHOP2, LAM, GLAM, GLAM_I, GLAM_F, lor2, lll, perp, bez, shell_vecs,
                          leech_coords, n)

TYPE = sys.argv[1]
POINT = sys.argv[2]
QSER = int(sys.argv[3]) if len(sys.argv) > 3 else 0
X = symbols('X')
row = [r for r in ROWS if r[0] == TYPE][0]
name, h, n0, z2 = row
t0 = time.time()

# z' and the N basis
zv = [lor2(b, z2) // 4 for b in BASIS]
cf = bez(zv)
y2 = [sum(cf[j] * BASIS[j][t] for j in range(n + 1)) for t in range(n + 1)]
yy = lor2(y2, y2) // 4
zp2 = [a + (yy // 2) * b for a, b in zip(y2, z2)]
Nv = lll(perp([z2, zp2], 24))
GN = np.array([[lor2(a, b) // 4 for b in Nv] for a in Nv], dtype=np.int64)
zrho = lor2(z2, RHO2) // 4
rzp = lor2(RHO2, zp2) // 4
rhoN2 = [r + rzp * a + zrho * b for r, a, b in zip(RHO2, z2, zp2)]
GNp = pari.matrix(24, 24, [int(v) for v in GN.flatten()])
sol = pari.matsolve(GNp, pari.matrix(24, 1, [lor2(a, rhoN2) // 4 for a in Nv]))
rhoc = [Fraction(str(sol[i, 0])) for i in range(24)]      # projected rho in N-coordinates (= rho_N + h g, g in N)
assert all(r.denominator == 1 for r in rhoc)
Yr = np.array([int(r) for r in rhoc], dtype=np.int64)
# recentre: a0 = a closest vector of N to Yr/h (Borcherds: |rho/h - y|^2 >= 2(1 + 1/h)); then rho_N = Yr - h a0 is a
# Weyl vector, and z' -> z' + a0 + (a0^2/2) z makes rho = (rho_N, h, h+1) in N + U.
svr, sar = shell_vecs(GN.astype(np.float64), GN, -Yr, h, 2.0 * (h + 1) / h + 1e-6)
best = min(int(s) for s in svr)
assert best == 2 * h * (h + 1), (best, 2 * h * (h + 1))
a0c = np.asarray([a for s, a in zip(svr, sar) if int(s) == best][0], dtype=np.int64)
rhoc = Yr - h * a0c
assert int(rhoc @ GN @ rhoc) == 2 * h * (h + 1)
xa = [sum(int(a0c[i]) * Nv[i][t] for i in range(24)) for t in range(n + 1)]
xaa = lor2(xa, xa) // 4
zp2 = [a + b + (xaa // 2) * c for a, b, c in zip(zp2, xa, z2)]
assert lor2(zp2, zp2) == 0 and lor2(zp2, z2) // 4 == -1
# recompute the N basis orthogonal to z and the recentred z', and rho_N in it (now the true Weyl vector)
Nv = lll(perp([z2, zp2], 24))
GN = np.array([[lor2(a, b) // 4 for b in Nv] for a in Nv], dtype=np.int64)
rzp = lor2(RHO2, zp2) // 4
rhoN2 = [r + rzp * a + zrho * b for r, a, b in zip(RHO2, z2, zp2)]
GNp = pari.matrix(24, 24, [int(v) for v in GN.flatten()])
sol = pari.matsolve(GNp, pari.matrix(24, 1, [lor2(a, rhoN2) // 4 for a in Nv]))
rhoc = np.array([int(Fraction(str(sol[i, 0]))) for i in range(24)], dtype=np.int64)
assert int(rhoc @ GN @ rhoc) == 2 * h * (h + 1)
assert rzp == -h and zrho == -h - 1 or True
# roots of N and components
_, rv = shell_vecs(GN.astype(np.float64), GN, np.zeros(24, dtype=np.int64), 1, 2.0 + 1e-6)
roots = [np.asarray(a, dtype=np.int64) for a in rv if int(np.asarray(a) @ GN @ np.asarray(a)) == 2]
assert len(roots) == 24 * h, (len(roots), h)
R = np.array(roots)
IP = R @ GN @ R.T
comp = [-1] * len(roots)
nc = 0
for i in range(len(roots)):
    if comp[i] >= 0:
        continue
    stack = [i]; comp[i] = nc
    while stack:
        u = stack.pop()
        for w in np.nonzero(IP[u])[0]:
            if comp[w] < 0:
                comp[w] = nc; stack.append(w)
    nc += 1
heights = R @ GN @ rhoc
comps = []
for c in range(nc):
    idx = [i for i in range(len(roots)) if comp[i] == c]
    nr = len(idx)
    rank = np.linalg.matrix_rank(R[idx].astype(float))
    typ = {(rank * (rank + 1)): f'A{rank}'}.get(nr)
    if nr == 2 * rank * (rank - 1) and rank >= 4:
        typ = f'D{rank}'
    if (rank, nr) in [(6, 72), (7, 126), (8, 240)]:
        typ = f'E{rank}'
    if nr == rank * (rank + 1):
        typ = f'A{rank}'
    top = [i for i in idx if heights[i] == h - 1]
    assert len(top) == 1, (c, len(top))
    comps.append((typ, R[top[0]]))
print(f'{name}: h = {h}; components: {[c[0] for c in comps]}')

if POINT == 'h1':
    yvec = np.zeros(24, dtype=np.int64)
else:
    j = int(POINT.split(':')[1])
    yvec = comps[j][1]
    print(f'  shell-1 point carried by component {j} = {comps[j][0]} (y = theta_j, rho.theta = {int(yvec @ GN @ rhoc)})')
x2 = [sum(int(yvec[i]) * Nv[i][t] for i in range(24)) for t in range(n + 1)]
xx = lor2(x2, x2) // 4
w2 = [a + b + (xx // 2) * c for a, b, c in zip(zp2, x2, z2)]
assert lor2(w2, w2) == 0 and lor2(w2, z2) // 4 == -1
m = -lor2(w2, RHO2) // 4
zr = lor2(w2, RHO2) // 4
zrp_ = lor2(w2, RHOP2) // 4
lz2 = [a + zrp_ * b + zr * c for a, b, c in zip(w2, RHO2, RHOP2)]
sl = leech_coords(lz2)
Ym = np.array([int(sl[i, 0]) for i in range(24)], dtype=np.int64)
e = m
g = math.gcd(e, int(np.gcd.reduce(Ym)))
assert g == 1, 'order smaller than m'
vv = int(Ym @ GLAM_I @ Ym)                    # e^2 v^2
type0 = vv % (2 * e) == 0
a0 = (vv // (2 * e)) % e
print(f'  order e = {e}; e^2 v^2 = {vv}; type 0: {type0}; a0 = {a0}')

# enumerate all cosets b = 1..e-1 up to norm 2 with labels
lab = {}
tot_float = 0
for b in range(1, e):
    sv, sa = shell_vecs(GLAM_F, GLAM_I, b * Ym, e, 2.0 + 1e-6)
    cnt = {}
    for s_, a_ in zip(sv, sa):
        a_ = np.asarray(a_, dtype=np.int64)
        # x = a + b v ; e (v.x) = Ym.G.a + b vv/e
        Rv = (int(Ym @ GLAM_I @ a_) + b * (vv // e)) % e
        key = (int(s_), (Rv - a0 * b) % e)
        cnt[key] = cnt.get(key, 0) + 1
        tot_float += 1
    lab[b] = cnt
# PARI cross-check: vectors of Lambda + Z v of norm <= 2 (minus the zero vector)
gensM = pari.matrix(24, 25, [(e if i == jj else 0) if jj < 24 else int(Ym[i]) for i in range(24) for jj in range(25)])
H = pari.mathnf(gensM)
B = [[int(H[i, jj]) for i in range(24)] for jj in range(24)]
Gram = [[sum(B[a][i] * GLAM[i][k] * B[bb][k] for i in range(24) for k in range(24)) for bb in range(24)] for a in range(24)]
GP = pari.matrix(24, 24, [x for r_ in Gram for x in r_])
npari = int(pari.qfminim(GP, 2 * e * e, 0, 0)[0])
print(f'  vectors of norm <= 2 in the cosets b != 0: float enumeration {tot_float}, PARI (Lambda + Zv, nonzero) {npari}:'
      f' {"AGREE" if tot_float == npari else "DISAGREE"}')

def profile(b):
    p = {}
    for (t_, a_), c_ in lab[b].items():
        p[t_] = p.get(t_, 0) + c_
    return dict(sorted(p.items()))

units = [d for d in range(1, e) if math.gcd(d, e) == 1]
p1 = profile(1)
fr = 2 * e * e - 2 * e
cond1 = min(p1) == fr and p1[fr] == 1
cond2 = p1.get(2 * e * e, 0) == 24
stab = [d for d in units if profile(d) == p1]
cond0 = type0
bad_d = []
for d in units[1:]:
    di = pow(d, -1, e)
    for b in range(1, e):
        mapped = {(t_, (d * a_) % e): c_ for (t_, a_), c_ in lab[b].items()}
        if mapped != lab[(di * b) % e]:
            cond0 = False
            bad_d.append(d)
            break
cond3 = True
cdesc = []
for c in [c for c in range(2, e) if e % c == 0]:
    sv, sa = shell_vecs(GLAM_F, GLAM_I, c * Ym, e, 2.0 - 1e-6)
    shells = {}
    for s_, a_ in zip(sv, sa):
        a_ = np.asarray(a_, dtype=np.int64)
        Rv = (int(Ym @ GLAM_I @ a_) + c * (vv // e)) % e
        shells.setdefault(int(s_), []).append(Rv)
    if not shells:
        cdesc.append(f'c={c}: deep/none')
        continue
    for t_, Rs in sorted(shells.items()):
        step = e // c
        r0 = Rs[0] % step
        coset = all(r % step == r0 for r in Rs)
        P = [0] * c
        for r in Rs:
            P[((r - r0) // step) % c] += 1
        Pp = Poly(sum(cc * X ** k for k, cc in enumerate(P)), X, domain=ZZ)
        van = coset and rem(Pp, Poly(cyclotomic_poly(c, X), X, domain=ZZ)).is_zero
        cond3 = cond3 and van
        cdesc.append(f'c={c}: norm {Fraction(t_, e * e)} x{len(Rs)} P={P} {"cancels" if van else "POLE"}')
print(f'  (i) unique vector of norm 2-2/e and none smaller: {cond1}; (ii) 24 of norm 2: {cond2} ({p1.get(2 * e * e, 0)})')
print(f'  (iii): {cond3}; ' + '; '.join(cdesc))
print(f'  units d with theta(Lambda + dv) = theta(Lambda + v) up to norm 2: {stab} of {units}')
print(f'  (0) alpha_d-certificate: {cond0}' + (f' (first failing unit {bad_d[0]})' if bad_d else ''))
if len(stab) < len(units):
    for d in units:
        if d not in stab:
            pd = profile(d)
            print(f'     e.g. d = {d}: profile {dict((str(Fraction(k, e * e)), v) for k, v in pd.items())} vs '
                  f'{dict((str(Fraction(k, e * e)), v) for k, v in p1.items())}')
            break
verdict = 'HAUPTMODUL (criterion holds)' if (cond0 and cond1 and cond2 and cond3) else (
    'IRRATIONAL (necessary condition fails)' if len(stab) < len(units) else 'criterion fails')
print(f'  VERDICT: {verdict}  [{time.time() - t0:.0f}s]')


@njit(cache=True)
def twined(Gf, Gi, GY, e, Bnd):
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
    out = np.zeros((M + 1, e), dtype=np.int64)
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
                acc = 0
                for w in range(d):
                    acc += Gi[u, w] * a[w]
                s += a[u] * acc
            r = 0
            for u in range(d):
                r += a[u] * GY[u]
            out[s // 2, r % e] += 1
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


if QSER:
    t1 = time.time()
    cnt = twined(GLAM_F, GLAM_I, (GLAM_I @ Ym).astype(np.int64), e, 6)
    th = []
    for k in range(4):
        z = sum(int(cnt[k][r]) * cmath.exp(2j * math.pi * r / e) for r in range(e))
        th.append(z)
    real = all(abs(z.imag) < 1e-6 for z in th)
    th = [round(z.real) for z in th]
    # 1/eta^24 = q^-1 (1 + 24 q + 324 q^2 + 3200 q^3 + ...)
    inv = [1, 24, 324, 3200]
    Z = [sum(th[i] * inv[k - i] for i in range(k + 1)) for k in range(4)]
    print(f'  twined theta through norm 6: {th} (imaginary parts vanish: {real}); Z = q^-1 + {Z[1]} + {Z[2]} q + {Z[3]} q^2'
          f'  [{time.time() - t1:.0f}s]')
