"""Identify the classes of the S53 audit's order-5 census (audit_S53/order5_invariants.out) with cusp points (session
opus-2026-10-22-S53; HYP-9368 (A)): the distribution of y.lambda over the 196560 minimal Leech vectors, y = 5x with x
the Fricke-norm vector of Lambda + v, for the A3^8 (h+1)-point and the A2^12 shell-1 point (both of order 5). Runs the
setup of fricke_cusp_certificate.py. Usage: python order5_compare.py [path to fricke_cusp_certificate.py]"""
import sys, os
from collections import Counter
cert = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fricke_cusp_certificate.py')
src = open(cert, encoding='utf-8').read()
setup = src.split('HITS = []')[0]
sys.argv = [cert, '1', '5', '0']
__file__ = cert
exec(compile(setup, cert, 'exec'))

sv4, sa4 = shell_vecs(GL_f, GL_i, np.zeros(24, dtype=np.int64), 1, 4.0 + 1e-6)
MIN = np.array([a_ for s_, a_ in zip(sv4, sa4) if int(s_) == 4], dtype=np.int64)
print('minimal vectors:', len(MIN))
for name, h, n0, z2 in rows:
    if name not in ('a2^12', 'a3^8'):
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
    sv, sa = shell_vecs(Gf, Gi, -Yr, h, 2.0 * (h + 1) / h + 2.0 * 1.5 / h + 1e-6)
    by_shell = {}
    for s_, a_ in zip(sv, sa):
        by_shell.setdefault(int(s_), []).append(a_)
    for sidx in (0, 1):
        a0 = by_shell[sorted(by_shell)[sidx]][0]
        x2 = [sum(int(a0[i]) * Nv[i][t] for i in range(24)) for t in range(n + 1)]
        xx = lor2(x2, x2) // 4
        zpp2 = [a + b + (xx // 2) * c for a, b, c in zip(zp2, x2, z2)]
        m = -lor2(zpp2, rho2) // 4
        if m != 5:
            continue
        zr = lor2(zpp2, rho2) // 4
        zrp_ = lor2(zpp2, rhop2) // 4
        lz2 = [a + zrp_ * b + zr * c for a, b, c in zip(zpp2, rho2, rhop2)]
        sol = pari.matsolve(GpL, pari.matrix(24, 1, [lor2(a, lz2) // 4 for a in lam]))
        Ym = np.array([int(sol[i, 0]) for i in range(24)], dtype=np.int64)
        # y = 5v = Ym (in Lambda), but reduce: choose the representative of minimal norm in Ym + 5 Lambda? use Ym as is
        # the audit's y has norm 40: find the norm-40 representative y = Ym + 5 lambda0 with |y/5 - lambda0'|... use
        # the unique Fricke-norm vector x of Lambda + v: y = 5x has norm 25 (2 - 2/5) = 40
        svL, saL = shell_vecs(GL_f, GL_i, Ym, m, 2.0 - 2.0 / m + 1e-6)
        xs = [a_ for s_, a_ in zip(svL, saL) if int(s_) == 2 * m * m - 2 * m]
        y = m * xs[0] + Ym
        ny = int(y @ GL_i @ y)
        dist = Counter(int(t) for t in (MIN @ GL_i @ y))
        print(f'{name} shell {sidx} (m = 5): |y|^2 = {ny}; y.lambda distribution {dict(sorted(dist.items()))}')
