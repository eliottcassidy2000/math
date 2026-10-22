"""Which component carries each shell-1 cusp class, and its THM-4760 verdict (session opus-2026-10-22-S53).

A shell-1 vector y of N + rho_N/h is a glue translate of the highest root theta_j of one component j. That component
is read off from the excess (h/2)|proj_c(y - rho_N/h)|^2 - (h/2)|proj_c(y_0 - rho_N/h)|^2, with y_0 a shell-0 vector:
it is 1 on the carrier and 0 elsewhere. Projections are exact (PARI matsolve on the Gram matrix of a root
basis of each component). The verdict is cusp_data of fricke_cusp_certificate.py, whose setup this script executes.
Every shell-1 vector is assigned its carrier (cheap); the certificate runs on REPS vectors per component.
Usage: python fricke_cusp_carriers.py REPS [types...]"""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
CERT = os.path.join(HERE, 'fricke_cusp_certificate.py')
REPS = int(sys.argv[1])
ONLY_ = sys.argv[2:]
src = open(CERT, encoding='utf-8').read()
setup = src.split('HITS = []')[0]
sys.argv = [CERT, '1', '5', '0']
__file__ = CERT
exec(compile(setup, CERT, 'exec'))
from collections import Counter

NAMES = {12: 'E6', 18: 'E7', 30: 'E8'}


def comp_name(size, h):
    rank = size // h
    if size == h * (h - 1):
        return f'A{h - 1}'
    if h % 2 == 0 and size == 2 * (h // 2 + 1) * (h // 2):
        return f'D{h // 2 + 1}'
    return NAMES.get(h, f'?{rank}')


for name, h, n0, z2 in rows:
    if ONLY_ and name not in ONLY_:
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
    svr, sar = shell_vecs(Gf, Gi, np.zeros(24, dtype=np.int64), 1, 2.0 + 1e-6)
    R = np.array([a_.copy() for s_, a_ in zip(svr, sar) if int(s_) == 2], dtype=np.int64)
    IP = R @ Gi @ R.T
    nr = len(R)
    comp = [-1] * nr
    cid = 0
    for s0 in range(nr):
        if comp[s0] >= 0:
            continue
        stack = [s0]; comp[s0] = cid
        while stack:
            u = stack.pop()
            for w_ in np.nonzero(IP[u])[0]:
                if comp[w_] < 0:
                    comp[w_] = cid; stack.append(int(w_))
        cid += 1
    members = [[i for i in range(nr) if comp[i] == c_] for c_ in range(cid)]
    bases = []
    for idx in members:
        # an exact basis of the component's span: greedy independent roots (PARI rank)
        chosen = []
        for i in idx:
            trial = chosen + [i]
            M_ = pari.matrix(len(trial), 24, [int(x) for k in trial for x in R[k]])
            if int(pari.matrank(M_)) == len(trial):
                chosen = trial
        B = R[chosen]
        Gc = B @ Gi @ B.T
        bases.append((B, pari.matrix(len(chosen), len(chosen), [int(x) for x in Gc.flatten()])))

    def parts(a0):
        H_ = h * a0 - Yr                                     # h (y - rho_N/h), integral
        out = []
        for B, Gcp in bases:
            rhs = B @ Gi @ H_
            t = pari.matsolve(Gcp, pari.matrix(len(rhs), 1, [int(x) for x in rhs]))
            pc2 = sum(pari(t[i, 0]) * int(rhs[i]) for i in range(len(rhs)))   # h^2 |proj|^2, exact
            out.append(pari(pc2) / (2 * h))                   # (h/2)|proj|^2
        return out

    sv, sa = shell_vecs(Gf, Gi, -Yr, h, 2.0 * (h + 1) / h + 2.0 * 1.5 / h + 1e-6)
    by_shell = {}
    for s_, a_ in zip(sv, sa):
        by_shell.setdefault(int(s_), []).append(a_)
    levels = sorted(by_shell)
    base = parts(by_shell[levels[0]][0])
    tally = Counter()
    groups = {}
    for a0 in by_shell[levels[1]]:
        ex = [p - b for p, b in zip(parts(a0), base)]              # exact rationals: 0 or 1
        carriers = [c_ for c_ in range(cid) if ex[c_] == 1]
        assert len(carriers) == 1 and all(ex[c_] == 0 for c_ in range(cid) if c_ != carriers[0]), ex
        groups.setdefault(carriers[0], []).append(a0)
    for c0, lst in sorted(groups.items()):
      for a0 in lst[:REPS]:
        carriers = [c0]
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
        tag = 'HAUPTMODUL' if verdict else ('IRRATIONAL' if key[6] else 'no')
        tally[(c0, comp_name(len(members[c0]), h), tag)] += 1
    print(f'  {name:9s} shell 1 (m = {h + 2}): ' + '; '.join(
        f'component {c0} = {cn} ({len(groups[c0])} vectors): {t} on {k}/{min(REPS, len(groups[c0]))} tested'
        for (c0, cn, t), k in sorted(tally.items())), flush=True)
