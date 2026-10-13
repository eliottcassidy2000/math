"""The 21-vertex Exoo-Ismailescu graph EI21 (HoG 51378) in Q(zeta55), by the direction-class method.

Run from the repository root:  python3 -B <script>            (about one minute)
                               python3 -B <script> --lscan    (adds the L-unit parameter scan, ~7 min)
Reuses the structure / census / LLL functions of ei19_fish_placement_q_zeta55_20261012.py.

EI21 has 21 vertices and 38 = 2*21 - 4 edges (one short of minimal rigidity), girth 4, |Aut| = 2,
10 unit 4-cycles, 18 direction classes, 8 relations, 9 free angles: a one-parameter family like the
fish, with no relation of support three (two of support four), hence no forced angle.  The L-unit
scan (117 parameter values, 524 faithful realizations) finds no further direction in Q(zeta55).
"""
import importlib.util, json, sys
from pathlib import Path
import numpy as np
import mpmath as mp
from scipy.optimize import least_squares
from sympy import Rational

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('base', ROOT/'04-computation/experiments/ei19_fish_placement_q_zeta55_20261012.py')
base = importlib.util.module_from_spec(spec); spec.loader.exec_module(base)
check, structure, small_support_relations, census, lll_first, positions = base.check, base.structure, base.small_support_relations, base.census, base.lll_first, base.positions

EI21 = [[4, 8, 19], [5, 9, 20], [7, 9, 20], [6, 8, 19], [0, 6, 14], [1, 7, 13], [3, 4, 11], [2, 5, 12], [0, 3, 11], [1, 2, 12],
        [13, 14, 18], [6, 8, 17, 18], [7, 9, 16, 18], [5, 10, 15, 17], [4, 10, 15, 16], [13, 14, 19, 20], [12, 14, 17, 20],
        [11, 13, 16, 19], [10, 11, 12, 19, 20], [0, 3, 15, 17, 18], [1, 2, 15, 16, 18]]


def automorphisms(ADJ):
    n = len(ADJ); adj = {v: set(a) for v, a in enumerate(ADJ)}
    order = sorted(range(n), key=lambda v: -len(adj[v])); found = []

    def ext(m):
        if len(m) == n:
            found.append(1); return
        v = order[len(m)]
        for h in range(n):
            if h in m.values() or len(adj[h]) != len(adj[v]):
                continue
            if all((m[w] in adj[h]) == (w in adj[v]) for w in m):
                m[v] = h; ext(m); del m[v]
    ext({})
    return len(found)


def lscan(S):
    m = S['m']; piv = S['pivots']; R = S['R']; free = [c for c in range(m) if c not in piv]
    gauge, param = free[0], free[1]; angles = [c for c in free if c not in (gauge, param)]
    rows = [[Rational(R[i, c]) for c in range(m)] for i in range(R.rows)]
    rowsf = [[float(R[i, c]) for c in range(m)] for i in range(R.rows)]

    def U_np(theta, pval):
        U = np.zeros(m, dtype=complex); U[gauge] = 1; U[param] = pval
        for c, t in zip(angles, theta):
            U[c] = np.exp(1j*t)
        for row, p in zip(rowsf, piv):
            U[p] = -sum(row[c]*U[c] for c in range(m) if c != p and row[c] != 0)
        return U

    def resid(theta, pval):
        U = U_np(theta, pval); return np.array([abs(U[p])**2-1 for p in piv])
    mp.mp.dps = 230
    zeta = mp.exp(2*mp.pi*1j/55); ZP = [zeta**k for k in range(40)]; lam = mp.pi
    basisL = [mp.re(w)+lam*mp.im(w) for w in ZP]

    def U_mp(theta, pval):
        U = [mp.mpc(0)]*m; U[gauge] = mp.mpc(1); U[param] = pval
        for c, t in zip(angles, theta):
            U[c] = mp.exp(1j*t)
        for row, p in zip(rows, piv):
            U[p] = -sum(mp.mpf(row[c].p)/row[c].q*U[c] for c in range(m) if c != p and row[c] != 0)
        return U

    def F_mp(theta, pval):
        U = U_mp(theta, pval); return mp.matrix([abs(U[p])**2-1 for p in piv])

    def J_mp(theta, pval):
        h = mp.mpf(10)**-100; M = mp.matrix(len(piv), len(theta))
        for k in range(len(theta)):
            tp = list(theta); tm = list(theta); tp[k] += h; tm[k] -= h; col = (F_mp(tp, pval)-F_mp(tm, pval))/(2*h)
            for r in range(len(piv)):
                M[r, k] = col[r]
        return M
    bank = {f'z110^{k}': mp.exp(2*mp.pi*1j*k/110) for k in range(110)}
    g = mp.exp(2*mp.pi*1j/11)
    for j in (1, 2, 3, 4, 5):
        bank[f'mers{j}'] = (2-g**j)/(2-g**(-j))
    b11 = mp.sqrt(-11); bank['thirds'] = (-5+b11)/6; bank['thirdsbar'] = (-5-b11)/6
    rng = np.random.default_rng(5); summary = {}; n = S['n']
    meth = 'lm' if len(piv) >= len(angles) else 'trf'
    for name, pval in bank.items():
        pv = complex(pval); found = {}
        for _ in range(60):
            r = least_squares(resid, rng.uniform(0, 2*np.pi, len(angles)), args=(pv,), method=meth, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=300)
            if r.cost < 1e-22:
                U = U_np(r.x, pv); P = positions(S, U); Z = np.array([[P[v].real, P[v].imag] for v in range(n)])
                D = np.sqrt(((Z[:, None, :]-Z[None, :, :])**2).sum(2)); iu = np.triu_indices(n, 1)
                if D[iu].min() < 1e-6:
                    continue
                sig = tuple(np.round(np.sort(D[iu]), 6))
                if sig not in found:
                    found[sig] = r.x
        heights = []
        for th in found.values():
            theta = [mp.mpf(t) for t in th]
            for _ in range(25):
                rr = F_mp(theta, pval)
                if mp.norm(rr) < mp.mpf(10)**-215:
                    break
                d = mp.lu_solve(J_mp(theta, pval), rr); theta = [t-dk for t, dk in zip(theta, d)]
            U = U_mp(theta, pval); hs = []
            for c in range(m):
                if c in (gauge, param):
                    continue
                v = lll_first(basisL+[mp.re(U[c])+lam*mp.im(U[c])], 200); hs.append(len(str(max(abs(x) for x in v[:41])))-1)
            heights.append(min(hs))
        summary[name] = dict(faithful=len(found), min_height_digits=heights)
    check(all(h >= 3 for s_ in summary.values() for h in s_['min_height_digits']), 'no further direction in L at any scanned parameter')
    return summary


def main(do_scan):
    S = structure(EI21)
    check((S['n'], len(S['E']), S['ncycles'], S['triangles'], S['m'], S['R'].rows) == (21, 38, 10, 0, 18, 8), 'EI21 structure')
    check(automorphisms(EI21) == 2, '|Aut(EI21)| = 2')
    small = small_support_relations(S, 4)
    check(all(len(r) == 4 for r in small) and len(small) == 2, 'no relation of support three; two of support four')
    found, nsol = census(S, 300, 1)
    faithful = {k: v for k, v in found.items() if v['mind'] > 1e-6}
    check(len(faithful) > 10 and all(v['dof'] == 1 for v in faithful.values()), 'one-dimensional family at every faithful sample')
    out = dict(vertices=21, edges=38, four_cycles=10, classes=18, relations=8, free_angles=9, automorphisms=2,
               support4_relations=[{str(k): str(v) for k, v in r.items()} for r in small],
               census=dict(starts=300, converged=nsol, faithful_samples=len(faithful), local_dof=1))
    if do_scan:
        out['lscan'] = lscan(S)
    out['checks'] = base.CHECKS
    return out


if __name__ == '__main__':
    print(json.dumps(main('--lscan' in sys.argv[1:]), indent=1, sort_keys=True, default=str))
