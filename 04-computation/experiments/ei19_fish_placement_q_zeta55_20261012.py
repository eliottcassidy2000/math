"""EI19 and the Hochberg-O'Donnell fish graph in Q(zeta55), by the direction-space method.

Run from the repository root:  python3 -B <script>           (about two minutes)
                               python3 -B <script> --lscan   (adds the fish L-unit parameter scan, ~40 min)
Dependencies: numpy, scipy, sympy (python-flint LLL), mpmath.

Method.  In a faithful unit-distance realization every unit 4-cycle is a parallelogram, so the
edges fall into direction classes (opposite sides of a 4-cycle have equal direction) and the
cycle space of the graph gives exact integer relations among the class vectors.  With one class
as gauge, realizations are the solutions of |dependent class| = 1 over the free angles.
A relation with three classes, n_a u_a + n_b u_b + n_c u_c = 0, forces the angles between the
three unit vectors (a rigid triangle of side lengths |n_a|, |n_b|, |n_c|), hence forces a
direction ratio of modulus one with a specific real part; such a ratio lies in a field F only
if F contains the corresponding square root.

Results.
EI19 (HoG 51376; 19 vertices, 35 edges, 11 unit 4-cycles): 13 classes, 6 relations, and two
  relations of support three,  u5 = -(3/2)(u12 + u13)  and  u5 = 3(u4 + u9)  (1-indexed classes).
  Hence in EVERY faithful realization  u12/u5 = (-1 +- 2 sqrt(-2))/3  and  u4/u5 = (1 +- sqrt(-35))/6.
  Neither sqrt(-2) nor sqrt(-35) lies in Q(zeta55) (quadratic subfields Q(sqrt5), Q(sqrt-11),
  Q(sqrt-55); witnesses: split primes at which -2, resp. -35, is a non-residue), so EI19 has no
  faithful placement in Q(zeta55) -- PROVED, independent of any realization census.  The census
  (direction space, random starts) finds 32 non-congruent real faithful realizations, all with the
  forced ratios; a vertex-space re-solve confirms the relation without the direction machinery.
FISH (HoG 51404; 23 vertices, 42 edges, 9 unit 4-cycles): 24 classes, 11 relations, no relation of
  support at most three; one of support four (the four spokes to vertices 0,1,2,3 form an inscribed
  parallelogram, hence a rectangle: spokes 0,2 and 1,3 are antipodal); the realization family is
  one-dimensional (local degree of freedom 1 at every faithful sample), so the field question is a
  Diophantine one and stays OPEN; the scan of 117 L-unit parameter values finds no further
  direction in L (LLL at 200 digits).
EI17 (control): 13 classes, 6 relations, no relation of support at most four; direction-space census
  with random starts finds 12 non-congruent real faithful realizations (every one hit many times).
"""
import itertools, json, sys, time
import numpy as np
import mpmath as mp
import sympy as sp
from sympy import ZZ, Rational
from sympy.polys.matrices import DomainMatrix
from scipy.optimize import least_squares

CHECKS = 0


def check(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)


GRAPHS = {
    'EI17': [[11, 12, 15], [10, 12, 16], [11, 13, 15], [10, 14, 16], [9, 13, 15], [8, 14, 16], [7, 9, 14, 15],
             [6, 8, 13, 16], [5, 7, 10, 12], [4, 6, 11, 12], [1, 3, 8, 11], [0, 2, 9, 10], [0, 1, 8, 9],
             [2, 4, 7, 14], [3, 5, 6, 13], [0, 2, 4, 6], [1, 3, 5, 7]],
    'EI19': [[4, 8, 15], [5, 8, 15], [3, 6, 18], [2, 16, 17], [0, 5, 14], [1, 4, 10], [2, 7, 13], [6, 11, 17],
             [0, 1, 10], [11, 12, 16], [5, 8, 16, 18], [7, 9, 14, 18], [9, 14, 15, 17], [6, 14, 16, 17],
             [4, 11, 12, 13], [0, 1, 12, 16, 18], [3, 9, 10, 13, 15], [3, 7, 12, 13, 18], [2, 10, 11, 15, 17]],
    'FISH': [[10, 20, 22], [11, 21, 22], [8, 19, 22], [9, 18, 22], [13, 17, 22], [12, 16, 22], [15, 17, 22],
             [14, 16, 22], [2, 10, 20], [3, 11, 21], [0, 8, 19], [1, 9, 18], [5, 14, 20], [4, 15, 21], [7, 12, 19],
             [6, 13, 18], [5, 7, 20], [4, 6, 21], [3, 11, 15, 19, 20], [2, 10, 14, 18, 21], [0, 8, 12, 16, 18, 21],
             [1, 9, 13, 17, 19, 20], [0, 1, 2, 3, 4, 5, 6, 7]]}


# ------------------------------------------------------------ direction structure
def structure(ADJ):
    n = len(ADJ); adj = {v: set(a) for v, a in enumerate(ADJ)}
    E = sorted({(min(i, j), max(i, j)) for i, nb in enumerate(ADJ) for j in nb}); ei = {e: k for k, e in enumerate(E)}
    cyc = set()
    for a in range(n):
        for b in adj[a]:
            for c in adj[b]:
                if c == a:
                    continue
                for d in adj[c]:
                    if d in (a, b) or a not in adj[d]:
                        continue
                    t = (a, b, c, d); rots = [t[i:]+t[:i] for i in range(4)]; rots += [tuple(reversed(r)) for r in rots]
                    cyc.add(min(rots))
    tri = sum(1 for a, b in E for c in adj[a] & adj[b])
    parent = list(range(len(E))); sign = [1]*len(E)

    def find(x):
        if parent[x] == x:
            return x, 1
        r, s = find(parent[x]); parent[x] = r; sign[x] *= s
        return r, sign[x]

    def dirsym(a, b):
        return (ei[(a, b)], 1) if a < b else (ei[(b, a)], -1)

    def union(x, sx, y, sy):
        rx, ax = find(x); ry, ay = find(y)
        if rx == ry:
            return
        parent[ry] = rx; sign[ry] = ax*sx*sy*ay
    for a, b, c, d in cyc:
        (x, sx), (y, sy) = dirsym(a, b), dirsym(d, c); union(x, sx, y, sy)
        (x, sx), (y, sy) = dirsym(b, c), dirsym(a, d); union(x, sx, y, sy)
    roots = sorted({find(k)[0] for k in range(len(E))}); cls = {r: i for i, r in enumerate(roots)}; m = len(roots)
    edge_dir = {}
    for e in E:
        r, sg = find(ei[e]); edge_dir[e] = (cls[r], sg)
    pos = {0: np.zeros(m, dtype=int)}; tree = set()

    def dvec(a, b):
        v = np.zeros(m, dtype=int)
        if a < b:
            c, s = edge_dir[(a, b)]; v[c] = s
        else:
            c, s = edge_dir[(b, a)]; v[c] = -s
        return v
    while len(pos) < n:
        for v in list(pos):
            for w in adj[v]:
                if w in pos:
                    continue
                pos[w] = pos[v]+dvec(v, w); tree.add((min(v, w), max(v, w)))
    rels = [pos[j]-pos[i]-dvec(i, j) for (i, j) in E if (i, j) not in tree]
    rels = [r for r in rels if r.any()]
    M = sp.Matrix([[int(x) for x in r] for r in rels])
    R, piv = M.rref()
    R = sp.Matrix(R[:len(piv), :])
    return dict(n=n, E=E, cycles=sorted(cyc), triangles=tri, m=m, edge_dir=edge_dir, pos=pos, R=R, pivots=list(piv), ncycles=len(cyc))


def small_support_relations(S, kmax=4):
    R = S['R']; m = S['m']; out = []
    for k in range(2, kmax+1):
        for sub in itertools.combinations(range(m), k):
            others = [c for c in range(m) if c not in sub]
            for y in R[:, others].T.nullspace():
                v = (y.T*R)
                if all(v[0, c] != 0 for c in sub):
                    out.append({c+1: v[0, c] for c in sub})
    return out


def solver(S):
    m = S['m']; R = S['R']; piv = S['pivots']; free = [c for c in range(m) if c not in piv]
    gauge = free[0]; free_angles = [c for c in free if c != gauge]
    rows = [[float(R[i, c]) for c in range(m)] for i in range(R.rows)]

    def classes_from(theta, extra=None):
        U = np.zeros(m, dtype=complex); U[gauge] = 1
        for c, t in zip(free_angles, theta):
            U[c] = np.exp(1j*t)
        if extra:
            for c, val in extra.items():
                U[c] = val
        for row, p in zip(rows, piv):
            U[p] = -sum(row[c]*U[c] for c in range(m) if c != p and row[c] != 0)
        return U

    def resid(theta):
        U = classes_from(theta); return np.array([abs(U[p])**2-1 for p in piv])
    return gauge, free_angles, classes_from, resid


def positions(S, U):
    return {v: complex(np.dot(S['pos'][v].astype(float), U)) for v in S['pos']}


def census(S, nstart, seed):
    gauge, free_angles, classes_from, resid = solver(S)
    rng = np.random.default_rng(seed); found = {}; n = S['n']; nsol = 0
    meth = 'lm' if len(S['pivots']) >= len(free_angles) else 'trf'
    for _ in range(nstart):
        th0 = rng.uniform(0, 2*np.pi, len(free_angles))
        r = least_squares(resid, th0, method=meth, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=300)
        if r.cost < 1e-22:
            nsol += 1
            U = classes_from(r.x); P = positions(S, U)
            Z = np.array([[P[v].real, P[v].imag] for v in range(n)])
            D = np.sqrt(((Z[:, None, :]-Z[None, :, :])**2).sum(2)); iu = np.triu_indices(n, 1)
            mind = float(D[iu].min())
            J = np.zeros((len(S['pivots']), len(free_angles))); h = 1e-7
            for k in range(len(free_angles)):
                e = np.zeros(len(free_angles)); e[k] = h; J[:, k] = (resid(r.x+e)-resid(r.x-e))/(2*h)
            dof = len(free_angles)-int((np.linalg.svd(J, compute_uv=False) > 1e-6).sum())
            check(max(abs(abs(P[j]-P[i])-1) for i, j in S['E']) < 1e-9, 'every edge unit in the reconstructed realization')
            sig = tuple(np.round(np.sort(D[iu]), 6))
            if sig not in found:
                found[sig] = dict(hits=0, theta=r.x.tolist(), mind=mind, dof=dof, U=[[u.real, u.imag] for u in U],
                                  z={v: [P[v].real, P[v].imag] for v in range(n)})
            found[sig]['hits'] += 1
    return found, nsol


def lll_first(values, digits):
    Sc = mp.mpf(10)**digits; k = len(values); rows = []
    for i in range(k):
        r = [0]*k; r[i] = 1; r.append(int(mp.nint(values[i]*Sc))); rows.append(r)
    M = DomainMatrix([[ZZ(v) for v in r] for r in rows], (k, k+1), ZZ)
    return [int(c) for c in M.lll().to_list()[0]][:k]


def nonresidue_witness(d):
    """A prime p = 1 mod 55 at which d is a quadratic non-residue: then sqrt(d) is not in Q(zeta55)."""
    p = 56
    while True:
        p += 55
        if sp.isprime(p) and sp.legendre_symbol(d % p, p) == -1:
            return p


# ------------------------------------------------------------------------ main
def ei19(S):
    out = {}
    check((S['n'], len(S['E']), S['ncycles'], S['triangles'], S['m'], S['R'].rows) == (19, 35, 11, 0, 13, 6), 'EI19 structure')
    small = small_support_relations(S, 3)
    out['support3_relations'] = [{str(k): str(v) for k, v in r.items()} for r in small]
    want1 = {5: 1, 12: Rational(3, 2), 13: Rational(3, 2)}
    want2 = {4: -3, 5: 1, 9: -3}
    norm = lambda r: {k: v/r[5] for k, v in r.items()} if 5 in r else None
    check(any(norm(r) == want1 for r in small) and any(norm(r) == want2 for r in small), 'the two forced relations')
    # forced ratios: u12/u5 and u13/u5 are unit vectors with sum -2/3; u4/u5 and u9/u5 are unit vectors with sum 1/3
    out['forced_ratios'] = {'u12/u5': '(-1 +- 2 sqrt(-2))/3', 'u4/u5': '(1 +- sqrt(-35))/6'}
    out['witness_primes'] = {'-2': nonresidue_witness(-2), '-35': nonresidue_witness(-35)}
    # which edges carry the classes (0-indexed HoG labels)
    cls_edges = {}
    for e, (c, s) in S['edge_dir'].items():
        cls_edges.setdefault(c+1, []).append((e, s))
    out['class_edges'] = {str(c): [[list(e), s] for e, s in cls_edges[c]] for c in (4, 5, 9, 12, 13)}
    # census in direction space
    found, nsol = census(S, 3000, 1)
    faithful = {k: v for k, v in found.items() if v['mind'] > 1e-6}
    out['census'] = dict(starts=3000, converged=nsol, distinct=len(found), faithful=len(faithful),
                         all_isolated=all(v['dof'] == 0 for v in faithful.values()), min_hits=min(v['hits'] for v in faithful.values()))
    target = (-1+2*np.sqrt(2)*1j)/3
    for v in faithful.values():
        U = [complex(*u) for u in v['U']]
        q = U[11]/U[4]; q2 = U[3]/U[4]
        check(abs(q-target) < 1e-9 or abs(q-target.conjugate()) < 1e-9, 'forced ratio u12/u5 in every realization')
        check(abs(q2.real-1/6) < 1e-9 and abs(abs(q2.imag)-np.sqrt(35)/6) < 1e-9, 'forced ratio u4/u5 in every realization')
    # independent vertex-space confirmation (35 unit-edge equations only), from perturbed starts
    E = S['E']; n = 19
    def vres(x):
        z = x.reshape(n, 2); d = z[[i for i, j in E]]-z[[j for i, j in E]]; return (d*d).sum(1)-1
    rng = np.random.default_rng(0); confirmed = 0
    for v in list(faithful.values())[:8]:
        z0 = np.array([v['z'][k] for k in range(n)])
        r = least_squares(vres, (z0+rng.normal(0, 3e-3, z0.shape)).ravel(), method='trf', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=3000)
        z = r.x.reshape(n, 2); P = {k: complex(*z[k]) for k in range(n)}
        U5 = P[6]-P[2]; U12 = P[13]-P[6]; U13 = P[2]-P[18]   # class 5 = edge (2,6); class 12 contains (6,13); class 13 contains (2,18) with sign -1
        if r.cost < 1e-20 and abs(U5+1.5*(U12+U13)) < 1e-9:
            confirmed += 1
    check(confirmed == 8, 'vertex-space re-solve confirms the forced relation')
    out['vertex_space_confirmations'] = confirmed
    return out


def fish(S, lscan):
    out = {}
    check((S['n'], len(S['E']), S['ncycles'], S['triangles'], S['m'], S['R'].rows) == (23, 42, 9, 0, 24, 11), 'fish structure')
    small = small_support_relations(S, 4)
    check(len(small) == 1 and set(small[0]) == {2, 4, 6, 8}, 'exactly one relation of support at most four: the spoke rectangle')
    out['support_le4_relations'] = [{str(k): str(v) for k, v in r.items()} for r in small]
    found, nsol = census(S, 300, 3)
    faithful = {k: v for k, v in found.items() if v['mind'] > 1e-6}
    check(len(faithful) > 10 and all(v['dof'] == 1 for v in faithful.values()), 'one-dimensional family at every faithful sample')
    # rectangle check: spokes to 0,2 antipodal and 1,3 antipodal in every faithful sample
    for v in faithful.values():
        P = {k: complex(*v['z'][k]) for k in range(23)}
        check(abs((P[0]-P[22])+(P[2]-P[22])) < 1e-9 and abs((P[1]-P[22])+(P[3]-P[22])) < 1e-9, 'inscribed parallelogram is a rectangle')
    out['census'] = dict(starts=300, converged=nsol, faithful_samples=len(faithful), local_dof=1)
    if lscan:
        out['lscan'] = fish_lscan(S)
    return out


def fish_lscan(S):
    """Fix one free direction to an L-unit from a bank and test the other classes for membership in L."""
    m = S['m']; piv = S['pivots']; R = S['R']; free = [c for c in range(m) if c not in piv]
    gauge, param = free[0], free[1]; angles = [c for c in free if c not in (gauge, param)]
    rowsf = [[float(R[i, c]) for c in range(m)] for i in range(R.rows)]
    rows = [[Rational(R[i, c]) for c in range(m)] for i in range(R.rows)]

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
    rng = np.random.default_rng(5); summary = {}
    for name, pval in bank.items():
        pv = complex(pval); found = {}
        for _ in range(60):
            r = least_squares(resid, rng.uniform(0, 2*np.pi, len(angles)), args=(pv,), method='lm', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=300)
            if r.cost < 1e-22:
                U = U_np(r.x, pv); P = positions(S, U); Z = np.array([[P[v].real, P[v].imag] for v in range(23)])
                D = np.sqrt(((Z[:, None, :]-Z[None, :, :])**2).sum(2)); iu = np.triu_indices(23, 1)
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
            U = U_mp(theta, pval)
            hs = []
            for c in range(m):
                if c in (gauge, param):
                    continue
                v = lll_first(basisL+[mp.re(U[c])+lam*mp.im(U[c])], 200)
                hs.append(len(str(max(abs(x) for x in v[:41])))-1)
            heights.append(min(hs))
        summary[name] = dict(faithful=len(found), min_height_digits=heights)
    check(all(h >= 3 for s in summary.values() for h in s['min_height_digits']), 'no further direction in L at any scanned parameter')
    return summary


def main(lscan):
    out = {}
    S17 = structure(GRAPHS['EI17'])
    check((S17['m'], S17['R'].rows, S17['ncycles']) == (13, 6, 9) and not small_support_relations(S17, 4), 'EI17: no forced relation of support <= 4')
    found, nsol = census(S17, 3000, 1)
    faithful = {k: v for k, v in found.items() if v['mind'] > 1e-6}
    check(len(faithful) == 12 and all(v['dof'] == 0 for v in faithful.values()), 'EI17: twelve isolated real faithful realizations in direction space')
    out['EI17'] = dict(classes=13, relations=6, census=dict(starts=3000, converged=nsol, faithful=12, min_hits=min(v['hits'] for v in faithful.values())))
    out['EI19'] = ei19(structure(GRAPHS['EI19']))
    out['FISH'] = fish(structure(GRAPHS['FISH']), lscan)
    out['checks'] = CHECKS
    return out


if __name__ == '__main__':
    print(json.dumps(main('--lscan' in sys.argv[1:]), indent=1, sort_keys=True, default=str))
