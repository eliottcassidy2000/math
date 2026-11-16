# Wave 39 (mac-mini-2026-11-15-frame). Run from a working directory; w39_lit_setup.py writes lit_base.pkl used by the others.
# Cone LP for the literal span-3 certificate, sound version.
# Generators of the cone at state s are the RAW vectors g = A_s e_j (j in the raw support R), deduplicated by exact
# primitive raw direction; constraints are imposed on raw generators directly:
#   phi_t(O A_l g) - rho phi_s(O g) >= eps |(O g, O A_l g)|   (every live edge s -l-> t, every generator g)
#   phi_s(O g) <= H(g) = (O g)[0]
#   phi_s(O X_w) >= 1 for the actual states of the bottleneck-free initial words w (length Linit).
# phi_s is parametrised on span{O g}.  Writes an exact integer certificate for w39_verify_literal_cert.py.
# usage: python3 conelp3.py AUTO.pkl Linit METHOD rho1 rho2 ...
import pickle, sys, time, itertools, numpy as np
from fractions import Fraction
from scipy.optimize import linprog
from scipy.sparse import csr_matrix
B = pickle.load(open('lit_base.pkl', 'rb'))
A, v0, h, O = B['A'], B['v0'], B['h'], B['O']
N = B['N']
Af = [a.astype(float) for a in A]
Of = O.astype(float)
AU = pickle.load(open(sys.argv[1], 'rb'))
order, edges = AU['order'], AU['edges']
Linit = int(sys.argv[2])
method = sys.argv[3]
rhos = [float(x) for x in sys.argv[4:]]
sidx = {s: i for i, s in enumerate(order)}
tag = sys.argv[1][:-4]
log = open(tag + '_lp3.log', 'a')


def say(*a):
    print(*a, flush=True)
    print(*a, file=log, flush=True)


R = set(np.nonzero(v0)[0])
fr = list(R)
while fr:
    nf = []
    for i in fr:
        for b in range(8):
            for j in np.nonzero(A[b][:, i])[0]:
                if j not in R:
                    R.add(j)
                    nf.append(j)
    fr = nf
R = np.array(sorted(R))
out = {}
for (a, l, b) in edges:
    out.setdefault(a, []).append((l, b))
init = []
for w in itertools.product(range(8), repeat=Linit):
    s = ()
    ok = True
    for l in w:
        nxt = [b for (ll, b) in out.get(sidx[s], []) if ll == l]
        if not nxt:
            ok = False
            break
        s = order[nxt[0]]
    if not ok:
        continue
    X = v0.copy()
    for l in w:
        X = A[l] @ X
    init.append((sidx[s], (O @ X).astype(float)))
live = set(i for i, _ in init)
fr = list(live)
while fr:
    nf = []
    for a in fr:
        for (l, b) in out.get(a, []):
            if b not in live:
                live.add(b)
                nf.append(b)
    fr = nf
live = sorted(live)
t0 = time.time()
G = {}
Q = {}
for a in live:
    M = np.identity(N, dtype=np.int64)
    for l in order[a]:
        M = A[l] @ M
    Gm = M[:, R]
    Gm = Gm[:, Gm.any(axis=0)]
    g = np.gcd.reduce(Gm, axis=0)
    Pm = Gm // g
    _, idx = np.unique(Pm.T, axis=0, return_index=True)
    Gm = Gm[:, np.sort(idx)].astype(float)
    G[a] = Gm
    Y = Of @ Gm
    nz = np.linalg.norm(Y, axis=0) > 0
    if nz.any():
        Yn = Y[:, nz] / np.linalg.norm(Y[:, nz], axis=0)
        U_, sv, _ = np.linalg.svd(Yn, full_matrices=False)
        Q[a] = U_[:, :int(np.sum(sv > 1e-10 * sv[0]))]
    else:
        Q[a] = np.zeros((Of.shape[0], 0))
off = {}
nv = 0
for a in live:
    off[a] = nv
    nv += Q[a].shape[1]
say('auto %s live %d initial %d; raw generators mean %.1f max %d; span dim mean %.1f max %d; variables %d (%.0fs)' % (
    tag, len(live), len(init), np.mean([G[a].shape[1] for a in live]), max(G[a].shape[1] for a in live),
    np.mean([Q[a].shape[1] for a in live]), max(Q[a].shape[1] for a in live), nv, time.time() - t0))


def assemble(rho, eps):
    rows, cols, vals, rhs = [], [], [], []
    r = 0
    for a in live:
        Gm = G[a]
        Y = Of @ Gm
        ka = Q[a].shape[1]
        Z = Q[a].T @ Y
        for (l, b) in out.get(a, []):
            Y2 = Of @ (Af[l] @ Gm)
            nrm = np.sqrt((Y ** 2).sum(axis=0) + (Y2 ** 2).sum(axis=0))
            keep = nrm > 0
            if not keep.any():
                continue
            Z2 = (Q[b].T @ Y2[:, keep]) / nrm[keep]
            Za = Z[:, keep] / nrm[keep]
            kb = Q[b].shape[1]
            gk = int(keep.sum())
            rr = np.arange(r, r + gk)
            rows.append(np.repeat(rr, kb)); cols.append(np.tile(off[b] + np.arange(kb), gk)); vals.append(-Z2.T.ravel())
            rows.append(np.repeat(rr, ka)); cols.append(np.tile(off[a] + np.arange(ka), gk)); vals.append(rho * Za.T.ravel())
            rhs.append(np.full(gk, -eps))
            r += gk
        n1 = np.linalg.norm(Y, axis=0)
        keep = n1 > 0
        gk = int(keep.sum())
        if gk:
            rr = np.arange(r, r + gk)
            rows.append(np.repeat(rr, ka)); cols.append(np.tile(off[a] + np.arange(ka), gk))
            vals.append((Z[:, keep] / n1[keep]).T.ravel())
            rhs.append(Y[0, keep] / n1[keep])
            r += gk
    for (a, y) in init:
        z = Q[a].T @ y
        ka = len(z)
        rows.append(np.full(ka, r)); cols.append(off[a] + np.arange(ka)); vals.append(-z)
        rhs.append(np.array([-1.0]))
        r += 1
    rows = np.concatenate(rows); cols = np.concatenate(cols); vals = np.concatenate(vals)
    m = vals != 0
    return csr_matrix((vals[m], (rows[m], cols[m])), shape=(r, nv)), np.concatenate(rhs), r


for rho in rhos:
    t0 = time.time()
    Am, bb, nr = assemble(rho, 1e-6)
    t1 = time.time()
    res = linprog(np.zeros(nv), A_ub=Am, b_ub=bb, bounds=[(None, None)] * nv, method=method)
    say('rho %.5f rows %d nnz %d assemble %.0fs solve %.0fs status %d %s' % (
        rho, nr, Am.nnz, t1 - t0, time.time() - t1, res.status, res.message[:70]))
    if res.status != 0:
        continue
    # exact integer certificate: phi_s = Q_s c_s scaled by D and rounded
    D = 2 ** 40
    Phi = {}
    for a in live:
        phi = Q[a] @ res.x[off[a]:off[a] + Q[a].shape[1]]
        Phi[a] = [int(round(x * D)) for x in phi]
    fr_rho = Fraction(rho).limit_denominator(100000)
    rr = Fraction(int(rho * 100000) - 1, 100000)     # certify slightly below the LP rate
    cert = {'p': rr.numerator, 'q': rr.denominator, 'D': D, 'Linit': Linit, 'L0': AU['L0'], 'LMAX': AU['LMAX'],
            'unit_ids': AU['unit_ids'], 'order': order, 'edges': edges, 'Phi': Phi}
    pickle.dump(cert, open(tag + '_cert_%.5f.pkl' % rho, 'wb'))
    say('certificate written for p/q = %s' % rr)
