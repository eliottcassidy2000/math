"""Dominance-pruned slow-factor tree for span-3 ribbon words (exploratory, floating-point state vectors).

A factor u is SLOW at rate rho if its operator has rank >= 2 and H(T[u']) < rho^|u'| for every prefix u'.
Finite tree of slow factors (depth L)  ==>  H(W) >= rho^(|W| - L - 1) for every bottleneck-free word W.

Node = (class, x):  class = the essential suffix s of u (the shortest suffix with rank A_s = rank A_u, so that
col A_u = col A_s);  x = raw state vector of the standalone tournament T[u].
Dominance (same depth): (s, x) dominates (s', x') if s is a suffix of s' (col A_s' is contained in col A_s) and
x <= x' coordinatewise.  Then every slow continuation of the second node is a slow continuation of the first,
so only the minimal vectors of each class are expanded.
usage: python3 domtree.py RHO MAXDEPTH NODECAP"""
import pickle, sys, time, numpy as np
B = pickle.load(open('lit_base.pkl', 'rb'))
A, O, Rb, v0, h = B['A'], B['O'], B['Rb'], B['v0'], B['h']
RHO = float(sys.argv[1]); MAXD = int(sys.argv[2]); CAP = int(sys.argv[3])
p = 2147483629
def inv_mod_matrix(M):
    n = M.shape[0]; M = [[int(x) % p for x in r] + [1 if i == j else 0 for j in range(n)] for i, r in enumerate(M)]
    for c in range(n):
        piv = next(i for i in range(c, n) if M[i][c]); M[c], M[piv] = M[piv], M[c]
        iv = pow(M[c][c], p - 2, p); M[c] = [x * iv % p for x in M[c]]
        for i in range(n):
            if i != c and M[i][c]:
                f = M[i][c]; M[i] = [(a - f * b) % p for a, b in zip(M[i], M[c])]
    return np.array([r[n:] for r in M], dtype=np.int64)
def pivots(M):
    M = [[int(x) % p for x in r] for r in M]; rows = len(M); cols = len(M[0]); piv = []; r = 0
    for c in range(cols):
        pr = next((i for i in range(r, rows) if M[i][c]), None)
        if pr is None: continue
        M[r], M[pr] = M[pr], M[r]; iv = pow(M[r][c], p - 2, p); M[r] = [x * iv % p for x in M[r]]
        for i in range(rows):
            if i != r and M[i][c]:
                f = M[i][c]; M[i] = [(a - f * b) % p for a, b in zip(M[i], M[r])]
        piv.append(c); r += 1
        if r == rows: break
    return piv
OR = O @ Rb; J = pivots(OR.tolist()); ORJinv = inv_mod_matrix(OR[:, J])
def modmul(X, Y):
    Xh, Xl = X >> 16, X & 0xFFFF
    return (((Xh @ Y) % p) * 65536 + (Xl @ Y)) % p
Mm = [modmul((O @ A[l] @ Rb)[:, J] % p, ORJinv) for l in range(8)]
d = Mm[0].shape[0]
def col_basis(L):
    L = L % p; basis = []; rows_used = []
    for c in range(L.shape[1]):
        v = L[:, c].copy()
        for b, r0 in zip(basis, rows_used):
            if v[r0]: v = (v - v[r0] * b) % p
        nz = np.nonzero(v)[0]
        if len(nz) == 0: continue
        r0 = int(nz[0]); v = (v * pow(int(v[r0]), p - 2, p)) % p
        basis.append(v); rows_used.append(r0)
    return np.array(basis, dtype=np.int64).T if basis else np.zeros((d, 0), dtype=np.int64)
BASIS = {(): np.identity(d, dtype=np.int64)}
def basis(w):
    b = BASIS.get(w)
    if b is None:
        b = col_basis(modmul(Mm[w[-1]], basis(w[:-1]))); BASIS[w] = b
    return b
def rank(w): return basis(w).shape[1]
ESS = {}
def essential(w):
    e = ESS.get(w)
    if e is None:
        R = rank(w)
        for k in range(1, len(w) + 1):
            if rank(w[-k:]) == R: e = w[-k:]; break
        ESS[w] = e
    return e
Af = [M.astype(np.float64) for M in A]; hf = h.astype(np.float64); v0f = v0.astype(np.float64)
def pareto(X):
    """rows of X that are not coordinatewise >= another row (duplicates removed)."""
    X = np.unique(X, axis=0)
    if len(X) <= 1: return X
    Hs = X @ hf; order = np.argsort(Hs, kind='stable'); X = X[order]
    keep = []
    for i in range(len(X)):
        if keep:
            K = X[keep]
            if np.any(np.all(K <= X[i] * (1 + 1e-12), axis=1)): continue
        keep.append(i)
    return X[keep]
level = {}
for l in range(8):
    s = essential((l,)); level.setdefault(s, []).append(v0f.copy())
level = {s: pareto(np.array(v)) for s, v in level.items()}
t0 = time.time(); total = 0
print("rate %.5f" % RHO, flush=True)
for depth in range(1, MAXD + 1):
    n = sum(len(X) for X in level.values()); total += n
    mn = min(float((X @ hf).min()) for X in level.values()) if level else float('nan')
    maxlen = max((len(s) for s in level), default=0)
    print("depth %3d: %7d nodes in %6d classes (longest class word %2d); least H %.6g ; rho^depth %.6g ; cached words %d  [%.0f s]" %
          (depth, n, len(level), maxlen, mn, RHO ** depth, len(BASIS), time.time() - t0), flush=True)
    if n == 0:
        print("TREE IS FINITE: no slow factor has %d vertices; total nodes expanded %d" % (depth, total)); break
    if n > CAP:
        print("node cap exceeded; stopping (tree not exhausted)"); break
    nxt = {}; m = depth + 1; thr = RHO ** m * (1 + 1e-9)
    for s, X in level.items():
        for l in range(8):
            t = essential(s + (l,))
            if rank(t) <= 1: continue
            le = l & ((1 << min(3, m - 1)) - 1)
            Y = X @ Af[le].T
            Y = Y[(Y @ hf) < thr]
            if len(Y): nxt.setdefault(t, []).append(Y)
    level = {t: pareto(np.vstack(v)) for t, v in nxt.items()}
    # dominance across classes: a class whose word is a suffix of another's contains the latter's column space
    keys = set(level)
    for t in sorted(level, key=len, reverse=True):
        X = level[t]
        for k in range(1, len(t)):
            s = t[-k:]
            if s in keys and s in level and len(level[s]) and len(X):
                K = level[s]
                mask = np.array([not np.any(np.all(K <= x * (1 + 1e-12), axis=1)) for x in X])
                X = X[mask]
        if len(X): level[t] = X
        else: del level[t]
pickle.dump({'ESS': {k: v for k, v in ESS.items()}}, open('domtree_%.4f.pkl' % RHO, 'wb'))
