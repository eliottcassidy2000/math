"""Two-sided multiplicativity of H on span-3 ribbon words.
For a standalone word u (first vertex, then letters) and an untruncated continuation v, compare H(uv) with
H(u) * H(T[v]), where T[v] is the tournament induced on the vertices of v (letters truncated at its start).
usage: python3 qm_ratio.py NU NV"""
import pickle, sys, itertools, numpy as np
B = pickle.load(open('lit_base.pkl', 'rb'))
A = [M.astype(np.float64) for M in B['A']]; v0 = B['v0'].astype(np.float64); h = B['h'].astype(np.float64)
NU, NV = int(sys.argv[1]), int(sys.argv[2])
def trunc(word):      # effective letters of a standalone word given as letters of vertices 2..n
    out = []
    for q, b in enumerate(word, start=2):
        out.append(b & ((1 << min(3, q - 1)) - 1))
    return tuple(out)
# does the engine truncate implicitly?
def H(word):
    x = v0.copy()
    for b in word: x = A[b] @ x
    return float(h @ x)
chk = all(H(w) == H(trunc(w)) for w in itertools.product(range(8), repeat=4))
print("engine ignores back arcs that point before the first vertex:", chk)
best = (0, None, None); worst = (1e99, None, None)
for nu in range(1, NU + 1):
    us = sorted({trunc(w) for w in itertools.product(range(8), repeat=nu - 1)})
    X = np.array([v0] * len(us))
    for j in range(nu - 1):
        Y = np.zeros_like(X)
        for b in range(8):
            idx = [i for i, w in enumerate(us) if w[j] == b]
            if idx: Y[idx] = X[idx] @ A[b].T
        X = Y
    Hu = X @ h
    for nv in range(1, NV + 1):
        mx = 0; mn = 1e99
        for v in itertools.product(range(8), repeat=nv):
            Z = X
            for b in v: Z = Z @ A[b].T
            Huv = Z @ h
            tv = trunc(v[1:])                       # standalone tournament on the vertices of v
            Hv = H(tv)
            r = Huv / (Hu * Hv)
            i1 = int(np.argmax(r)); i0 = int(np.argmin(r))
            if r[i1] > mx: mx = r[i1]
            if r[i0] < mn: mn = r[i0]
            if r[i1] > best[0]: best = (float(r[i1]), us[i1], v, Hu[i1], Hv, Huv[i1])
            if r[i0] < worst[0]: worst = (float(r[i0]), us[i0], v, Hu[i0], Hv, Huv[i0])
        print("|u| = %d (%d words), |v| = %d: min ratio %.4f, max ratio %.4f" % (nu, len(us), nv, mn, mx), flush=True)
print("overall max H(uv)/(H(u)H(T[v])) =", best)
print("overall min =", worst)
