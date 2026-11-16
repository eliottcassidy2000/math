"""Beam search for large H(uv)/(H(u) H(T[v])): does the two-sided multiplicativity constant saturate?"""
import pickle, sys, itertools, numpy as np
B = pickle.load(open('lit_base.pkl', 'rb'))
A = [M.astype(object) for M in B['A']]; v0 = B['v0'].astype(object); h = B['h'].astype(object)
Af = [M.astype(np.float64) for M in B['A']]
NU, NVMAX, BEAM = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
def trunc(word):
    return tuple(b & ((1 << min(3, q - 1)) - 1) for q, b in enumerate(word, start=2))
def Hx(word):
    x = v0.copy()
    for b in word: x = A[b].dot(x)
    return int(h.dot(x))
us = sorted({trunc(w) for w in itertools.product(range(8), repeat=NU - 1)})
X = np.array([B['v0'].astype(np.float64)] * len(us))
for j in range(NU - 1):
    Y = np.zeros_like(X)
    for b in range(8):
        idx = [i for i, w in enumerate(us) if w[j] == b]
        if idx: Y[idx] = X[idx] @ Af[b].T
    X = Y
hf = B['h'].astype(np.float64)
Hu = X @ hf
beam = [((), X)]
for nv in range(1, NVMAX + 1):
    cand = []
    for v, Z in beam:
        for b in range(8):
            Z2 = Z @ Af[b].T; v2 = v + (b,)
            Hv = Hx(trunc(v2[1:]))
            r = (Z2 @ hf) / (Hu * Hv)
            i1 = int(np.argmax(r))
            cand.append((float(r[i1]), v2, Z2, i1, Hv))
    cand.sort(key=lambda t: -t[0])
    r, v2, Z2, i1, Hv = cand[0]
    exact = Hx(us[i1] + v2)
    print("|v| = %2d: best ratio %.1f  u = %s v = %s  H(u) = %d  H(T[v]) = %d  H(uv) = %d (exact)" %
          (nv, r, us[i1], v2, Hx(us[i1]), Hv, exact), flush=True)
    beam = [(c[1], c[2]) for c in cand[:BEAM]]
