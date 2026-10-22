#!/usr/bin/env python3
"""Targeted search for short vectors of tensor rank 2 in L_raw = E (x) J:
z = x1(x)y1 + x2(x)y2, x_i E8-roots (h_E = 1), y_i minimal vectors of J (h_J = 14):
trace norm(z) = 2 h(z,z) = 56 + 4 Re(h_E(x1,x2) h_J(y1,y2))."""
import json, subprocess, sys, time
import numpy as np
from fractions import Fraction as Fr
cand = int(sys.argv[1]); nsample = int(sys.argv[2]) if len(sys.argv) > 2 else 400
s7 = np.sqrt(7.0)
def cx(p, q): return (float(Fr(p)) + float(Fr(q))*1j*s7)/2
HJ = np.array([[cx(p, q) for p, q in row] for row in json.load(open(f'HJ_cand{cand}.json'))])
i7 = 1j/s7
HE = np.array([[2, 3, -3*i7, -2*i7], [3, 5, -5*i7, -3*i7], [3*i7, 5*i7, 1, 0], [2*i7, 3*i7, 0, 1]])
w = (1 + 1j*s7)/2
def traceform(H):
    m = H.shape[0]
    G = np.zeros((2*m, 2*m))
    G[:m, :m] = 2*H.real; G[m:, :m] = 2*(w*H).real; G[:m, m:] = 2*(np.conj(w)*H).real; G[m:, m:] = 2*(abs(w)**2*H).real
    return np.rint(G).astype(np.int64)
GE, GJ = traceform(HE), traceform(HJ)
def minvecs(G, B):
    m = G.shape[0]
    s = "G=[" + ";".join(",".join(map(str, r)) for r in G) + "];\nv=qfminim(G," + str(B) + ",,2)[3];\n" \
        "for(j=1,#v,print(concat(\"\",Vec(v[,j]))));\nquit\n"
    out = subprocess.run(["gp", "-q", "-s", "1000000000"], input=s, capture_output=True, text=True).stdout
    V = np.array([[int(t) for t in line.strip()[1:-1].split(",")] for line in out.strip().splitlines() if line.startswith("[")])
    return np.vstack([V, -V])
t0 = time.time()
XE = minvecs(GE, 2); YJ = minvecs(GJ, 28)
print(f"E8 roots: {len(XE)}, J minimal vectors: {len(YJ)}  ({time.time()-t0:.1f}s)")
def tocomplex(V, m): return V[:, :m] + w*V[:, m:]
X = tocomplex(XE, 4); Y = tocomplex(YJ, 12)
nE = np.einsum('ij,jk,ik->i', X, HE, np.conj(X)).real; nJ = np.einsum('ij,jk,ik->i', Y, HJ, np.conj(Y)).real
print("h_E on roots:", set(np.round(nE, 6)), " h_J on min vectors:", set(np.round(nJ, 6)))
Bm = X @ HE @ np.conj(X).T          # h_E(x1,x2)
B = np.unique(np.round(Bm.flatten(), 9))
B = B[np.abs(np.abs(B) - 1) > 1e-6]   # drop x2 = unit*x1 (degenerate: pure tensor or zero)
# drop b = 1 (x2 = x1) handled implicitly (pure tensors); keep all
print("distinct h_E(x1,x2) values over root pairs:", len(B), " max |b| over non-collinear pairs:",
      max(abs(b) for b in B if abs(abs(b)-1) > 1e-9))
rng = np.random.default_rng(0)
idx = rng.choice(len(Y), size=min(nsample, len(Y)), replace=False)
best = (1e9, None)
hist = {}
YHc = HJ @ np.conj(Y).T            # 12 x N
for i in idx:
    c = Y[i] @ YHc                  # h_J(y1, y2) for all y2
    vals = (B[:, None] * c[None, :]).real          # Re(b c)
    mn = vals.min(axis=0)
    norms = np.rint(56 + 4*mn).astype(int)
    for nv in np.unique(norms):
        hist[nv] = hist.get(nv, 0) + int((norms == nv).sum())
    j = int(np.argmin(mn))
    if 56 + 4*mn[j] < best[0]:
        bi = int(np.argmin(vals[:, j]))
        best = (56 + 4*mn[j], (int(i), j, complex(B[bi]), complex(c[j])))
print(f"sampled {len(idx)} y1's: smallest trace norm of x1(x)y1 + x2(x)y2 found = {best[0]:.6f}  (b, c) = {best[1][2:]}")
print("histogram of min-over-b norms (norm: count) for small norms:", {k: v for k, v in sorted(hist.items()) if k <= 28})
print(f"time {time.time()-t0:.1f}s")
