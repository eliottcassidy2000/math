#!/usr/bin/env python3
"""image_orbit2.py: the orbit of image subspaces col(A_W) of literal span-3 prefix products under the eight letters,
with vectorised row reduction modulo a prime (31-bit) in numpy.  Canonical form = reduced row-echelon basis of the
image.  usage: image_orbit2.py CAP [MAXDEPTH]"""
import sys, numpy as np
A = np.load('A.npy').astype(np.int64); N = A.shape[1]; P = 2147483629
CAP = int(sys.argv[1]); MAXD = int(sys.argv[2]) if len(sys.argv) > 2 else 60
def rref(M):
    M = M % P; M = M.copy(); r = 0; rows, cols = M.shape
    for c in range(cols):
        if r == rows: break
        nz = np.nonzero(M[r:, c])[0]
        if len(nz) == 0: continue
        piv = r + nz[0]
        if piv != r: M[[r, piv]] = M[[piv, r]]
        inv = pow(int(M[r, c]), P - 2, P); M[r] = (M[r] * inv) % P
        col = M[:, c].copy(); col[r] = 0
        M = (M - np.outer(col, M[r])) % P
        r += 1
    return M[:r]
def key(B): return B.tobytes() + bytes([B.shape[0]])
# start: the image of the one-vertex state is the ray v0; the orbit of images of A_W (as a semigroup) starts from the identity
start = rref(np.identity(N, dtype=np.int64)); seen = {key(start): 0}; frontier = [start]
from collections import Counter
for depth in range(1, MAXD + 1):
    nf = []
    for B in frontier:
        for b in range(8):
            img = rref((A[b] @ B.T).T % P)     # images A_b * (basis vectors) as rows
            if img.shape[0] == 0: continue
            k = key(img)
            if k not in seen: seen[k] = depth; nf.append(img)
    frontier = nf
    ranks = Counter(len(bytes_) for bytes_ in [])  # placeholder
    print(f"depth {depth}: distinct image subspaces {len(seen)} (new {len(nf)})", flush=True)
    if not nf: print("CLOSED: the image orbit is finite with", len(seen), "subspaces"); break
    if len(seen) > CAP: print(f"cap {CAP} exceeded: the orbit has more than {CAP} subspaces (not finite within the cap)"); break
rk = Counter(); 
for k, d in seen.items(): rk[k[-1]] += 1
print("rank distribution of the subspaces found:", sorted(rk.items()))
