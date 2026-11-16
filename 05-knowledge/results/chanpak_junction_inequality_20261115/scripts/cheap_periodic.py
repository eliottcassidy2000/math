"""Exhaustive search for cheap periodic span-3 words.
For a period P: all primitive necklaces u over the 8 letters whose periodic word u^Z has growth below RHO.
Pruning (exact supermultiplicativity, H(T) >= H(T[S]) H(T - S)): growth(u)^P >= H(T[v]) for every prefix v of u
(v standalone, letters truncated), and a periodic word containing three consecutive empty letters has a cut.
For every survivor: growth = spectral radius of the minimal-representation product; stable rank of the minimal
operator (exact, python-flint): rank of O M_u^64 Rb.  Bottleneck-free periodic <=> stable rank >= 2.
usage: python3 cheap_periodic.py P RHO"""
import pickle, sys, itertools, numpy as np, time
import flint
B = pickle.load(open('lit_base.pkl', 'rb'))
A, O, Rb, v0, h = B['A'], B['O'], B['Rb'], B['v0'], B['h']
Af = [M.astype(np.float64) for M in A]
P, RHO = int(sys.argv[1]), float(sys.argv[2])
OR = O @ Rb
pinv = np.linalg.pinv(OR.astype(float))
Mm = [(O @ A[l] @ Rb).astype(float) @ pinv for l in range(8)]
thr = RHO ** P
FA = [flint.fmpz_mat(M.tolist()) for M in A]; FO = flint.fmpz_mat(O.tolist()); FR = flint.fmpz_mat(Rb.tolist())
FI = flint.fmpz_mat([[1 if i == j else 0 for j in range(A[0].shape[0])] for i in range(A[0].shape[0])])
def stable_rank(w):
    """exact: the rank of O M^k Rb is non-increasing in k, and it is the stable rank from k = 61 on.
    Stage 1: k = 8; if the rank there is <= 1 the stable rank is <= 1.  Stage 2: k = 64."""
    M = FI
    for l in w: M = FA[l] * M
    for _ in range(3): M = M * M            # M^8
    r = (FO * M * FR).rank()
    if r <= 1: return r
    for _ in range(3): M = M * M            # M^64
    return (FO * M * FR).rank()
t0 = time.time()
# level-by-level enumeration of prefixes with H(T[prefix]) <= thr, chunked by the first two letters
found = []; nsurv = 0
h_f = h.astype(np.float64); v0f = v0.astype(np.float64)
for a1 in range(8):
  for a2 in range(8):
    words = [(a1, a2)]
    x = v0f.copy(); x = Af[a2 & 1] @ x          # vertex 1 (letter truncated to nothing), vertex 2 (truncated to {1})
    X = np.array([x])
    for m in range(3, P + 1):
        nw = []; nX = []
        for b in range(8):
            be = b & 3 if m == 3 else b
            Y = X @ Af[be].T
            Hv = Y @ h_f
            keep = np.nonzero(Hv <= thr + 1e-9)[0]
            for i in keep:
                w = words[i] + (b,)
                if b == 0 and w[-2] == 0 and w[-3] == 0: continue
                nw.append(w)
            nX.append(Y[[i for i in keep if not (b == 0 and words[i][-1] == 0 and words[i][-2] == 0)]])
        words = nw; X = np.vstack(nX) if nX else np.zeros((0, len(v0)))
        if not words: break
    for w in words:
        # cyclic conditions: canonical rotation, primitive, no 000 across the seam
        rots = [w[i:] + w[:i] for i in range(P)]
        if min(rots) != w: continue
        if any(w == rots[d] for d in range(1, P) if P % d == 0): continue
        ww = w + w[:2]
        if any(ww[i] == 0 and ww[i+1] == 0 and ww[i+2] == 0 for i in range(P)): continue
        nsurv += 1
        Pm = np.identity(Mm[0].shape[0])
        for l in w: Pm = Mm[l] @ Pm
        g = max(abs(np.linalg.eigvals(Pm))) ** (1.0 / P)
        if g < RHO:
            found.append((g, w))
found.sort()
print("period %d, growth < %.4f: %d necklaces survive the prefix pruning, %d cheap cycles (%.0f s)" % (P, RHO, nsurv, len(found), time.time() - t0), flush=True)
pickle.dump({'found': found}, open('cheap_P%d_found.pkl' % P, 'wb'))
bf = []
for g, w in found:
    r = stable_rank(w)
    if r >= 2: bf.append((g, w, r))
print("bottleneck-free among them (stable rank >= 2):", len(bf))
for g, w, r in bf[:40]: print("  growth %.6f  stable rank %d  word %s" % (g, r, w))
pickle.dump({'found': found, 'bf': bf}, open('cheap_P%d.pkl' % P, 'wb'))
