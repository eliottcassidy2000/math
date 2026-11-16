"""The tree of 'bad' factors for the parsing certificate.

A factor u (letters of consecutive vertices, untruncated) is BAD at rate rho if
  (i)  rank of its operator in the minimal representation is >= 2 (no bottleneck inside u), and
  (ii) H(T[u[1..i]]) < rho^i for every prefix (T[.] = induced tournament, i = number of vertices).
If the tree of bad factors is finite, of depth L, then every bottleneck-free word W satisfies
H(W) >= rho^(|W| - L - 1): cut W greedily into prefixes that are not bad; H(T) >= H(T[S]) H(T - S).

Ranks are computed modulo a prime (a lower bound for the true rank), so a node is discarded as a bottleneck only
when its rank mod p is <= 1; such discards are re-checked exactly by verify_badtree.py.
usage: python3 badtree.py RHO LEVELCAP [MAXDEPTH]"""
import pickle, sys, numpy as np, time
B = pickle.load(open('lit_base.pkl', 'rb'))
A, O, Rb, v0, h = B['A'], B['O'], B['Rb'], B['v0'], B['h']
RHO = float(sys.argv[1]); CAP = int(sys.argv[2]); MAXD = int(sys.argv[3]) if len(sys.argv) > 3 else 200
p = 2147483629
def inv_mod_matrix(M, p):
    n = M.shape[0]; M = [[int(x) % p for x in r] + [1 if i == j else 0 for j in range(n)] for i, r in enumerate(M)]
    for c in range(n):
        piv = next(i for i in range(c, n) if M[i][c]); M[c], M[piv] = M[piv], M[c]
        iv = pow(M[c][c], p - 2, p); M[c] = [x * iv % p for x in M[c]]
        for i in range(n):
            if i != c and M[i][c]:
                f = M[i][c]; M[i] = [(a - f * b) % p for a, b in zip(M[i], M[c])]
    return np.array([r[n:] for r in M], dtype=np.int64)
OR = (O @ Rb)
# pivot columns of O Rb mod p
def pivots(M, p):
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
J = pivots(OR.tolist(), p); assert len(J) == OR.shape[0]
ORJinv = inv_mod_matrix(OR[:, J], p)
def modmul(X, Y):                      # exact product mod p of int64 matrices with entries < 2^31
    Xh, Xl = X >> 16, X & 0xFFFF
    return (((Xh @ Y) % p) * 65536 + (Xl @ Y)) % p
Mm = [modmul((O @ A[l] @ Rb)[:, J] % p, ORJinv) for l in range(8)]       # 61 x 61 minimal representation mod p
d = Mm[0].shape[0]
def col_basis(L):
    """column echelon basis of the column space of L mod p (returns d x r matrix)."""
    L = L.copy() % p; basis = []; rows_used = []
    for c in range(L.shape[1]):
        v = L[:, c].copy()
        for b, r0 in zip(basis, rows_used):
            if v[r0]: v = (v - v[r0] * b) % p
        nz = np.nonzero(v)[0]
        if len(nz) == 0: continue
        r0 = int(nz[0]); v = (v * pow(int(v[r0]), p - 2, p)) % p
        basis.append(v); rows_used.append(r0)
    return np.array(basis, dtype=np.int64).T if basis else np.zeros((d, 0), dtype=np.int64)
Af = [M.astype(np.float64) for M in A]; hf = h.astype(np.float64); v0f = v0.astype(np.float64)
# level 1: one vertex with any untruncated letter; operator = Mm[l]; standalone state = v0
level = []
for l in range(8):
    level.append(((l,), col_basis(Mm[l]), v0f.copy()))
t0 = time.time(); total = 0
print("rate %.5f; level cap %d" % (RHO, CAP), flush=True)
stats = []
for depth in range(1, MAXD + 1):
    # nodes in 'level' have 'depth' vertices; keep those that are bad
    bad = []
    nrank = 0; ngood = 0
    for w, L, x in level:
        if L.shape[1] <= 1: nrank += 1; continue
        if hf @ x >= RHO ** depth * (1 - 1e-9): ngood += 1; continue
        bad.append((w, L, x))
    total += len(bad)
    stats.append((depth, len(level), nrank, ngood, len(bad)))
    print("depth %3d: %7d candidates, %7d bottlenecks (rank <= 1 mod p), %7d good prefixes, %7d bad  [%.0f s]" %
          (depth, len(level), nrank, ngood, len(bad), time.time() - t0), flush=True)
    if not bad:
        print("TREE IS FINITE: depth of the bad tree = %d, total bad nodes %d" % (depth - 1, total)); break
    if len(bad) > CAP:
        print("level cap exceeded; stopping (tree not exhausted)"); break
    nxt = []
    for w, L, x in bad:
        m = depth + 1
        for l in range(8):
            le = l & ((1 << min(3, m - 1)) - 1)          # truncated letter for the standalone tournament
            L2 = col_basis(modmul(Mm[l], L))
            nxt.append((w + (l,), L2, Af[le] @ x))
    level = nxt
pickle.dump(stats, open('badtree_stats_%.4f.pkl' % RHO, 'wb'))
