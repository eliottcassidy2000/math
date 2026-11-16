# Wave 41 (mac-mini-2026-11-16-dials). Run in a directory holding lit_base.pkl from w39_lit_setup.py (slowtree) / masks files from gentourng (arc_affine).
# Wave 41: the slow-factor tree of Corollary D (chanpak_junction_inequality_20261115) with DOMINANCE pruning.
# Node = literal factor u (m letters = m vertices).  Kept iff
#   rank M_u >= 2 (minimal 61-dim representation; the column space Im M_u is tracked exactly mod p),
#   H(T[u_<=i]) < rho^i for every prefix length i (slow), and
#   not dominated: no other kept node at the same level has the same column space and a coordinatewise smaller
#   standalone raw state (then every continuation of it has smaller counts and the same future ranks).
# The tree is finite at rate rho > 5^(1/3)  <=>  HYP-9742 (literal span 3) holds with that rate (Corollary D (iii)).
# This exploration uses float64 normalised states Y = X / rho^m (exact recheck needed for a certificate).
# usage: python3 slowtree.py RHO MAXLEVEL [TIME_LIMIT_S]
import pickle, sys, time, numpy as np
from collections import defaultdict
B = pickle.load(open('lit_base.pkl', 'rb'))
A, v0, h, O, Rb = B['A'], B['v0'], B['h'], B['O'], B['Rb']
N = B['N']
rho = float(sys.argv[1]); MAXL = int(sys.argv[2]); TLIM = float(sys.argv[3]) if len(sys.argv) > 3 else 1e9
p = 33554393
# minimal representation mod p: pick 61 columns C of Rb with O Rb_C invertible mod p
OR = (O @ Rb) % p


def rref_mod(M):
    """reduced row echelon form mod p of integer matrix M (rows), returns (R, pivots)."""
    M = M.copy() % p
    r = 0
    piv = []
    rows, cols = M.shape
    for c in range(cols):
        nz = np.nonzero(M[r:, c])[0]
        if len(nz) == 0:
            continue
        i = r + nz[0]
        if i != r:
            M[[r, i]] = M[[i, r]]
        inv = pow(int(M[r, c]), p - 2, p)
        M[r] = (M[r] * inv) % p
        f = M[:, c].copy()
        f[r] = 0
        nzr = np.nonzero(f)[0]
        if len(nzr):
            M[nzr] = (M[nzr] - np.outer(f[nzr], M[r])) % p
        piv.append(c)
        r += 1
        if r == rows:
            break
    return M[:r], piv


# choose independent columns of O Rb
_, piv = rref_mod(OR.T.copy())          # row space of (O Rb)^T = column space of O Rb; pivots index columns of O? no
# columns of OR that are independent: rref of OR gives pivot columns
R_, pc = rref_mod(OR.copy())
assert len(pc) == 61, len(pc)
C = pc
ORC = OR[:, C] % p                      # 61 x 61 invertible mod p
# inverse mod p via rref of [ORC | I]
aug = np.concatenate([ORC, np.identity(61, dtype=np.int64)], axis=1)
Raug, pv = rref_mod(aug)
assert pv[:61] == list(range(61))
ORCinv = Raug[:, 61:]
Mm = []
for l in range(8):
    Mm.append((((O @ A[l]) % p @ (Rb[:, C] % p)) % p @ ORCinv) % p)
# sanity: M_l O x = O A_l x on a reachable vector
x = Rb[:, 5]
for l in range(8):
    assert (((Mm[l] @ ((O @ x) % p)) % p) == ((O @ (A[l] @ x)) % p)).all()


def colspace_key(Bm):
    """canonical basis (RREF of the transpose) of the column space of Bm (61 x k) mod p."""
    R, pv = rref_mod(Bm.T.copy())
    return R, R.tobytes()


Af = [a.astype(np.float64) for a in A]
hf = h.astype(np.float64)
e_empty = np.zeros(N); e_empty[0] = 1.0
t0 = time.time()
level = []
for l in range(8):
    Y = Af[l] @ e_empty / rho                 # one vertex, normalised by rho^1
    Bm, key = colspace_key(Mm[l] % p)
    level.append((Y, Bm, (l,)))
stats = []
for m in range(1, MAXL + 1):
    # prune: rank, slowness, dominance
    groups = defaultdict(list)
    n_good = n_bad = n_dom = 0
    for (Y, Bm, w) in level:
        r = Bm.shape[0]
        if r <= 1:
            n_bad += 1
            continue
        if hf @ Y >= 1.0:
            n_good += 1
            continue
        key = Bm.tobytes()
        lst = groups[key]
        dominated = False
        for (Y2, _, _) in lst:
            if np.all(Y2 <= Y + 1e-12):
                dominated = True
                break
        if dominated:
            n_dom += 1
            continue
        keep = [(Y2, B2, w2) for (Y2, B2, w2) in lst if not np.all(Y <= Y2 + 1e-12)]
        n_dom += len(lst) - len(keep)
        keep.append((Y, Bm, w))
        groups[key] = keep
    kept = [t for lst in groups.values() for t in lst]
    hmax = max((hf @ Y for (Y, _, _) in kept), default=0)
    rk = np.bincount([Bm.shape[0] for (_, Bm, _) in kept], minlength=62) if kept else []
    print('level %d: candidates %d bottleneck %d good %d dominated %d kept %d keys %d max H/rho^m %.4f ranks %s (%.0fs)' % (
        m, len(level), n_bad, n_good, n_dom, len(kept), len(groups), hmax,
        {i: int(c) for i, c in enumerate(rk) if c}, time.time() - t0), flush=True)
    stats.append((m, len(kept)))
    if not kept:
        print('TREE FINITE: every factor of %d letters is good or a bottleneck at rate %.4f' % (m, rho))
        break
    if time.time() - t0 > TLIM:
        print('time limit')
        break
    if m == MAXL:
        break
    nxt = []
    for (Y, Bm, w) in kept:
        for l in range(8):
            Y2 = Af[l] @ Y / rho
            B2 = (Mm[l] @ Bm.T) % p          # images of basis vectors (61 x r)
            R2, _ = rref_mod(B2.T.copy())
            nxt.append((Y2, R2, w + (l,)))
    level = nxt
pickle.dump({'rho': rho, 'stats': stats, 'last': [(w, float(hf @ Y)) for (Y, _, w) in level[:2000]]},
            open('slowtree_%.4f.pkl' % rho, 'wb'))
