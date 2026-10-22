# Tensorized loop engine: one einsum per core.  The decomposition weights of a core (loop_engine.decompositions_fast)
# form a tensor W[d_1, ..., d_E] over per-edge propagator descriptors; stacking the descriptor matrices of edge e as
# P_e[d, :, :] (diagonals for loops), the core sum is einsum(W, P_1, ..., P_E) / |G|.  mac-mini-2026-10-22-sierp.
import sys
import numpy as np
sys.path.insert(0, '.')
from loop_engine import cores, decompositions_fast, automorphisms, propagator2

DESCS = [('m1', 1), ('m1', -1), ('m1', '*')] + [('m2', a, b) for a in (1, -1, '*') for b in (1, -1, '*')]
DIDX = {d: i for i, d in enumerate(DESCS)}

class Core:
    def __init__(self, V, edges):
        self.V, self.edges = V, edges
        self.G = automorphisms(V, edges)
        dec = decompositions_fast(V, edges)
        E = len(edges)
        used = [sorted({DIDX[p[e]] for p in dec}) for e in range(E)]
        self.used = used
        pos = [{d: k for k, d in enumerate(u)} for u in used]
        W = np.zeros([len(u) for u in used])
        for p, w in dec.items():
            W[tuple(pos[e][DIDX[p[e]]] for e in range(E))] += w
        self.W = W
        letters = 'abcdefgh'; dl = 'ijklmnopqrstuvwxyz'
        terms = []
        for e, (a, b) in enumerate(edges):
            terms.append(dl[e] + (letters[a] if a == b else letters[a] + letters[b]))
        self.spec = dl[:E] + ',' + ','.join(terms) + '->'
        self.path = None

    def value(self, Pstack):
        """Pstack[d] = propagator matrix for descriptor d (all 12)"""
        ops = []
        for e, (a, b) in enumerate(self.edges):
            M = Pstack[self.used[e]]
            ops.append(np.einsum('kii->ki', M) if a == b else M)
        if self.path is None:
            self.path = np.einsum_path(self.spec, self.W, *ops, optimize='optimal')[0]
        return np.einsum(self.spec, self.W, *ops, optimize=self.path) / self.G

_CORES = {}
def get_cores(excess):
    if excess not in _CORES: _CORES[excess] = [Core(V, e) for V, e in cores(excess)]
    return _CORES[excess]

def stack(Sh, z=1.0):
    return np.array([propagator2(d, Sh.astype(complex) if np.iscomplexobj(z) else Sh, z) for d in DESCS])

def X(Sh, excess, z=1.0):
    P = stack(Sh, z)
    return sum(c.value(P) for c in get_cores(excess))

if __name__ == '__main__':
    import time
    from two_loop_regular import two_loop
    rng = np.random.default_rng(5); n = 11; d = 5
    A = np.array([[1 if 0 < (j - i) % n <= d else 0 for j in range(n)] for i in range(n)])
    for _ in range(4000):
        u, v, w = rng.choice(n, 3, replace=False)
        if A[u, v] and A[v, w] and A[w, u]: A[u, v] = A[v, w] = A[w, u] = 0; A[v, u] = A[w, v] = A[u, w] = 1
    Sh = (A - A.T) / n
    I = np.eye(n); S2 = Sh @ Sh
    G = -S2 @ np.linalg.inv(I + S2); H = S2 @ np.linalg.inv(I - S2); O = Sh @ np.linalg.inv(I - S2)
    gd, hd = np.diag(G), np.diag(H)
    ref = np.sum(2 * gd * hd - 1.5 * gd ** 2) + np.sum(G * (H * H - O * O - (2 / 3) * G * G)) - (hd - gd) @ G @ (hd - gd)
    t0 = time.time(); x1 = X(Sh, 1); t1 = time.time()
    print(f"X_1: tensor engine {x1:+.12e}, closed form {ref:+.12e} ({t1 - t0:.2f}s)")
    t0 = time.time(); cs = get_cores(2); t1 = time.time()
    print(f"built {len(cs)} excess-2 cores in {t1 - t0:.1f}s; W sizes {[c.W.size for c in cs]}")
    t0 = time.time(); x2 = X(Sh, 2); t1 = time.time()
    print(f"X_2 = {x2:+.12e} ({t1 - t0:.2f}s)")

# ---------------------------------------------------------------------------------------------------------------
# Basis reduction: every descriptor matrix is a combination of five matrices
#   B0 = R_+, B1 = R_-, B2 = Q_(++), B3 = Q_(+-) = Q_(-+), B4 = Q_(--),   Q_(s1 s2) = -R_s1 R_s2 (I - Shat^2)(I + Shat^2)^-1,
# so the weight tensor is transformed axis by axis into a 5^E tensor, and the contraction uses the 5 basis matrices.
def desc_coeffs():
    C = np.zeros((len(DESCS), 5))
    def sgn_list(s): return [1, -1] if s == '*' else [s]
    for i, d in enumerate(DESCS):
        if d[0] == 'm1':
            for s in sgn_list(d[1]): C[i, 0 if s == 1 else 1] += 1
        else:
            for s1 in sgn_list(d[1]):
                for s2 in sgn_list(d[2]):
                    C[i, {(1, 1): 2, (1, -1): 3, (-1, 1): 3, (-1, -1): 4}[(s1, s2)]] += 1
    return C
CMAT = desc_coeffs()

class CoreB(Core):
    def __init__(self, V, edges):
        super().__init__(V, edges)
        Wb = self.W
        for e in range(len(edges)):
            Ce = CMAT[self.used[e]]                      # (len(used_e), 5)
            Wb = np.tensordot(Wb, Ce, axes=([0], [0]))   # contract the leading axis, append a basis axis at the end
        self.Wb = Wb
        self.pathb = None
    def value_b(self, B):
        ops = [np.einsum('kii->ki', B) if a == b else B for (a, b) in self.edges]
        if self.pathb is None:
            self.pathb = np.einsum_path(self.spec, self.Wb, *ops, optimize='greedy')[0]
        return np.einsum(self.spec, self.Wb, *ops, optimize=self.pathb) / self.G

def basis(Sh, z=1.0):
    n = Sh.shape[0]; I = np.eye(n); X = z * Sh
    Rp = X @ np.linalg.inv(I - X); Rm = -X @ np.linalg.inv(I + X)
    Cm = (I - X @ X) @ np.linalg.inv(I + X @ X)
    return np.array([Rp, Rm, -Rp @ Rp @ Cm, -Rp @ Rm @ Cm, -Rm @ Rm @ Cm])

_CORESB = {}
def XB(Sh, excess, z=1.0):
    if excess not in _CORESB: _CORESB[excess] = [CoreB(V, e) for V, e in cores(excess)]
    Bm = basis(Sh.astype(complex) if isinstance(z, complex) else Sh, z)
    return sum(c.value_b(Bm) for c in _CORESB[excess])

# Hand-coded contraction for the tetrahedron K_4 (edges (01),(02),(03),(12),(13),(23)):
#   T1[i,j,l,a,b,c] = B_i[a,b] B_j[a,c] B_l[b,c],  T2[k,m,n,a,b,c] = sum_d B_k[a,d] B_m[b,d] B_n[c,d],
#   value = sum W[i,j,k,l,m,n] sum_abc T1 T2.
K4_EDGES = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
def _k4_value(core, B):
    n = B.shape[1]; r = B.shape[0]
    T1 = np.einsum('iab,jac,lbc->ijlabc', B, B, B).reshape(r ** 3, n ** 3)
    T2 = np.einsum('kad,mbd,ncd->kmnabc', B, B, B, optimize=True).reshape(r ** 3, n ** 3)
    M = T1 @ T2.T                                         # [(i,j,l), (k,m,n)]
    Wt = core.Wb.transpose(0, 1, 3, 2, 4, 5).reshape(r ** 3, r ** 3)
    return np.sum(Wt * M) / core.G
_orig_value_b = CoreB.value_b
def _value_b(self, B):
    if self.V == 4 and list(self.edges) == K4_EDGES: return _k4_value(self, B)
    return _orig_value_b(self, B)
CoreB.value_b = _value_b

# Hand-coded contraction for the doubled square (edges (01),(01),(02),(13),(23),(23)):
#   Y[i0,i1,i2,b,c] = sum_a B_i0[a,b] B_i1[a,b] B_i2[a,c],  Z[i3,i4,i5,b,c] = sum_d B_i3[b,d] B_i4[c,d] B_i5[c,d].
SQ_EDGES = [(0, 1), (0, 1), (0, 2), (1, 3), (2, 3), (2, 3)]
def _sq_value(core, B):
    n = B.shape[1]; r = B.shape[0]
    Y = np.einsum('iab,jab,kac->ijkbc', B, B, B, optimize=True).reshape(r ** 3, n * n)
    Z = np.einsum('lbd,mcd,ncd->lmnbc', B, B, B, optimize=True).reshape(r ** 3, n * n)
    M = Y @ Z.T
    return np.sum(core.Wb.reshape(r ** 3, r ** 3) * M) / core.G
_prev_value_b = CoreB.value_b
def _value_b2(self, B):
    if self.V == 4 and list(self.edges) == SQ_EDGES: return _sq_value(self, B)
    return _prev_value_b(self, B)
CoreB.value_b = _value_b2
