# Four-loop (excess-3) extension of the linear-forest loop engine (THM-4757, HYP-9367).
#
# Cores: cores_enum.py (111 at excess 3, checked against the configuration-model count).  Weight tensors in the
# 5-matrix basis of loop_engine_fast (B0 = R_+, B1 = R_-, B2 = Q_++, B3 = Q_+-, B4 = Q_--): decomp.c, a port of
# loop_engine.decompositions_fast checked key by key against it (decomp_e3/coreNNN.bin, gen_decomp_e3.py).
# A core contributes <W, Gamma> / |G|, Gamma[b_1..b_E] = sum_(f: V -> [n]) prod_e B_(b_e)[f(tail_e), f(head_e)].
#
# Contraction: an exact dynamic program over subsets of the operands (W and the E edge operands; at most 10)
# minimises the total multiply-add count of a pairwise contraction tree under a memory cap; each pairwise step is a
# batched matrix product (numpy matmul).  numpy's own 'greedy'/'optimal' paths were unusable for 6-edge cores.
#
# Doubly regular tournaments (S S^T = nI - J): every propagator f(Shat) with f(0) = 0 equals alpha P + gamma Shat,
# P = I - J/n, alpha = (f(iw) + f(-iw))/2, gamma = (f(iw) - f(-iw))/(2iw), w = n^(-1/2).  So a core is
# <W x_e C2(z), Gamma2> / |G| with the z-independent 2^E tensor Gamma2 of {P, Shat} homomorphism numbers.
# mac-mini-2026-10-22-loop4.
import os, sys, json, subprocess
from math import prod
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

# ---------------------------------------------------------------------------------------------------------------
# cores and weight tensors
_CORES = {}
def load_cores(excess):
    if excess not in _CORES:
        fn = os.path.join(HERE, f'cores_e{excess}.json')
        if not os.path.exists(fn):
            subprocess.run([sys.executable, os.path.join(HERE, 'cores_enum.py'), str(excess)], cwd=HERE, check=True)
        _CORES[excess] = [(c['V'], [tuple(e) for e in c['edges']], c['aut']) for c in json.load(open(fn))]
    return _CORES[excess]

_W = {}
def weight_tensor(excess, i):
    """5-basis weight tensor of core i (int64, shape (5,)*E)"""
    key = (excess, i)
    if key in _W: return _W[key]
    V, edges, _ = load_cores(excess)[i]
    E = len(edges)
    fn = os.path.join(HERE, f'decomp_e{excess}', f'core{i:03d}.bin')
    if not os.path.exists(fn):
        os.makedirs(os.path.dirname(fn), exist_ok=True)
        inp = f"{V} {E}\n" + ''.join(f"{a} {b}\n" for a, b in edges)
        subprocess.run([os.path.join(HERE, 'decomp'), fn], input=inp, capture_output=True, text=True, check=True)
    W = np.fromfile(fn, dtype=np.int64).reshape((5,) * E)
    if excess <= 2: _W[key] = W           # small; the excess-3 tensors (up to 16 MB each) are not cached
    return W

# ---------------------------------------------------------------------------------------------------------------
# DP-optimal pairwise contraction
def dp_plan(op_inds, dims, out_inds=(), memcap=2 ** 25):
    """Exact optimum over all pairwise contraction trees (outer products allowed) of sum of prod(dims of the union of
    the two operands' open indices), with every intermediate (not the output) at most memcap entries.
    op_inds: list of index tuples; returns (cost, tree, opn) with tree = operand number or (left, right)."""
    N = len(op_inds); full = (1 << N) - 1
    occ = {}
    for k, inds in enumerate(op_inds):
        for i in inds: occ[i] = occ.get(i, 0) | (1 << k)
    outset = set(out_inds)
    opn = {}; size = {}
    for S in range(1, full + 1):
        o = frozenset(i for i, m in occ.items() if (m & S) and ((m & ~S) or i in outset))
        opn[S] = o; size[S] = prod(dims[i] for i in o)
    best = {1 << k: (0, k) for k in range(N)}
    for S in sorted(range(1, full + 1), key=lambda s: bin(s).count('1')):
        if bin(S).count('1') < 2 or (S != full and size[S] > memcap): continue
        bc = None; sub = (S - 1) & S
        while sub:
            oth = S ^ sub
            if sub < oth and sub in best and oth in best:
                c = best[sub][0] + best[oth][0] + prod(dims[i] for i in opn[sub] | opn[oth])
                if bc is None or c < bc[0]: bc = (c, (best[sub][1], best[oth][1]))
            sub = (sub - 1) & S
        if bc is not None: best[S] = bc
    if full not in best: raise ValueError('no contraction tree within the memory cap')
    return best[full][0], best[full][1], opn

def _pair(A, ia, B, ib, keep):
    """contract A (indices ia) with B (indices ib); keep = indices that stay open; returns (C, ic)"""
    # sum out indices private to one operand first
    for X, ix, other in ((0, ia, ib), (1, ib, ia)):
        priv = [i for i in ix if i not in other and i not in keep]
        if priv:
            T = A if X == 0 else B
            T = T.sum(axis=tuple(ix.index(i) for i in priv)); ix = tuple(i for i in ix if i not in priv)
            if X == 0: A, ia = T, ix
            else: B, ib = T, ix
    batch = [i for i in ia if i in ib and i in keep]
    contr = [i for i in ia if i in ib and i not in keep]
    fa = [i for i in ia if i not in ib]
    fb = [i for i in ib if i not in ia]
    At = A.transpose([ia.index(i) for i in batch + fa + contr])
    Bt = B.transpose([ib.index(i) for i in batch + contr + fb])
    sb = [A.shape[ia.index(i)] for i in batch]
    sa = [A.shape[ia.index(i)] for i in fa]; sc = [A.shape[ia.index(i)] for i in contr]
    sf = [B.shape[ib.index(i)] for i in fb]
    C = np.matmul(At.reshape(prod(sb), prod(sa), prod(sc)), Bt.reshape(prod(sb), prod(sc), prod(sf)))
    return C.reshape(sb + sa + sf), tuple(batch + fa + fb)

def execute(tree, ops, op_inds, opn, out_inds=()):
    def rec(t):
        if isinstance(t, int): return ops[t], tuple(op_inds[t]), 1 << t
        A, ia, SA = rec(t[0]); B, ib, SB = rec(t[1])
        S = SA | SB
        C, ic = _pair(A, ia, B, ib, opn[S])
        return C, ic, S
    C, ic, _ = rec(tree)
    if out_inds:
        C = C.transpose([ic.index(i) for i in out_inds])
    return C

class Network:
    """graph tensor network of a core: edge operands B[b_e, f(tail), f(head)] (loops: diagonals), optional weight
    tensor W over all b_e; plans cached per (r, n, with_W)."""
    def __init__(self, V, edges):
        self.V, self.edges = V, edges
        self.E = len(edges); self.plans = {}
    def inds(self):
        ops = []
        for e, (a, b) in enumerate(self.edges):
            ops.append((('b', e), ('v', a)) if a == b else (('b', e), ('v', a), ('v', b)))
        return ops
    def plan(self, r, n, with_W):
        key = (r, n, with_W)
        if key not in self.plans:
            ops = self.inds()
            dims = {('b', e): r for e in range(self.E)}
            dims.update({('v', v): n for v in range(self.V)})
            if with_W:
                ops = [tuple(('b', e) for e in range(self.E))] + ops; out = ()
            else:
                out = tuple(('b', e) for e in range(self.E))
            cost, tree, opn = dp_plan(ops, dims, out)
            self.plans[key] = (cost, tree, opn, ops, out)
        return self.plans[key]
    def contract(self, Bstack, W=None):
        """Bstack: (r, n, n) basis matrices.  With W: the scalar <W, Gamma>; without: Gamma (r,)*E."""
        r, n = Bstack.shape[0], Bstack.shape[1]
        diag = np.einsum('kii->ki', Bstack)
        cost, tree, opn, op_inds, out = self.plan(r, n, W is not None)
        edge_ops = [diag if a == b else Bstack for (a, b) in self.edges]
        ops = ([W.astype(Bstack.dtype)] if W is not None else []) + edge_ops
        return execute(tree, ops, op_inds, opn, out)

_NETS = {}
def network(excess, i):
    if (excess, i) not in _NETS:
        V, edges, _ = load_cores(excess)[i]
        _NETS[(excess, i)] = Network(V, edges)
    return _NETS[(excess, i)]

# ---------------------------------------------------------------------------------------------------------------
# general regular tournaments
def basis(Sh, z=1.0):
    n = Sh.shape[0]; I = np.eye(n); X = z * Sh
    Rp = X @ np.linalg.inv(I - X); Rm = -X @ np.linalg.inv(I + X)
    Cm = (I - X @ X) @ np.linalg.inv(I + X @ X)
    return np.array([Rp, Rm, -Rp @ Rp @ Cm, -Rp @ Rm @ Cm, -Rm @ Rm @ Cm])

def X_e(Sh, excess, z=1.0, cores=None, per_core=False):
    """X_excess(z) = sum over cores of <W, Gamma(B(z Shat))> / |G| (complex z allowed)"""
    Bm = basis(Sh.astype(complex) if np.iscomplexobj(z) or isinstance(z, complex) else Sh, z)
    vals = []
    for i, (V, edges, aut) in enumerate(load_cores(excess)):
        if cores is not None and i not in cores: continue
        W = weight_tensor(excess, i)
        vals.append(network(excess, i).contract(Bm, W) / aut)
    return vals if per_core else sum(vals)

def series(Sh, excess, Lmax, K=32, r=0.5, cores=None, per_core=False):
    """[z^L] X_excess for L <= Lmax by a Cauchy sum on |z| = r; uses X(conj z) = conj X(z) and X(-z) = X(z)
    (only even L survive), so K/4 + 1 evaluations."""
    assert K % 4 == 0
    pts = {}
    for k in range(K // 4 + 1):
        z = r * np.exp(2j * np.pi * k / K)
        v = X_e(Sh, excess, z, cores=cores, per_core=True)
        pts[k] = np.array(v)
    vals = []
    for k in range(K):
        kk = k % (K // 2)                                   # X(-z) = X(z)
        if kk <= K // 4: vals.append(pts[kk])
        else: vals.append(np.conj(pts[K // 2 - kk]))         # z_(K/2 - kk) = -conj(z_kk)... conj symmetry
    vals = np.array(vals)                                   # (K, ncores)
    co = np.array([(np.mean(vals * np.exp(-2j * np.pi * np.arange(K) * L / K)[:, None], axis=0) / r ** L).real
                   for L in range(Lmax + 1)])               # (Lmax+1, ncores)
    return co if per_core else co.sum(axis=1)

# ---------------------------------------------------------------------------------------------------------------
# doubly regular tournaments: 2-matrix basis {P, Shat}
def drt_coeffs(n, z=1.0):
    """C2[k] = (alpha_k, gamma_k) for the five basis functions at z (complex allowed)"""
    w = n ** -0.5
    def funcs(lam):
        X = z * lam
        Rp = X / (1 - X); Rm = -X / (1 + X); Cm = (1 - X * X) / (1 + X * X)
        return np.array([Rp, Rm, -Rp * Rp * Cm, -Rp * Rm * Cm, -Rm * Rm * Cm])
    fp, fm = funcs(1j * w), funcs(-1j * w)
    al = (fp + fm) / 2; ga = (fp - fm) / (2j * w)
    if not np.iscomplexobj(z) and not isinstance(z, complex): al, ga = al.real, ga.real
    return np.stack([al, ga], axis=1)                       # (5, 2)

def drt_gamma2(S, excess, i):
    """Gamma2 (2,)*E for core i: {P, Shat} homomorphism numbers of the core"""
    n = S.shape[0]; Sh = S / n
    P = np.eye(n) - np.ones((n, n)) / n
    return network(excess, i).contract(np.array([P, Sh]))

def reduce_W(W, C):
    """W x_e C along every axis (C: 5 x r)"""
    T = W.astype(C.dtype)
    for _ in range(W.ndim):
        T = np.tensordot(T, C, axes=([0], [0]))             # leading axis contracted, new axis appended
    return T

def check_drt(S):
    n = S.shape[0]
    return np.array_equal(S @ S.T, n * np.eye(n, dtype=S.dtype) - np.ones((n, n), dtype=S.dtype))
