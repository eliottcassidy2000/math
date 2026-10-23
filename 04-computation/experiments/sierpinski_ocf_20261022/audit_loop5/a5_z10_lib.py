# Audit (HYP-9371 section 8): exact z^10 coefficients of the loop terms X_e from first principles.
#   [z^L] X_e = sum over cores of excess e (1/|Aut|) sum_(lambda_e >= 1, sum lambda = L) K(core, lambda) hom(core; S^lambda) / n^L
# K = first-principles cluster kernel (a5_kernel.c, written from the THM-4757 section 4 definition), hom = exact integer
# homomorphism number of the core with S^(lambda_e) on edge e (reference orientation; loops: diagonal), S = A - A^T.
# Everything exact (python ints / Fractions); no code from loop_engine*.py or decomp.c.
import os, sys, json, itertools, subprocess
from fractions import Fraction
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.dirname(HERE)
LET = 'abcdefgh'

def load(fn):
    return [(c['V'], [tuple(e) for e in c['edges']], c['aut']) for c in json.load(open(os.path.join(D, fn)))]

CORE_FILES = {1: 'cores_e1.json', 2: 'cores_e2.json', 3: 'cores_e3.json', 4: 'cores_e4.json', 5: 'cores_e5_v5.json'}

def compositions(L, E, loops):
    """lambda in N^E, lambda_e >= 1, sum = L, loops get even lambda (odd powers of a skew matrix have zero diagonal)"""
    def rec(i, rem):
        if i == E:
            if rem == 0: yield ()
            return
        lo = 2 if i in loops else 1
        for l in range(lo, rem - (E - i - 1) + 1):
            if i in loops and l % 2: continue
            for rest in rec(i + 1, rem - l): yield (l,) + rest
    yield from rec(0, L)

def kernels(binary, V, edges, lams):
    if not lams: return []
    inp = f"{V} {len(edges)}\n" + ''.join(f"{a} {b}\n" for a, b in edges) + f"{len(lams)}\n" + \
          ''.join(' '.join(map(str, l)) + '\n' for l in lams)
    out = subprocess.run([binary], input=inp, capture_output=True, text=True, check=True).stdout.split()
    assert len(out) == len(lams)
    return [int(x) for x in out]

def kernel_table(binary, L, excesses, cache_fn=None, vmax=None):
    """{excess: [(core index, V, edges, aut, [(lambda, K), ...]), ...]} for cores with E <= L"""
    cache = json.load(open(cache_fn)) if cache_fn and os.path.exists(cache_fn) else {}
    tab = {}
    for e in excesses:
        rows = []
        for i, (V, edges, aut) in enumerate(load(CORE_FILES[e])):
            if len(edges) > L or (vmax and V > vmax): continue
            key = f"{L}:{e}:{i}"
            loops = {k for k, (a, b) in enumerate(edges) if a == b}
            lams = list(compositions(L, len(edges), loops))
            if key in cache: Ks = cache[key]
            else:
                Ks = kernels(binary, V, edges, lams); cache[key] = Ks
            rows.append((i, V, edges, aut, list(zip(lams, Ks))))
        tab[e] = rows
    if cache_fn: json.dump(cache, open(cache_fn, 'w'))
    return tab

class Homs:
    """exact integer homomorphism numbers hom(core; S^lambda) for one tournament"""
    def __init__(self, S, Lmax):
        S = np.asarray(S, dtype=np.int64)
        self.n = S.shape[0]
        self.P = [np.eye(self.n, dtype=np.int64)]
        for _ in range(Lmax): self.P.append(self.P[-1] @ S)
        self.mx = [int(np.abs(p).max()) for p in self.P]
        assert max(self.mx) < 2 ** 40, self.mx
        self.paths = {}
    def hom(self, V, edges, lam):
        # overflow guard: every intermediate of any contraction order is bounded by n^V prod_e max|S^lambda_e|
        bound = self.n ** V
        for l in lam: bound *= self.mx[l]
        terms, ops = [], []
        for (a, b), l in zip(edges, lam):
            if a == b: terms.append(LET[a]); ops.append(np.diagonal(self.P[l]).copy())
            else: terms.append(LET[a] + LET[b]); ops.append(self.P[l])
        spec = ','.join(terms) + '->'
        if bound >= 2 ** 62:                                  # exact python-int fallback
            return int(np.einsum(spec, *[o.astype(object) for o in ops]))
        key = (V, tuple(edges), self.n)
        if key not in _PATHS:
            _PATHS[key] = np.einsum_path(spec, *[o.astype(float) for o in ops], optimize='greedy')[0]
        return int(np.einsum(spec, *ops, optimize=_PATHS[key]))

_PATHS = {}                                                   # contraction paths depend only on (core, n)

def X_coeff(tab_e, homs, L):
    """exact [z^L] X_e as a Fraction"""
    tot = Fraction(0)
    for i, V, edges, aut, lk in tab_e:
        s = 0
        for lam, K in lk:
            if K: s += K * homs.hom(V, edges, lam)
        tot += Fraction(s, aut)
    return tot / homs.n ** L

def logphi_coeffs(N, n, LM):
    """exact [z^L] log Phi_T(z), Phi_T(z) = sum_L N_L (z/n)^L"""
    phi = [Fraction(int(N[L]), n ** L) if L < len(N) else Fraction(0) for L in range(LM + 1)]
    c = [Fraction(0)] * (LM + 1)
    for L in range(1, LM + 1):
        s = L * phi[L]
        for j in range(1, L): s -= j * c[j] * phi[L - j]
        c[L] = s / L
    return c
