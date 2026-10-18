#!/usr/bin/env python3
"""S41 audit A13 (independent referee code, floating point):
1. The edge counts of every design-loop recipe of golden_forcing_search_20261017.py (base, dM, contact1, dM_C2, MM,
   contact2), rebuilt from the recipe definitions: |A'|, |C'|, unit edges in A' and C', the number of mixed edges
   (1/2) sum_q m_q n_q^2 over unit differences q of C' (m_q = #{a : a, a + dq in A'}, n_q = #{b : b, b + q in C'}),
   and the golden two-distance edge count |E(A')| + (1/2) sum_q m_q.
2. Bonus (THM-4668 / note section 3b): float census of dA + u d1 + v d2 + e d3, e = Z s/h, h = sqrt(4 - N) > 0,
   d3 in D(C) u D(L) minus 0, L = L+ u (-L+), L+ = (0, 1, 1-a, 1/4, 5/4, 5/4-a, 1/2), a = P_5; candidate count.
Usage: python a13_golden_recipes_and_e.py
"""
import os, sys, json, time
sys.dont_write_bytecode = True
os.environ.setdefault('NUMBA_CACHE_DIR', os.path.join(os.environ.get('TEMP', '/tmp'), 's41_audit_numba_cache'))
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from a9_golden_numeric import diffset, census

TOL = 1e-9


def key(x):
    return (int(round(x.real * 1e7)), int(round(x.imag * 1e7)))


class PointSet:
    def __init__(self, pts=()):
        self.pts = []
        self.idx = {}
        for p in pts:
            self.add(p)

    def find(self, x):
        kx, ky = key(x)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                i = self.idx.get((kx + dx, ky + dy))
                if i is not None and abs(self.pts[i] - x) < 1e-8:
                    return i
        return None

    def add(self, x):
        if self.find(x) is None:
            self.idx[key(x)] = len(self.pts)
            self.pts.append(complex(x))


def unit_edges(X):
    X = np.array(X)
    n = 0
    for i in range(len(X) - 1):
        n += int(np.sum(np.abs(np.abs(X[i + 1:] - X[i]) - 1) < TOL))
    return n


def main():
    t0 = time.time()
    z = np.exp(2j * np.pi / 15)
    al = z ** 3
    j = 1j * np.sqrt(11)
    wv = float([r.real for r in np.roots([1, 2, -4, -4]) if abs(r.imag) < 1e-12 and 1 < r.real < 2][0])
    phi = (1 + al + al ** 4).real
    r5 = 2 * phi - 1
    P = [complex(p) for p in [0, 1, 1 + al, 1 + al + al ** 2, -al ** 4, (3 - j * r5) / 8, 3 / 4, (3 + j * r5) / 8]]
    eta = 1 + z ** 5
    rho = (5 + j) / 6
    M = [complex(m) for m in [0, 1, eta, 1 + eta, rho, rho * eta, rho * (1 + eta)]]
    d = 1 + al
    A0 = PointSet(P + [p + m for p in P for m in M]).pts
    C0 = PointSet(M + P[:5]).pts
    UC = PointSet([b - a for a in C0 for b in C0 if abs(a - b) > TOL and abs(abs(b - a) - 1) < TOL]).pts

    def grow(base, steps):
        S = PointSet(base)
        for a in list(base):
            for s in steps:
                S.add(a + s)
        return S.pts
    rec = {}
    rec['base'] = (A0, C0)
    rec['contact1'] = (grow(A0, [d * q for q in UC]), C0)
    AdM = grow(A0, [d * m for m in M])
    rec['dM'] = (AdM, C0)
    rec['dM_C2'] = (AdM, grow(C0, UC))
    rec['MM'] = (AdM, PointSet([a + b for a in M for b in M] + P[:5]).pts)
    rec['contact2'] = (grow(rec['contact1'][0], [d * q for q in UC]), C0)
    claimed = dict(base=(54, 10, 146, 15, 56, 181), dM=(364, 10, 1139, 15, 2142, 1943), contact1=(1099, 10, 3794, 15, 5302, 6222),
                   dM_C2=(364, 174, 1139, 532, 2297178, 1943), MM=(364, 29, 1139, 73, 72358, 1943),
                   contact2=(10597, 10, 44745, 15, 87846, None))
    out = {}
    for name, (A, C) in rec.items():
        SA = PointSet(A)
        SC = PointSet(C)
        Uq = PointSet([b - a for a in C for b in C if abs(a - b) > TOL and abs(abs(b - a) - 1) < TOL]).pts
        mixed2 = 0
        tdx = 0
        for q in Uq:
            m_q = sum(1 for a in A if SA.find(a + d * q) is not None)
            n_q = sum(1 for b in C if SC.find(b + q) is not None)
            mixed2 += m_q * n_q * n_q
            tdx += m_q
        eA = unit_edges(A)
        eC = unit_edges(C)
        got = (len(A), len(C), eA, eC, mixed2 // 2, eA + tdx // 2)
        out[name] = dict(got=got, claimed=claimed[name])
        print('%-9s |A| %5d |C| %3d  A edges %6d  C edges %4d  mixed %8d  two-distance %6d   claimed %s  [%.0f s]'
              % ((name,) + got + (claimed[name], time.time() - t0)), flush=True)
    # e-direction census
    Cc = C0
    a5 = P[5]
    Lp = [0, 1, 1 - a5, 0.25, 1.25, 1.25 - a5, 0.5]
    L = PointSet([complex(x) for x in Lp] + [-complex(x) for x in Lp]).pts
    u = (wv + 3 - j) / (wv + 3 + j)
    Z = d + u
    N = abs(Z) ** 2
    s = 1j * np.sqrt(4 / N - 1)
    v = Z * (-1 + s) / 2
    h = np.sqrt(4 - N)
    e = Z * s / h
    DA = diffset(A0)
    DC = diffset(Cc)
    D3 = PointSet(list(DC) + list(diffset(L))).pts
    a_ = np.concatenate([[0], DA]).astype(np.complex128)
    b_ = (u * np.concatenate([[0], DC])).astype(np.complex128)
    c_ = (v * np.concatenate([[0], DC])).astype(np.complex128)
    e_ = (e * np.array(D3)).astype(np.complex128)
    ncand = len(a_) * len(b_) * len(c_) * len(e_)
    hits, best = census(a_, b_, c_, e_, 1e-6, 0.0)
    pure = [hh for hh in hits if hh[0] == 0 and hh[1] == 0 and hh[2] == 0]
    other = sorted(hh[4] for hh in hits if not (hh[0] == 0 and hh[1] == 0 and hh[2] == 0))
    exact_like = sum(1 for hh in hits if hh[4] < 1e-11)
    unitD3 = sum(1 for x in D3 if abs(abs(x) - 1) < TOL)
    out['e_census'] = dict(abs_e=abs(e), D3=len(D3), candidates=ncand, below_1e11=exact_like, pure=len(pure),
                           unit_D3=unitD3, smallest_other=(other[0] if other else None))
    print('e-direction: |e| = %.15f, |D(C) u D(L)| = %d, candidates %d (= %.3g); ||x|^2-1| < 1e-11: %d (pure e-axis %d, unit '
          'elements of D3: %d); smallest deviation among other candidates %.2e [%.0f s]'
          % (abs(e), len(D3), ncand, ncand, exact_like, len(pure), unitD3, other[0] if other else -1, time.time() - t0))
    json.dump(out, open(os.path.join(HERE, 'a13_golden_recipes_and_e.json'), 'w'), indent=1, default=str)


if __name__ == '__main__':
    main()
