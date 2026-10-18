#!/usr/bin/env python3
"""Independent affine sl3 check (audit2): V^k(sl3) vacuum module, maximal submodule = radical of the contravariant form
(all weight blocks), simple-quotient character vs Dhat(q^u,z)/Dhat(q,z) = PE[sum_{u not| n} q^n chi_adj] (u = 2 for
k = -3/2, u = 4 for k = -9/4), and the singular-support / arc-space analysis (ker pi, minimal generators).
usage: python run_sl3_audit.py k Dmax [u]"""
import sys, time, itertools
from fractions import Fraction
from vamod import ModeVA, GradedAnalysis

# basis of sl3: E_ij (i != j) and H1 = E11-E22, H2 = E22-E33 ; types 0..7
PAIRS = [(0, 1), (1, 2), (0, 2), (1, 0), (2, 1), (2, 0)]
NT = 8  # 6 root vectors + H1 + H2


def mat_of(t):
    M = [[0] * 3 for _ in range(3)]
    if t < 6:
        i, j = PAIRS[t]; M[i][j] = 1
    elif t == 6:
        M[0][0] = 1; M[1][1] = -1
    else:
        M[1][1] = 1; M[2][2] = -1
    return M


def mm(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


def decompose(M):
    """express traceless 3x3 matrix in the basis; returns dict t -> coeff"""
    out = {}
    for t, (i, j) in enumerate(PAIRS):
        if M[i][j]:
            out[t] = Fraction(M[i][j])
    # diagonal: d = diag(a,b,c), a+b+c = 0 -> a H1 + (a+b) H2
    a, b = M[0][0], M[1][1]
    if a:
        out[6] = Fraction(a)
    if a + b:
        out[7] = Fraction(a + b)
    return out


BR = {}
TR = {}
for s, t in itertools.product(range(NT), repeat=2):
    A, B = mat_of(s), mat_of(t)
    C = [[mm(A, B)[i][j] - mm(B, A)[i][j] for j in range(3)] for i in range(3)]
    BR[(s, t)] = decompose(C)
    TR[(s, t)] = sum(mm(A, B)[i][i] for i in range(3))


def weight_code(t):
    if t >= 6:
        return 0
    i, j = PAIRS[t]
    v = [0, 0, 0]; v[i] += 1; v[j] -= 1
    return v[0] * 100 + v[1]


TRANSPOSE = {t: (PAIRS.index((PAIRS[t][1], PAIRS[t][0])) if t < 6 else t) for t in range(NT)}


class AffSl3(ModeVA):
    def __init__(self, k, P=2147483647):
        super().__init__(P)
        self.k = Fraction(k)
        self.gens = [(f't{t}', 2, 0, weight_code(t)) for t in range(NT)]

    def omega(self, X):
        return (TRANSPOSE[X[0]], -X[1])

    def bracket_apply(self, X, Y, rest):
        (s, a2), (t, b2) = X, Y
        m = a2 // 2
        terms = [((u, a2 + b2), self.F(co)) for u, co in BR[(s, t)].items()]
        cen = self.F(self.k * m * TR[(s, t)]) if a2 + b2 == 0 else 0
        return self.lie_bracket_apply(terms, cen, rest)


if __name__ == '__main__':
    k = Fraction(sys.argv[1]); Dmax = int(sys.argv[2]); u = int(sys.argv[3]) if len(sys.argv) > 3 else None
    va = AffSl3(k)
    t0 = time.time()
    G = va.gram_all(2 * Dmax)
    I = va.radical(G)
    print(f'sl3 k={k}: gram+radical to weight {Dmax} in {time.time()-t0:.1f}s', flush=True)
    dimL = {}
    for (W2, ch), (B, Gm) in G.items():
        dimL[(W2 // 2, ch)] = len(B) - len(I[(W2, ch)])
    print('dim L by weight:', [sum(v for (w, ch), v in dimL.items() if w == d) for d in range(Dmax + 1)])
    print('dim L zero-weight:', [dimL.get((d, 0), 0) for d in range(Dmax + 1)])
    if u:
        # PE[sum_{u not| n} q^n chi_adj] charge-resolved
        from collections import defaultdict
        roots = [weight_code(t) for t in range(6)]
        ser = defaultdict(int); ser[(0, 0)] = 1
        letters = []
        for n in range(1, Dmax + 1):
            if n % u == 0:
                continue
            for r in roots + [0, 0]:
                letters.append((n, r))
        for (n, r) in letters:
            new = defaultdict(int)
            for (w, ch), v in ser.items():
                kk = 0
                while w + n * kk <= Dmax:
                    new[(w + n * kk, ch + r * kk)] += v
                    kk += 1
            ser = new
        bad = [(key, dimL.get(key, 0), ser.get(key, 0)) for key in set(dimL) | set(ser) if key[0] <= Dmax and dimL.get(key, 0) != ser.get(key, 0)]
        print(f'simple quotient vs Dhat(q^{u},z)/Dhat(q,z) (all weight blocks):', 'MATCH' if not bad else f'MISMATCH {sorted(bad)[:6]}')
    t0 = time.time()
    GA = GradedAnalysis(va, I, 2 * Dmax)
    GA.run()
    summ = GA.summary_by_weight()
    print(f'graded analysis in {time.time()-t0:.1f}s')
    for W2 in sorted(summ):
        o = summ[W2]
        gens = o['gens']
        li0 = sum(g[2] for g in gens if g[0] == 0)
        nonarc = [g for g in gens if g[0] > 0]
        print(f" weight {W2//2}: dimV {o['V']} dimI {o['I']} dimL {o['V']-o['I']} arcHS {o['V']-o['J']} ker {o['ker']} Li0-gens {li0} non-arc {nonarc}")
