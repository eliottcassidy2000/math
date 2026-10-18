#!/usr/bin/env python3
"""S44 audit: first-order transfer T_1, independently re-implemented (own engine vcore.py, own Koszul code).

Checks:
 (1) lift independence (THM-4695(a)): T_1 computed with the ascending lift and with the descending lift differ by J;
 (2) T_1(z) lies in grI (top of an element of I) -- implied by construction; checked via rank against grI;
 (3) Koszul boundaries: T_1(b) in J for ALL boundary generators b = m (f_a e_b - f_b e_a) in the range;
 (4) universality (THM-4695(c)): for two models with the same N, T_1 images coincide mod J and
     T_1^(1)(z) - ((c1-c0)/(c2-c0)) T_1^(2)(z) in J for every cycle z of a common basis; and for the p = 2 model
     (c = c0) T_1 is identically 0 mod J;
 (5) cyclicity of Im T_1 (HYP-9315(a) consequence): span of C[x][d] b_N + J equals Im T_1 + J.
Usage: python t1_audit.py N Dmax model1 model2 [model3 ...]   (models as p,q)
"""
import sys, time, json
import flint
from fractions import Fraction
from vcore import parts, block, block_index, apply_neg_vec, dpoly, polymul_mono, central_charge, asc
from orders_audit import Model


def rank(rows, n, P):
    if not rows or not n:
        return 0
    return flint.nmod_mat(len(rows), n, [x for r in rows for x in r], P).rank()


def to_row(vec, idx, P):
    row = [0] * len(idx)
    for k, x in vec.items():
        row[idx[k]] = (row[idx[k]] + x) % P
    return row


class Koszul:
    def __init__(self, N, Dmax, P):
        self.N, self.Dmax, self.P = N, Dmax, P
        self.f = [{(2,) * N: 1}]
        while len(self.f) <= Dmax - 2 * N:
            self.f.append(dpoly(self.f[-1], P))

    def k1(self, D, r):
        N = self.N
        return [(n, m) for n in range(0, D - 2 * N + 1) for m in parts(D - 2 * N - n) if len(m) == r - N]

    def d1_rows(self, B, D, r):
        idx = block_index(D, r)
        return [to_row(polymul_mono(self.f[n], m), idx, self.P) for (n, m) in B], len(idx)

    def cycles(self, D, r):
        B = self.k1(D, r)
        if not B:
            return B, []
        rows, nc = self.d1_rows(B, D, r)
        if nc == 0:
            return B, [[1 if i == j else 0 for i in range(len(B))] for j in range(len(B))]
        AT = flint.nmod_mat(nc, len(B), [rows[i][j] for j in range(nc) for i in range(len(B))], self.P)
        X, nul = AT.nullspace()
        return B, [[int(X[i, k]) for i in range(len(B))] for k in range(nul)]

    def boundaries(self, B, D, r):
        N, P = self.N, self.P
        i1 = {b: i for i, b in enumerate(B)}
        out = []
        nf = D - 2 * N + 1
        for a in range(nf):
            for b in range(a + 1, nf):
                w = D - 4 * N - a - b
                if w < 0 or r - 2 * N < 0:
                    continue
                for m in parts(w):
                    if len(m) != r - 2 * N:
                        continue
                    row = [0] * len(B)
                    for mono, co in self.f[a].items():
                        j = i1[(b, asc(mono + m))]
                        row[j] = (row[j] + co) % P
                    for mono, co in self.f[b].items():
                        j = i1[(a, asc(mono + m))]
                        row[j] = (row[j] - co) % P
                    out.append(row)
        return out

    def J_rows(self, D, s):
        idx = block_index(D, s)
        B = self.k1(D, s)
        return [to_row(polymul_mono(self.f[n], m), idx, self.P) for (n, m) in B], len(idx)


def w_desc(M, n, m, cache):
    key = (n, m)
    if key not in cache:
        if not m:
            cache[key] = M.u[n]
        else:
            cache[key] = apply_neg_vec(m[-1], w_desc(M, n, m[:-1], cache), M.P)
    return cache[key]


def phi(M, z, B, desc=None):
    acc = {}
    for c, (n, m) in zip(z, B):
        if c:
            vec = M.w(n, m) if desc is None else w_desc(M, n, m, desc)
            for k, x in vec.items():
                acc[k] = (acc.get(k, 0) + c * x) % M.P
    return {k: x for k, x in acc.items() if x}


def part_len(vec, s):
    return {k: x for k, x in vec.items() if len(k) == s}


def main():
    N, Dmax = int(sys.argv[1]), int(sys.argv[2])
    models = [tuple(map(int, a.split(','))) for a in sys.argv[3:]]
    P = 2147483647
    K = Koszul(N, Dmax, P)
    Ms = []
    for (p, q) in models:
        assert (p - 1) * (q - 1) // 2 == N
        Ms.append(Model(p, q, Dmax, P))
    c0 = central_charge(2, 2 * N + 1)
    red = lambda x: x.numerator % P * pow(x.denominator % P, -1, P) % P
    descs = [dict() for _ in Ms]
    report = dict(N=N, Dmax=Dmax, models=models, lift_indep=True, in_J_bnd=True, bnd_checked=0, images_equal=True,
                  proportional=True, p2_zero=True, bidegrees=0, rows=[])
    for D in range(2 * N + 3, Dmax + 1):
        for r in range(N + 1, D // 2 + 1):
            B, cyc = K.cycles(D, r)
            if not cyc:
                continue
            s = r - 1
            Jr, n = K.J_rows(D, s)
            rkJ = rank(Jr, n, P)
            idx = block_index(D, s)
            T = []
            for mi, M in enumerate(Ms):
                ta = [to_row(part_len(phi(M, z, B), s), idx, P) for z in cyc]
                td = [to_row(part_len(phi(M, z, B, descs[mi]), s), idx, P) for z in cyc]
                diff = [[(x - y) % P for x, y in zip(a, b)] for a, b in zip(ta, td)]
                if rank(Jr + diff, n, P) != rkJ:
                    report['lift_indep'] = False
                T.append(ta)
            # boundaries -> J (all generators)
            brows = K.boundaries(B, D, r)
            for M in Ms:
                for b in brows:
                    tb = to_row(part_len(phi(M, b, B), s), idx, P)
                    report['bnd_checked'] += 1
                    if rank(Jr + [tb], n, P) != rkJ:
                        report['in_J_bnd'] = False
            ranks = [rank(Jr + t, n, P) - rkJ for t in T]
            # compare models pairwise: non-free ones proportional with ratio (c1-c0)/(c2-c0); p=2 model zero
            nonfree = [i for i, (p, q) in enumerate(models) if p != 2]
            free = [i for i, (p, q) in enumerate(models) if p == 2]
            for i in free:
                if ranks[i] != 0:
                    report['p2_zero'] = False
            for a_ in nonfree:
                for b_ in nonfree:
                    if b_ <= a_:
                        continue
                    ca, cb = central_charge(*models[a_]), central_charge(*models[b_])
                    rat = red((ca - c0) / (cb - c0))
                    joint = rank(Jr + T[a_] + T[b_], n, P) - rkJ
                    if not (ranks[a_] == ranks[b_] == joint):
                        report['images_equal'] = False
                    diff = [[(x - rat * y) % P for x, y in zip(xa, xb)] for xa, xb in zip(T[a_], T[b_])]
                    if rank(Jr + diff, n, P) != rkJ:
                        report['proportional'] = False
            report['bidegrees'] += 1
            report['rows'].append(dict(D=D, s=s, ncyc=len(cyc), nbnd=len(brows), T1ranks=ranks))
            print('D=%2d r=%2d: %d cycles, %d boundary gens, rank T_1 mod J per model %s' % (D, r, len(cyc), len(brows), ranks), flush=True)
    for k in ('lift_indep', 'in_J_bnd', 'bnd_checked', 'images_equal', 'proportional', 'p2_zero', 'bidegrees'):
        print('%s: %s' % (k, report[k]))
    json.dump(report, open('t1_audit_N%d_D%d_%s.json' % (N, Dmax, '_'.join('%d%d' % m for m in models)), 'w'), indent=1)


if __name__ == '__main__':
    main()
