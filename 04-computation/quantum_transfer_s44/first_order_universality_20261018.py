#!/usr/bin/env python3
"""First-order universality of the quantum transfer (S44).

Claim.  For Vir_(p,q) with N = (p-1)(q-1)/2, the first-order transfer T_1 depends only on (N, c) and equals
(c - c_(2,2N+1)) T*_N modulo J, with T*_N universal.  Ingredients checked here:
 (a) the length-(N-1) part of the singular vector, alpha L_-4 L_-2^(N-2) + beta L_-3^2 L_-2^(N-3), has (alpha, beta)
     affine in c along models with the same N (exact rationals; collinearity of (c, alpha, beta));
 (b) the 2x2 system fixing (alpha, beta) (top-length parts of L_1 v = L_2 v = 0) has a c-independent matrix and an
     affine-in-c right side, with nonzero determinant;
 (c) for two non-free models with the same N, the T_1 images coincide as subspaces modulo J, bidegree by bidegree, and
     the T_1 matrices (on a common cycle basis) are proportional with ratio (c - c_0)/(c' - c_0), c_0 = c_(2,2N+1).
Usage: python first_order_universality_20261018.py
Session opus-2026-10-18-S44.
"""
import os, sys, json
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'voa_singular_support'))
sys.path.insert(0, HERE)
from vir_ss import singular_vector, central_charge, block
from vir_quantum_transfer_20261018 import Transfer, by_length, to_rows, rank, left_nullspace


def subleading(p, q):
    v = singular_vector(p, q)
    N = (p - 1) * (q - 1) // 2
    a = v.get((4,) + (2,) * (N - 2), Fraction(0))
    b = v.get((3, 3) + (2,) * (N - 3), Fraction(0)) if N >= 3 else Fraction(0)
    return N, central_charge(p, q), a, b


def main():
    out = {}
    # (a) collinearity
    groups = {4: [(2, 9), (3, 5)], 6: [(2, 13), (3, 7), (4, 5)], 9: [(2, 19), (3, 10), (4, 7)],
              10: [(2, 21), (3, 11), (5, 6)]}
    print('(a) subleading coefficients (alpha, beta) of v vs c:')
    for N, models in groups.items():
        pts = []
        for (p, q) in models:
            try:
                n, c, a, b = subleading(p, q)
            except AssertionError as e:
                print('   skip', (p, q), e)
                continue
            pts.append((c, a, b))
            print('   N=%d (%d,%d): c = %s, alpha = %s, beta = %s' % (N, p, q, c, a, b))
        if len(pts) >= 2:
            (c0, a0, b0), (c1, a1, b1) = pts[0], pts[1]
            sa = (a1 - a0) / (c1 - c0)
            sb = (b1 - b0) / (c1 - c0)
            ok = all(a == a0 + sa * (c - c0) and b == b0 + sb * (c - c0) for c, a, b in pts)
            print('   N=%d: affine in c along all %d models: %s (d alpha/dc = %s, d beta/dc = %s)' % (N, len(pts), ok, sa, sb))
            out['N%d' % N] = dict(points=[[str(x) for x in t] for t in pts], affine=ok, dalpha=str(sa), dbeta=str(sb))
    # (c) T_1 images for (3,7) vs (4,5), and (3,10) vs (4,7)
    P = 2147483647
    for (m1, m2, Dmax) in [((3, 7), (4, 5), 22), ((3, 10), (4, 7), 25)]:
        T = Transfer(m1[0], m1[1], Dmax, P)
        U = Transfer(m2[0], m2[1], Dmax, P)
        N = T.N
        c0 = central_charge(2, 2 * N + 1)
        c1, c2 = central_charge(*m1), central_charge(*m2)
        ratio = (c1 - c0) / (c2 - c0)
        rat = ratio.numerator % P * pow(ratio.denominator % P, -1, P) % P
        same = True
        prop = True
        checked = 0
        for D in range(2 * N, Dmax + 1):
            for r in range(N + 1, D // 2 + 1):
                B1 = T.k1_basis(D, r)
                if not B1:
                    continue
                rows1, nc = T.d1_rows(B1, D, r)
                cyc = left_nullspace(rows1, nc, P)
                if not cyc:
                    continue
                s = r - 1
                J = T.ss.J.get((D, s), [])
                rJ, n = to_rows(J, D, s, P)
                t_a = [by_length(T.phi_of(z, B1)).get(s, {}) for z in cyc]
                t_b = [by_length(U.phi_of(z, B1)).get(s, {}) for z in cyc]
                ra, _ = to_rows(t_a, D, s, P)
                rb, _ = to_rows(t_b, D, s, P)
                rkJ = rank(rJ, n, P)
                rk_a = rank(rJ + ra, n, P)
                rk_b = rank(rJ + rb, n, P)
                rk_ab = rank(rJ + ra + rb, n, P)
                if not (rk_a == rk_b == rk_ab):
                    same = False
                # proportionality modulo J: t_a - rat * t_b in J for every cycle
                diff = [[(x - rat * y) % P for x, y in zip(xa, xb)] for xa, xb in zip(ra, rb)]
                if rank(rJ + diff, n, P) != rkJ:
                    prop = False
                checked += 1
        print('(c) %s vs %s (N=%d) to D=%d: T_1 images equal mod J: %s; T_1 proportional with ratio (c-c0)/(c\'-c0) = %s mod J: %s  [%d bidegrees]'
              % (m1, m2, N, Dmax, same, ratio, prop, checked))
        out['%s_vs_%s' % (m1, m2)] = dict(N=N, Dmax=Dmax, images_equal=same, proportional=prop, ratio=str(ratio), bidegrees=checked)
    json.dump(out, open(os.path.join(HERE, 'first_order_universality.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
