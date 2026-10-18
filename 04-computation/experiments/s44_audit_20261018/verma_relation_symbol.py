#!/usr/bin/env python3
"""S44 audit, note on HYP-9315(b): the E_infinity^1 page is gr(ker Phi), ker Phi = kernel of the Verma map
M(c,2N) -> I (the Verma relations at relative levels 2p-1 and 2q-1, by the BGG/Rocha-Caridi structure).
Claim tested: at relative level 2p-1 (p < q) ker Phi is one-dimensional and its top (source-length) symbol is the
Koszul cycle x_0^(p-2) s_1 = x_0^(p-1) e_1 - N x_0^(p-2) x_1 e_0 -- the unique cycle of length N+p-1 and weight
2N+2p-1 -- which is why that cycle never transfers and the degree-N generator list stops at Li 2p-3.
Own engine.  Usage: python verma_relation_symbol.py
"""
import flint
from orders_audit import Model
from vcore import basis_by_length, index_of


def main():
    P = 2147483647
    for (p, q) in [(3, 4), (3, 5), (3, 7), (4, 5), (4, 7), (5, 6), (5, 7), (6, 7)]:
        N = (p - 1) * (q - 1) // 2
        D = 2 * N + 2 * p - 1
        M = Model(p, q, D, P)
        prs = M.pairs(D)
        basis = basis_by_length(D)
        idx = index_of(D)
        data = [0] * (len(prs) * len(basis))
        for i, (n, m) in enumerate(prs):
            for mono, co in M.w(n, m).items():
                data[i * len(basis) + idx[mono]] = co
        A = flint.nmod_mat(len(prs), len(basis), data, P)
        # left kernel: x A = 0  <=>  A^T x^T = 0
        X, nul = A.transpose().nullspace()
        vecs = [[int(X[i, k]) for i in range(len(prs))] for k in range(nul)]
        out = []
        for v in vecs:
            L = max(len(prs[i][1]) for i in range(len(prs)) if v[i])
            top = {prs[i]: v[i] for i in range(len(prs)) if v[i] and len(prs[i][1]) == L}
            a = top.get((1, (2,) * (p - 1)), 0)
            b = top.get((0, (2,) * (p - 2) + (3,)), 0)
            is_euler = (set(top) <= {(1, (2,) * (p - 1)), (0, (2,) * (p - 2) + (3,))}) and a and (b == (-N * a) % P)
            out.append((N + L, len(top), bool(is_euler)))
        print('(%d,%d) N=%d: relative level 2p-1=%d, dim ker Phi = %d; (source length of top, #terms, top == x0^(p-2) s1): %s'
              % (p, q, N, 2 * p - 1, nul, out))


if __name__ == '__main__':
    main()
