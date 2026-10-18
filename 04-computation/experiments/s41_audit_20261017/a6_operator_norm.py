#!/usr/bin/env python3
"""S41 audit A6 (independent referee code): the operator-norm claim of THM-4667 statement 7 (ideal transport on
last-digit-balanced vectors has norm 1.283, 1.324, 1.335, 1.346 for d = 1..4).

In the energy norm ||K||_d^2 = sum_b 2^-b sum_r K(r,b)^2/(4 2^d), put y(s,b) = 2^(-b/2) K(s,b)/sqrt(4 2^d). The one-step
ideal transport becomes the causal block-Toeplitz operator y -> 2^(-1/2) sum_o (-1)^(o+1) 2^(-o/2) S^o (x) Pi_o, Pi_o the
injection s -> F_o(s), Z/3^d -> Z/3^(d+1). Its norm is sup_|z|=1 ||T(z)|| with T(z) = 2^(-1/2) sum_o (-1)^(o+1)
(z/sqrt2)^o Pi_o (period 2 3^(d+1) in o folded geometrically).

1. Own symbol (forward maps, own fold), sup over a fine circle grid with local refinement, on the balanced subspace
   {x : sum_i x(s + i 3^(d-1)) = 0} and on the full space.
2. A check that avoids the symbol: the L-band truncation of the operator itself (block lower-triangular matrix acting on
   balanced vectors), largest singular value for growing L.
3. An explicit finite real balanced input sequence K_(d-1) whose energy grows under one ideal step.
4. Invariance: T maps balanced vectors to balanced vectors (R8.2).
Usage: python a6_operator_norm.py
"""
import os, sys, json
sys.dont_write_bytecode = True
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def Fmap(Q, A):
    a = 1 if A % 2 == 1 else 2
    cA = (3 + 2 ** a - 2 ** A) // 3
    inv = pow(2, -A, Q)
    return lambda s: ((3 * s + cA) * inv) % Q


def branch_targets(d):
    Q = 3 ** (d + 1)
    q = Q // 3
    L = 2 * Q
    tg = np.zeros((L, q), dtype=np.int64)
    for o in range(L):
        f = Fmap(Q, o + 1)
        tg[o] = [f(s) for s in range(q)]
    # period check: o -> o + L gives the same map
    f2 = Fmap(Q, L + 1)
    assert all(f2(s) == tg[0][s] for s in range(q))
    return tg, Q, q, L


def symbol(tg, Q, q, L, z):
    wv = z / np.sqrt(2)
    M = np.zeros((Q, q), dtype=complex)
    cols = np.arange(q)
    for o in range(L):
        M[tg[o], cols] += (-1) ** (o + 1) * wv ** o
    return M / (1 - wv ** L) / np.sqrt(2)


def balanced_basis(q):
    lo = q // 3
    C = np.zeros((lo, q))
    for s in range(lo):
        C[s, [s, s + lo, s + 2 * lo]] = 1
    _, sv, Vt = np.linalg.svd(C)
    return Vt[lo:].T            # orthonormal basis of the null space, q x (q - lo)


def is_balanced_out(Q, M):
    lo = Q // 3
    return np.abs(M[:lo] + M[lo:2 * lo] + M[2 * lo:]).max()


def main():
    res = {}
    for d in range(1, 5):
        tg, Q, q, L = branch_targets(d)
        Bq = balanced_basis(q)
        NZ = 4096
        best = (0, 0)
        bestfull = 0
        inv_err = 0
        for th in 2 * np.pi * np.arange(NZ) / NZ:
            M = symbol(tg, Q, q, L, np.exp(1j * th))
            s = np.linalg.svd(M @ Bq, compute_uv=False)[0]
            sf = np.linalg.svd(M, compute_uv=False)[0]
            inv_err = max(inv_err, is_balanced_out(Q, M @ Bq))
            if s > best[0]:
                best = (s, th)
            bestfull = max(bestfull, sf)
        lo, hi = best[1] - 2 * np.pi / NZ, best[1] + 2 * np.pi / NZ
        for _ in range(3):
            ths = np.linspace(lo, hi, 401)
            vals = [np.linalg.svd(symbol(tg, Q, q, L, np.exp(1j * t)) @ Bq, compute_uv=False)[0] for t in ths]
            k = int(np.argmax(vals))
            best = max(best, (vals[k], ths[k]))
            lo, hi = ths[max(k - 2, 0)], ths[min(k + 2, 400)]
        res[d] = dict(norm_balanced=float(best[0]), theta_over_2pi=float(best[1] / (2 * np.pi)), norm_full=float(bestfull),
                      max_imbalance_of_image=float(inv_err))
        print('1. d = %d: sup_z ||T(z)|_balanced|| = %.6f at theta/2pi = %.5f; full space %.6f; image of balanced is balanced '
              '(max defect %.1e)' % (d, best[0], best[1] / (2 * np.pi), bestfull, inv_err), flush=True)
        # 2. truncated operator, no symbol
        if d <= 3:
            for Lb in ((50, 100, 200, 400) if d == 1 else ((50, 100, 200) if d == 2 else (40, 80))):
                nb = Bq.shape[1]
                Mt = np.zeros((Q * Lb, nb * Lb))
                blocks = []
                for o in range(Lb):
                    Pi = np.zeros((Q, q))
                    Pi[tg[o % L], np.arange(q)] = 1
                    blocks.append(((-1) ** (o + 1)) * 2 ** (-o / 2) / np.sqrt(2) * (Pi @ Bq))
                for b in range(Lb):
                    for o in range(b + 1):
                        Mt[b * Q:(b + 1) * Q, (b - o) * nb:(b - o + 1) * nb] = blocks[o]
                s = np.linalg.svd(Mt, compute_uv=False)[0]
                res[d]['truncated_%d' % Lb] = float(s)
                print('2.        L = %3d bands: largest singular value of the truncated operator on balanced inputs = %.6f'
                      % (Lb, s), flush=True)
        # 3. explicit real witness from the symbol's top singular vector
        th = best[1]
        M = symbol(tg, Q, q, L, np.exp(1j * th)) @ Bq
        U, S, Vh = np.linalg.svd(M)
        v = Bq @ Vh[0].conj()                       # complex balanced input direction
        Lw = 600
        win = np.sin(np.pi * (np.arange(Lw) + 0.5) / Lw)
        # the causal convolution maps v z^(-b) to z^(-b) T(z) v, so the input frequency is z^(-1) = exp(-i th)
        y = np.array([np.real(v * np.exp(-1j * th * b)) * win[b] for b in range(Lw)])   # [b][s]
        yo = np.zeros((Lw, Q))
        for o in range(min(Lw, 220)):           # 2^(-o/2) < 1e-33 beyond 220
            c = ((-1) ** (o + 1)) * 2 ** (-o / 2) / np.sqrt(2)
            yo[o:, tg[o % L]] += c * y[:Lw - o, :]          # injective targets: no repeated indices
        ratio = np.sqrt((yo ** 2).sum() / (y ** 2).sum())
        # back to K coordinates: K(s,b) = 2^(b/2) sqrt(4 2^(d-1)) y(s,b); the ratio of energies is ratio^2
        bal = np.abs(y[:, :q // 3] + y[:, q // 3:2 * q // 3] + y[:, 2 * q // 3:]).max()
        res[d]['real_witness_energy_ratio'] = float(ratio ** 2)
        print('3.        explicit real balanced input on %d bands (max lift defect %.1e): ||T K||^2/||K||^2 = %.4f'
              % (Lw, bal, ratio ** 2), flush=True)
    json.dump(res, open(os.path.join(HERE, 'a6_operator_norm.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
