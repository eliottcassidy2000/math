#!/usr/bin/env python3
"""Is the ideal ROOT transport a contraction? Operator norm via its band symbol (S41).

clock_duality_oct15's ideal recursion (R8) maps K_(d-1) (residues mod 3^d, bands b >= 0) to K_d (residues mod 3^(d+1)):
    K_d(r, b) = sum_(j >= 0, a) (-1)^a K_(d-1)(s_(j,a)(r), b - o_(j,a)),   o_(j,a) = a - 1 + 2j,
where a = 2 if B^j(r) = 1 (mod 3), a = 1 if B^j(r) = 2 (mod 3), s = (2^a B^j(r) - 1)/3 mod 3^d. The energy
E_d = sum_b 2^-b sum_r K_d^2 / (4 2^d) is a weighted l2 norm. With y(b) = 2^(-b/2) x(b) the bands become an unweighted
one-sided sequence space and the transport is the analytic Toeplitz operator with matrix symbol

    T(z) = 2^(-1/2) sum_(j, a) (-1)^a (z / sqrt 2)^(o_(j,a)) Pi_(j,a),        |z| = 1,

Pi_(j,a) the partial permutation r <- s_(j,a)(r). B has order 3^(d+1) modulo 3^(d+1), so the j-sum is a finite matrix
with geometric factors. The operator norm on vector-valued H^2 equals sup_(|z|=1) ||T(z)||, and the last-digit balanced
vectors form an invariant subspace (R8.2). The diagonal of T(z)*T(z) is the identity (the energy isometry of THM-4667),
so ||T|| <= 1 on the balanced subspace would make the ideal energy non-increasing at EVERY depth.
Usage: python root_transport_symbol_norm_20261017.py DMAX [NZ]
Session opus-2026-10-17-S41.
"""
import sys, os, json
import numpy as np


def symbol_parts(d):
    """Return list of (o, sign, rows, cols) with the j-sum folded mod period q = 3^(d+1)."""
    q = 3 ** (d + 1); qo = 3 ** d
    parts = []
    for r in range(q):
        rho = r
        for j in range(q):
            if rho % 3 != 0:
                a = 2 if rho % 3 == 1 else 1
                s = ((2 ** a * rho - 1) // 3) % qo
                parts.append((a - 1 + 2 * j, (-1) ** a, r, s))
            rho = (4 * rho + 1) % q
    parts = tuple(np.array(c) for c in zip(*parts))
    return parts, q, qo


def T_of_z(parts, q, qo, z):
    M = np.zeros((q, qo), dtype=complex)
    w = z / np.sqrt(2)
    fold = 1.0 / (1.0 - w ** (2 * q))           # j -> j + q multiplies by w^(2q)
    o, sg, r, s = parts
    np.add.at(M, (r, s), sg * w ** o)
    return M * fold / np.sqrt(2)


def balanced_basis(d):
    """Orthonormal basis of vectors on Z/3^d whose three last-digit lifts sum to zero (d >= 1)."""
    qo = 3 ** d; lo = qo // 3
    V = []
    for s in range(lo):
        for v in ([1, -1, 0], [1, 1, -2]):
            x = np.zeros(qo)
            for i in range(3):
                x[s + i * lo] = v[i]
            V.append(x / np.linalg.norm(x))
    return np.array(V).T


def main():
    DMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    NZ = int(sys.argv[2]) if len(sys.argv) > 2 else 720
    rows = []
    for d in range(1, DMAX + 1):
        parts, q, qo = symbol_parts(d)
        Q = balanced_basis(d)
        best = (0, None); bestfull = (0, None)
        zs = np.exp(2j * np.pi * np.arange(NZ) / NZ)
        for z in zs:
            M = T_of_z(parts, q, qo, z)
            sv = np.linalg.svd(M @ Q, compute_uv=False)[0]
            svf = np.linalg.svd(M, compute_uv=False)[0]
            if sv > best[0]:
                best = (sv, z)
            if svf > bestfull[0]:
                bestfull = (svf, z)
        # refine near the best z
        th0 = np.angle(best[1])
        for th in th0 + np.linspace(-2 * np.pi / NZ, 2 * np.pi / NZ, 201):
            z = np.exp(1j * th)
            sv = np.linalg.svd(T_of_z(parts, q, qo, z) @ Q, compute_uv=False)[0]
            if sv > best[0]:
                best = (sv, z)
        M1 = T_of_z(parts, q, qo, 1.0)
        row = dict(d=d, norm_balanced=float(best[0]), z_angle_over_2pi=float(np.angle(best[1]) / (2 * np.pi)),
                   norm_full=float(bestfull[0]), norm_balanced_at_1=float(np.linalg.svd(M1 @ Q, compute_uv=False)[0]))
        rows.append(row)
        print('d = %d: sup_z ||T(z)|| on balanced = %.6f (at angle %.5f * 2pi); at z = 1: %.6f; full space sup %.6f'
              % (d, row['norm_balanced'], row['z_angle_over_2pi'], row['norm_balanced_at_1'], row['norm_full']), flush=True)
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(rows, open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'root_transport_symbol_norm_20261017.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
