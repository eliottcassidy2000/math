#!/usr/bin/env python3
"""exploration: the m-dimensional invariant subspaces behind the factorisation of the characteristic polynomial."""
import sys

import numpy as np
from flint import fmpz_mat

from explore_charpoly import big_matrix
from explore_invariant import twisted_matrix

np.set_printoptions(linewidth=250, precision=4, suppress=True)


def rref(A, tol=1e-9):
    A = A.copy()
    rows, cols = A.shape
    piv = []
    r = 0
    for c in range(cols):
        if r >= rows:
            break
        p = np.argmax(np.abs(A[r:, c])) + r
        if abs(A[p, c]) < tol:
            continue
        A[[r, p]] = A[[p, r]]
        A[r] /= A[r, c]
        for i in range(rows):
            if i != r:
                A[i] -= A[i, c] * A[r]
        piv.append(c)
        r += 1
    return A[:r], piv


def main():
    m = int(sys.argv[1])
    k = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    B = big_matrix(m, 1)
    cp = fmpz_mat(B).charpoly()
    c, facs = cp.factor()
    T = twisted_matrix(m, k)
    ev, R = np.linalg.eig(T)          # right eigenvectors: T R = R diag(ev)
    evl, L = np.linalg.eig(T.T)       # left eigenvectors
    zeta = np.exp(2j * np.pi / m)
    for f, mu in facs:
        co = [int(x) for x in f.coeffs()]
        val = np.array([abs(np.polyval(co[::-1], e)) for e in ev])
        sel = val < 1e-6 * max(1.0, val.max())
        print(f"factor of degree {f.degree()}: {sel.sum()} eigenvalues of T_{k}; moduli {np.round(np.sort(np.abs(ev[sel]))[::-1], 4)}")
        if 0 < sel.sum() < 2 * m:
            V = R[:, sel].T               # rows span the right-invariant subspace
            Rr, piv = rref(V)
            print("  right invariant subspace (column vectors u(eta,beta), index 2 eta + beta), RREF; pivots", piv)
            # express the entries as powers of zeta when possible
            for row in Rr:
                desc = []
                for v in row:
                    if abs(v) < 1e-8:
                        desc.append("   0   ")
                    else:
                        desc.append(f"{abs(v):.3f}@{(np.angle(v) / (2 * np.pi) * m) % m:5.2f}")
                print("   ", " ".join(desc))
            val2 = np.array([abs(np.polyval(co[::-1], e)) for e in evl])
            sel2 = val2 < 1e-6 * max(1.0, val2.max())
            V2 = L[:, sel2].T
            Rl, piv2 = rref(V2)
            print("  left invariant subspace (row vectors), RREF; pivots", piv2)
            for row in Rl:
                desc = []
                for v in row:
                    if abs(v) < 1e-8:
                        desc.append("   0   ")
                    else:
                        desc.append(f"{abs(v):.3f}@{(np.angle(v) / (2 * np.pi) * m) % m:5.2f}")
                print("   ", " ".join(desc))


if __name__ == "__main__":
    main()
