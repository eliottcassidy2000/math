#!/usr/bin/env python3
"""Lane twoapex, deliverable (c): leading eigenvalues of the twisted operator T_k on the base orbit modulo m
(sparse, exact structure as in c4_twisted_base.py), for small orbits (m = 3, 4, 8).

usage: python c6_twisted_eigs.py m k [neig]
Prints the leading eigenvalues of T_k (so that sum_{|w| = n} e(k H(w)/m) grows like lambda^n), their moduli, and
for each whether the corresponding left eigenfunction is seen by the initial state / the sum (weights).
"""
import json
import sys

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


def main():
    m, k = int(sys.argv[1]), int(sys.argv[2])
    neig = int(sys.argv[3]) if len(sys.argv) > 3 else 12
    d = json.load(open("minrep17.json"))
    r = d["dimension"]
    A = [np.array(M, dtype=np.int64) for M in d["A"]]
    x0 = np.array(d["x0"], dtype=np.int64)
    C = [A[l][1:, :] for l in range(4)]
    y0 = (x0[1:] - x0[0]) % m
    pw = m ** np.arange(r - 1, dtype=np.int64)

    def enc(Y):
        return Y @ pw

    def dec(codes):
        Y = np.empty((len(codes), r - 1), np.int64)
        c = codes.copy()
        for j in range(r - 1):
            Y[:, j] = c % m
            c //= m
        return Y

    def step(Y, l):
        Z = Y @ C[l]
        return (Z[:, 1:] - Z[:, :1]) % m, Z[:, 0] % m

    seen = np.array([int(enc(y0[None, :])[0])], dtype=np.int64)
    frontier = seen.copy()
    while True:
        Y = dec(frontier)
        new = np.unique(np.concatenate([enc(step(Y, l)[0]) for l in range(4)]))
        new = new[~np.isin(new, seen, assume_unique=True)]
        if len(new) == 0:
            break
        seen = np.union1d(seen, new)
        frontier = new
    n = len(seen)
    Y = dec(seen)
    rows, cols, vals = [], [], []
    for l in range(4):
        Yn, delta = step(Y, l)
        pos = np.searchsorted(seen, enc(Yn))
        rows.append(pos)
        cols.append(np.arange(n))
        vals.append(np.exp(2j * np.pi * k * delta / m))
    T = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n))
    print(f"m = {m}, k = {k}: base orbit {n} states")
    ev, vec = spla.eigs(T, k=neig, which="LM", tol=1e-13, maxiter=50000)
    order = np.argsort(-np.abs(ev))
    start = int(np.searchsorted(seen, enc(y0[None, :])[0]))
    # left eigenvectors are needed to know whether nu_0 = delta_start excites the mode; use the transpose problem
    evl, vecl = spla.eigs(T.T.tocsr(), k=neig, which="LM", tol=1e-13, maxiter=50000)
    print("eigenvalue (of T_k, so the sum over 4^n words grows like lambda^n): value, modulus, modulus / 4, arg / 2 pi,")
    print("  |sum of the right eigenvector| (does the mode reach the total sum), normalised to l1 norm 1")
    for i in order:
        v = vec[:, i]
        v = v / np.abs(v).sum()
        print(f"  {ev[i].real:+.12f} {ev[i].imag:+.12f} i   |.| = {abs(ev[i]):.12f}   /4 = {abs(ev[i]) / 4:.10f}   "
              f"arg/2pi = {np.angle(ev[i]) / (2 * np.pi):+.6f}   |sum v| = {abs(v.sum()):.3e}")
    print("left eigenfunctions at the initial state (|f(start)| / max |f|):")
    for i in np.argsort(-np.abs(evl)):
        f = vecl[:, i]
        print(f"  lambda = {evl[i].real:+.10f} {evl[i].imag:+.10f} i : {abs(f[start]) / np.abs(f).max():.3e}")


if __name__ == "__main__":
    main()
