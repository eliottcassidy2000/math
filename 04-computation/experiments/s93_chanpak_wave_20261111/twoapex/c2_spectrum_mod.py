#!/usr/bin/env python3
"""Lane twoapex, deliverable (c): spectrum of the Markov chain of prefix states modulo m (for a closed orbit) and
the modes seen by the residue of H.

usage: python c2_spectrum_mod.py m [neig]
1. orbit of x(empty) modulo m (as c1_orbit_mod.py), transition operator P on measures: (P mu)(s') =
   (1/4) sum over (s, l) with s A_l = s' of mu(s);
2. the common left eigenvector x_0 over Z (x_0 A_l = x_0 for the four letters) with x_0 . readout != 0: translation
   by x_0 commutes with the letters and shifts H by x_0 . readout; is the orbit invariant under it modulo m?
3. the twisted sums c_n(k) = 4^-n sum_{|w| = n} e(k H(w) / m) for n <= NMAX from the chain (exact up to float
   rounding), and their decay: least-squares fit of log |c_n| on a window, and the leading eigenvalues of P whose
   eigenvectors have a nonzero k-th Fourier weight  sum_s v(s) e(k H(s) / m).
"""
import json
import sys
import time
from fractions import Fraction

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from ta_lattice import nullspace_frac


def main():
    m = int(sys.argv[1])
    neig = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    d = json.load(open("minrep17.json"))
    r = d["dimension"]
    Aint = d["A"]
    A = [np.array(M, dtype=np.int64) % m for M in Aint]
    x0 = np.array(d["x0"], dtype=np.int64) % m
    pw = m ** np.arange(r, dtype=np.int64)

    def enc(X):
        return X @ pw

    def dec(codes):
        X = np.empty((len(codes), r), np.int64)
        c = codes.copy()
        for j in range(r):
            X[:, j] = c % m
            c //= m
        return X

    seen = np.array([int(enc(x0[None, :])[0])], dtype=np.int64)
    frontier = seen.copy()
    while True:
        X = dec(frontier)
        new = np.unique(np.concatenate([enc((X @ A[l]) % m) for l in range(4)]))
        new = new[~np.isin(new, seen, assume_unique=True)]
        if len(new) == 0:
            break
        seen = np.union1d(seen, new)
        frontier = new
    n = len(seen)
    X = dec(seen)
    hval = X[:, 0]
    nxt = np.empty((4, n), np.int64)
    for l in range(4):
        codes = enc((X @ A[l]) % m)
        nxt[l] = np.searchsorted(seen, codes)
    print(f"m = {m}: {n} states")

    # ---- 2. common left eigenvectors with eigenvalue 1: x (A_l - I) = 0 for all l
    rows = []
    for l in range(4):
        for j in range(r):
            rows.append([Fraction(Aint[l][i][j]) - (1 if i == j else 0) for i in range(r)])
    ns = nullspace_frac(rows)
    print(f"common left fixed vectors (x A_l = x for all l): dimension {len(ns)}")
    for v in ns:
        den = 1
        for q in v:
            den = den * q.denominator // np.gcd(den, q.denominator)
        vi = [int(q * den) for q in v]
        g = 0
        for q in vi:
            g = np.gcd(g, abs(q))
        vi = [q // g for q in vi]
        print("  x_0 =", vi, "; x_0 . readout =", vi[0])
        shift = np.array(vi, dtype=np.int64) % m
        codes = enc((X + shift) % m)
        pos = np.searchsorted(seen, codes)
        pos[pos >= n] = 0
        inv = bool((seen[pos] == codes).all())
        print(f"  orbit invariant under translation by x_0 modulo {m}: {inv}")
    # ---- recurrent classes: strongly connected components
    G = sp.csr_matrix((np.ones(4 * n), (np.tile(np.arange(n), 4), nxt.ravel())), shape=(n, n))
    ncomp, lab = sp.csgraph.connected_components(G, directed=True, connection="strong")
    sizes = np.bincount(lab)
    # a component is recurrent iff no edge leaves it
    leaves = np.zeros(ncomp, bool)
    for l in range(4):
        leaves |= np.bincount(lab, weights=(lab[nxt[l]] != lab).astype(float), minlength=ncomp) > 0
    rec = np.nonzero(~leaves)[0]
    print(f"strongly connected components: {ncomp}; recurrent classes: {len(rec)} with sizes {sorted(sizes[rec].tolist())}")

    # ---- 3. twisted sums from the chain
    start = int(np.searchsorted(seen, enc(x0[None, :])[0]))
    dist = np.zeros(n)
    dist[start] = 1.0
    NMAX = 120
    om = np.exp(2j * np.pi * hval / m)
    cn = {k: [] for k in range(1, m)}
    for step in range(1, NMAX + 1):
        new = np.zeros(n)
        for l in range(4):
            np.add.at(new, nxt[l], dist)
        dist = new / 4.0
        for k in range(1, m):
            cn[k].append(np.sum(dist * om ** k))
    for k in range(1, (m + 1) // 2 + (1 if m % 2 == 0 else 0)):
        if m % 2 == 0 and k % 2 == 0 and False:
            continue
        a = np.abs(np.array(cn[k]))
        print(f"k = {k}: |c_n| at n = 8, 12, 16, 20, 24, 32, 40, 48: " +
              " ".join(f"{a[t - 1]:.3e}" for t in (8, 12, 16, 20, 24, 32, 40, 48)))
        lo, hi = 20, 60
        sel = [t for t in range(lo, hi + 1) if a[t - 1] > 1e-13]
        if len(sel) > 5:
            sl = np.polyfit(sel, np.log(a[np.array(sel) - 1]), 1)[0]
            print(f"   least-squares rate of |c_n| on n = {sel[0]}..{sel[-1]}: {np.exp(sl):.5f} per letter "
                  f"(spectral radius {4 * np.exp(sl):.4f} of 4)")

    # ---- eigenvalues of P seen by the residue of H
    P = sp.csr_matrix((np.full(4 * n, 0.25), (nxt.ravel(), np.tile(np.arange(n), 4))), shape=(n, n))
    t0 = time.time()
    vals, vecs = spla.eigs(P, k=neig, which="LM", tol=1e-10, maxiter=20000)
    order = np.argsort(-np.abs(vals))
    print(f"leading {neig} eigenvalues of P ({time.time() - t0:.1f} s): modulus, argument/2pi, Fourier weights |w_k|")
    for i in order:
        v = vecs[:, i]
        v = v / np.abs(v).sum()
        w = [abs(np.sum(v * om ** k)) for k in range(1, m)]
        tag = "  <- seen by H" if max(w) > 1e-6 else ""
        print(f"  |lambda| = {abs(vals[i]):.6f}  arg/2pi = {np.angle(vals[i]) / (2 * np.pi):+.4f}  "
              + " ".join(f"{x:.2e}" for x in w[: (m - 1)]) + tag)


if __name__ == "__main__":
    main()
