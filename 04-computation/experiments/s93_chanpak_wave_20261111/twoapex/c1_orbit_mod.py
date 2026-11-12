#!/usr/bin/env python3
"""Lane twoapex, deliverable (c): the prefix states modulo m as a finite Markov chain.

usage: python c1_orbit_mod.py m [maxstates]
State of a prefix w: x(w) = (H(w s_j))_j modulo m (minimal representation, minrep17.json), x(w l) = x(w) A_l.
Breadth-first search of the orbit of x(empty) under the four letters (encoded base m in an int64), the number of
states reached by words of each length, and - when the orbit closes - the exact distribution of H = x_1 modulo m
over all 4^n words for n up to NMAX (integer counts), the limiting shares, and the rate of convergence.
"""
import json
import sys
import time

import numpy as np


def main():
    m = int(sys.argv[1])
    maxstates = int(sys.argv[2]) if len(sys.argv) > 2 else 30_000_000
    nmax = int(sys.argv[3]) if len(sys.argv) > 3 else 64
    d = json.load(open("minrep17.json"))
    r = d["dimension"]
    A = [np.array(M, dtype=np.int64) % m for M in d["A"]]
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

    t0 = time.time()
    seen = np.array([int(enc(x0[None, :])[0])], dtype=np.int64)
    frontier = seen.copy()
    depth = 0
    sizes = [1]
    closed = False
    while True:
        X = dec(frontier)
        new = np.concatenate([enc((X @ A[l]) % m) for l in range(4)])
        new = np.unique(new)
        new = new[~np.isin(new, seen, assume_unique=True)]
        depth += 1
        if len(new) == 0:
            closed = True
            break
        seen = np.union1d(seen, new)
        frontier = new
        sizes.append(len(seen))
        if len(seen) > maxstates:
            break
    print(f"m = {m}: orbit of the prefix state: {len(seen)} states after depth {depth} "
          f"({'closed' if closed else 'NOT closed (stopped)'}), {time.time() - t0:.1f} s")
    print("states reached within depth 0, 1, 2, ...:", sizes[:40])
    if not closed:
        return
    # transition table
    n = len(seen)
    X = dec(seen)
    nxt = np.empty((4, n), np.int64)
    for l in range(4):
        codes = enc((X @ A[l]) % m)
        nxt[l] = np.searchsorted(seen, codes)
        assert (seen[nxt[l]] == codes).all()
    hval = X[:, 0]                       # H mod m (readout = first coordinate)
    start = int(np.searchsorted(seen, enc(x0[None, :])[0]))
    # exact counts with Python integers are too slow for large orbits; use float64 with exact integer range check
    # (counts < 2^53 for n <= 26) and, beyond, normalised float distributions.
    dist = np.zeros(n)
    dist[start] = 1.0
    print("n, shares of the classes of H mod m over all 4^n words (times m), and the largest deviation from 1:")
    hist = []
    for step in range(1, nmax + 1):
        new = np.zeros(n)
        for l in range(4):
            np.add.at(new, nxt[l], dist)
        dist = new / 4.0
        sh = np.bincount(hval, weights=dist, minlength=m) * m
        hist.append(sh)
        if step <= 20 or step % 4 == 0:
            cls = [c for c in range(m) if sh[c] > 0 or m % 2 == 1]
            if m % 2 == 0:
                # only odd classes occur; normalise by the number of odd classes
                txt = " ".join(f"{c}:{sh[c] / 2:.6f}" for c in range(1, m, 2))
            else:
                txt = " ".join(f"{c}:{sh[c]:.6f}" for c in range(m))
            print(f"  n = {step:2d}: {txt}")
    hist = np.array(hist)
    lim = hist[-1]
    print("limit shares (n =", nmax, "):", " ".join(f"{v:.9f}" for v in lim))
    # rate: deviation from the limit
    dev = np.abs(hist - lim).max(axis=1)
    for a, b in ((8, 12), (12, 16), (16, 20), (20, 28), (28, 40)):
        if b < nmax and dev[b - 1] > 0 and dev[a - 1] > 0:
            print(f"  max deviation from the limit at n = {a}: {dev[a - 1]:.3e}, at n = {b}: {dev[b - 1]:.3e}, "
                  f"rate per letter {(dev[b - 1] / dev[a - 1]) ** (1 / (b - a)):.5f}")
    # recurrent structure: stationary distribution support
    sup = dist > 1e-15
    print(f"states carrying mass at n = {nmax}: {int(sup.sum())} of {n}")


if __name__ == "__main__":
    main()
