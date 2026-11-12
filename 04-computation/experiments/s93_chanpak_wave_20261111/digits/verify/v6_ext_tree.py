#!/usr/bin/env python3
"""Verifier (lane digits), check V6: the unrestricted extension tree H -> 2H + 1, 2H - 1 (one-vertex steps).

K_1 = { T : some one-vertex extension has H = 2H(T) + 1 and some has H = 2H(T) - 1 },
K_2 = { T : some extension with 2H + 1 lies in K_1 and some extension with 2H - 1 lies in K_1 }.

Own counting code (own subset DP and own insertion formula).  Class representatives of orders <= 8 are taken from the
lane's tour.iso_classes (their numbers 1, 1, 2, 4, 12, 56, 456, 6880 are the known numbers of tournaments).
Order 9 is covered WITHOUT canonical forms: every tournament of order 9 is a one-vertex extension of a tournament of
order 8, so all 6880 * 256 labelled extensions are examined (each class at least once).

Usage: v6_ext_tree.py [9]     (with the argument 9 the order-9 slice is run as well)
"""
import sys
import time

import numpy as np
from numba import njit, prange

sys.path.insert(0, sys.path[0] + "/..")
from tour import iso_classes  # noqa: E402


@njit(cache=True)
def ext_vals(n, adj, out):
    """out[O] = H(T + v), v -> O, V - O -> v.  Returns H(T).  Insertion formula with own tables."""
    full = (1 << n) - 1
    E = np.zeros((full + 1, n), np.int64)      # paths of T[mask] ending at v
    S = np.zeros((full + 1, n), np.int64)      # paths of T[mask] starting at v
    for v in range(n):
        E[1 << v, v] = 1
        S[1 << v, v] = 1
    for m in range(1, full + 1):
        for v in range(n):
            c = E[m, v]
            if c != 0:
                for u in range(n):
                    if (adj[v] >> u) & 1 and not (m >> u) & 1:
                        E[m | (1 << u), u] += c
            c = S[m, v]
            if c != 0:
                for u in range(n):
                    if (adj[u] >> v) & 1 and not (m >> u) & 1:
                        S[m | (1 << u), u] += c
    h = 0
    for v in range(n):
        h += E[full, v]
    Q = np.zeros((n, n), np.int64)             # Q[e, f]: P1 ends at e, P2 starts at f, both nonempty, covering
    for m in range(1, full):
        comp = full ^ m
        for e in range(n):
            ce = E[m, e]
            if ce == 0:
                continue
            for f in range(n):
                cf = S[comp, f]
                if cf != 0:
                    Q[e, f] += ce * cf
    for O in range(full + 1):
        t = 0
        for f in range(n):
            if (O >> f) & 1:
                t += S[full, f]                # v first
                for e in range(n):
                    if not (O >> e) & 1:
                        t += Q[e, f]           # ... e v f ...
            else:
                t += E[full, f]                # v last
        out[O] = t
    return h


@njit(cache=True)
def extend(n, adj, O, a):
    for i in range(n):
        a[i] = adj[i]
        if not (O >> i) & 1:
            a[i] |= 1 << n
    a[n] = O


@njit(cache=True)
def in_K1(n, adj, buf):
    h = ext_vals(n, adj, buf)
    p = False
    q = False
    for O in range(1 << n):
        if buf[O] == 2 * h + 1:
            p = True
        elif buf[O] == 2 * h - 1:
            q = True
    return p and q


@njit(cache=True)
def depth2(n, adj):
    """0: not in K_1; 1: in K_1 but not in K_2; 2: in K_2."""
    buf = np.zeros(1 << n, np.int64)
    buf2 = np.zeros(1 << (n + 1), np.int64)
    a = np.zeros(n + 1, np.int64)
    h = ext_vals(n, adj, buf)
    hasp = False
    hasm = False
    for O in range(1 << n):
        if buf[O] == 2 * h + 1:
            hasp = True
        elif buf[O] == 2 * h - 1:
            hasm = True
    if not (hasp and hasm):
        return 0
    okp = False
    for O in range(1 << n):
        if buf[O] == 2 * h + 1:
            extend(n, adj, O, a)
            if in_K1(n + 1, a, buf2):
                okp = True
                break
    if not okp:
        return 1
    for O in range(1 << n):
        if buf[O] == 2 * h - 1:
            extend(n, adj, O, a)
            if in_K1(n + 1, a, buf2):
                return 2
    return 1


@njit(parallel=True, cache=True)
def order9_slice(reps):
    """reps: (count, 8) adjacency rows of the order-8 classes.  Every one-vertex extension (order 9)."""
    nc = reps.shape[0]
    k1 = np.zeros(nc, np.int64)
    k2 = np.zeros(nc, np.int64)
    for c in prange(nc):
        a = np.zeros(9, np.int64)
        for O in range(256):
            extend(8, reps[c], O, a)
            d = depth2(9, a)
            if d >= 1:
                k1[c] += 1
            if d >= 2:
                k2[c] += 1
    return k1.sum(), k2.sum()


def main():
    do9 = len(sys.argv) > 1 and sys.argv[1] == "9"
    t0 = time.time()
    classes = iso_classes(8)
    print(f"class representatives (lane enumeration): {[len(classes[n]) for n in range(1, 9)]}   [{time.time() - t0:.1f} s]", flush=True)
    for n in range(1, 9):
        c1 = 0
        c2 = 0
        for key, adj in classes[n].items():
            d = depth2(n, np.array(adj, np.int64))
            c1 += d >= 1
            c2 += d >= 2
        print(f"   order {n}: classes {len(classes[n])}, in K_1: {c1}, in K_2: {c2}   [{time.time() - t0:.1f} s]", flush=True)
    if do9:
        reps = np.array([[int(x) for x in adj] for adj in classes[8].values()], dtype=np.int64)
        k1, k2 = order9_slice(reps)
        print(f"   order 9: {reps.shape[0] * 256} labelled one-vertex extensions of the order-8 classes (every class of order 9 "
              f"occurs): in K_1: {k1} (labelled instances), in K_2: {k2}   [{time.time() - t0:.1f} s]")
    print(f"time {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
