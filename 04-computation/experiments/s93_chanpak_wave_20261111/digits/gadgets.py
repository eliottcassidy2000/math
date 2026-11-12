#!/usr/bin/env python3
"""Lane "digits", part (2): the gadget matrices of THM-5130 (A), recomputed from the definitions.

Pair (B, O): state s = (r, M, h, E), h = H(B), r = h - W_B(O), M = U_B(O), E = End_B(O).
Gadget (X, P): B+ = B => X => v, v -> O, B - O -> v; new pair (B+, P); s+ = G(X, P) s with

    r+ = r_X r + kappa_X M
    M+ = u (r + h) + (4u - 2 sigma + 2 u3) M - sigma E
    h+ = H_X r + Pi M
    E+ = eps (r - h) + nu M.

This script (1) computes the nine constants of every gadget with X of order <= NMAX (one representative per
isomorphism class of X, every nonempty P), (2) checks the law against a direct count of the state of (B+, P) on
random bases, (3) writes the distinct matrices to gadget_matrices_NMAX.json.

Usage: gadgets.py NMAX [NCHECK]
"""
import json
import random
import sys

import numpy as np

from tour import ham_tables, iso_classes, pair_stats


def pair_state(n, adj, O):
    """(r, M, h, E) of the pair (B, O); O a bitmask."""
    h, start, end, Q, hsub = pair_stats(n, adj)
    W = 0
    M = 0
    E = 0
    for f in range(n):
        if (O >> f) & 1:
            M += int(start[f])
            E += int(end[f])
            for e in range(n):
                if e != f:
                    M += int(Q[e, f])
                    if (O >> e) & 1:
                        W += int(Q[e, f])
    return (int(h) - W, M, int(h), E)


def gadget_constants(n, adj, P):
    """the constants (H_X, Pi, w, w3, sigma, u, eps, nu, u3) of the gadget (X, P)."""
    Et, St = ham_tables(n, adj)
    full = (1 << n) - 1
    hsub = [1] + [int(Et[m].sum()) for m in range(1, full + 1)]
    H = hsub[full]
    Pi = sum(hsub[m] * hsub[full ^ m] for m in range(full + 1))
    # pi2n[R] = ordered pairs of disjoint NONEMPTY paths covering R
    pi2n = [0] * (full + 1)
    for R in range(full + 1):
        s = R
        t = 0
        while s:
            if s != R:
                t += hsub[s] * hsub[R ^ s]
            s = (s - 1) & R
        pi2n[R] = t
    Pl = [i for i in range(n) if (P >> i) & 1]
    sigma = sum(int(St[full, o]) for o in Pl)
    eps = sum(int(Et[full, o]) for o in Pl)
    w = 0
    w3 = 0
    u = sigma
    nu = eps
    u3 = 0
    for m in range(1, full):
        comp = full ^ m
        for o in Pl:
            if (m >> o) & 1:
                # a path on m starting at o (resp. ending at o), the rest one nonempty path / two nonempty paths
                u += int(St[m, o]) * hsub[comp]
                nu += int(Et[m, o]) * hsub[comp]
                u3 += int(St[m, o]) * pi2n[comp]
        for e in Pl:
            if not (m >> e) & 1:
                continue
            ce = int(Et[m, e])
            if ce == 0:
                continue
            # second path on a nonempty subset s2 of comp starting at f; third path (possibly none) on the rest
            s2 = comp
            while s2:
                rest = comp ^ s2
                for f in Pl:
                    if f != e and (s2 >> f) & 1:
                        cf = int(St[s2, f])
                        if cf:
                            if rest == 0:
                                w += ce * cf
                            else:
                                w3 += ce * cf * hsub[rest]
                s2 = (s2 - 1) & comp
    return (H, Pi, w, w3, sigma, u, eps, nu, u3)


def gadget_matrix(c):
    H, Pi, w, w3, sigma, u, eps, nu, u3 = c
    rX = H - w
    kappa = Pi - 2 * w - 2 * w3
    return ((rX, kappa, 0, 0),
            (u, 4 * u - 2 * sigma + 2 * u3, u, -sigma),
            (H, Pi, 0, 0),
            (eps, nu, -eps, 0))


def closing(nB, adjB, O, nX, adjX):
    """B+ = B => X => v with v -> O, B - O -> v; vertices: B, then X, then v."""
    n = nB + nX + 1
    a = np.zeros(n, np.int64)
    xmask = ((1 << nX) - 1) << nB
    v = n - 1
    for i in range(nB):
        a[i] = int(adjB[i]) | xmask
        if not (O >> i) & 1:
            a[i] |= 1 << v
    for i in range(nX):
        a[nB + i] = (int(adjX[i]) << nB) | (1 << v)
    a[v] = O
    return n, a


def random_tournament(n, rng):
    a = np.zeros(n, np.int64)
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < 0.5:
                a[i] |= 1 << j
            else:
                a[j] |= 1 << i
    return a


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    ncheck = int(sys.argv[2]) if len(sys.argv) > 2 else 400
    classes = iso_classes(nmax)
    rng = random.Random(20261111)
    mats = {}
    ngad = 0
    for n in range(1, nmax + 1):
        for key, adj in classes[n].items():
            for P in range(1, 1 << n):
                c = gadget_constants(n, adj, P)
                G = gadget_matrix(c)
                ngad += 1
                if G not in mats:
                    mats[G] = {"n": n, "adj": [int(x) for x in adj], "P": P, "const": list(c), "count": 1}
                else:
                    mats[G]["count"] += 1
        print(f"order {n}: gadgets so far {ngad}, distinct matrices {len(mats)}", flush=True)
    # check of the law on random bases
    keys = list(mats.keys())
    bad = 0
    for t in range(ncheck):
        G = rng.choice(keys)
        info = mats[G]
        nB = rng.randint(1, 6)
        adjB = random_tournament(nB, rng)
        O = rng.randint(1, (1 << nB) - 1)
        s = pair_state(nB, adjB, O)
        nn, a = closing(nB, adjB, O, info["n"], np.array(info["adj"], np.int64))
        s2 = pair_state(nn, a, info["P"] << nB)
        pred = tuple(sum(G[i][j] * s[j] for j in range(4)) for i in range(4))
        if pred != s2:
            bad += 1
            if bad < 5:
                print("LAW FAILS", info, nB, [int(x) for x in adjB], O, s, pred, s2)
    print(f"gadget law checked on {ncheck} random (base, gadget) pairs with bases of order <= 6: {bad} failures")
    out = [{"G": [list(r) for r in G], **info} for G, info in mats.items()]
    with open(f"gadget_matrices_{nmax}.json", "w") as f:
        json.dump(out, f)
    print(f"{ngad} gadgets (iso classes of X of order <= {nmax}, all nonempty P), {len(mats)} distinct matrices -> gadget_matrices_{nmax}.json")


if __name__ == "__main__":
    main()
