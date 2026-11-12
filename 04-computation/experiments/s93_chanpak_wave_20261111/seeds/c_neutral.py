#!/usr/bin/env python3
"""S93 lane "seeds", deliverable (c): neutral gadgets.

Gadget (X, P) of THM-5130:  r+ = r_X r + kappa_X M,  r_X = H(X) - W_X(P),  kappa_X(P) = Pi_2(X) - 2 W_X(P) - 2 W3_X(P).
Both W and W3 are sums over the pairs of P, so  kappa_X(P) = Pi_2(X) - sum over pairs {e, f} of P of k_ef,
k_ef = 2 (QQ(e, f) + QQ(f, e)),  QQ(e, f) = number of triples (D1, D2, D3) of disjoint paths covering X, D1, D2 nonempty,
D3 possibly empty, D1 ending at e, D2 starting at f.
For every class X of order <= NMAX and EVERY nonempty P: the gadgets with kappa = 0, by (order, |P|, r_X, H(X));
violations of HYP-9662 (D) (|P| != 2 or W != 2 H(X)); neutral gadgets with |r_X| = 1; sign statistics of
kappa_X({e, f}) over all pairs.
Usage: python3 c_neutral.py NMAX   (NMAX <= 9)
"""
import sys
import time
from collections import Counter

import numpy as np
from numba import njit, prange

from seedlib import dp_tables

NCH = 48
MAXREC = 20000


@njit(cache=True)
def decode(n, code, out):
    for v in range(n):
        out[v] = 0
    b = 0
    for j in range(1, n):
        for i in range(j):
            if (code >> b) & 1:
                out[i] |= 1 << j
            else:
                out[j] |= 1 << i
            b += 1


@njit(cache=True)
def gadget_tables(n, fe, gs, hS, UR, K, Wm):
    """K[e, f] = k_ef (symmetric), Wm[e, f] = W_ef (symmetric); returns Pi_2."""
    full = (1 << n) - 1
    for R in range(full + 1):
        for f in range(n):
            UR[R, f] = 0
    for R in range(1, full + 1):
        S2 = R
        while S2:
            h = hS[R ^ S2]
            for f in range(n):
                if (S2 >> f) & 1:
                    g = gs[S2, f]
                    if g:
                        UR[R, f] += g * h
            S2 = (S2 - 1) & R
    for e in range(n):
        for f in range(n):
            K[e, f] = 0
            Wm[e, f] = 0
    for S1 in range(1, full):
        C = full ^ S1
        for e in range(n):
            if (S1 >> e) & 1:
                a = fe[S1, e]
                if a:
                    for f in range(n):
                        if (C >> f) & 1:
                            K[e, f] += a * UR[C, f]
                            Wm[e, f] += a * gs[C, f]
    for e in range(n):
        for f in range(e + 1, n):
            s = 2 * (K[e, f] + K[f, e])
            K[e, f] = s
            K[f, e] = s
            w = Wm[e, f] + Wm[f, e]
            Wm[e, f] = w
            Wm[f, e] = w
    pi2 = 0
    for S in range(full + 1):
        pi2 += hS[S] * hS[full ^ S]
    return pi2


@njit(parallel=True, cache=True)
def scan(n, reps, rec, nrec, signs, extremes):
    nrep = reps.shape[0]
    for ch in prange(NCH):
        fe = np.zeros((1 << n, n), np.int64)
        gs = np.zeros((1 << n, n), np.int64)
        hS = np.zeros(1 << n, np.int64)
        UR = np.zeros((1 << n, n), np.int64)
        K = np.zeros((n, n), np.int64)
        Wm = np.zeros((n, n), np.int64)
        Ks = np.zeros(1 << n, np.int64)
        Ws = np.zeros(1 << n, np.int64)
        out = np.zeros(n, np.int64)
        full = (1 << n) - 1
        for idx in range(ch, nrep, NCH):
            decode(n, reps[idx], out)
            dp_tables(n, out, fe, gs, hS)
            pi2 = gadget_tables(n, fe, gs, hS, UR, K, Wm)
            H = hS[full]
            for e in range(n):
                for f in range(e + 1, n):
                    kap = pi2 - K[e, f]
                    if kap > 0:
                        signs[ch, 2] += 1
                    elif kap == 0:
                        signs[ch, 1] += 1
                    else:
                        signs[ch, 0] += 1
                    # extremes of kappa / Pi_2 (as fractions kept by cross-multiplication with floats)
                    x = kap / pi2
                    if x > extremes[ch, 1]:
                        extremes[ch, 1] = x
                    if x < extremes[ch, 0]:
                        extremes[ch, 0] = x
            Ks[0] = 0
            Ws[0] = 0
            for P in range(1, full + 1):
                low = P & (-P)
                l = 0
                while (low >> l) != 1:
                    l += 1
                R = P ^ low
                k = Ks[R]
                w = Ws[R]
                for e in range(n):
                    if (R >> e) & 1:
                        k += K[l, e]
                        w += Wm[l, e]
                Ks[P] = k
                Ws[P] = w
                if k == pi2:
                    c = nrec[ch]
                    if c < MAXREC:
                        rec[ch, c, 0] = reps[idx]
                        rec[ch, c, 1] = P
                        rec[ch, c, 2] = H - w
                        rec[ch, c, 3] = H
                    nrec[ch] += 1


def popcount(x):
    return bin(int(x)).count("1")


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    lines = []

    def P(s=""):
        print(s, flush=True)
        lines.append(s)

    P("S93 lane seeds: neutral gadgets, kappa_X(P) = 0, all classes X of each order, every nonempty P")
    for n in range(2, nmax + 1):
        t0 = time.time()
        reps = np.atleast_1d(np.loadtxt(f"cls_{n}.txt", dtype=np.int64))
        rec = np.zeros((NCH, MAXREC, 4), np.int64)
        nrec = np.zeros(NCH, np.int64)
        signs = np.zeros((NCH, 3), np.int64)
        extremes = np.zeros((NCH, 2), np.float64)
        extremes[:, 0] = 1e18
        extremes[:, 1] = -1e18
        scan(n, reps, rec, nrec, signs, extremes)
        if (nrec > MAXREC).any():
            P("  WARNING: record buffer overflow")
        stats = Counter()
        wit = {}
        viol = []
        unit = []
        for ch in range(NCH):
            for c in range(min(int(nrec[ch]), MAXREC)):
                code, Pm, rX, H = (int(x) for x in rec[ch, c])
                sz = popcount(Pm)
                key = (sz, rX, H)
                stats[key] += 1
                wit.setdefault(key, (code, Pm))
                if sz != 2 or rX != -H:
                    viol.append((code, Pm, rX, H))
                if abs(rX) == 1:
                    unit.append((code, Pm, rX, H))
        sg = signs.sum(axis=0)
        P(f"order {n}: {reps.shape[0]} classes; kappa = 0 for {int(nrec.sum())} gadgets (labelled P of the representatives)  [{time.time() - t0:.1f} s]")
        P(f"   pairs {{e, f}}: kappa < 0: {int(sg[0])}, = 0: {int(sg[1])}, > 0: {int(sg[2])};  range of kappa / Pi_2 over pairs: [{extremes[:, 0].min():.4f}, {extremes[:, 1].max():.4f}]")
        P(f"   by (|P|, r_X, H(X)): " + ", ".join(f"{k}: {stats[k]}" for k in sorted(stats, key=lambda t: (t[0], t[2]))))
        P(f"   violations of HYP-9662 (D) (|P| != 2 or r_X != -H(X)): {len(viol)}" + (f"; first {viol[:6]}" if viol else ""))
        hs1 = sum(v for k, v in stats.items() if k[2] == 1)
        P(f"   neutral with |r_X| = 1: {len(unit)}; all of them have H(X) = 1: {all(u[3] == 1 for u in unit)}; number with H(X) = 1: {hs1} (2n - 3 = {2 * n - 3})")
    with open(f"c_neutral_{nmax}.out", "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
