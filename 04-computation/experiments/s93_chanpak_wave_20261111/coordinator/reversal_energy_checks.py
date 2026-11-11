#!/usr/bin/env python3
"""S93 coordinator: exact finite checks of the lemmas about the reversal graph
    Gamma_w = {(x, rev_w(x)) : 0 <= x < 2^w}  (rev_w = reversal of the w-bit block).

(L1) The map eps in {-1,0,1}^w -> (D(eps), D(rev eps)), D(eps) = sum eps_i 2^i, is injective.
     Equivalently: no nonzero delta in {-2,...,2}^w has delta(2) = delta(1/2) = 0.
(L2) Additive energy  #{(x1,x2,x3,x4) : p(x1) - p(x2) = p(x3) - p(x4)} = 6^w  (p(x) = (x, rev x)).
(L3) Q_w(k, k') = #{(x, x') : x != x', (x - x', rev x - rev x') in Z (k, k')}, (k, k') primitive:
     max over all primitive directions, compared with 6^(w/2).
(L4) S_w = sum over ordered pairs x != x' of gcd(X, Xi)/max(|X|, |Xi|), all pairs and pairs with X, Xi odd.
(L5) N_d = #{(x, x') : x != x', d | X, d | Xi} for odd d, compared with 6^(w/2), for d >= 2^(w/2 + 1).
Usage: python3 reversal_energy_checks.py [WMAX]   (default 12)
"""
import sys
from collections import Counter

import numpy as np


def rev_table(w):
    r = np.zeros(1 << w, np.int64)
    for i in range(w):
        r |= ((np.arange(1 << w, dtype=np.int64) >> i) & 1) << (w - 1 - i)
    return r


def check_injective(w):
    """all eps in {-1,0,1}^w: images (D, Drev) distinct."""
    n = 3 ** w
    idx = np.arange(n, dtype=np.int64)
    D = np.zeros(n, np.int64)
    R = np.zeros(n, np.int64)
    t = idx.copy()
    for i in range(w):
        e = (t % 3) - 1
        t //= 3
        D += e << i
        R += e << (w - 1 - i)
    key = (D + (1 << w)) * (1 << (w + 2)) + (R + (1 << w))
    return np.unique(key).shape[0] == n


def main():
    wmax = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    print("(L1) injectivity of eps -> (D(eps), D(rev eps)) on {-1,0,1}^w")
    for w in range(1, 15):
        print(f"  w = {w:2d}: {check_injective(w)}", flush=True)
    print("(L2)-(L5)")
    for w in range(2, wmax + 1):
        size = 1 << w
        rv = rev_table(w)
        x = np.arange(size, dtype=np.int64)
        X = (x[:, None] - x[None, :]).ravel()
        Xi = (rv[:, None] - rv[None, :]).ravel()
        keep = X != 0
        X = X[keep]
        Xi = Xi[keep]
        # (L2) energy through difference multiplicities
        key = (X + size) * (4 * size) + (Xi + size)
        _, cnt = np.unique(key, return_counts=True)
        energy = int((cnt.astype(np.int64) ** 2).sum()) + size * size   # zero difference: 2^w pairs, squared
        g = np.gcd(X, Xi)
        mx = np.maximum(np.abs(X), np.abs(Xi))
        S_all = float((g / mx).sum())
        odd = ((X & 1) == 1) & ((Xi & 1) == 1)
        S_odd = float((g[odd] / mx[odd]).sum())
        # (L3) directions
        k = X // g
        kp = Xi // g
        # identify (k, k') with (-k, -k')
        sgn = np.where((k < 0) | ((k == 0) & (kp < 0)), -1, 1)
        k *= sgn
        kp *= sgn
        dkey = k * (4 * size) + (kp + size)
        dk, dcnt = np.unique(dkey, return_counts=True)
        imax = int(np.argmax(dcnt))
        kbest = int(dk[imax] // (4 * size))
        kpbest = int(dk[imax] % (4 * size)) - size
        qmax = int(dcnt[imax])
        # the same restricted to odd directions
        dko, dcnto = np.unique(dkey[odd], return_counts=True)
        qmax_odd = int(dcnto.max()) if dcnto.size else 0
        six = 6 ** (w / 2)
        # (L5)
        worst = 0.0
        worst_d = 0
        lo = 1 << (w // 2 + 1)
        go = g[odd]
        for d in range(lo | 1, size, 2):
            nd = int(np.count_nonzero(go % d == 0))
            if nd / six > worst:
                worst, worst_d = nd / six, d
        print(f"  w = {w:2d}: energy = {energy} (6^w = {6 ** w}, equal: {energy == 6 ** w}); "
              f"S_all/6^(w/2) = {S_all / six:.4f}, S_odd/6^(w/2) = {S_odd / six:.4f}; "
              f"max Q = {qmax} at ({kbest},{kpbest}) = {qmax / six:.3f} * 6^(w/2), odd dirs: {qmax_odd / six:.3f}; "
              f"max over odd d >= 2^(w/2+1) of N_d/6^(w/2) = {worst:.4f} (d = {worst_d})", flush=True)


if __name__ == "__main__":
    main()
