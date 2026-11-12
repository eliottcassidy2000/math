#!/usr/bin/env python3
"""Independent verifier (part D/E) - the gcd sum S_w and the per-pair lattice-line bounds, exactly.

Odd pairs (x, x') of w-bit words are parametrised by their signed digit strings eps = x - x' (digitwise) with
eps_0 != 0 and eps_(w-1) != 0; a string with z zero digits stands for 2^z pairs (immediate: at a zero digit the two
words carry the same bit, 0 or 1).  X = D(eps) = sum eps_i 2^i, Xi = D*(eps) = sum eps_i 2^(w-1-i).

For each w this script computes
  * the exact histogram hist[M] = number of odd pairs of height M = max(|X|,|Xi|)/gcd(X,Xi)   (int64),
    hence S_w = sum hist[M]/M (correctly rounded with math.fsum; also as an exact Fraction for small w);
  * for every (t', s) in a list of combos, the per-pair lattice bounds summed over all odd pairs:
      L1 = sum 4 (1 + g min(2^(t'-1)/|X|, 2^s/|Xi|))                         (the bound in the proof of (C))
      L2 = sum 4 (1 + min(floor((2^(t'-1)-1) g/|X|), floor((2^s-1) g/|Xi|)))   (integer form of the same argument)
      L3 = sum 4 (1 + min(floor((2^(t'-1)-1) g/(2|X|)), floor((2^s-1) g/(2|Xi|))))  (parity-aware: j must be even)
  * a brute-force cross-check over all pairs (x, x') for w <= wbrute.
Usage: python de_gcdsum.py wmax nmin nmax out.json [wbrute]
The combos for width w are all (t', s) with t' >= 1, s >= 0 and nmin <= w + t' + s <= nmax.
"""
import sys, json, time, math
from fractions import Fraction
import numpy as np
from numba import njit, prange


@njit(cache=True)
def gcd(a, b):
    while b:
        a, b = b, a % b
    return a


def half_tables(L, first_plus=False, last_nonzero=False):
    codes = np.arange(3 ** L, dtype=np.int64)
    D = np.zeros(3 ** L, dtype=np.int64)
    Ds = np.zeros(3 ** L, dtype=np.int64)
    z = np.zeros(3 ** L, dtype=np.int64)
    mask = np.ones(3 ** L, dtype=bool)
    for i in range(L):
        dig = (codes // 3 ** i) % 3 - 1
        D += dig * (1 << i)
        Ds += dig * (1 << (L - 1 - i))
        z += (dig == 0)
        if i == 0 and first_plus:
            mask &= (dig == 1)
        if i == L - 1 and last_nonzero:
            mask &= (dig != 0)
    return D[mask].copy(), Ds[mask].copy(), z[mask].copy()


@njit(parallel=True, cache=True)
def kernel(w, wl, wh, Dlo, Dslo, zlo, Dhi, Dshi, zhi, ctm1, cs, nchunks):
    nlo = Dlo.shape[0]
    nhi = Dhi.shape[0]
    ncomb = ctm1.shape[0]
    hist = np.zeros((nchunks, 1 << w), dtype=np.int64)
    L1 = np.zeros((nchunks, ncomb), dtype=np.float64)
    L2 = np.zeros((nchunks, ncomb), dtype=np.int64)
    L3 = np.zeros((nchunks, ncomb), dtype=np.int64)
    tot = np.zeros(nchunks, dtype=np.int64)
    twl = 1 << wl
    twh = 1 << wh
    for c in prange(nchunks):
        acc = np.zeros(ncomb, dtype=np.float64)
        for i in range(c, nlo, nchunks):
            for k in range(ncomb):
                acc[k] = 0.0
            for j in range(nhi):
                X = Dlo[i] + Dhi[j] * twl
                Xi = Dslo[i] * twh + Dshi[j]
                wt = 1 << (1 + zlo[i] + zhi[j])      # factor 2: sign symmetry eps -> -eps
                aX = X if X > 0 else -X
                aXi = Xi if Xi > 0 else -Xi
                g = gcd(aX, aXi)
                mx = aX if aX > aXi else aXi
                hist[c, mx // g] += wt
                tot[c] += wt
                for k in range(ncomb):
                    A = 1 << ctm1[k]
                    B = 1 << cs[k]
                    v1 = (A * g) / aX
                    v2 = (B * g) / aXi
                    acc[k] += wt * (v1 if v1 < v2 else v2)
                    f1 = ((A - 1) * g) // aX
                    f2 = ((B - 1) * g) // aXi
                    L2[c, k] += wt * (1 + (f1 if f1 < f2 else f2))
                    e1 = ((A - 1) * g) // (2 * aX)
                    e2 = ((B - 1) * g) // (2 * aXi)
                    L3[c, k] += wt * (1 + (e1 if e1 < e2 else e2))
            for k in range(ncomb):
                L1[c, k] += acc[k]
    return hist.sum(axis=0), L1.sum(axis=0), L2.sum(axis=0), L3.sum(axis=0), tot.sum()


@njit(cache=True)
def brute(w, ctm1, cs):
    """All ordered pairs (x, x') of w-bit words, directly from the definition of odd pair."""
    N = 1 << w
    rev = np.zeros(N, dtype=np.int64)
    for x in range(N):
        r = 0
        for i in range(w):
            r = (r << 1) | ((x >> i) & 1)
        rev[x] = r
    ncomb = ctm1.shape[0]
    hist = np.zeros(N, dtype=np.int64)
    L2 = np.zeros(ncomb, dtype=np.int64)
    tot = 0
    for x in range(N):
        for y in range(N):
            X = x - y
            Xi = rev[x] - rev[y]
            if (X & 1) == 0 or (Xi & 1) == 0:
                continue
            aX = abs(X)
            aXi = abs(Xi)
            g = gcd(aX, aXi)
            mx = aX if aX > aXi else aXi
            hist[mx // g] += 1
            tot += 1
            for k in range(ncomb):
                A = 1 << ctm1[k]
                B = 1 << cs[k]
                f1 = ((A - 1) * g) // aX
                f2 = ((B - 1) * g) // aXi
                L2[k] += 1 + (f1 if f1 < f2 else f2)
    return hist, L2, tot


def main():
    wmax = int(sys.argv[1])
    nmin = int(sys.argv[2])
    nmax = int(sys.argv[3])
    out_json = sys.argv[4]
    wbrute = int(sys.argv[5]) if len(sys.argv) > 5 else 11
    res = {}
    for w in range(1, wmax + 1):
        t0 = time.time()
        combos = [(tp, s) for tp in range(1, nmax) for s in range(0, nmax) if nmin <= w + tp + s <= nmax]
        ctm1 = np.array([tp - 1 for tp, s in combos], dtype=np.int64)
        cs = np.array([s for tp, s in combos], dtype=np.int64)
        if w == 1:
            # the two odd pairs (0,1), (1,0): X = Xi = +-1, g = 1
            hist = np.zeros(2, dtype=np.int64)
            hist[1] = 2
            tot = 2
            L1 = np.array([2.0 * min(2 ** a, 2 ** b) for a, b in zip(ctm1, cs)])
            L2 = np.array([2 * (1 + min(2 ** a - 1, 2 ** b - 1)) for a, b in zip(ctm1, cs)], dtype=np.int64)
            L3 = np.array([2 * (1 + min((2 ** a - 1) // 2, (2 ** b - 1) // 2)) for a, b in zip(ctm1, cs)], dtype=np.int64)
        else:
            wl = w // 2
            wh = w - wl
            Dlo, Dslo, zlo = half_tables(wl, first_plus=True)
            Dhi, Dshi, zhi = half_tables(wh, last_nonzero=True)
            if len(combos) == 0:
                ctm1 = np.zeros(0, dtype=np.int64)
                cs = np.zeros(0, dtype=np.int64)
            hist, L1, L2, L3, tot = kernel(w, wl, wh, Dlo, Dslo, zlo, Dhi, Dshi, zhi, ctm1, cs, 8)
        tot = int(tot)
        nodd = 2 if w == 1 else 4 ** (w - 1)
        assert tot == nodd, (w, tot, nodd)
        assert int(hist.sum()) == nodd
        Ms = np.nonzero(hist)[0]
        S = math.fsum(int(hist[M]) / int(M) for M in Ms)
        Sexact = None
        if w <= 12:
            Sexact = sum(Fraction(int(hist[M]), int(M)) for M in Ms)
            assert abs(float(Sexact) - S) < 1e-9 * S
        if w <= wbrute:
            bh, bL2, btot = brute(w, ctm1, cs)
            assert btot == nodd
            assert np.array_equal(bh, hist[:1 << w]) if w > 1 else (bh[1] == 2)
            assert np.array_equal(bL2, L2), (w, bL2[:5], L2[:5])
            bchk = 'brute-force cross-check OK'
        else:
            bchk = ''
        # dyadic histogram of heights and partial sums
        dy = []
        for i in range(w + 1):
            lo, hi = 1 << i, min(1 << (i + 1), hist.shape[0])
            if lo >= hist.shape[0]:
                break
            dy.append(int(hist[lo:hi].sum()))
        entry = {'w': w, 'odd_pairs': nodd, 'S': S, 'S_over_6h': S / 6 ** (w / 2), 'height1': int(hist[1]),
                 'dyadic_hist': dy,
                 'combos': {'%d,%d' % (tp, s): {'L1': 4.0 * nodd + 4.0 * float(L1[k]), 'L2': 4 * int(L2[k]),
                                                'L3': 4 * int(L3[k])} for k, (tp, s) in enumerate(combos)}}
        if Sexact is not None:
            entry['S_exact'] = [str(Sexact.numerator), str(Sexact.denominator)]
        res[str(w)] = entry
        print("w=%2d odd pairs=%d  S_w=%.6f  S_w/6^(w/2)=%.4f  S_w/2^(3w/2)=%.5f  height-1 pairs=%d (=%.4f*6^(w/2))  "
              "max height hit=%d  combos=%d  %s  %.1fs" % (w, nodd, S, S / 6 ** (w / 2), S / 2 ** (1.5 * w), int(hist[1]),
                                                           int(hist[1]) / 6 ** (w / 2), int(Ms[-1]), len(combos), bchk,
                                                           time.time() - t0))
        sys.stdout.flush()
        with open(out_json, 'w') as f:
            json.dump(res, f)
    print("written", out_json)


if __name__ == '__main__':
    main()
