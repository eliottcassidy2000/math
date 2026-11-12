#!/usr/bin/env python3
"""Adversarial verifier, lane twoapex: INDEPENDENT census of the two-apex words of length <= LMAX.

usage: python v2_census.py LMAX THREADS STATE_DIR [resume]

Differences from the lane's b1_census.py (so that agreement is informative):
  * own 24-state transfer (v1_core.py, validated there against a literal Hamiltonian-path count);
  * only PRIMITIVE words are enumerated (first letter not 0, last letter not 3): 9/16 of the words;
    P_n(N) = number of primitive words of length exactly n with H = N.  Since H(0 w) = H(w) and H(w 3) = H(w),
        cumprim_n = sum_{m <= n} P_m,   c_n = c_(n-1) + cumprim_n   (c_n(N) = number of ALL words of length n with H = N),
        cum_L = sum_{n <= L} c_n;
    the lane counted all 4^n words directly, and derived its primitive counts as c_L - c_(L-1);
  * different split of the word (suffix of length <= 6), different histogram layout, no modular bookkeeping
    except the class of H modulo 3 and modulo 8 (all words, including those above the bound).
Bound: odd N < 2^26.  Words with H >= 2^26 are counted (and listed when there are few).
The state after each length >= 15 is kept in STATE_DIR (outside the repository) so that the last length can be run
in a separate call ("resume").
"""
import os
import pickle
import sys
import time

import numpy as np
from numba import get_num_threads, njit, prange, set_num_threads

from v1_core import FIN, INIT, MAT

BEXP = 26
HALF = 1 << (BEXP - 1)


def prefixes(n1):
    """alpha vectors of the words of length n1 >= 1 whose first letter is not 0."""
    cur = np.stack([INIT @ MAT[l] for l in (1, 2, 3)], axis=0)
    for _ in range(n1 - 1):
        cur = np.stack([cur @ MAT[l] for l in range(4)], axis=1).reshape(-1, 24)
    return np.ascontiguousarray(cur.astype(np.int32))


def suffixes(n2):
    """beta vectors of the words of length n2 >= 1 whose last letter is not 3."""
    cur = np.stack([MAT[l] @ FIN for l in (0, 1, 2)], axis=0)
    for _ in range(n2 - 1):
        cur = np.stack([cur @ MAT[l].T for l in range(4)], axis=0).reshape(-1, 24)
    return np.ascontiguousarray(cur.astype(np.int32))


@njit(parallel=True, cache=True)
def census(al, be, half, hist, cls, bigbuf):
    n1 = al.shape[0]
    n2 = be.shape[0]
    T = hist.shape[0]
    cap = bigbuf.shape[1]
    nbig = np.zeros(T, np.int64)
    for c in prange(T):
        lo = (n1 * c) // T
        hi = (n1 * (c + 1)) // T
        cnt = 0
        loc = np.zeros(24, np.int64)
        buf = np.empty(n2, np.int64)
        for i in range(lo, hi):
            for j in range(n2):
                h = np.int64(0)
                for s in range(24):
                    h += np.int64(al[i, s]) * np.int64(be[j, s])
                buf[j] = h
                loc[(h % 3) * 8 + (h & 7)] += 1
            for j in range(n2):
                h = buf[j]
                k = h >> 1
                if k < half:
                    hist[c, k] += 1
                else:
                    if cnt < cap:
                        bigbuf[c, cnt] = h
                    cnt += 1
        nbig[c] = cnt
        for t in range(24):
            cls[c, t] = loc[t]
    return nbig


def main():
    LMAX = int(sys.argv[1])
    T = int(sys.argv[2])
    sdir = sys.argv[3]
    resume = len(sys.argv) > 4 and sys.argv[4] == "resume"
    os.makedirs(sdir, exist_ok=True)
    set_num_threads(T)
    t00 = time.time()
    hist = np.zeros((T, HALF), np.uint16)
    cls = np.zeros((T, 24), np.int64)
    CAP = 4096
    bigbuf = np.zeros((T, CAP), np.int64)

    if resume:
        st = pickle.load(open(os.path.join(sdir, "small.pkl"), "rb"))
        cumprim = np.load(os.path.join(sdir, "cumprim.npy"))
        c_n = np.load(os.path.join(sdir, "c_n.npy"))
        cum = np.load(os.path.join(sdir, "cum.npy"))
        minlen = np.load(os.path.join(sdir, "minlen.npy"))
        cumcls, R, cumbig, cbig = st["cumcls"], st["R"], st["cumbig"], st["cbig"]
        bigvals_prev, cumbigvals, NL, nstart = st["bigvals_prev"], st["cumbigvals"], st["NL"], st["n"] + 1
        print(f"resumed after length {nstart - 1}")
    else:
        cumprim = np.zeros(HALF, np.uint32)       # sum_{m <= n} P_m
        c_n = np.zeros(HALF, np.uint32)           # all words of exact length n
        cum = np.zeros(HALF, np.uint32)           # all words of length <= n
        minlen = np.full(HALF, 255, np.uint8)
        cumprim[0] = 1                            # n = 0: the empty word, H = 1
        c_n[0] = 1
        cum[0] = 1
        minlen[0] = 0
        cumcls = np.zeros(24, np.int64)           # sum_{m <= n} classes of the primitive words
        cumcls[(1 % 3) * 8 + 1] = 1
        R = cumcls.copy()                         # classes of all words of exact length n
        cumbig = 0                                # sum_{m <= n} number of primitive words with H >= 2^26
        cbig = 0                                  # number of all words of exact length n with H >= 2^26
        bigvals_prev = {}                         # multiplicities of the values >= 2^26, all words of length n - 1
        cumbigvals = {}                           # the same among the primitive words of length <= n
        NL = {}
        nstart = 1
    print(f"threads {get_num_threads()}, bound 2^{BEXP}", flush=True)
    for n in range(nstart, LMAX + 1):
        t0 = time.time()
        lim = min(HALF, 3 ** n + 1)               # H <= 2 * 3^n, index H >> 1 <= 3^n
        if n == 1:
            nprim = 2
            P = np.zeros(lim, np.uint32)
            pc = np.zeros(24, np.int64)
            for l in (1, 2):
                h = int(INIT @ MAT[l] @ FIN)
                P[h >> 1] += 1
                pc[(h % 3) * 8 + (h & 7)] += 1
            nb = 0
            bigs = np.zeros(0, np.int64)
            tc = 0.0
        else:
            n2 = min(n // 2, 6)
            n1 = n - n2
            al = prefixes(n1)
            be = suffixes(n2)
            nprim = al.shape[0] * be.shape[0]
            assert nprim == 9 * 4 ** (n - 2)
            hist[:, :lim] = 0
            tc = time.time()
            nbig = census(al, be, HALF, hist, cls, bigbuf)
            tc = time.time() - tc
            P = hist[:, :lim].sum(axis=0, dtype=np.uint32)
            pc = cls.sum(axis=0)
            nb = int(nbig.sum())
            bigs = np.concatenate([bigbuf[c, : min(int(nbig[c]), CAP)] for c in range(T)]) if nb <= CAP else None
            del al, be
        assert int(P.sum(dtype=np.int64)) + nb == nprim, "histogram overflow or lost word"
        assert int(pc.sum()) == nprim
        cumprim[:lim] += P
        c_n[:lim] += cumprim[:lim]
        newly = (cum[:lim] == 0) & (c_n[:lim] != 0)
        minlen[:lim][newly] = n
        cum[:lim] += c_n[:lim]
        cumcls += pc
        R = R + cumcls
        assert int(R.sum()) == 4 ** n, "class recursion does not account for all 4^n words"
        cumbig += nb
        cbig += cumbig
        assert int(c_n[:lim].sum(dtype=np.int64)) + cbig == 4 ** n, "primitive decomposition does not give 4^n words"
        exact_big = True
        if nb:
            if bigs is None:
                exact_big = False
            else:
                for h in bigs.tolist():
                    cumbigvals[h] = cumbigvals.get(h, 0) + 1
        bigvals = dict(bigvals_prev)
        for h, mlt in cumbigvals.items():
            bigvals[h] = bigvals.get(h, 0) + mlt
        bigvals_prev = bigvals
        c64 = c_n[:lim].astype(np.int64)
        E_low = int((c64 * c64).sum())
        E = E_low + sum(v * v for v in bigvals.values())
        distinct = int(np.count_nonzero(c64))
        mxc = int(c64.max())
        amx = 2 * int(np.argmax(c64)) + 1
        del c64
        miss = np.flatnonzero(cum[:lim] == 0)
        missN = 2 * miss + 1
        later = missN[missN > 21]
        NL[n] = int(later[0]) if len(later) else (2 * lim + 1 if lim < HALF else None)
        fac3 = [3 * float(R.reshape(3, 8)[c].sum()) / 4 ** n for c in range(3)]
        r8 = R.reshape(3, 8).sum(axis=0)
        fac8 = [4 * float(r8[c]) / 4 ** n for c in (1, 3, 5, 7)]
        print(f"n = {n:2d}: primitive words {nprim}, with H >= 2^26: {nb}; all words with H >= 2^26: {cbig}; "
              f"distinct odd values < 2^26 at exact length: {distinct} (+{len(bigvals)} above, "
              f"{'complete' if exact_big else 'not listed'}); cumulative distinct: {int(np.count_nonzero(cum[:lim]))}; "
              f"misses <= 21: {missN[missN <= 21].tolist()}; N_{n} = {NL[n]}", flush=True)
        print(f"        E_n = {E}{'' if exact_big else ' (lower bound: values >= 2^26 not counted)'}; "
              f"E_n (3/16)^n = {E * (3 / 16) ** n:.5f}; max_N c_n(N) = {mxc} at N = {amx}"
              f"{'' if not bigvals else '; max multiplicity above the bound ' + str(max(bigvals.values()))}; "
              f"max/(4/3)^n = {mxc / (4 / 3) ** n:.2f}; class factors mod 3: "
              + " ".join(f"{t:.6f}" for t in fac3) + "; odd classes mod 8: " + " ".join(f"{t:.6f}" for t in fac8)
              + f"; census {tc:.1f} s, total {time.time() - t0:.1f} s", flush=True)
        if not exact_big:
            cumbigvals = {}
            bigvals_prev = {}
        if n in (14, 16, 17) or n == LMAX:
            nmiss = int((cum == 0).sum())
            allmiss = 2 * np.flatnonzero(cum == 0)[:6] + 1
            print(f"  --- words of length <= {n}: total multiplicity below the bound {int(cum.sum(dtype=np.int64))}; "
                  f"primitive words recorded {int(cumprim.sum(dtype=np.int64))}; odd N < 2^26 not attained: {nmiss}"
                  f" (first: {allmiss.tolist()})")
            print("      window          zeros      min   median       mean | prim: min  median     mean | max minlen")
            for k in range(19, BEXP):
                lo, hi = (1 << k) >> 1, (1 << (k + 1)) >> 1
                seg, ps, ml = cum[lo:hi], cumprim[lo:hi], minlen[lo:hi]
                att = seg != 0
                print(f"      [2^{k}, 2^{k + 1}) {int((seg == 0).sum()):9d} {int(seg.min()):8d} {int(np.median(seg)):8d} "
                      f"{seg.mean():10.1f} | {int(ps.min()):9d} {int(np.median(ps)):7d} {ps.mean():8.1f} | "
                      f"{int(ml[att].max()) if att.any() else -1}")
            prof = [(m, int(((minlen == m) & (cum != 0)).sum())) for m in range(max(0, n - 4), n + 1)]
            print("      number of odd N of minimal length m:", prof, flush=True)
        if n == LMAX:
            print(f"  --- dispersion of the primitive multiplicities (length <= {n}), blocks of 256 consecutive odd N")
            for k in range(20, BEXP):
                lo, hi = (1 << k) >> 1, (1 << (k + 1)) >> 1
                seg = cumprim[lo:hi].astype(np.float64)
                blk = seg.reshape(-1, 256)
                res = blk - blk.mean(axis=1, keepdims=True)
                var = float((res ** 2).sum() / (blk.shape[0] * 255))
                lag = float((res[:, 1:] * res[:, :-1]).sum() / (res ** 2).sum())
                odd = float((cumprim[lo:hi] & 1).mean())
                blk2 = seg.reshape(-1, 4096)
                var2 = float(((blk2 - blk2.mean(axis=1, keepdims=True)) ** 2).sum() / (blk2.shape[0] * 4095))
                # global variance in the window around the window mean, and by class of N modulo 3
                N = 2 * np.arange(lo, hi, dtype=np.int64) + 1
                m3 = [float(seg[N % 3 == c].mean()) / float(seg.mean()) for c in range(3)]
                print(f"      [2^{k}, 2^{k + 1}): mean {seg.mean():9.2f}  var/mean {var / seg.mean():.4f}  "
                      f"var/mean of prim/2 {var / seg.mean() / 2:.4f}  (blocks of 4096: {var2 / seg.mean() / 2:.4f})  "
                      f"lag-1 {lag:+.4f}  share of odd prim {odd:.6f}  min {int(seg.min())}  "
                      f"by-value class factors mod 3 (primitive): " + " ".join(f"{t:.5f}" for t in m3))
        if n >= 15 and n < 17:
            np.save(os.path.join(sdir, "cumprim.npy"), cumprim)
            np.save(os.path.join(sdir, "c_n.npy"), c_n)
            np.save(os.path.join(sdir, "cum.npy"), cum)
            np.save(os.path.join(sdir, "minlen.npy"), minlen)
            pickle.dump({"cumcls": cumcls, "R": R, "cumbig": cumbig, "cbig": cbig, "bigvals_prev": bigvals_prev,
                         "cumbigvals": cumbigvals, "NL": NL, "n": n}, open(os.path.join(sdir, "small.pkl"), "wb"))
    print("N_L:", NL)
    print(f"total time {time.time() - t00:.1f} s")


if __name__ == "__main__":
    main()
