#!/usr/bin/env python3
"""S93 lane "seeds", deliverable (d): unit pairs (B, O) with B a one-apex or a two-apex tournament, EVERY vertex set O.

One-apex B(y): spine z_1 -> ... -> z_n (vertices 0..n-1), apex a (vertex n), a -> z_i iff y_i = 1; the word is the
integer w with y_1 its most significant bit (so U_B(a) = 1 + w).
Two-apex B(a, b): spine as above, apexes u (vertex n) -> v (vertex n+1), u -> z_i iff a_i = 1, v -> z_i iff b_i = 1.
For every word and every O: r = H(B) - W_B(O), M = U_B(O); records with r = +1 or -1.
Output: d_<family>_<n>.npz with the records (word(s), O mask, sign, M, H(B)), and a summary.
Usage: python3 d_apex_families.py one NMAX | two NMAX
"""
import sys
import time

import numpy as np
from numba import njit, prange

from fastdp import dp_range, wu_tables, subset_sums_w, popcount

NCH = 48


@njit(cache=True)
def build_one(n, w, out):
    for i in range(n):
        m = 0
        for j in range(i + 1, n):
            m |= 1 << j
        out[i] = m
    out[n] = 0
    for i in range(n):
        if (w >> (n - 1 - i)) & 1:
            out[n] |= 1 << i
        else:
            out[i] |= 1 << n


@njit(cache=True)
def build_two(n, wa, wb, out):
    for i in range(n):
        m = 0
        for j in range(i + 1, n):
            m |= 1 << j
        out[i] = m
    out[n] = 1 << (n + 1)
    out[n + 1] = 0
    for i in range(n):
        if (wa >> (n - 1 - i)) & 1:
            out[n] |= 1 << i
        else:
            out[i] |= 1 << n
        if (wb >> (n - 1 - i)) & 1:
            out[n + 1] |= 1 << i
        else:
            out[i] |= 1 << (n + 1)


@njit(parallel=True, cache=True)
def scan(family, n, rec, nrec, maxrec):
    N = n + 1 if family == 1 else n + 2
    nwords = (1 << n) if family == 1 else (1 << (2 * n))
    for ch in prange(NCH):
        fe = np.zeros((1 << N, N), np.int64)
        gs = np.zeros((1 << N, N), np.int64)
        hS = np.zeros(1 << N, np.int64)
        Wm = np.zeros((N, N), np.int64)
        U = np.zeros(N, np.int64)
        Wsub = np.zeros(1 << N, np.int64)
        Msub = np.zeros(1 << N, np.int64)
        out = np.zeros(N, np.int64)
        full = (1 << N) - 1
        for w in range(ch, nwords, NCH):
            if family == 1:
                build_one(n, w, out)
            else:
                build_two(n, w >> n, w & ((1 << n) - 1), out)
            dp_range(N, out, fe, gs, hS, 1)
            wu_tables(N, fe, gs, hS, Wm, U)
            subset_sums_w(N, Wm, U, Wsub, Msub)
            H = hS[full]
            for O in range(1, full + 1):
                x = Wsub[O]
                if x == H - 1 or x == H + 1:
                    c = nrec[ch]
                    if c < maxrec:
                        rec[ch, c, 0] = w
                        rec[ch, c, 1] = O
                        rec[ch, c, 2] = 1 if x == H - 1 else -1
                        rec[ch, c, 3] = Msub[O]
                        rec[ch, c, 4] = H
                    nrec[ch] += 1


def main():
    fam = sys.argv[1]
    nmax = int(sys.argv[2])
    family = 1 if fam == "one" else 2
    lines = []

    def P(s=""):
        print(s, flush=True)
        lines.append(s)

    P(f"S93 lane seeds: unit pairs over {fam}-apex bases, every vertex set O")
    allrec = {}
    nmin = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    for n in range(nmin, nmax + 1):
        t0 = time.time()
        maxrec = 60000
        rec = np.zeros((NCH, maxrec, 5), np.int64)
        nrec = np.zeros(NCH, np.int64)
        scan(family, n, rec, nrec, maxrec)
        if (nrec > maxrec).any():
            P("  WARNING: record buffer overflow")
        R = np.concatenate([rec[ch, :min(int(nrec[ch]), maxrec)] for ch in range(NCH)])
        allrec[n] = R
        for sg, name in ((1, "seeds"), (-1, "anti-seeds")):
            sub = R[R[:, 2] == sg]
            sizes = np.array([popcount(int(o)) for o in sub[:, 1]], np.int64) if len(sub) else np.zeros(0, np.int64)
            napex = n if family == 1 else n
            apexmask = (1 << n) if family == 1 else (3 << n)
            with_apex = int(((sub[:, 1] & apexmask) != 0).sum()) if len(sub) else 0
            bysize = sorted((int(s), int((sizes == s).sum())) for s in set(sizes.tolist()))
            P(f"n = {n:2d}: {name}: {len(sub)} pairs; by |O|: {bysize}; with an apex in O: {with_apex}; distinct slopes {len(set(sub[:, 3].tolist()))}; "
              f"largest slope {int(sub[:, 3].max()) if len(sub) else 0}  [{time.time() - t0:.1f} s]")
    tag = f"{nmax}" if nmin == 1 else f"{nmin}_{nmax}"
    np.savez_compressed(f"d_{fam}_{tag}.npz", **{f"rec_{n}": allrec[n] for n in allrec})
    # coverage of odd parts and slopes by the whole family (words of length <= nmax), closure under doubling
    for sg, name in ((1, "seeds"), (-1, "anti-seeds")):
        slopes = set()
        for n in allrec:
            R = allrec[n]
            slopes |= set(R[R[:, 2] == sg][:, 3].tolist())
        best = {}
        for M in slopes:
            c, v = M, 0
            while c % 2 == 0:
                c //= 2
                v += 1
            if c not in best or v < best[c]:
                best[c] = v
        P(f"{name}: distinct slopes over lengths <= {nmax}: {len(slopes)}; distinct odd parts: {len(best)}")
        for k in range(4, 17, 2):
            X = 1 << k
            cov = sum(1 for c in best if c <= X)
            odd = sum(1 for c, v in best.items() if c <= X and v == 0)
            dens = sum(1 for c, v in best.items() for t in range(v, 40) if X <= c * (1 << t) < 2 * X)
            P(f"   odd c <= 2^{k}: odd parts covered {cov} of {X // 2}; odd slopes {odd}; integers M in [2^{k}, 2^{k + 1}) in the doubling closure: {dens} ({dens / X:.3f})")
    with open(f"d_{fam}_{tag}.out", "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
