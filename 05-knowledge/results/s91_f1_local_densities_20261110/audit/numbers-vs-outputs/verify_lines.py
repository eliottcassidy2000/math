#!/usr/bin/env python3
"""Auditor numbers-vs-outputs: independent recomputation of the line counts of THM-5129 (F) / HYP-9661 (A).

Own kernel (tables st/en written for the audit).  Rooted tournaments of order n + 1 are (base of order n, cut), root =
the new vertex z, which beats exactly the vertices of the cut.  Bases: reps_le8.npz (orders 1..8) and codes9.npy.

Usage: python3 verify_lines.py REPS_LE8.npz CODES9.npy
"""
import sys
import time

import numpy as np
from numba import njit, prange

NMAX = 1 << 15


@njit(cache=False)
def tables(n, out, hS, st, en):
    full = (1 << n) - 1
    hS[0] = 1
    for S in range(1, full + 1):
        tot = 0
        for v in range(n):
            st[S, v] = 0
            en[S, v] = 0
            if (S >> v) & 1:
                R = S ^ (1 << v)
                if R == 0:
                    st[S, v] = 1
                    en[S, v] = 1
                else:
                    a = 0
                    b = 0
                    for w in range(n):
                        if (R >> w) & 1:
                            if (out[v] >> w) & 1:
                                a += st[R, w]
                            else:
                                b += en[R, w]
                    st[S, v] = a
                    en[S, v] = b
                tot += st[S, v]
        hS[S] = tot


@njit(parallel=True, cache=False)
def scan(n, reps, M, nmax, small, dmax):
    """marks M[N >> 1, D >> 1] for N = D + 2 U 2^j < nmax; small[D, U] = 1 for D <= dmax, U < small.shape[1]."""
    full = (1 << n) - 1
    for r in prange(reps.shape[0]):
        hS = np.zeros(1 << n, np.int64)
        st = np.zeros((1 << n, n), np.int64)
        en = np.zeros((1 << n, n), np.int64)
        out = np.zeros(n, np.int64)
        Q = np.zeros((n, n), np.int64)
        ym = np.zeros(n, np.int64)
        code = reps[r]
        for j in range(1, n):
            base = j * (j - 1) // 2
            for i in range(j):
                if (code >> (base + i)) & 1:
                    out[i] |= (1 << j)
                else:
                    out[j] |= (1 << i)
        tables(n, out, hS, st, en)
        for S in range(1, full + 1):
            C = full ^ S
            for a in range(n):
                if (S >> a) & 1:
                    ym[a] += st[S, a] * hS[C]
                    if en[S, a]:
                        for b in range(n):
                            if (C >> b) & 1:
                                Q[a, b] += en[S, a] * st[C, b]
        HB = hS[full]
        for cut in range(1 << n):
            H = 0
            U = HB
            for b in range(n):
                if (cut >> b) & 1:
                    H += st[full, b]
                    U += ym[b]
                else:
                    H += en[full, b]
                    for c in range(n):
                        if (cut >> c) & 1:
                            H += Q[b, c]
            if H <= dmax and U < small.shape[1]:
                small[H, U] = 1
            N = H + 2 * U
            while N < nmax:
                M[N >> 1, H >> 1] = 1
                U *= 2
                N = H + 2 * U


@njit(cache=False)
def brute_HU(n, out, y):
    """literal: H by DP on the full tournament, U from the definition by DP over subsets (for validation only)."""
    full = (1 << n) - 1
    hS = np.zeros(1 << n, np.int64)
    st = np.zeros((1 << n, n), np.int64)
    en = np.zeros((1 << n, n), np.int64)
    tables(n, out, hS, st, en)
    u = 0
    for S in range(1, full + 1):
        if (S >> y) & 1:
            u += st[S, y] * hS[full ^ S]
    return hS[full], u


@njit(cache=False)
def validate(n, reps, nsamp, seed):
    """compare the cut formulas with the direct computation on the (n+1)-vertex tournament for random (rep, cut)."""
    np.random.seed(seed)
    full = (1 << n) - 1
    hS = np.zeros(1 << n, np.int64)
    st = np.zeros((1 << n, n), np.int64)
    en = np.zeros((1 << n, n), np.int64)
    out = np.zeros(n, np.int64)
    out2 = np.zeros(n + 1, np.int64)
    bad = 0
    for t in range(nsamp):
        r = np.random.randint(0, reps.shape[0])
        cut = np.random.randint(0, 1 << n)
        code = reps[r]
        for i in range(n):
            out[i] = 0
        for j in range(1, n):
            base = j * (j - 1) // 2
            for i in range(j):
                if (code >> (base + i)) & 1:
                    out[i] |= (1 << j)
                else:
                    out[j] |= (1 << i)
        tables(n, out, hS, st, en)
        H = 0
        U = hS[full]
        for b in range(n):
            if (cut >> b) & 1:
                H += st[full, b]
                for S in range(1, full + 1):
                    if (S >> b) & 1:
                        U += st[S, b] * hS[full ^ S]
            else:
                H += en[full, b]
                for c in range(n):
                    if (cut >> c) & 1:
                        for S in range(1, full + 1):
                            if ((S >> b) & 1) and not ((S >> c) & 1):
                                H += en[S, b] * st[full ^ S, c]
        for v in range(n):
            out2[v] = out[v]
            if not ((cut >> v) & 1):
                out2[v] |= (1 << n)
        out2[n] = cut
        h2, u2 = brute_HU(n + 1, out2, n)
        if h2 != H or u2 != U:
            bad += 1
    return bad


def stats(M, label):
    r = M.sum(axis=1).astype(np.int64)
    odd = np.arange(1, NMAX, 2)
    rr = r[odd >> 1]
    print(f"{label}: odd N < {NMAX} with r(N) = 0: {odd[rr == 0].tolist()}")
    shares = []
    for k in range(6, 15):
        sel = (odd >= (1 << k)) & (odd < (1 << (k + 1)))
        c = rr[sel]
        nn = odd[sel]
        amin = nn[c == c.min()].tolist()
        sh = float(np.mean(c / (nn / 6.0)))
        exact_pts = (nn // 3 + 1) // 2        # number of odd D <= N/3
        sh_exact = float(c.sum() / exact_pts.sum())
        shares.append(sh)
        print(f"  [2^{k}, 2^{k + 1}): min {c.min()} at N = {amin[:6]}{'...' if len(amin) > 6 else ''}, median {np.median(c):.1f}, mean {c.mean():.2f}, "
              f"share (mean of r/(N/6)) {sh:.4f}; pooled share of the exact lattice points {sh_exact:.4f}; min of r/(N/64) {float(np.min(c * 64.0 / nn)):.3f}")
    sel = (odd >= 129) & (odd < (1 << 14))
    viol = odd[sel][rr[sel] * 64 < odd[sel]]
    print(f"  r(N) >= N/64 on [129, 2^14): violations {len(viol)}" + (f", first {viol[:8].tolist()}, least {viol.min()}" if len(viol) else ""))
    ratio = rr[sel] * 64.0 / odd[sel]
    print(f"  least value of r(N) / (N/64) on [129, 2^14): {ratio.min():.4f} at N = {odd[sel][ratio.argmin()]}")
    sel2 = (odd >= 129) & (odd < (1 << 15))
    viol2 = odd[sel2][rr[sel2] * 64 < odd[sel2]]
    print(f"  r(N) >= N/64 on [129, 2^15): violations {len(viol2)}" + (f", least {viol2.min()}" if len(viol2) else ""))
    print("  shares increasing over k = 6..11:", all(shares[i] < shares[i + 1] for i in range(5)), [round(s, 4) for s in shares])
    for lo, hi in ((1 << 10, 1 << 11), (1 << 12, 1 << 13)):
        num = np.zeros((2, 10))
        den = np.zeros((2, 10))
        for N in range(lo + 1, hi, 2):
            D = np.arange(1, N // 3 + 1, 2)
            Uv = (N - D) // 2
            k = np.minimum(9, (30 * D) // N)
            par = Uv & 1
            hit = M[N >> 1, D >> 1]
            np.add.at(den, (par, k), 1)
            np.add.at(num, (par, k), hit)
        print(f"  N in [{lo}, {hi}): share of lattice points that are sources by decile of 3D/N")
        print("     U odd :", " ".join(f"{num[1, k] / den[1, k]:.4f}" for k in range(10)), f" overall {num[1].sum() / den[1].sum():.4f}")
        print("     U even:", " ".join(f"{num[0, k] / den[0, k]:.4f}" for k in range(10)), f" overall {num[0].sum() / den[0].sum():.4f}")
        print("     both  :", " ".join(f"{(num[0, k] + num[1, k]) / (den[0, k] + den[1, k]):.4f}" for k in range(10)),
              f" overall {num.sum() / den.sum():.4f}")
    return rr


def main():
    z = np.load(sys.argv[1])
    codes9 = np.load(sys.argv[2]).astype(np.int64)
    reps = {n: z["n%d" % n].astype(np.int64) for n in range(1, 9)}
    t0 = time.time()
    for n in (3, 5, 7, 8):
        print(f"validation of the cut formulas against the direct computation, bases of order {n}: failures", int(validate(n, reps[n], 300, 12345 + n)))
    print("validation, bases of order 9 (codes9):", int(validate(9, codes9, 200, 999)))
    # codes9 = all (rep of order 8, cut) such that the new vertex has least out-degree ?
    mine = set()
    for code in reps[8].tolist():
        outdeg = [0] * 8
        for j in range(1, 8):
            base = j * (j - 1) // 2
            for i in range(j):
                if (code >> (base + i)) & 1:
                    outdeg[i] += 1
                else:
                    outdeg[j] += 1
        for cut in range(256):
            k = bin(cut).count("1")
            ok = True
            c9 = code
            for b in range(8):
                if (cut >> b) & 1:
                    if outdeg[b] < k:
                        ok = False
                        break
                else:
                    if outdeg[b] + 1 < k:
                        ok = False
                        break
                    c9 |= 1 << (28 + b)
            if ok:
                mine.add(c9)
    theirs = set(codes9.tolist())
    print(f"codes9.npy: {len(codes9)} codes, {len(theirs)} distinct; (class rep of order 8, cut) with the new vertex of least out-degree: {len(mine)}; "
          f"equal as sets: {mine == theirs}; codes9 subset of mine: {theirs <= mine}; mine subset of codes9: {mine <= theirs}")
    M = np.zeros((NMAX >> 1, (NMAX // 3 >> 1) + 2), dtype=np.uint8)
    small = np.zeros((74, 1 << 18), dtype=np.uint8)
    U = 1
    while 1 + 2 * U < NMAX:
        M[(1 + 2 * U) >> 1, 0] = 1
        U *= 2
    small[1, 1] = 1
    for n in range(1, 9):
        scan(n, reps[n], M, NMAX, small, 73)
    print(f"rooted order <= 9: number of pairs (D, U) with D <= 73: {int(small.sum())}  (largest U among them: {int(np.nonzero(small.sum(axis=0))[0].max())})  [{time.time() - t0:.0f} s]")
    r9 = stats(M, "rooted order <= 9")
    scan(9, codes9, M, NMAX, small, 73)
    print(f"rooted order <= 10: number of pairs (D, U) with D <= 73, U < 2^18: {int(small.sum())}  [{time.time() - t0:.0f} s]")
    r10 = stats(M, "rooted order <= 10")
    odd = np.arange(1, NMAX, 2)
    print("r_10 >= r_9 everywhere:", bool(np.all(r10 >= r9)))
    for N in (131, 259, 521, 1027, 2101, 4109, 8203, 16447, 129, 4117, 15427, 32615):
        print(f"  N = {N}: r_9 = {int(r9[N >> 1])}, r_10 = {int(r10[N >> 1])}, N/64 = {N / 64:.2f}")
    np.save("r10_odd.npy", r10)
    np.save("r9_odd.npy", r9)


if __name__ == "__main__":
    main()
