#!/usr/bin/env python3
"""S91 (THM-5129 (B')): the complete supply of one-vertex sources at the offsets 129 = 2^7 + 1 and 257 = 2^8 + 1.

By THM-5129 (A) the slopes of offset D are Slopes_0(D) * 2^j, built from the strong tournaments whose path counts
divide into D.  129 = 3 * 43 needs the strong tournaments with H = 3, 43, 129; 257 is prime and needs those with
H = 257.  The least path count of a strong tournament of order n is f(n) = 3, 5, 9, 15, 25, 45, 75, 125, 225, 375
for n = 3..12 (Busch's formula; orders <= 10 exhaustive in this repository, and f(11) = 225, f(12) > 257 are
verified below, so nothing is taken from the literature), and f is nondecreasing (THM-5129 (A2)).  So H = 129
forces order <= 10 and H = 257 forces order <= 11.

Enumeration (exhaustive):
  * orders <= 8: class representatives reps_le8.npz; order 9: the list codes9.npy (every class);
  * order 10: every tournament of order 10 is (a tournament of the list codes9, a cut of the new vertex); all pairs
    with H <= 257 are kept (labelled, with repetitions);
  * order 11: a strong tournament T of order 11 with H(T) <= 257 has a vertex v with T - v strong, and then
    H(T - v) <= H(T) - 2; so T is (a strong tournament of order 10 with H <= 255 from the previous list, a cut);
    all such extensions with H <= 257 are kept: their least strong path count is f(11), none has H = 129;
  * order 12: the same step from the strong tournaments of order 11 with H <= 255 finds no strong tournament of
    order 12 with H <= 257, so f(12) > 257 and no strong tournament of order >= 12 has H <= 257.
For every strong tournament found with H in {43, 129, 257}: H, Pi_2 and U at every root (source_supply.strong_data).

Output: the numbers of labelled tournaments found, the sets Root and Pi_2, the odd parts of the slopes at the
offsets 129 and 257 with their least 2-exponents (complete at every order of the base), the rows of the 9-grid, and
the cross-check that every pair of the order <= 9 census with these offsets is in the list.
Exit status 1 on failure.  Needs codes9.npy (written by census.py; not kept in the repository); about 3 minutes.
"""
import os
import sys
import time
from collections import defaultdict

import numpy as np
from numba import njit, prange

import source_supply as ss

HERE = os.path.dirname(os.path.abspath(__file__))


@njit(cache=True)
def base_tables(n, code, f, g, hS, out, Q):
    full = (1 << n) - 1
    for i in range(n):
        out[i] = 0
    for j in range(1, n):
        base = j * (j - 1) // 2
        for i in range(j):
            if (code >> (base + i)) & 1:
                out[i] |= (1 << j)
            else:
                out[j] |= (1 << i)
    hS[0] = 1
    for S in range(1, full + 1):
        tot = 0
        for v in range(n):
            f[S, v] = 0
            g[S, v] = 0
            if ((S >> v) & 1) == 0:
                continue
            R = S ^ (1 << v)
            if R == 0:
                f[S, v] = 1
                g[S, v] = 1
                tot += 1
                continue
            a = 0
            b = 0
            for w in range(n):
                if (R >> w) & 1:
                    if (out[w] >> v) & 1:
                        a += f[R, w]
                    else:
                        b += g[R, w]
            f[S, v] = a
            g[S, v] = b
            tot += a
        hS[S] = tot
    for a in range(n):
        for b in range(n):
            Q[a, b] = 0
    for S in range(1, full + 1):
        C = full ^ S
        for a in range(n):
            if (S >> a) & 1 and f[S, a]:
                for b in range(n):
                    if (C >> b) & 1:
                        Q[a, b] += f[S, a] * g[C, b]


@njit(parallel=True, cache=True)
def extend(n, codes, hbase_max, hlo, hhi, counts, outbuf, fill):
    """For every base with H <= hbase_max and every cut: if hlo <= H(base + new vertex) <= hhi, count it and, when
    fill is set, store the code of the extension in outbuf[r, :]."""
    full = (1 << n) - 1
    shift = n * (n - 1) // 2
    for r in prange(codes.shape[0]):
        f = np.zeros((1 << n, n), np.int64)
        g = np.zeros((1 << n, n), np.int64)
        hS = np.zeros(1 << n, np.int64)
        out = np.zeros(n, np.int64)
        Q = np.zeros((n, n), np.int64)
        base_tables(n, codes[r], f, g, hS, out, Q)
        c = 0
        if hS[full] <= hbase_max:
            for cut in range(1 << n):
                H = 0
                for b in range(n):
                    if (cut >> b) & 1:
                        H += g[full, b]
                    else:
                        H += f[full, b]
                        for cc in range(n):
                            if (cut >> cc) & 1:
                                H += Q[b, cc]
                if hlo <= H <= hhi:
                    if fill:
                        outbuf[r, c] = codes[r] | (((~cut) & full) << shift)
                    c += 1
        counts[r] = c


def extensions(n, codes, hbase_max, hlo, hhi):
    counts = np.zeros(len(codes), np.int64)
    dummy = np.zeros((1, 1), np.int64)
    extend(n, codes, hbase_max, hlo, hhi, counts, dummy, False)
    width = int(counts.max()) if len(counts) else 0
    if width == 0:
        return np.zeros(0, np.int64)
    sel = np.nonzero(counts)[0]
    sub = codes[sel]
    buf = np.zeros((len(sub), width), np.int64)
    c2 = np.zeros(len(sub), np.int64)
    extend(n, sub, hbase_max, hlo, hhi, c2, buf, True)
    parts = [buf[i, :c2[i]] for i in range(len(sub))]
    return np.unique(np.concatenate(parts))


def odd_part(u):
    j = (u & -u).bit_length() - 1
    return u >> j, j


def main():
    t0 = time.time()
    ok = True
    root = defaultdict(set)
    pi2 = defaultdict(set)
    nfound = defaultdict(lambda: defaultdict(int))
    targets = (3, 43, 129, 257)

    def absorb(n, codes):
        res = np.zeros((len(codes), 3 + n), np.int64)
        ss.strong_data(n, codes, res)
        for row in res:
            h = int(row[0])
            if row[2] == 1 and h in targets:
                nfound[h][n] += 1
                pi2[h].add(int(row[1]))
                for y in range(n):
                    root[h].add(int(row[3 + y]))
        return res

    z = np.load(os.path.join(HERE, "reps_le8.npz"))
    for n in range(3, 9):
        absorb(n, z["n%d" % n].astype(np.int64))
    codes9 = np.load(os.path.join(HERE, "codes9.npy")).astype(np.int64)
    res9 = absorb(9, codes9)
    print(f"orders <= 9 done ({time.time() - t0:.0f} s)", flush=True)
    # self-test of the extension step: the strong tournaments of order 9 with H <= 150 obtained as extensions of the
    # class representatives of order 8 must be the same isomorphism-invariant signatures (H, Pi_2, sorted U) as
    # those of the list codes9
    def signatures(res, hmax):
        out = set()
        for row in res:
            if row[2] == 1 and row[0] <= hmax:
                out.add((int(row[0]), int(row[1]), tuple(sorted(int(v) for v in row[3:]))))
        return out
    t9 = extensions(8, z["n8"].astype(np.int64), 150, 1, 150)
    rt9 = np.zeros((len(t9), 3 + 9), np.int64)
    ss.strong_data(9, t9, rt9)
    sig_a = signatures(rt9, 150)
    sig_b = signatures(res9, 150)
    print(f"self-test (order 9, H <= 150): {len(sig_a)} signatures of strong tournaments by extension of order 8, "
          f"{len(sig_b)} from the list codes9, equal: {sig_a == sig_b}")
    ok &= sig_a == sig_b and len(sig_a) > 0
    # order 10: all (base of order 9, cut) with H <= 257
    c10 = extensions(9, codes9, 257, 1, 257)
    res10 = absorb(10, c10)
    strong10 = c10[(res10[:, 2] == 1) & (res10[:, 0] <= 255)]
    print(f"order 10: {len(c10)} labelled tournaments with H <= 257, of which {int((res10[:, 2] == 1).sum())} strong; "
          f"least H of a strong one: {int(res10[res10[:, 2] == 1][:, 0].min())}  ({time.time() - t0:.0f} s)", flush=True)
    ok &= int(res10[res10[:, 2] == 1][:, 0].min()) == 125
    # order 11: every extension with H <= 257 of the strong tournaments of order 10 with H <= 255
    c11 = extensions(10, strong10, 255, 1, 257)
    res11 = absorb(11, c11)
    st11 = res11[:, 2] == 1
    f11 = int(res11[st11][:, 0].min())
    print(f"order 11: {len(c11)} labelled extensions with H <= 257 of the {len(strong10)} strong tournaments of order 10 "
          f"with H <= 255, of which {int(st11.sum())} strong; least H of a strong one: {f11}; strong with H = 129: "
          f"{int((st11 & (res11[:, 0] == 129)).sum())}; strong with H = 257: {int((st11 & (res11[:, 0] == 257)).sum())}  "
          f"({time.time() - t0:.0f} s)", flush=True)
    ok &= f11 == 225 and int((st11 & (res11[:, 0] == 129)).sum()) == 0
    # order 12: extensions with H <= 257 of the strong tournaments of order 11 with H <= 255: none is strong
    # (an order-12 tournament has 66 arcs and does not fit a 64-bit code: the cuts are counted, not stored.  The
    # extension of a strong base is strong iff the cut is neither empty nor full, and those two cuts keep H.)
    strong11 = c11[st11 & (res11[:, 0] <= 255)]
    cnt12 = np.zeros(len(strong11), np.int64)
    extend(11, strong11, 255, 1, 257, cnt12, np.zeros((1, 1), np.int64), False)
    n12s = int((cnt12 - 2).sum())
    print(f"order 12: the {len(strong11)} strong tournaments of order 11 with H <= 255 have {int(cnt12.sum())} cuts "
          f"with H <= 257, two per base (the empty and the full cut): {bool(np.all(cnt12 == 2))}; strong extensions "
          f"with H <= 257: {n12s} (so f(12) > 257)  ({time.time() - t0:.0f} s)", flush=True)
    ok &= n12s == 0 and bool(np.all(cnt12 == 2))
    for h in targets:
        print(f"strong tournaments with H = {h}: labelled copies found by order {dict(nfound[h])}; "
              f"{len(root[h])} values of U, {len(pi2[h])} values of Pi_2")
    ok &= root[3] == {3} and pi2[3] == {12}
    # slopes: offsets 129 = 3 * 43 and 257
    slopes = {}
    s129 = set(root[129]) | set(pi2[129]) | {129}
    for u in root[43]:
        s129.add(3 * u)            # C_3 before, root in the 43-component
        s129.add(12 * u)           # root in the 43-component, C_3 after
    for p in pi2[43]:
        s129.add(3 * p)            # root in C_3 (U = 3), 43-component after; also trivial root between
        s129.add(12 * p)           # trivial root first, both orders of the components (Pi_2 products)
    s129.add(3 * 43)               # both components before a trivial root (= 129, already there)
    s129.add(43 * 3)               # 43-component before, root in C_3
    s129.add(43 * 12)              # 43-component before, trivial root, C_3 after
    slopes[129] = s129
    slopes[257] = set(root[257]) | set(pi2[257]) | {257}
    prim = {}
    for D, s in slopes.items():
        m = {}
        for u in s:
            dl, j = odd_part(u)
            if dl not in m or j < m[dl]:
                m[dl] = j
        prim[D] = m
    for D, lev in ((129, 7), (257, 8)):
        m = prim[D]
        print(f"offset {D} (level s = {lev}): {len(m)} odd parts, {sum(1 for j in m.values() if j <= lev)} with j <= s; "
              f"d:j ('*' = j > s)")
        line = "    "
        for dl, j in sorted(m.items()):
            c = f"{dl}:{j}{'*' if j > lev else ''} "
            if len(line) + len(c) > 118:
                print(line)
                line = "    "
            line += c
        print(line)
        print(f"    9-grid at level {lev}: least exponent j = {m.get(9)} (every order of the base); "
              f"dyadic slope (d = 1): {m.get(1)}")
        hi = (3 * D) // 2
        c = 0
        for u in range(D, hi + 1):
            dl, j = odd_part(u)
            if dl in m and j >= m[dl]:
                c += 1
        print(f"    window [h, 3h/2] at h = {D}: {c} slopes among {hi - D + 1} integers (share {c / (hi - D + 1):.3f})")
    # cross-check with the census of order <= 9
    cz = os.path.join(HERE, "rooted_le9.npz")
    if os.path.exists(cz):
        zz = np.load(cz)
        rows = np.concatenate([zz[f"m{m}"] for m in range(1, 10)], axis=0)
        for D in (129, 257):
            us = set(rows[rows[:, 0] == D][:, 4].tolist())
            bad = []
            for u in us:
                dl, j = odd_part(u)
                if dl not in prim[D] or j < prim[D][dl]:
                    bad.append(u)
            print(f"census of order <= 9 at offset {D}: {len(us)} slopes, not in the list: {len(bad)}")
            ok &= len(bad) == 0
    ok &= prim[129].get(9) == 4 and prim[257].get(9) == 5
    print(f"({time.time() - t0:.0f} s)")
    print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
