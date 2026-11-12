#!/usr/bin/env python3
"""Independent check (verifier C) of the lattice-line count of THM-5141 (C) by DIRECT solution of the equation.

usage: python c3_lattice_direct.py NMAX      (default 15)

For every n <= NMAX, every split n = t' + w + s (t' >= 1, w >= 1, s >= 0), EVERY odd pair (x, x') of w-bit windows
(solutions or not), every P in [2^(t'-1), 2^t') and each parity cp of c:
    gamma = 2(q(x) - q(x')) + (P mod 2)(x_(w-1) - x'_(w-1)) + cp (x_0 - x'_0),   c = (P Xi + gamma)/X
is a solution iff c is an integer of (2^s, 2^(s+1)] with c = cp (mod 2).
 (1) the set of solutions, mapped back to (h, h') = (P|x|Z, P|x'|Z) with Z = rev_s(2^(s+1) - c), is compared EXACTLY
     with the set of ordered colliding pairs x(h) = x(h') of n-bit words obtained by sorting the values;
 (2) for every (t', s, x, x') the number of solutions is compared with 4 (1 + 2^max(t'-1,s) gcd/max(|X|,|Xi|));
 (3) for every parity class (P mod 2, c mod 2) the number is compared with 1 + g min(2^(t'-1)/|X|, 2^s/|Xi|)
     and with the exact progression count 1 + floor((2^(t'-1) - 1) g / (2|X|)), 1 + floor((2^s - 1) g / (2|Xi|));
 (4) for each parity class the solutions are checked to be equally spaced with step a multiple of (Xi, X)/g.
"""
import sys
import time
import numpy as np
from numba import njit

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 15


@njit(cache=False)
def build_q(nmax):
    L = 1 << nmax
    q = np.zeros(L, dtype=np.int64)
    for h in range(1, L):
        if h & 1:
            q[h] = q[h >> 1] + (((h >> 1) + 1) >> 1)
        else:
            q[h] = q[h >> 1]
    return q


@njit(cache=False)
def revw(x, w):
    r = 0
    for i in range(w):
        if (x >> i) & 1:
            r |= 1 << (w - 1 - i)
    return r


@njit(cache=False)
def gcd(a, b):
    while b:
        a, b = b, a % b
    return a


@njit(cache=False)
def collisions(n, q):
    lo = 1 << (n - 1)
    m = lo
    vals = np.empty(m, dtype=np.int64)
    for i in range(m):
        h = lo + i
        vals[i] = 2 * h + 1 - 2 * q[h]
    order = np.argsort(vals)
    total = 0
    i = 0
    while i < m:
        j = i + 1
        while j < m and vals[order[j]] == vals[order[i]]:
            j += 1
        total += (j - i) * (j - i - 1)
        i = j
    codes = np.empty(total, dtype=np.int64)
    k = 0
    i = 0
    while i < m:
        j = i + 1
        while j < m and vals[order[j]] == vals[order[i]]:
            j += 1
        for a in range(i, j):
            for b in range(i, j):
                if a != b:
                    codes[k] = ((lo + order[a]) << 24) | (lo + order[b])
                    k += 1
        i = j
    return np.sort(codes)


@njit(cache=False)
def direct(n, q, codes, st, fl):
    """codes: output array (size = expected number); st: int stats; fl: float stats.
    st[0] solutions, st[1] (t',s,x,x') groups examined, st[2] groups with a solution, st[3] violations of the
    stated bound, st[4] parity classes examined, st[5] violations of 1 + g min, st[6] violations of the exact
    progression count, st[7] overflow of codes, st[8] largest group, st[9] spacing violations,
    st[10] classes where count == exact progression count (bound attained) , st[11] classes with a solution"""
    k = 0
    cap = codes.shape[0]
    lastP = np.full(4, -1, dtype=np.int64)
    lastc = np.full(4, -1, dtype=np.int64)
    for tp in range(1, n):
        for w in range(1, n - tp + 1):
            s = n - tp - w
            top = 1 << (w - 1)
            Plo = 1 << (tp - 1)
            Phi = 1 << tp
            clo = 1 << s
            chi = 1 << (s + 1)
            e = max(tp - 1, s)
            for x in range(1 << w):
                rx = revw(x, w)
                for xp in range(1 << w):
                    if ((x ^ xp) & 1) == 0 or ((x ^ xp) & top) == 0:
                        continue
                    X = x - xp
                    Xi = rx - revw(xp, w)
                    aX = abs(X)
                    aXi = abs(Xi)
                    g = gcd(aX, aXi)
                    M = max(aX, aXi)
                    base = 2 * (q[x] - q[xp])
                    dt = ((x >> (w - 1)) & 1) - ((xp >> (w - 1)) & 1)
                    db = (x & 1) - (xp & 1)
                    tot = 0
                    c00 = 0
                    c01 = 0
                    c10 = 0
                    c11 = 0
                    # last solution seen in each class, for the spacing check: index = 2*pp + cp
                    for ii in range(4):
                        lastP[ii] = -1
                        lastc[ii] = -1
                    for P in range(Plo, Phi):
                        pp = P & 1
                        for cp in range(2):
                            gam = base + pp * dt + cp * db
                            num = P * Xi + gam
                            if num % X != 0:
                                continue
                            c = num // X
                            if c <= clo or c > chi:
                                continue
                            if (c & 1) != cp:
                                continue
                            tot += 1
                            cls = 2 * pp + cp
                            if cls == 0:
                                c00 += 1
                            elif cls == 1:
                                c01 += 1
                            elif cls == 2:
                                c10 += 1
                            else:
                                c11 += 1
                            if lastP[cls] >= 0:
                                dP = P - lastP[cls]
                                dc = c - lastc[cls]
                                # (dc, dP) = j (Xi, X)/g with j an integer (j even because both differences are even)
                                if dc * X != dP * Xi or (dP * g) % aX != 0 or ((dP * g) // aX) % 2 != 0:
                                    st[9] += 1
                            lastP[cls] = P
                            lastc[cls] = c
                            Z = revw(chi - c, s)
                            h = (P << (w + s)) | (x << s) | Z
                            hp = (P << (w + s)) | (xp << s) | Z
                            if k < cap:
                                codes[k] = (h << 24) | hp
                            else:
                                st[7] += 1
                            k += 1
                    st[1] += 1
                    if tot > 0:
                        st[2] += 1
                    if tot * M > 4 * (M + (1 << e) * g):
                        st[3] += 1
                    r = tot * M / float(4 * (M + (1 << e) * g))
                    if r > fl[0]:
                        fl[0] = r
                    if tot > st[8]:
                        st[8] = tot
                    b1 = 1 + ((1 << (tp - 1)) - 1) * g // (2 * aX)
                    b2 = 1 + ((1 << s) - 1) * g // (2 * aXi)
                    for cls in range(4):
                        if cls == 0:
                            cc = c00
                        elif cls == 1:
                            cc = c01
                        elif cls == 2:
                            cc = c10
                        else:
                            cc = c11
                        st[4] += 1
                        if cc > 0:
                            st[11] += 1
                            if (cc - 1) * aX > g * (1 << (tp - 1)) or (cc - 1) * aXi > g * (1 << s):
                                st[5] += 1
                            if cc > b1 or cc > b2:
                                st[6] += 1
                            if cc == min(b1, b2):
                                st[10] += 1
    st[0] = k


def main():
    t0 = time.time()
    q = build_q(max(NMAX, 2))
    print("THM-5141 (C), lattice line by direct solution of c X - P Xi = gamma, n = 2 ..", NMAX, flush=True)
    for n in range(2, NMAX + 1):
        ref = collisions(n, q)
        codes = np.empty(ref.shape[0] + 16, dtype=np.int64)
        st = np.zeros(16, dtype=np.int64)
        fl = np.zeros(4, dtype=np.float64)
        direct(n, q, codes, st, fl)
        k = int(st[0])
        same = (k == ref.shape[0]) and bool(np.array_equal(np.sort(codes[:k]), ref))
        print("n = %2d: solutions of the equation in the boxes: %8d; ordered colliding pairs: %8d; the two sets are equal: %s"
              % (n, k, ref.shape[0], same), flush=True)
        print("        groups (t', s, x, x') examined (all odd pairs): %d, with a solution: %d; violations of 4(1 + 2^max(t'-1,s) g/max): %d;"
              % (st[1], st[2], st[3]),
              "max count/bound %.4f; largest group %d" % (fl[0], st[8]), flush=True)
        print("        parity classes examined: %d, with a solution: %d; violations of 1 + g min(...): %d; of the exact progression count: %d;"
              % (st[4], st[11], st[5], st[6]),
              "classes attaining the exact progression count: %d; spacing violations: %d" % (st[10], st[9]), flush=True)
    print("time %.1f s" % (time.time() - t0))


main()
