#!/usr/bin/env python3
"""Verifier (lane digits), check V5: the two-apex family G(2,0) (claims D7, D8), own code.

(1) the 36-state transfer against a direct count (own subset DP on the explicit tournament): all words of length <= 7
    and random words of lengths 8..12.
(2) exhaustive search for exact affine pairs: base words of length <= LEN, suffix words p of length 1 or 2 (20 words).
    Unlike the lane, the multiplier is NOT taken from a list: it is read off the data,
        a_p = (H(wpp) - H(wp)) / (H(wp) - H(w))   (must be an integer; H(wp) != H(w)),
    single test: H(w p^(j+1)) = a_p H(w p^j) + b_p for j = 0, 1, 2;
    pair test: both single tests, the two maps do not commute, and the rules hold at every node of w {p, q}^j, j <= 2.

Usage: v5_twoapex.py LEN
"""
import random
import sys
import time

import numpy as np
from numba import njit, prange


@njit(cache=True)
def step(vec, a, b, out):
    for k in range(36):
        out[k] = 0
    for arr in range(2):
        bx = a if arr == 0 else b
        by = b if arr == 0 else a
        for r1 in range(3):
            for r2 in range(3):
                for r3 in range(2):
                    c = vec[arr * 18 + r1 * 6 + r2 * 2 + r3]
                    if c == 0:
                        continue
                    n1 = 1 if bx == 0 else 2
                    out[arr * 18 + n1 * 6 + r2 * 2 + r3] += c
                    if r2 != 0 or bx == 1:
                        n2 = 1 if by == 0 else 2
                        out[arr * 18 + r1 * 6 + n2 * 2 + r3] += c
                    if r3 != 0 or by == 1:
                        out[arr * 18 + r1 * 6 + r2 * 2 + 1] += c


@njit(cache=True)
def final(vec):
    tot = 0
    for arr in range(2):
        for r1 in range(2):
            for r2 in range(2):
                if r2 == 0 and arr == 1:
                    continue
                for r3 in range(2):
                    tot += vec[arr * 18 + r1 * 6 + r2 * 2 + r3]
    return tot


@njit(cache=True)
def H_transfer(word, n):
    vec = np.zeros(36, np.int64)
    out = np.zeros(36, np.int64)
    vec[0] = 1
    vec[18] = 1
    for i in range(n):
        l = word[i]
        step(vec, l & 1, (l >> 1) & 1, out)
        for k in range(36):
            vec[k] = out[k]
    return final(vec)


@njit(cache=True)
def H_direct(word, n):
    """own count: spine 0..n-1 transitive, u = n, v = n + 1, u -> v; letter l: bit 0 = (u -> s_i), bit 1 = (v -> s_i)."""
    V = n + 2
    u = n
    v = n + 1
    adj = np.zeros(V, np.int64)
    for i in range(n):
        for j in range(i + 1, n):
            adj[i] |= 1 << j
        l = word[i]
        if l & 1:
            adj[u] |= 1 << i
        else:
            adj[i] |= 1 << u
        if (l >> 1) & 1:
            adj[v] |= 1 << i
        else:
            adj[i] |= 1 << v
    adj[u] |= 1 << v
    full = 1 << V
    E = np.zeros((full, V), np.int64)
    for x in range(V):
        E[1 << x, x] = 1
    for m in range(1, full):
        for x in range(V):
            c = E[m, x]
            if c == 0:
                continue
            o = adj[x] & ~m
            for y in range(V):
                if (o >> y) & 1:
                    E[m | (1 << y), y] += c
    t = 0
    for x in range(V):
        t += E[full - 1, x]
    return t


@njit(cache=True)
def check_transfer_all(nmax):
    bad = 0
    cnt = 0
    word = np.zeros(nmax, np.int64)
    for n in range(0, nmax + 1):
        for code in range(4 ** n):
            c = code
            for i in range(n):
                word[i] = c & 3
                c >>= 2
            cnt += 1
            if H_transfer(word, n) != H_direct(word, n):
                bad += 1
    return cnt, bad


@njit(cache=True)
def apply_word(vec, pw, pl, out, tmp):
    for k in range(36):
        out[k] = vec[k]
    for i in range(pl):
        l = pw[i]
        step(out, l & 1, (l >> 1) & 1, tmp)
        for k in range(36):
            out[k] = tmp[k]


@njit(cache=True)
def test_base(vec, cw, cl, stats, chain, tree, th, tmp, okc, ac, bc):
    """stats[c, 0] = single-exact & moving count, stats[c, 1 + min(a, 12)] histogram of multipliers;
    returns (prefilter pairs, passing pairs)."""
    nc = cl.shape[0]
    h = final(vec)
    nok = 0
    for c in range(nc):
        for k in range(36):
            chain[0, k] = vec[k]
        x = np.zeros(4, np.int64)
        x[0] = h
        for j in range(3):
            apply_word(chain[j], cw[c], cl[c], chain[j + 1], tmp)
            x[j + 1] = final(chain[j + 1])
        okc[c] = 0
        d0 = x[1] - x[0]
        if d0 == 0:
            continue
        d1 = x[2] - x[1]
        if d1 % d0 != 0:
            continue
        a = d1 // d0
        if x[3] - x[2] != a * d1:
            continue
        okc[c] = 1
        ac[c] = a
        bc[c] = x[1] - a * x[0]
        nok += 1
        stats[c, 0] += 1
        if a >= 0:
            stats[c, 1 + min(a, 12)] += 1
    pre = 0
    npass = 0
    if nok >= 2:
        for c1 in range(nc):
            if okc[c1] == 0:
                continue
            for c2 in range(c1 + 1, nc):
                if okc[c2] == 0:
                    continue
                if bc[c1] * (ac[c2] - 1) == bc[c2] * (ac[c1] - 1):
                    continue
                pre += 1
                for k in range(36):
                    tree[0, k] = vec[k]
                th[0] = h
                good = True
                for i in range(7):
                    apply_word(tree[i], cw[c1], cl[c1], tree[2 * i + 1], tmp)
                    apply_word(tree[i], cw[c2], cl[c2], tree[2 * i + 2], tmp)
                    th[2 * i + 1] = final(tree[2 * i + 1])
                    th[2 * i + 2] = final(tree[2 * i + 2])
                    if th[2 * i + 1] != ac[c1] * th[i] + bc[c1] or th[2 * i + 2] != ac[c2] * th[i] + bc[c2]:
                        good = False
                        break
                if good:
                    npass += 1
    return pre, npass


@njit(parallel=True, cache=True)
def run(n, cw, cl, nthreads_stats):
    """all base words of length exactly n (n >= 2: split over the first two letters)."""
    nc = cl.shape[0]
    stats = np.zeros((16, nc, 14), np.int64)
    pres = np.zeros(16, np.int64)
    passes = np.zeros(16, np.int64)
    words = np.zeros(16, np.int64)
    for pre2 in prange(16):
        vec = np.zeros((n + 1, 36), np.int64)
        tmp = np.zeros(36, np.int64)
        out = np.zeros(36, np.int64)
        chain = np.zeros((4, 36), np.int64)
        tree = np.zeros((15, 36), np.int64)
        th = np.zeros(15, np.int64)
        okc = np.zeros(nc, np.int64)
        ac = np.zeros(nc, np.int64)
        bc = np.zeros(nc, np.int64)
        vec[0, 0] = 1
        vec[0, 18] = 1
        l0 = pre2 & 3
        l1 = pre2 >> 2
        step(vec[0], l0 & 1, (l0 >> 1) & 1, out)
        for k in range(36):
            vec[1, k] = out[k]
        step(vec[1], l1 & 1, (l1 >> 1) & 1, out)
        for k in range(36):
            vec[2, k] = out[k]
        letters = np.zeros(n + 2, np.int64)
        d = 2
        letters[2] = -1
        while d >= 2:
            if d == n:
                words[pre2] += 1
                p, q = test_base(vec[n], cw, cl, stats[pre2], chain, tree, th, tmp, okc, ac, bc)
                pres[pre2] += p
                passes[pre2] += q
                d -= 1
                continue
            letters[d] += 1
            if letters[d] == 4:
                d -= 1
                continue
            l = letters[d]
            step(vec[d], l & 1, (l >> 1) & 1, out)
            for k in range(36):
                vec[d + 1, k] = out[k]
            d += 1
            if d < n:
                letters[d] = -1
    return words.sum(), pres.sum(), passes.sum(), stats.sum(axis=0)


def main():
    maxlen = int(sys.argv[1]) if len(sys.argv) > 1 else 9
    t0 = time.time()
    cnt, bad = check_transfer_all(7)
    print(f"(1) 36-state transfer vs own direct count: all {cnt} words of length <= 7, mismatches {bad}   [{time.time() - t0:.1f} s]")
    rng = random.Random(3)
    rb = 0
    for _ in range(60):
        n = rng.randint(8, 12)
        w = np.array([rng.randrange(4) for _ in range(n)], np.int64)
        if H_transfer(w, n) != H_direct(w, n):
            rb += 1
    print(f"    random words of lengths 8..12: 60 tests, mismatches {rb}   [{time.time() - t0:.1f} s]")
    cands = [(l,) for l in range(4)] + [(l1, l2) for l1 in range(4) for l2 in range(4)]
    nc = len(cands)
    cw = np.zeros((nc, 2), np.int64)
    cl = np.zeros(nc, np.int64)
    for i, p in enumerate(cands):
        cl[i] = len(p)
        for j, l in enumerate(p):
            cw[i, j] = l
    tot_stats = np.zeros((nc, 14), np.int64)
    tw = tp = tq = 0
    # lengths 0, 1 by the same routine applied directly
    for n in range(0, maxlen + 1):
        if n < 2:
            # tiny: do it with a one-thread helper on explicit vectors
            vecs = []
            for code in range(4 ** n):
                v = np.zeros(36, np.int64)
                o = np.zeros(36, np.int64)
                v[0] = 1
                v[18] = 1
                c = code
                for i in range(n):
                    step(v, c & 1, (c >> 1) & 1, o)
                    v = o.copy()
                    c >>= 2
                vecs.append(v)
            w = len(vecs)
            p = q = 0
            st = np.zeros((nc, 14), np.int64)
            for v in vecs:
                pp, qq = test_base(v, cw, cl, st, np.zeros((4, 36), np.int64), np.zeros((15, 36), np.int64), np.zeros(15, np.int64),
                                   np.zeros(36, np.int64), np.zeros(nc, np.int64), np.zeros(nc, np.int64), np.zeros(nc, np.int64))
                p += pp
                q += qq
        else:
            w, p, q, st = run(n, cw, cl, 0)
        tw += w
        tp += p
        tq += q
        tot_stats += st
        print(f"    length {n}: {w} base words; non-commuting pairs passing both single tests {p}; passing the tree test to depth 2: {q}   [{time.time() - t0:.1f} s]", flush=True)
    print(f"(2) base words of length <= {maxlen}: {tw}; non-commuting (p, q) with both single rays exact to depth 2: {tp}; "
          f"EXACT NON-COMMUTING PAIRS to depth 2: {tq}")
    print("    single maps exact to depth 2 with H moving: suffix p: count; multipliers seen (a: count)")
    for i, p in enumerate(cands):
        if tot_stats[i, 0]:
            hist = {a: int(tot_stats[i, 1 + a]) for a in range(13) if tot_stats[i, 1 + a]}
            print(f"       p = {''.join(map(str, p))}: {int(tot_stats[i, 0])}; {hist}")
    print(f"time {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
