#!/usr/bin/env python3
"""Lane "digits": direct search for exact binary digit pairs inside the two-apex family G(2,0).

For a word w (letters 0 = (0,0), 1 = (1,0), 2 = (0,1), 3 = (1,1); 36-state transfer) and two suffix words p, q:
    b_p(w) = H(w p) - a_p H(w),   b_q(w) = H(w q) - a_q H(w)      (a_p, a_q given multipliers).
w generates an exact affine pair iff b_p and b_q are constant on the whole sub-monoid orbit w {p, q}^*.
This script tests the orbit to depth DEPTH for every word w of length <= LEN (exhaustive), for the pairs
(p, q) listed in PAIRS, and reports the words that pass with b_p != b_q-type non-commuting maps.

Usage: twoapex_tree.py LEN DEPTH
"""
import sys
import time

import numpy as np
from numba import njit


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
def apply_word(vec, word, wlen, out, tmp):
    for k in range(36):
        out[k] = vec[k]
    for i in range(wlen):
        l = word[i]
        step(out, l & 1, (l >> 1) & 1, tmp)
        for k in range(36):
            out[k] = tmp[k]


@njit(cache=True)
def orbit_ok(vec, p, plen, ap, q, qlen, aq, depth, bp, bq, bufs, level):
    """is b_p = bp and b_q = bq on the orbit of vec under {p, q} to the given depth?"""
    h = final(vec)
    tmp = np.zeros(36, np.int64)
    vp = bufs[level, 0]
    vq = bufs[level, 1]
    apply_word(vec, p, plen, vp, tmp)
    apply_word(vec, q, qlen, vq, tmp)
    if final(vp) - ap * h != bp:
        return False
    if final(vq) - aq * h != bq:
        return False
    if depth == 0:
        return True
    if not orbit_ok(vp, p, plen, ap, q, qlen, aq, depth - 1, bp, bq, bufs, level + 1):
        return False
    if not orbit_ok(vq, p, plen, ap, q, qlen, aq, depth - 1, bp, bq, bufs, level + 1):
        return False
    return True


@njit(cache=True)
def search(n, p, plen, ap, q, qlen, aq, depth, res, maxres):
    """all words of length exactly n; returns the number of words passing; res rows: code, H, bp, bq."""
    stack = np.zeros((n + 1, 36), np.int64)
    stack[0, 0] = 1
    stack[0, 18] = 1
    letters = np.zeros(n + 1, np.int64)
    bufs = np.zeros((depth + 3, 2, 36), np.int64)
    tmp = np.zeros(36, np.int64)
    vp = np.zeros(36, np.int64)
    vq = np.zeros(36, np.int64)
    out = np.zeros(36, np.int64)
    cnt = 0
    cnt_nc = 0
    if n == 0:
        d = 0
    d = 0
    letters[0] = -1
    while d >= 0:
        if d == n:
            vec = stack[n]
            h = final(vec)
            apply_word(vec, p, plen, vp, tmp)
            apply_word(vec, q, qlen, vq, tmp)
            bp = final(vp) - ap * h
            bq = final(vq) - aq * h
            # non-commuting maps and h moving
            if bp * (aq - 1) != bq * (ap - 1):
                if orbit_ok(vec, p, plen, ap, q, qlen, aq, depth, bp, bq, bufs, 0):
                    code = 0
                    for i in range(n):
                        code = code * 4 + letters[i]
                    if cnt_nc < maxres:
                        res[cnt_nc, 0] = code
                        res[cnt_nc, 1] = h
                        res[cnt_nc, 2] = bp
                        res[cnt_nc, 3] = bq
                    cnt_nc += 1
            cnt += 1
            d -= 1
            continue
        letters[d] += 1
        if letters[d] == 4:
            d -= 1
            continue
        l = letters[d]
        step(stack[d], l & 1, (l >> 1) & 1, out)
        for k in range(36):
            stack[d + 1, k] = out[k]
        d += 1
        if d < n:
            letters[d] = -1
    return cnt, cnt_nc


def wstr(code, n):
    s = []
    for i in range(n):
        s.append("0123"[code % 4])
        code //= 4
    return "".join(reversed(s))


def main():
    maxlen = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    depth = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    t0 = time.time()
    # suffix words with their multipliers: single letters 0, 1, 2 (multipliers 2 and 3), and two-letter words (4, 9)
    cands = []
    for l in (0, 1, 2):
        cands.append(((l,), 2))
        cands.append(((l,), 3))
    for l1 in range(4):
        for l2 in range(4):
            if (l1, l2) != (3, 3):
                cands.append(((l1, l2), 4))
    pairs = []
    for i in range(len(cands)):
        for j in range(i + 1, len(cands)):
            if cands[i][0] == cands[j][0]:
                continue
            pairs.append((cands[i], cands[j]))
    print(f"suffix candidates: {len(cands)}; pairs tested: {len(pairs)}; base words: all of length <= {maxlen}; orbit depth {depth}")
    res = np.zeros((200, 4), np.int64)
    total_hits = 0
    for (p, ap), (q, aq) in pairs:
        pa = np.array(p, np.int64)
        qa = np.array(q, np.int64)
        hits = []
        for n in range(0, maxlen + 1):
            cnt, nc = search(n, pa, len(p), ap, qa, len(q), aq, depth, res, 200)
            for i in range(min(nc, 200)):
                hits.append((n, wstr(int(res[i, 0]), n), int(res[i, 1]), int(res[i, 2]), int(res[i, 3])))
            if nc > 200:
                hits.append((n, "...", nc, 0, 0))
        if hits:
            total_hits += len(hits)
            print(f"PAIR p = {''.join(map(str, p))} (a = {ap}), q = {''.join(map(str, q))} (a = {aq}): {len(hits)} base words pass to depth {depth} with non-commuting maps: {hits[:12]}")
    print(f"total passing (pair, base word): {total_hits}   [{time.time() - t0:.1f} s]")


if __name__ == "__main__":
    main()
