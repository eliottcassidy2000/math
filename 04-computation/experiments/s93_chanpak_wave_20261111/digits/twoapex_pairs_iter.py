#!/usr/bin/env python3
"""Lane "digits", finishing pass: exact affine pairs inside the two-apex family G(2,0), iterative version.

Same question as twoapex_tree.py (whose recursive JIT code crashes for LEN >= 8), larger candidate set, no recursion.

Word w over the letters 0 = (0,0), 1 = (1,0), 2 = (0,1), 3 = (1,1) (36-state transfer of seed/two_apex_probe.py).
Candidate = (suffix word p, multiplier a):  the four letters with a in {1, 2, 3}; the sixteen two-letter words with
a in {1, 4, 9}  (these are all the integer eigenvalues >= 1 of words of length <= 2, twoapex_classify_2.out).
b_p(w) = H(w p) - a H(w).
  single test: b_p constant on w, w p, ..., w p^DEPTH  (exact class of the single map; "moving" iff H(w p) != H(w)).
  pair test:   b_p and b_q constant on every node of w {p, q}^j, j <= DEPTH, and the two maps do not commute:
               b_p (a_q - 1) != b_q (a_p - 1).
Exhaustive over all base words w of length <= LEN.

Usage: twoapex_pairs_iter.py LEN DEPTH
"""
import sys
import time

import numpy as np
from numba import njit


@njit(cache=False)
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


@njit(cache=False)
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


@njit(cache=False)
def apply_word(vec, cw, c, cl, out, tmp):
    for k in range(36):
        out[k] = vec[k]
    for i in range(cl[c]):
        l = cw[c, i]
        step(out, l & 1, (l >> 1) & 1, tmp)
        for k in range(36):
            out[k] = tmp[k]


@njit(cache=False)
def run(n, cw, cl, ca, depth, single_pass, single_moving, res, maxres):
    """all base words of length exactly n.  Returns (words, pair passes)."""
    nc = cl.shape[0]
    stack = np.zeros((n + 1, 36), np.int64)
    stack[0, 0] = 1
    stack[0, 18] = 1
    letters = np.zeros(n + 2, np.int64)
    out = np.zeros(36, np.int64)
    tmp = np.zeros(36, np.int64)
    chain = np.zeros((depth + 2, 36), np.int64)
    nnodes = (1 << (depth + 2)) - 1
    tree = np.zeros((nnodes, 36), np.int64)
    th = np.zeros(nnodes, np.int64)
    okc = np.zeros(nc, np.int64)
    bc = np.zeros(nc, np.int64)
    words = 0
    npass = 0
    d = 0
    letters[0] = -1
    while d >= 0:
        if d == n:
            words += 1
            vec = stack[n]
            h = final(vec)
            nok = 0
            for c in range(nc):
                # single test along w p^j
                for k in range(36):
                    chain[0, k] = vec[k]
                hprev = h
                good = True
                b = 0
                for j in range(depth + 1):
                    apply_word(chain[j], cw, c, cl, chain[j + 1], tmp)
                    hj = final(chain[j + 1])
                    bj = hj - ca[c] * hprev
                    if j == 0:
                        b = bj
                    elif bj != b:
                        good = False
                        break
                    hprev = hj
                okc[c] = 1 if good else 0
                bc[c] = b
                if good:
                    nok += 1
                    single_pass[c] += 1
                    if b != (1 - ca[c]) * h:          # H(w p) != H(w)
                        single_moving[c] += 1
            if nok >= 2:
                for c1 in range(nc):
                    if okc[c1] == 0:
                        continue
                    for c2 in range(c1 + 1, nc):
                        if okc[c2] == 0:
                            continue
                        same = cl[c1] == cl[c2]
                        if same:
                            for i in range(cl[c1]):
                                if cw[c1, i] != cw[c2, i]:
                                    same = False
                        if same:
                            continue
                        if bc[c1] * (ca[c2] - 1) == bc[c2] * (ca[c1] - 1):
                            continue
                        # pair test on the binary tree of depth `depth` (+1 level of children for the offsets)
                        for k in range(36):
                            tree[0, k] = vec[k]
                        th[0] = h
                        good = True
                        for i in range((1 << (depth + 1)) - 1):
                            apply_word(tree[i], cw, c1, cl, tree[2 * i + 1], tmp)
                            apply_word(tree[i], cw, c2, cl, tree[2 * i + 2], tmp)
                            th[2 * i + 1] = final(tree[2 * i + 1])
                            th[2 * i + 2] = final(tree[2 * i + 2])
                            if th[2 * i + 1] - ca[c1] * th[i] != bc[c1] or th[2 * i + 2] - ca[c2] * th[i] != bc[c2]:
                                good = False
                                break
                        if good:
                            if npass < maxres:
                                code = 0
                                for i in range(n):
                                    code = code * 4 + letters[i]
                                res[npass, 0] = n
                                res[npass, 1] = code
                                res[npass, 2] = c1
                                res[npass, 3] = c2
                                res[npass, 4] = h
                                res[npass, 5] = bc[c1]
                                res[npass, 6] = bc[c2]
                            npass += 1
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
    return words, npass


def wstr(code, n):
    s = []
    for i in range(n):
        s.append("0123"[code % 4])
        code //= 4
    return "".join(reversed(s))


def main():
    maxlen = int(sys.argv[1]) if len(sys.argv) > 1 else 9
    depth = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    t0 = time.time()
    cands = []
    for l in range(4):
        for a in (1, 2, 3):
            cands.append(((l,), a))
    for l1 in range(4):
        for l2 in range(4):
            for a in (1, 4, 9):
                cands.append(((l1, l2), a))
    nc = len(cands)
    cw = np.zeros((nc, 2), np.int64)
    cl = np.zeros(nc, np.int64)
    ca = np.zeros(nc, np.int64)
    for i, (p, a) in enumerate(cands):
        cl[i] = len(p)
        ca[i] = a
        for j, l in enumerate(p):
            cw[i, j] = l
    print(f"candidates: {nc} (letters with a in 1, 2, 3; two-letter words with a in 1, 4, 9); base words: all of length <= {maxlen}; depth {depth}")
    single_pass = np.zeros(nc, np.int64)
    single_moving = np.zeros(nc, np.int64)
    res = np.zeros((500, 7), np.int64)
    tot = 0
    for n in range(0, maxlen + 1):
        words, npass = run(n, cw, cl, ca, depth, single_pass, single_moving, res, 500)
        tot += npass
        print(f"   length {n}: {words} base words, non-commuting exact pairs passing to depth {depth}: {npass}   [{time.time() - t0:.1f} s]", flush=True)
        for i in range(min(npass, 20)):
            r = res[i]
            p, a1 = cands[int(r[2])]
            q, a2 = cands[int(r[3])]
            print(f"      PASS: base {wstr(int(r[1]), n) or 'empty'} (H = {int(r[4])}): p = {''.join(map(str, p))} (a = {a1}, b = {int(r[5])}), q = {''.join(map(str, q))} (a = {a2}, b = {int(r[6])})")
    print(f"total non-commuting exact pairs (base word, pair) passing to depth {depth}: {tot}")
    print("single maps: candidate (word, a): base words whose p-ray is exact to the depth / of these with H moving")
    for i, (p, a) in enumerate(cands):
        if single_pass[i]:
            print(f"   p = {''.join(map(str, p))}, a = {a}: {int(single_pass[i])} / {int(single_moving[i])}")
    print(f"time {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
