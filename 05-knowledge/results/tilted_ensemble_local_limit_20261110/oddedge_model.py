#!/usr/bin/env python3
"""Rank law of a uniformly random alternating m x m matrix over F_2 conditioned on the parity of its number of
edges (entries above the diagonal).  Exact Markov chain on (rank, edge parity, flag), flag = [degree vector is zero]
= [all-ones vector lies in the radical]; checked against brute force for m <= 6."""
from fractions import Fraction
from itertools import combinations
import sys

def chain(mmax):
    # state (r2, e, f): r2 = rank, e = edge parity, f = 1 iff beta*1 = 0 ; counts of matrices
    cur = {(0, 0, 1): 1}        # m = 1: the zero 1x1 matrix
    out = {1: dict(cur)}
    for m in range(1, mmax):
        nxt = {}
        def add(k, v):
            if v: nxt[k] = nxt.get(k, 0) + v
        for (r2, e, f), cnt in cur.items():
            img = 1 << r2
            if f:   # every u in the image has even weight; d = 0 lies in the image (u = 0)
                # u in image: all even weight. new flag = [u == d = 0]
                add((r2, e, 1), cnt * 1)               # u = 0 = d
                add((r2, e, 0), cnt * (img - 1))       # other image vectors, even weight
                ev_out = (1 << (m - 1)) - img          # even-weight vectors outside the image
                od_out = (1 << (m - 1))                # odd-weight vectors (none in the image)
            else:   # half of the image has odd weight; d is a non-zero image vector of even weight
                add((r2, e, 1), cnt * 1)                       # u = d (even weight)
                add((r2, e, 0), cnt * (img // 2 - 1))          # other even-weight image vectors
                add((r2, e ^ 1, 0), cnt * (img // 2))          # odd-weight image vectors
                ev_out = (1 << (m - 1)) - img // 2
                od_out = (1 << (m - 1)) - img // 2
            add((r2 + 2, e, 0), cnt * ev_out)
            add((r2 + 2, e ^ 1, 0), cnt * od_out)
        cur = nxt; out[m + 1] = dict(cur)
    return out

def brute(m):
    prs = list(combinations(range(m), 2)); res = {}
    for code in range(1 << len(prs)):
        rows = [0] * m
        for k, (i, j) in enumerate(prs):
            if code >> k & 1: rows[i] |= 1 << j; rows[j] |= 1 << i
        e = bin(code).count('1') & 1
        d = 0
        for i in range(m): d |= (bin(rows[i]).count('1') & 1) << i
        f = 1 if d == 0 else 0
        rr = rows[:]; rk = 0
        for c in range(m):
            p = next((i for i in range(rk, m) if rr[i] >> c & 1), None)
            if p is None: continue
            rr[p], rr[rk] = rr[rk], rr[p]
            for i in range(m):
                if i != rk and rr[i] >> c & 1: rr[i] ^= rr[rk]
            rk += 1
        res[(rk, e, f)] = res.get((rk, e, f), 0) + 1
    return res

mmax = int(sys.argv[1]) if len(sys.argv) > 1 else 11
ch = chain(mmax)
for m in range(2, 7):
    assert ch[m] == brute(m), (m, ch[m], brute(m))
print("chain = brute force for m = 2..6")
for m in range(2, mmax + 1):
    tot = {0: 0, 1: 0}; byr = {0: {}, 1: {}}
    for (r2, e, f), c in ch[m].items():
        tot[e] += c; byr[e][r2] = byr[e].get(r2, 0) + c
    assert tot[0] + tot[1] == 1 << (m * (m - 1) // 2)
    line = "m = %2d odd-edge:" % m
    for r2 in sorted(byr[1], reverse=True):
        line += "  rank %d: %.7f" % (r2, byr[1][r2] / tot[1])
    print(line)
    line = "        uncond.: "
    allr = {}
    for e in (0, 1):
        for r2, c in byr[e].items(): allr[r2] = allr.get(r2, 0) + c
    for r2 in sorted(allr, reverse=True):
        line += "  rank %d: %.7f" % (r2, allr[r2] / (tot[0] + tot[1]))
    print(line)
    # E[2^(-rank/2)] and E[2^-rank] under odd-edge conditioning
    e1 = sum(Fraction(c, tot[1]) * Fraction(1, 2 ** (r2 // 2)) for r2, c in byr[1].items())
    print("        odd-edge E[2^(-rank/2)] = %.6f ; 2^(m/2) * that = %.4f" % (float(e1), float(e1) * 2 ** (m / 2)))
