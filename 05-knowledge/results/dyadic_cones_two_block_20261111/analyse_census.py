#!/usr/bin/env python3
"""Summaries of the root censuses.
usage: analyse_census.py cones FILE [FILE ...]     (outputs of cone_census; merged, least j kept)
       analyse_census.py blocks FILE [FILE ...]    (outputs of block_census; merged)"""
import sys, collections
def load_cones(files):
    rs = {}
    for f in files:
        for line in open(f):
            if line.startswith('r='):
                p = line.split(); r = int(p[0][2:]); j = int(p[1][2:])
                if r not in rs or j < rs[r]: rs[r] = j
    return rs
def load_blocks(files):
    T = collections.defaultdict(dict)
    for f in files:
        for line in open(f):
            if line[0] == '#': continue
            c, r, j = map(int, line.split())
            if r not in T[c] or j < T[c][r]: T[c][r] = j
    return T
def complete_to(S, sign, limit):
    r = 1
    while r <= limit and sign * r in S: r += 2
    return r - 2
if __name__ == "__main__":
    mode = sys.argv[1]; files = sys.argv[2:]
    if mode == "cones":
        rs = load_cones(files)
        pos = sorted(r for r in rs if r > 0); neg = sorted(-r for r in rs if r < 0)
        print("cone roots: %d positive offsets, %d negative offsets" % (len(pos), len(neg)))
        print("positive offsets complete (all odd r) up to", complete_to(rs, 1, 19999), "; negative up to", complete_to(rs, -1, 19999))
        for lim in (199, 499, 999, 1999):
            mp = [r for r in range(1, lim + 1, 2) if r not in rs]; mn = [r for r in range(1, lim + 1, 2) if -r not in rs]
            print("odd r <= %d without a cone root: positive %d %s ; negative %d %s" % (lim, len(mp), mp if len(mp) <= 40 else str(mp[:40]) + '...', len(mn), mn if len(mn) <= 40 else str(mn[:40]) + '...'))
        print("offsets 1 + 2^e (rows of weight three): ", [(e, rs.get(1 + 2**e)) for e in range(1, 15)])
        print("offsets 2^e - 1: ", [(e, rs.get(2**e - 1)) for e in range(1, 15)])
        print("largest exponent j among cone roots with |r| <= 999:", max(j for r, j in rs.items() if abs(r) <= 999))
        print("window check (offset <= response for positive offsets):", all(r <= 2**j for r, j in rs.items() if r > 0))
    else:
        T = load_blocks(files)
        print("types (c, r), c odd <= 255, |r| < 2047: %d" % sum(len(v) for v in T.values()))
        for C in (1, 3, 7, 15, 31, 63, 127, 255):
            R = 10**9; worst = None
            for c in range(1, C + 1, 2):
                k = min(complete_to(T[c], 1, 2045), complete_to(T[c], -1, 2045))
                if k < R: R = k; worst = c
            print("box: every odd c <= %3d and every odd r with |r| <= %4d is a root type (limiting c = %s)" % (C, R, worst))
        for c in (1, 3, 5, 7, 9, 11, 13, 15):
            print("  c=%2d: positive complete to %4d, negative complete to %4d, largest j in |r| <= 63: %s" % (c, complete_to(T[c], 1, 2045), complete_to(T[c], -1, 2045), max((T[c][r] for r in T[c] if abs(r) <= 63), default=None)))
