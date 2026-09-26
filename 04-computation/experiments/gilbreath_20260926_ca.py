#!/usr/bin/env python3
"""gilbreath_20260926_ca.py -- Gilbreath's difference triangle as a cellular automaton, measured
(session collatz-oscillation-20260926, opus, 2026-09-26).

Row 0 = the primes 2, 3, 5, 7, ...; row r+1 = absolute differences of consecutive entries of row r.
Gilbreath's conjecture: every row starts with 1. On the sublattice {0, 2} the rule |a - b| is XOR (in
units of 2), i.e. Rule 90 = Pascal's triangle mod 2, and |1 - a| = 1 for a in {0, 2}: so once a row is
(1, then only 0s and 2s) the leading 1 persists forever. The only threat is a 'defect' (an entry >= 4)
travelling leftward one cell per row through the 0/2 sea, shrinking by 2 each time it meets a 2.
This script computes the triangle for the primes below P, and records per row: the leading entry, the
frontier F(r) = position of the first entry >= 4 (the 0/2 prefix length), the largest entry, and the
distribution of the sizes of maximal inverted zero triangles inside the 0/2 region (Rule 90 predicts
sizes 2^k - 1; sizes 3 and 7 are the cyclic triangle and the Paley heptagon of the repo's tournament
thread, and 2^k - 1 prime = Mersenne prime = a Paley size). It also measures the defect survival: for each
defect of size 2j at distance F, the number of 2s on its leftward diagonal path, to compare with the
random model (a defect of size 2j survives distance L iff it meets fewer than j twos: probability
sum_(i<j) C(L,i) 2^-L, the same one-sided binomial tail as the Collatz no-descent count but with the
threshold j/L -> 0 instead of log_3 2).
Usage: python3 gilbreath_20260926_ca.py [P=200000]
"""
import sys
import numpy as np


def primes_below(P):
    s = np.ones(P, dtype=bool); s[:2] = False
    for i in range(2, int(P ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0].astype(np.int64)


def main():
    P = int(sys.argv[1]) if len(sys.argv) > 1 else 200000
    row = primes_below(P)
    n0 = len(row)
    print("primes below %d: %d; rows computed: %d" % (P, n0, n0 - 1))
    lead_ok = True; frontier = []; maxes = []
    KEEP = 160          # keep the first KEEP cells of every row for the triangle analysis
    diagram = []
    first_bad = None
    r = 0
    while len(row) > 1:
        row = np.abs(np.diff(row))
        r += 1
        if row[0] != 1:
            lead_ok = False
            if first_bad is None:
                first_bad = r
        big = np.nonzero(row >= 4)[0]
        F = int(big[0]) if len(big) else len(row)
        frontier.append(F); maxes.append(int(row.max()))
        diagram.append(row[:KEEP].copy())
    print("leading entry == 1 in every row: %s%s" % (lead_ok, "" if lead_ok else " (first failure at row %d)" % first_bad))
    fr = np.array(frontier)
    for r in (1, 2, 5, 10, 50, 100, 500, 1000, 5000, 10000, min(len(fr), 17000)):
        if r <= len(fr):
            print("   row %6d: frontier F = %7d  (first entry >= 4), max entry = %d" % (r, fr[r - 1], maxes[r - 1]))
    print("   min frontier over rows 100..%d: %d at row %d; the leading 1 is protected for F more rows by Odlyzko's argument" % (len(fr), fr[99:].min(), 100 + int(np.argmin(fr[99:]))))
    # zero-triangle sizes in the 0/2 region: an inverted triangle of side t at (row r, col c) means rows r..r+t-1
    # have zeros at columns c..c+t-1-(row offset) ... we detect maximal runs of zeros and check that a run of length
    # t at row r is followed by runs t-1, t-2, ..., 1 in the rows below (the Rule-90 triangle), then histogram t.
    D = np.array([np.pad(d, (0, KEEP - len(d)), constant_values=-1) for d in diagram[:6000]])
    sizes = {}
    R, C = D.shape
    for r0 in range(1, R - 1):
        rowv = D[r0]
        c = 1
        while c < C:
            if rowv[c] == 0:
                c1 = c
                while c1 < C and rowv[c1] == 0:
                    c1 += 1
                t = c1 - c
                # maximal run bordered by 2s (or the leading 1)?
                left_ok = rowv[c - 1] in (2, 1); right_ok = c1 < C and rowv[c1] == 2
                # triangle check: rows below have zero runs at [c .. c1-1-s] for s = 1..t-1
                tri = left_ok and right_ok
                for s in range(1, t):
                    if r0 + s >= R or not np.all(D[r0 + s][c:c1 - s] == 0) or D[r0 + s][c1 - s] != 2:
                        tri = False; break
                if tri and t >= 1:
                    sizes[t] = sizes.get(t, 0) + 1
                c = c1
            else:
                c += 1
    ks = sorted(sizes)
    print("maximal inverted zero triangles (bordered by 2s) in rows 1..%d, columns < %d: size -> count" % (min(R, 6000), KEEP))
    print("   " + "  ".join("%d:%d" % (t, sizes[t]) for t in ks[:20]))
    print("   sizes that occur: %s; of the form 2^k - 1: %s" % (ks[:20], [t for t in ks if (t + 1) & t == 0]))
    # defect survival statistics: for rows r with a defect 2j at the frontier F(r), how far does it travel?
    surv = []
    for r in range(1, min(len(diagram), 3000)):
        d = diagram[r - 1]
        big = np.nonzero(d >= 4)[0]
        if len(big) == 0 or big[0] >= KEEP - 1:
            continue
        F = int(big[0]); size = int(d[F]); j = size // 2
        # follow the leftward diagonal: row r+s, column F-s
        twos = 0; s = 0; alive = True
        while s < F and r - 1 + s + 1 < len(diagram):
            nxt = diagram[r - 1 + s + 1]
            c = F - s - 1
            if c < 0 or c >= len(nxt):
                break
            if nxt[c] < 4:
                alive = False; break
            s += 1
        surv.append((size, F, s))
    if surv:
        import collections
        bysize = collections.defaultdict(list)
        for size, F, s in surv:
            bysize[size].append((F, s))
        print("defects at the frontier (rows < 3000, columns < %d): size -> (count, mean start distance F, mean travel before dropping below 4, max travel)" % KEEP)
        for size in sorted(bysize)[:8]:
            L = bysize[size]
            print("   %2d: (%4d, %6.1f, %5.2f, %3d)" % (size, len(L), np.mean([f for f, _ in L]), np.mean([t for _, t in L]), max(t for _, t in L)))


if __name__ == '__main__':
    main()
