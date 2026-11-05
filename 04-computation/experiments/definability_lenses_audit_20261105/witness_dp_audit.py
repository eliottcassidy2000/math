#!/usr/bin/env python3
"""Independent subset-DP check of the minimal-arc witnesses quoted in THM-4884 (3) / note section 5.
Chain 0..m-1 transitive (i -> j for i < j) with the listed arcs reversed (j -> i).  Reports H, the shared heads
(vertices that are the head of >= 2 backward arcs, i.e. beaten by two later vertices) and shared tails (tails of
>= 2 backward arcs, i.e. beating two earlier vertices), the cutwidth, and the double-event test (mode-1 violation)."""
from collections import Counter

def H(m, arcs):
    adj = [0] * m
    for i in range(m):
        for j in range(i + 1, m):
            adj[i] |= 1 << j
    for (t, h) in arcs:  # t -> h, t > h
        adj[h] &= ~(1 << t); adj[t] |= 1 << h
    full = (1 << m) - 1
    dp = [[0] * m for _ in range(1 << m)]
    for v in range(m): dp[1 << v][v] = 1
    for S in range(1, 1 << m):
        for v in range(m):
            c = dp[S][v]
            if not c: continue
            nb = adj[v] & ~S
            while nb:
                w = (nb & -nb).bit_length() - 1; nb &= nb - 1
                dp[S | (1 << w)][w] += c
    return sum(dp[full][v] for v in range(m))

claims = {  # value: (m, arcs as (tail, head)) transcribed from the note, section 5
    43: (7, [(2, 0), (4, 2), (6, 2), (6, 4)]),
    47: (7, [(6, 0), (6, 1)]),
    63: (8, [(3, 0), (6, 3), (7, 3)]),
    79: (8, [(7, 0), (7, 2)]),
    103: (9, [(3, 0), (6, 3), (8, 3)]),
    113: (9, [(2, 0), (6, 2), (6, 3), (8, 6)]),
    235: (10, [(6, 0), (6, 4), (8, 6), (9, 6)]),
}
for val, (m, arcs) in claims.items():
    h = H(m, arcs)
    heads = Counter(hd for (_, hd) in arcs); tails = Counter(t for (t, _) in arcs)
    sharedH = sorted(v for v, c in heads.items() if c >= 2); sharedT = sorted(v for v, c in tails.items() if c >= 2)
    cw = max(sum(1 for (t, hd) in arcs if hd <= k < t) for k in range(m - 1))
    print(f"claimed {val}: m={m} arcs={arcs} -> H={h} {'OK' if h == val else 'MISMATCH'}; shared heads {sharedH}, shared tails {sharedT}; cutwidth {cw}; double event present: {bool(sharedH or sharedT)}")

# cross-check against the S79 decoded file
import re, sys
path = r"C:/Users/Eliott/AppData/Local/Temp/s79/minfa_decoded.txt"
try:
    for line in open(path):
        mm = re.match(r"\s*(\d+): m=(\d+) arcs=(\d+) spans=\[[^\]]*\] arcs=\[(.*)\]", line)
        if not mm: continue
        val = int(mm.group(1))
        if val in claims:
            arcs = [tuple(int(x) for x in p.split(',')) for p in re.findall(r"\((\d+, \d+)\)", mm.group(4))]
            print(f"S79 file: {val}: m={mm.group(2)} arcs={arcs}  matches note: {arcs == claims[val][1] and int(mm.group(2)) == claims[val][0]}")
except FileNotFoundError:
    print("S79 file not found")
