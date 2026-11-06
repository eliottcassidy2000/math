#!/usr/bin/env python3
"""merge_census.py (klein-2026-11-06): merge the parts of a split hspec census (census_<n>_<r>.txt) into one file
census_<n>.txt with the same line types (T h c3 str count / W h witness / S h witness / C h hc count / # summary).
usage: merge_census.py n part1 part2 ... > census_n.txt"""
import sys
from collections import defaultdict
n = int(sys.argv[1]); parts = sys.argv[2:]
T = defaultdict(int); C = defaultdict(int); W = {}; S = {}
total = strong = 0; maxh = 0; maxhc = 0
for p in parts:
    with open(p) as fh:
        for line in fh:
            a = line.split()
            if not a: continue
            if a[0] == 'T': T[(int(a[1]), int(a[2]), int(a[3]))] += int(a[4])
            elif a[0] == 'C': C[(int(a[1]), int(a[2]))] += int(a[3]); maxhc = max(maxhc, int(a[2]))
            elif a[0] == 'W':
                h = int(a[1]); W[h] = min(W.get(h, a[2]), a[2])
            elif a[0] == 'S':
                h = int(a[1]); S[h] = min(S.get(h, a[2]), a[2])
            elif a[0] == '#' and a[1] == 'order':
                total += int(a[3]); strong += int(a[6]); maxh = max(maxh, int(a[11]))
print(f"# order {n}: {total} tournaments read, {strong} strong, max h = {maxh}")
if C: print(f"# max hc = {maxhc}")
for k in sorted(T): print(f"T {k[0]} {k[1]} {k[2]} {T[k]}")
for h in sorted(W): print(f"W {h} {W[h]}")
for h in sorted(S): print(f"S {h} {S[h]}")
for k in sorted(C): print(f"C {k[0]} {k[1]} {C[k]}")
hs = {k[0] for k in T}; hss = {k[0] for k in T if k[2] == 1}
print(f"# distinct h = {len(hs)}, distinct strong h = {len(hss)}")
