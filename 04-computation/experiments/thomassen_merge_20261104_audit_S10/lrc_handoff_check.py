#!/usr/bin/env python3
"""What the LRC radius-7 handoff file actually lists (note section 0, 'Coincidence in the repository').
Run from the repository root:  python3 <this folder>/lrc_handoff_check.py"""
import itertools, re
F = "05-knowledge/results/radius7_multicluster_handoff_opus_S314.txt"
pref = {}
for line in open(F):
    m = re.match(r"^\[([\d, ]+)\] \+ \[(\d+)\]", line)
    if m:
        P = tuple(int(x) for x in m.group(1).split(","))
        pref.setdefault(P, set()).add(int(m.group(2)))
allsub = set(itertools.combinations(range(1, 13), 5))
print("distinct 5-prefixes:", len(pref), " = all 5-subsets of {1..12}:", set(pref) == allsub)
print("prefixes having the lift 105:", sum(1 for P in pref if 105 in pref[P]))
print("lifts of [4,6,8,9,12]:", len(pref[(4, 6, 8, 9, 12)]), sorted(pref[(4, 6, 8, 9, 12)]))
print("residues mod 13 of these lifts:", sorted({x % 13 for x in pref[(4, 6, 8, 9, 12)]}))
print("is every lift of [4,6,8,9,12] = 14 mod 13?", all(x % 13 == 14 % 13 for x in pref[(4, 6, 8, 9, 12)]))
