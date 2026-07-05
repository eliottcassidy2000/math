#!/usr/bin/env python3
"""
mac-mini-2026-07-05-S53 -- HYP-4109: l=4 post-pass.
Reads the 6 parallel workers' outputs (l4_part_*.out), merges stats, and runs
exact M on any CELL lines (sets the C pyramid could not certify).  Usage:
  python3 lrc_l4_postpass_macmini_S53.py <dir-with-l4_part_files>
"""
from fractions import Fraction as F
import sys, glob, re

sys.path.insert(0, '04-computation')
from lonely_profile import profile

BETA = F(2, 25)

def M_exact(S):
    for cap in (11, 8, 6, 4, 3, 2):
        p = profile(sorted(S), F(1, cap))
        m = p.M()
        if m is not None:
            return m
    return None

d = sys.argv[1] if len(sys.argv) > 1 else '.'
tot = {}
cells = []
for f in sorted(glob.glob(d + '/l4_part_*.out')):
    for line in open(f):
        if line.startswith('CELL'):
            cells.append([int(x) for x in line.split()[1:]])
        m = re.match(r'l=4 \[(\d+),(\d+)\) done: (.*)', line)
        if m:
            for kv in m.group(3).split():
                k, v = kv.split('=')
                tot[k] = tot.get(k, 0) + int(v)
            print(f"  {f.split('/')[-1]}: {m.group(3)}")
print(f"\nl=4 MERGED: {tot}")
print(f"unresolved CELL sets: {len(cells)}")
sub = []
for W in cells:
    if max(W) > 3000:
        print(f"  BIG (needs witness extension, not profile): {W}")
        continue
    M = M_exact(W)
    tag = " << SUB-2/25" if M < BETA else ""
    print(f"  M = {M}  W = {W}{tag}")
    if M < BETA:
        sub.append((M, W))
print(f"\nl=4 VERDICT: {'FLOOR >= 2/25 on the FULL chain domain' if not sub and not any(max(W) > 3000 for W in cells) else 'see above'}")
