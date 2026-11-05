"""Exploratory (after HYP-9471, not pre-registered): plane maps of order <= K that map a slope band [lo, hi] of
t = R/H into itself (Mobius maps are monotone, so the band endpoints decide), their covering density on the band,
and an input file for plane_bfs.c with only the maps preserving [2, 3].
Usage: python3 slope_preserving_subalphabet.py K mats1.txt ... matsK.txt  (writes sel_2_3.txt)"""
import sys
K = int(sys.argv[1]); files = sys.argv[2:2 + K]
mats = set()
for f in files:
    for line in open(f):
        m = tuple(map(int, line.split()))
        if m[3] == m[4] == m[5] == 0:
            A = (m[0], m[2], m[6], m[8])
            if A[0] * A[3] - A[1] * A[2]: mats.add(A)
mats = sorted(mats)
def img(A, t): a, b, c, d = A; return (c + d * t) / (a + b * t)
for lo, hi in ((2.0, 3.0), (2.0, 2.75), (2.2, 3.2), (2.0, 4.0)):
    sel = [A for A in mats if lo <= img(A, lo) <= hi and lo <= img(A, hi) <= hi]
    Dm = [sum((A[0] + A[1] * t) ** -2 for A in sel) for t in (lo, (lo + hi) / 2, hi)]
    print(f'slope band [{lo},{hi}]: maps preserving it {len(sel)} of {len(mats)}; D at lo/mid/hi = {[round(x, 4) for x in Dm]}')
    if (lo, hi) == (2.0, 3.0):
        with open('sel_2_3.txt', 'w') as o:
            for a, b, c, d in sel: o.write(f'{a} 0 {b} 0 0 0 {c} 0 {d}\n')
