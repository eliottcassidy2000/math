# Lane architecture (S94).  Necessary condition for a FULL-CONE covering certificate of the port-chain orbit
# (HYP-9640 (B)): sum over the invertible block matrices G of 1/|det G| >= 1  (PROOFS.txt, P3).
# Input: the base-2 port-chain block matrices recomputed by the S93 lane thinorbit (distinct 3x3 matrices of marked
# blocks of order <= k), read-only from the S93 worktree.  Exact integer determinants, exact rational sums.
from fractions import Fraction
import sys
D = "C:/tmp/wt-s93/04-computation/experiments/s93_chanpak_wave_20261111/thinorbit/"
def det3(m):
    a, b, c, d, e, f, g, h, i = m
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
for k in (3, 4, 5, 6):
    mats = set()
    for line in open(D + "blocks2_k%d.txt" % k):
        t = line.split()
        if len(t) >= 9:
            mats.add(tuple(int(x) for x in t[:9]))
    dets = [abs(det3(m)) for m in mats]
    inv = [d for d in dets if d != 0]
    s = sum(Fraction(1, d) for d in inv)
    hist = {}
    for d in inv:
        hist[d] = hist.get(d, 0) + 1
    small = sorted(hist.items())[:12]
    # top-left entry (multiplier of h) statistics, to compare with the 1-dimensional readout
    print("order <= %d: %d distinct matrices, %d invertible, sum 1/|det| = %s = %.6f ; least |det| histogram %s"
          % (k, len(mats), len(inv), s, float(s), small))
print("A full-cone certificate (PROOFS.txt P3) needs sum 1/|det| >= 1 over the blocks it uses.")
