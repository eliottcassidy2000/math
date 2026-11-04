#!/usr/bin/env python3
"""Partial check of HYP-9462's 24-hole list for the product closure of UH path counts (input: uh_classes.txt from
./ta uh 11, i.e. UH classes of orders 4..11 only; orders 12-14 have H >= 375 by Busch 2006)."""
vals = {}
for l in open('uh_classes.txt'):
    n, H, aut, r, k, code = l.split()
    vals.setdefault(int(n), set()).add(int(H))
for n in sorted(vals):
    print("UH order %d: %d distinct H values, min %d max %d" % (n, len(vals[n]), min(vals[n]), max(vals[n])))
allv = sorted(set().union(*vals.values()))
LIM = 1125
reach = {1}
ch = True
while ch:
    ch = False
    for x in list(reach):
        for v in allv:
            if v > 1 and x * v <= LIM and x * v not in reach:
                reach.add(x * v); ch = True
holes = sorted(set(range(1, LIM + 1, 2)) - reach)
claimed = [3, 7, 13, 21, 41, 43, 67, 77, 83, 163, 173, 231, 233, 247, 251, 277, 299, 421, 449, 503, 533, 569, 637, 821]
print("holes of the UH-product closure using UH orders 4..11 only:", len(holes))
print("claimed 24 holes (UH orders 4..14) that my closure COVERS (would be a contradiction):", sorted(set(claimed) - set(holes)))
extra = sorted(set(holes) - set(claimed))
print("my extra holes (must be covered by UH orders 12..14 with H in [375,1125]):", len(extra), extra[:80])
print("smallest extra hole:", extra[0] if extra else None)
