"""
Enumerate the Farey-76 cells whose roof exceeds 1/7, their exact roof-superlevel
intervals, and the minimal subset summing to >= m_P.  Feeds the Lean AP76 certificate
(instantiate FareyRoofBridge.muGood_ge_sum_intervals).  boxeph-2026-07-07-S2.

For AP{1..k}, k=76: on Farey-76 cell (p/q, p'/q') the roof(x)=(q x - p)+(p' - q' x) is
affine, = 1/q at p/q and 1/q' at p'/q'.  {roof>1/7} within the cell is a sub-interval.
Only cells with min(q,q')<=6 contribute (roof node values 1/q).
"""
from fractions import Fraction as F
from math import gcd

k = 76
thr = F(1, 7)
mP = F(14249, 252252)

# Farey_k sequence in [0,1]
def farey(n):
    fr = set()
    for q in range(1, n+1):
        for p in range(0, q+1):
            fr.add(F(p, q))
    return sorted(fr)

fs = farey(k)
print(f"Farey_{k}: {len(fs)} fractions")

# for each consecutive pair (cell), compute roof-superlevel interval and length
cells = []
for a, b in zip(fs, fs[1:]):
    p, q = a.numerator, a.denominator
    p2, q2 = b.numerator, b.denominator
    # determinant check q*p2 - p*q2 == 1
    assert q*p2 - p*q2 == 1, (a, b)
    # roof(x) = (q x - p) + (p2 - q2 x) = (q-q2) x + (p2 - p)
    m = q - q2         # slope
    c = F(p2 - p)      # intercept
    # roof at endpoints
    rL = m*a + c       # = 1/q
    rR = m*b + c       # = 1/q2
    assert rL == F(1, q) and rR == F(1, q2), (rL, rR, q, q2)
    # superlevel {roof>thr} intersect (a,b)
    # roof(x) > thr
    lo, hi = a, b
    if m == 0:
        # roof constant = c ; contributes fully if c>thr
        if c > thr:
            seg = (a, b)
        else:
            seg = None
    elif m > 0:  # increasing: roof>thr for x > x*
        xstar = (thr - c)/m
        loo = max(a, xstar)
        seg = (loo, b) if loo < b else None
    else:  # decreasing
        xstar = (thr - c)/m
        hii = min(b, xstar)
        seg = (a, hii) if a < hii else None
    if seg and seg[1] > seg[0]:
        length = seg[1]-seg[0]
        cells.append({
            'p':p,'q':q,'p2':p2,'q2':q2,'a':a,'b':b,
            'segL':seg[0],'segR':seg[1],'len':length,
            'minq':min(q,q2)
        })

total = sum(cl['len'] for cl in cells)
print(f"contributing cells: {len(cells)}")
print(f"total roof-superlevel measure = {total} = {float(total):.8f}")
print(f"m_P = {mP} = {float(mP):.8f}")
print(f"published AP76 value 2314528732/40290957525 = {float(F(2314528732,40290957525)):.8f}")
print(f"total >= m_P ? {total >= mP}   (margin {float(total-mP):.6f})")
print()
# sort cells by length desc; greedily pick until sum >= mP
cells_sorted = sorted(cells, key=lambda c: -c['len'])
acc = F(0); picked = []
for cl in cells_sorted:
    if acc >= mP: break
    acc += cl['len']; picked.append(cl)
print(f"MINIMAL greedy subset to clear m_P: {len(picked)} intervals, sum={float(acc):.8f} >= m_P ({acc>=mP})")
print(f"  (margin {float(acc-mP):.6f})")
print()
print("The picked intervals (for Lean instantiation):")
for cl in sorted(picked, key=lambda c: c['segL']):
    print(f"  cell (p/q,p'/q')=({cl['p']}/{cl['q']},{cl['p2']}/{cl['q2']})  "
          f"interval ({cl['segL']}, {cl['segR']})  len={cl['len']}={float(cl['len']):.6f}")
print()
print("All contributing cells by node denominator:")
from collections import defaultdict
byq = defaultdict(lambda: [0, F(0)])
for cl in cells:
    byq[cl['minq']][0]+=1; byq[cl['minq']][1]+=cl['len']
for q in sorted(byq):
    print(f"  minq={q}: {byq[q][0]} cells, total len {float(byq[q][1]):.6f}")
