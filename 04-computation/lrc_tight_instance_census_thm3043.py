"""THE REAL TRANSFER FROM FC(2): its proof SHAPE is
   (i) exclude non-rigid configurations by a monodromy/inverse argument
   (ii) reducing to a RIGID family (flat top f_D = a(x+y)^D)
   (iii) kill the rigid family by a separate ARITHMETIC argument (transcendence).
LRC's Route B has the same architecture: inverse theorem -> AP-uniqueness -> kill APs.
FC therefore tells us WHERE the difficulty must sit, and (i)/(ii) is the analogue of the
LRC wall.  Test the (ii) half empirically: is the COVERING case exactly the AP family?

Census: all speed multisets of size n from [1,B], compute mu(Safe) exactly at threshold
1/(n+1), and classify the covering ones (mu = 0).
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations

def safe_measure(vs, n1):
    bad = []
    for v in vs:
        w = F(1, n1*v)
        for k in range(0, v+1):
            c = F(k, v)
            lo, hi = max(F(0), c-w), min(F(1), c+w)
            if lo < hi: bad.append((lo, hi))
    bad.sort(); merged = []
    for s, e in bad:
        if merged and s <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else: merged.append((s, e))
    return F(1) - sum(e-s for s, e in merged)

def is_dilated_AP(vs):
    """is vs a dilated AP {d, 2d, ..., nd} (as a SET, after sorting)?"""
    s = sorted(vs); d = s[0]
    return all(s[i] == (i+1)*d for i in range(len(s)))

def normalise(vs):
    g = reduce(gcd, vs)
    return tuple(v//g for v in vs)

for n, B in [(3, 14), (4, 14), (5, 12)]:
    n1 = n+1
    cover, total = [], 0
    for vs in combinations(range(1, B+1), n):
        total += 1
        if safe_measure(vs, n1) == 0:
            cover.append(vs)
    aps = [v for v in cover if is_dilated_AP(normalise(v))]
    non = [v for v in cover if not is_dilated_AP(normalise(v))]
    print(f"n={n}, speeds from [1,{B}], threshold 1/{n1}:  {total} sets, {len(cover)} COVER")
    print(f"   dilated-AP after normalising by gcd: {len(aps)}   NOT AP: {len(non)}")
    if non:
        print(f"   *** NON-AP COVERING SETS (these are the interesting ones) ***")
        for v in non[:14]:
            print(f"       {v}   normalised {normalise(v)}")
    else:
        print(f"   -> every covering set is a dilated AP: AP-uniqueness holds in this box")
    print()
