"""
lrc14_cap_minimizer_wide_search_kps.py  (kind-pasteur-2026-06-27-S31ag)
Double-check the canonical cap_8, cap_9 against a WIDER minimizer search.
THM-576 found (over P subset {1..16}): j=4 min {1,11,12,13}=1979/4004=cap_9,
j=5 min {1,5,7,8,9}=2243/5880=cap_8. Verify no smaller meas(lonely(P)) exists
over a larger range (P subset {1..28}, fixing 1 in P by dilation-normalization).
"""
import sys, itertools
from fractions import Fraction as F
from math import comb

def lonely_measure(P, thr=F(1, 14)):
    P = [p for p in P if p != 0]
    if not P: return F(1)
    bps = set([F(0), F(1)])
    for p in P:
        for k in range(p + 1):
            for v in (F(k)/p - thr/p, F(k)/p + thr/p):
                if 0 <= v <= 1: bps.add(v)
    bps = sorted(bps); total = F(0)
    for a, b in zip(bps, bps[1:]):
        if b <= a: continue
        mid = (a + b)/2
        if all(min((p*mid)%1, 1-(p*mid)%1) >= thr for p in P): total += b - a
    return total

if __name__ == "__main__":
    sys.stdout.reconfigure(line_buffering=True)
    C142 = comb(14,2)
    canon = {4: F(1979,4004), 5: F(2243,5880)}
    RANGE = 28
    for j in (4, 5):
        best = None; bestP = None; cnt = 0
        # WLOG 1 in P (primitive/normalized); choose other j-1 from {2..RANGE}
        for rest in itertools.combinations(range(2, RANGE+1), j-1):
            P = (1,) + rest; cnt += 1
            m = lonely_measure(P)
            if best is None or m < best: best = m; bestP = P
        smooth = F(comb(14-j,2), C142)
        print(f"j={j}: searched {cnt} sets (1 in P, rest in 2..{RANGE})")
        print(f"   min meas = {best} = {float(best):.6f}  at P={bestP}")
        print(f"   canonical cap = {canon[j]}  {'MATCH' if best==canon[j] else 'MISMATCH!!'}")
        print(f"   smooth C(14-j,2)/91 = {smooth}, dip = {best-smooth}")
