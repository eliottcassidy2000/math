#!/usr/bin/env python3
"""
lrc14_apbase_rigidity_macmini_S59  (mac-mini-2026-07-05-S59)

Verifies the AP-BASE case of klein-S140's loose-branch rigidity (HYP-4151):
  V = (dilated 11-AP c*{1..11}) u {X}, primitive => M(V) in {1/13} U [2/25, inf).
  i.e. NO such family lands in the gap (1/13, 2/25).

PROOF (elementary; this script confirms it):
  The 11-AP c*{1..11} has M=1/12 at the witnesses t = j/(12c), gcd(j,12)=1.
  At such t, min over the 11 runners = 1/12; X is safe iff ||X j/(12c)|| >= 1/12
  iff X*j mod 12c NOT in (-c, c).
  - c=1: X unsafe at ALL witnesses iff X == 0 mod 12 => the ladder {1..11,12k},
    M = k/(12k+1) (= 1/13 at k=1, >= 2/25 for k>=2). X != 0 mod 12 => safe => M >= 1/12.
  - c>=2: take the c witnesses j = 1,13,...,1+12(c-1) (all == 1 mod 12). Their
    residues X*j mod 12c = {e + 12t mod 12c : t=0..c-1} (e = X mod 12c), spaced 12,
    spanning [e, e+12(c-1)]. At most floor(c/6)+1 < c fall in the width-2c danger band
    (-c, c); primitivity gcd(c,X)=1 keeps e != 0. So >= 1 witness is SAFE => M >= 1/12 > 2/25.
Hence M(V) is never in the gap. (Reduces klein's open half to "all 11-subfamilies loose".)
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import sys
def out(*a): print(*a); sys.stdout.flush()
def gcd_all(xs): return reduce(gcd, xs)
def nd(x):
    x = x % 1
    return x if x <= 1-x else 1-x
def M_exact(sp):
    vs = sorted(set(sp)); Q = set()
    for i in range(len(vs)):
        Q.add(2*vs[i])
        for j in range(i+1, len(vs)): Q.add(vs[i]+vs[j]); Q.add(vs[j]-vs[i])
    best = F(0)
    for q in sorted(Q):
        if q < 2: continue
        for a in range(1, q):
            m = min(nd(v*F(a, q)) for v in sp)
            if m > best: best = m
    return best
if __name__ == "__main__":
    c113, c225, c112 = F(1,13), F(2,25), F(1,12)
    ingap = at13 = tested = 0
    for c in range(1, 13):
        base = [c*i for i in range(1, 12)]
        for X in range(1, 12*c*3 + 30):
            if X in base: continue
            S = sorted(base + [X])
            if len(S) != 12 or gcd_all(S) != 1: continue
            tested += 1
            M = M_exact(S)
            if c113 < M < c225:
                ingap += 1
                if ingap <= 8: out(f"  !! IN GAP c={c} X={X}: {S} M={M}")
            if M == c113: at13 += 1
    out(f"tested {tested} primitive dilated-11-AP-base families (c=1..12)")
    out(f"  IN GAP (1/13,2/25): {ingap};  at 1/13: {at13}")
    out(f"  ladder {{1..11,12k}} k=1..6: " + ", ".join(str(M_exact(list(range(1,12))+[12*k])) for k in range(1,7)))
    out(f"  c>=2 free-rider (M=1/12): 2*{{1..11}}+1 -> {M_exact([2*i for i in range(1,12)]+[1])}, "
        f"3*{{1..11}}+2 -> {M_exact([3*i for i in range(1,12)]+[2])}")
    out("  VERDICT: " + ("AP-base gap-emptiness HOLDS (0 in gap)" if ingap == 0 else "FAILS"))
