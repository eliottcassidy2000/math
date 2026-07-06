#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
The compressed extremizer is the consecutive comb: M(S) >= v_min/(v_min+v_max), witness t*=1/(v_min+v_max),
SATURATED exactly by full consecutive combs {v_min,...,v_max}. opus-2026-07-05-S85.

Honest scope: the LOWER BOUND and the ratio-13 threshold are THM-526 (my rational witness t*=1/(v_min+v_max)
lies in THM-526's safe interval J=(1/(14 v_min),13/(14 v_max))). The SHARPENING is the EXACT value + the
extremizer: among v_max<=13 v_min families the min of M is 1/14, attained UNIQUELY by the AP {1..13}; and the
midrange witness gives the base's EXACT slack v_min/(v_min+v_max)-1/14 for THM-608 peeling.

Also: the single-killer denominator is Q=w+1 (killer w == -1 mod (w+1) acts as a reflected unit runner),
so t*=a/(w+1); this is why the deep well's Eisenstein denominator 183=Phi6(14) equals killer+1=182+1.
"""
import sys
from fractions import Fraction as Fr
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
def norm(x):
    x=x-int(x)
    if x<0:x+=1
    return min(x,1-x)
def exact_M(S):
    S=sorted(set(S));cands=set()
    for v in S:
        for k in range(v):cands.add(Fr(2*k+1,2*v))
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            for den in (S[i]+S[j],abs(S[i]-S[j])):
                if den:
                    for s in range(den):cands.add(Fr(s,den))
    b=Fr(0);arg=None
    for t in cands:
        v=min(norm(x*t) for x in S)
        if v>b:b=v;arg=t
    return b,arg

print("(A) EXACT VALUE for full consecutive combs C={a..a+r-1}: M = a/(2a+r-1) = v_min/(v_min+v_max)")
allsat=True
for (a,r) in [(1,13),(2,13),(5,13),(1,12),(2,12),(3,5),(2,3),(7,10),(4,9)]:
    C=list(range(a,a+r)); M,t=exact_M(C); b=Fr(min(C),min(C)+max(C))
    sat=(M==b); allsat=allsat and sat
    print("   {%2d..%2d} (r=%2d): M=%-8s v_min/(v_min+v_max)=%-8s witness t*=%-8s saturates=%s"%(a,a+r-1,r,str(M),str(b),str(t),sat))
print("   all full consecutive combs saturate:",allsat)
print()
print("(B) EXTREMIZER: among v_max<=13 v_min families, min M = 1/14 UNIQUELY at the AP {1..13}")
print("    (the only full comb hitting the ratio-13 boundary with v_min=1). Sparse/non-full are strictly above:")
for S,tag in [([1,3,5,7,9,11,13],"odd comb (sparse)"),([1,2,4,5,6,7,8,9,10,11,12,13,14],"missing 3"),(list(range(1,14)),"AP {1..13}")]:
    M,_=exact_M(S); b=Fr(min(S),min(S)+max(S))
    print("   %-22s M=%-8s bound=%-8s  gap above bound=%s"%(tag,str(M),str(b),str(M-b)))
print()
print("(C) EXACT SLACK for THM-608 peeling: a compressed base B has delta = M(B)-1/14 >= v_min/(v_min+v_max)-1/14.")
for B in [[1,2,3,4,5,6,7,8,9,10,11,12],[2,3,4,5,6,7,8,9,10,11,12,13],[1,3,5,7,9,11,13]]:
    M,_=exact_M(B); b=Fr(min(B),min(B)+max(B))
    print("   B=%-34s M(B)=%-7s bound=%-7s slack delta=%s"%(str(B),str(M),str(b),str(M-Fr(1,14))))
print()
print("(D) SINGLE-KILLER denominator Q=w+1 (killer==-1 mod Q); explains deep-well 183=killer+1=Phi6(14):")
for w in [13,26,182]:
    S=list(range(1,13))+[w]; M,t=exact_M(S)
    print("   {1..12,%d}: M=%-8s t*=%s Q=%d  = w+1=%d ? %s  (w+1=Phi6(14)? %s)"%(w,str(M),str(t),t.denominator,w+1,t.denominator==w+1, w+1==14*14-14+1))
print("DONE.")
