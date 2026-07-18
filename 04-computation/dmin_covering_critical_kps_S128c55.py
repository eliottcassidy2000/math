#!/usr/bin/env python3
"""dmin_covering_critical_kps_S128c55.py -- the CRITICAL cases: COVERING families in the
d_min <= 5 stratum.  Covering forces 13|some speed and 14|some speed; with core <= 12 those
must be killers, so the near-equal killer pair is essentially (13a, 14b) with |14b-13a| <= 5.
Construct these systematically and compute exact M.  Data only."""
import sys
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations
sys.stdout.reconfigure(line_buffering=True)
def nd(x):
    fx=x-int(x); fx=fx+1 if fx<0 else fx; return min(fx,1-fx)
def M_exact(V):
    cand=set()
    for v in V:
        for a in range(1,2*v): cand.add(F(a,2*v))
    for i in range(len(V)):
        for j in range(i+1,len(V)):
            for s in (V[i]+V[j],abs(V[i]-V[j])):
                if s:
                    for a in range(1,s): cand.add(F(a,s))
    b=F(0)
    for t in cand:
        if 0<t<1:
            m=min(nd(v*t) for v in V)
            if m>b: b=m
    return b
def covering(V): return all(any(v%q==0 for v in V) for q in range(2,15))
# find near-equal killer pairs (13a, 14b), |14b-13a| <= 5
pairs=[]
for a in range(10,30):
    for b in range(10,30):
        x,y=13*a,14*b
        if x!=y and abs(y-x)<=5 and min(x,y)>156:
            pairs.append((min(x,y),max(x,y),abs(y-x)))
pairs=sorted(set(pairs))[:8]
print("near-equal killer pairs (13a,14b) with gap <= 5 and > 156:")
for p in pairs: print("   %s (d=%d)"%(str(p[:2]),p[2]))
print()
print("COVERING families = core (11 of {1..12}) + such a pair; exact M:")
print("  %-12s %-3s %-42s %-11s %-8s %-6s"%("killers","d","core (11 of 1..12)","M exact","M float","M*14"))
rows=[]
for (k1,k2,d) in pairs:
    for drop in range(1,13):
        core=[x for x in range(1,13) if x!=drop]
        V=sorted(core+[k1,k2])
        if len(V)!=13: continue
        if reduce(gcd,V)!=1: continue
        if not covering(V): continue
        Mx=M_exact(V)
        rows.append((float(Mx),d,drop,tuple(V),Mx))
        print("  %-12s %-3d drop %-2d %-33s %-11s %-8.6f %-6.3f"%(
            "%d,%d"%(k1,k2),d,drop,"{1..12}\{%d}"%drop,Mx,float(Mx),float(Mx)*14))
        break   # one core per pair is enough for the table
print()
if rows:
    rows.sort()
    print("  covering d_min<=5 families found: %d"%len(rows))
    print("  MIN M = %s = %.6f  (%.4f x threshold 1/14)"%(rows[0][4],rows[0][0],rows[0][0]*14))
    print("  MAX M = %s = %.6f  (%.4f x)"%(rows[-1][4],rows[-1][0],rows[-1][0]*14))
    print("  any at or below 1/14: %d"%len([r for r in rows if r[0]<=1/14+1e-12]))
else:
    print("  NO covering family constructed in this shape")
print("DONE")
