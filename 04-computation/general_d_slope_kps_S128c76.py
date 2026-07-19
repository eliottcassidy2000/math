#!/usr/bin/env python3
"""general_d_slope_kps_S128c76.py -- kind-pasteur S128 cont.76.
CHECKING THE GENERAL-d SLOPE.
In the continuum model the teeth sit at centres (7/6)*frac(-d_i u) - 1/12, width 1/6, in
[0,1], and BAD means longest surviving piece <= 1/6.
For small u, frac(-d_i u) = 1 - d_i u, so tooth i has LEFT EDGE 1 - (7/6) d_i u, and the
tooth with the LARGEST d_i enters first.  So the linear branch is
        F(u) = 1 - (7/6) d_max u ,      entry at   u = 5/(7 d_max) .
At d_max = 3 that is slope 7/2 and entry 5/21 -- matching THM-1173 exactly.
PREDICTION TO TEST: per-run width 1/(7 d_max), and total bad growing with d_max until it
exceeds the safe measure 0.164.  Compute the exact total bad measure for many (d2,d3,d4).
PRINT DATA ONLY."""
import sys, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
def Finf(u,DS):
    cuts=[]
    for d in DS:
        g=(-d*u)%1
        c=F(7,6)*g-F(1,12)
        a=max(c-F(1,12),F(0)); b=min(c+F(1,12),F(1))
        if b>a: cuts.append((a,b))
    cuts.sort(); cur=F(0); L=F(0)
    for a,b in cuts:
        if a>cur and a-cur>L: L=a-cur
        cur=max(cur,b)
    if 1-cur>L: L=1-cur
    return L
def badstats(DS,N=8400):
    thr=F(1,6)
    bad=[t for t in range(N) if Finf(F(t,N),DS)<=thr]
    if not bad: return 0.0,0,0.0
    runs=[]; cur=[bad[0]]
    for x in bad[1:]:
        if x==cur[-1]+1: cur.append(x)
        else: runs.append(cur); cur=[x]
    runs.append(cur)
    if len(runs)>1 and bad[0]==0 and bad[-1]==N-1:
        runs[0]=runs[-1]+runs[0]; runs.pop()
    return len(bad)/N, len(runs), max(len(r) for r in runs)/N
print("### verify the linear branch slope for general d_max ###")
print("  d-triple        predicted slope (7/6)d_max   entry 5/(7 d_max)   F(entry)   matches 1/6")
for DS in [(1,2,3),(1,2,4),(1,3,5),(2,4,6),(1,2,7),(3,6,9),(1,5,10)]:
    dm=max(DS); ent=F(5,7*dm)
    print("  %-15s %-28.5f %-19s %-10s %s"%(str(DS),7*dm/6,str(ent),str(Finf(ent,DS)),Finf(ent,DS)==F(1,6)))
print()
print("### total bad measure vs d_max ###")
print("  d-triple        d_max   total bad   #runs   per-run   predicted 1/(7 d_max)   vs 0.164")
rows=[]
for DS in [(1,2,3),(1,2,4),(1,3,4),(2,3,4),(1,2,5),(1,3,5),(2,4,5),(1,2,6),(2,4,6),
           (1,2,7),(3,5,7),(1,2,8),(2,5,8),(1,2,9),(3,6,9),(1,2,10),(1,5,10),(2,7,12)]:
    tb,nr,pr=badstats(DS)
    dm=max(DS)
    rows.append((tb,DS,dm,nr,pr))
    print("  %-15s %-7d %-11.6f %-7d %-9.6f %-23.6f %s"%(
        str(DS),dm,tb,nr,pr,1/(7*dm),"OK" if tb<0.164 else "*** EXCEEDS"))
print()
bad=[r for r in rows if r[0]>=0.164]
print("  triples whose total bad exceeds the 0.164 safe measure: %d"%len(bad))
if bad:
    for tb,DS,dm,nr,pr in bad[:8]: print("    %s  d_max=%d  total=%.6f"%(str(DS),dm,tb))
else:
    print("    NONE -- total bad stays below 0.164 on every triple tested")
rows.sort(reverse=True)
print("  largest total bad seen: %.6f at %s (d_max=%d, %d runs)"%(rows[0][0],str(rows[0][1]),rows[0][2],rows[0][3]))
print()
print("### does total bad grow like 2*d_max/21 as I predicted? ###")
print("  d_max   observed max total bad   my prediction 2*d_max/21")
for dm in range(3,13):
    obs=max((r[0] for r in rows if r[2]==dm), default=None)
    if obs is None: continue
    print("  %-7d %-24.6f %.6f"%(dm,obs,2*dm/21))
print("DONE")
