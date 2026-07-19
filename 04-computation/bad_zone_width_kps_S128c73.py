#!/usr/bin/env python3
"""bad_zone_width_kps_S128c73.py -- kind-pasteur S128 cont.73.
COMPONENT-AWARE CONDITIONING AT j ~ k1/4.

THM-1144 showed some k1-gaps fail (7*k4*L < 1), clustered near index j ~ k1/4, so no
gap-uniform proof works.  The component-aware repair: a core-safe component I of length ell
covers a j-range of length ell*k1.  If the BAD-j set has width w < ell*k1, then I must
contain a GOOD gap and the four-comb theorem follows.

The 495-core atlas gives ell >= 1/70, so the sufficient condition is

        (bad fraction) := w / k1  <  1/70 = 0.014286 .

Measure the bad fraction exactly.  PRINT DATA ONLY."""
import sys, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
def gap_profile(k1,k2,k3,k4):
    """for each k1-gap index j, the longest surviving piece; exact"""
    out=[]
    for j in range(k1):
        A=F(j,k1)+F(1,14*k1); B=F(j+1,k1)-F(1,14*k1)
        cuts=[]
        for k in (k2,k3,k4):
            i=int(A*k)-1
            while F(i,k)-F(1,14*k) < B:
                x=F(i,k)-F(1,14*k); y=F(i,k)+F(1,14*k)
                if y>A and x<B: cuts.append((max(x,A),min(y,B)))
                i+=1
        cuts.sort(); cur=A; L=F(0)
        for x,y in cuts:
            if x>cur and x-cur>L: L=x-cur
            cur=max(cur,y)
        if B>cur and B-cur>L: L=B-cur
        out.append(L)
    return out
def bad_stats(ks):
    k1,k2,k3,k4=ks
    prof=gap_profile(*ks)
    thr=F(1,7*k4)
    bad=[j for j,L in enumerate(prof) if L<=thr]
    if not bad: return 0,0,0.0,None
    # longest RUN of consecutive bad indices (cyclically)
    runs=[]; cur=[bad[0]]
    for x in bad[1:]:
        if x==cur[-1]+1: cur.append(x)
        else: runs.append(cur); cur=[x]
    runs.append(cur)
    if len(runs)>1 and bad[0]==0 and bad[-1]==k1-1:
        runs[0]=runs[-1]+runs[0]; runs.pop()
    w=max(len(r) for r in runs)
    return len(bad),w,w/k1,max(runs,key=len)[0]
print("### bad-j statistics ; sufficient condition: max-run/k1 < 1/70 = 0.014286 ###")
print("  killers                   #bad   longest run   run/k1     run start   j/k1")
worst=None
CASES=[(157,158,159,160),(197,198,199,200),(237,238,239,240),(277,278,279,280),
       (317,318,319,320),(300,301,302,303),(371,374,377,379),(157,163,169,175),
       (211,212,214,215),(163,167,173,179)]
for ks in CASES:
    nb,w,fr,st=bad_stats(ks)
    print("  %-25s %-6d %-13d %-10.6f %-11s %s"%(str(ks),nb,w,fr,st,
        ("%.3f"%(st/ks[0])) if st is not None else "-"))
    if worst is None or fr>worst[0]: worst=(fr,ks,w,st)
print()
print("  worst bad-fraction among these: %.6f at %s (run %d starting j=%d)"%(
    worst[0],str(worst[1]),worst[2],worst[3]))
print("  sufficient condition needs < 0.014286 -> %s"%("HOLDS" if worst[0]<1/70 else "FAILS"))
print()
print("### systematic: consecutive quadruples, k1 in [157,400] ###")
gw=None
for k1 in range(157,401,3):
    ks=(k1,k1+1,k1+2,k1+3)
    nb,w,fr,st=bad_stats(ks)
    if gw is None or fr>gw[0]: gw=(fr,ks,w,st,nb)
    if k1%60<3: print("  ... k1=%d  running worst bad-fraction %.6f at %s"%(k1,gw[0],str(gw[1])))
print()
print("  WORST bad-fraction over the scan: %.6f"%gw[0])
print("    at %s : %d bad indices, longest run %d starting at j=%d (j/k1 = %.4f)"%(
    str(gw[1]),gw[4],gw[2],gw[3],gw[3]/gw[1][0]))
print("  atlas gives ell >= 1/70 = %.6f"%(1/70))
print("  VERDICT: %s"%("bad run fits inside a component -> conditioning FAILS"
    if gw[0]>=1/70 else "bad run is narrower than every component -> conditioning HOLDS"))
print("DONE")
