#!/usr/bin/env python3
"""single_run_proof_kps_S128c74.py -- kind-pasteur S128 cont.74.
PROVING THE SINGLE-RUN BOUND FROM THE LINEAR DESCENT.

THE STRUCTURE.  Inside the k1-gap of index j, the tooth of comb k_i sits at an offset that
drifts LINEARLY in j.  Per unit j the offset moves by -d_i/(k1 k_i) in absolute terms, i.e.
by -d_i/k1 when normalised by k_i's own period 1/k_i.  So for CONSECUTIVE killers
(d_2,d_3,d_4) = (1,2,3) the three normalised offsets are

        u,  2u,  3u      with      u = j/k1      (mod 1).

Hence the whole configuration -- and therefore the longest surviving piece -- is a function
F(u) OF THE SINGLE PARAMETER u.  The point set {0,u,2u,3u} is maximally spread at u = 1/4,
which is exactly where THM-1145 measured the bad zone (j/k1 = 0.238-0.242).
IF F is unimodal on [0,1/2] with its minimum at u = 1/4, then {u : F(u) <= threshold} is an
INTERVAL, i.e. the bad set is a single run -- which is what THM-1145 needs.
TEST: unimodality, the location of the minimum, and the width of the bad interval.
PRINT DATA ONLY."""
import sys
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
def profile(k1):
    """longest surviving piece for each k1-gap index j, consecutive quadruple"""
    k2,k3,k4=k1+1,k1+2,k1+3
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
    return out,k4
print("### (1) is F unimodal in j, with minimum at j ~ k1/4 ? ###")
print("  k1     argmin j   argmin/k1   #local minima   #sign changes of the difference")
for k1 in [157,197,237,277,317,357,397]:
    prof,k4=profile(k1)
    am=min(range(k1),key=lambda j:prof[j])
    # count local minima on the circle
    lm=sum(1 for j in range(k1) if prof[j]<=prof[(j-1)%k1] and prof[j]<=prof[(j+1)%k1])
    sc=sum(1 for j in range(k1) if (prof[(j+1)%k1]-prof[j]>0) != (prof[j]-prof[(j-1)%k1]>0))
    print("  %-6d %-10d %-11.4f %-15d %d"%(k1,am,am/k1,lm,sc))
print()
print("### (2) the bad set: is it ONE run, and how wide? ###")
print("  k1     #bad   #runs   longest run   run/k1     run centre/k1   threshold-normalised")
worst=None
for k1 in range(157,420,10):
    prof,k4=profile(k1)
    thr=F(1,7*k4)
    bad=[j for j in range(k1) if prof[j]<=thr]
    if not bad:
        print("  %-6d %-6d %-7d %-13d %-10s %-15s -"%(k1,0,0,0,"-","-")); continue
    runs=[]; cur=[bad[0]]
    for x in bad[1:]:
        if x==cur[-1]+1: cur.append(x)
        else: runs.append(cur); cur=[x]
    runs.append(cur)
    if len(runs)>1 and bad[0]==0 and bad[-1]==k1-1:
        runs[0]=runs[-1]+runs[0]; runs.pop()
    big=max(runs,key=len)
    fr=len(big)/k1; ctr=(big[0]+big[-1])/2/k1
    if worst is None or fr>worst[0]: worst=(fr,k1,len(runs),len(big))
    if k1%50<10:
        print("  %-6d %-6d %-7d %-13d %-10.6f %-15.4f %s"%(k1,len(bad),len(runs),len(big),fr,ctr,"ok"))
print()
print("  WORST: bad fraction %.6f at k1=%d (%d runs, longest %d)"%(worst[0],worst[1],worst[2],worst[3]))
print("  THM-1145 asserted a single run of width <= 0.0457*k1 : %s"%(
    "CONSISTENT" if worst[0]<=0.0458 else "VIOLATED"))
print()
print("### (3) multi-run check: does the bad set EVER split? ###")
split=[]
for k1 in range(157,420):
    prof,k4=profile(k1)
    thr=F(1,7*k4)
    bad=[j for j in range(k1) if prof[j]<=thr]
    if not bad: continue
    runs=1
    for a,b in zip(bad,bad[1:]):
        if b!=a+1: runs+=1
    if bad[0]==0 and bad[-1]==k1-1 and runs>1: runs-=1
    if runs>1: split.append((k1,runs,len(bad)))
print("  k1 values in [157,420) whose bad set is NOT a single run: %d"%len(split))
if split: print("   ",split[:10])
else: print("    -> the bad set is a SINGLE RUN for every k1 tested")
print("DONE")
