#!/usr/bin/env python3
"""
lrc14_netmax_search_kps-S6-wf.py
Search for the E (0 in E, |E|=k) that MAXIMIZES meas(N(E)) = 1 - mu_{1/7}(E),
i.e. the worst case for the spread bound.  EXACT engine for verification;
grid for the descent.  We want to confirm max meas(N) <= cap_k with margin, and
identify the maximizer (expected: consecutive).
"""
import sys, itertools, random
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
random.seed(424242)
ONE7=F(1,7)

def net_intervals(E):
    E=sorted(set(E)); n=len(E)
    bp=set([F(0),F(1)])
    for i in range(n):
        for j in range(i+1,n):
            d=E[j]-E[i]
            for m in range(0,d+1): bp.add(F(m,d))
    bp=sorted(b for b in bp if 0<=b<=1)
    good=[]
    for a,b in zip(bp,bp[1:]):
        if b<=a: continue
        mid=(a+b)/2
        order=sorted(range(n),key=lambda i:(E[i]*mid)%1)
        floors=[(E[order[t]]*mid).__floor__() for t in range(n)]
        lo=a; hi=b; feasible=True
        for t in range(n):
            o1=order[t]; o2=order[(t+1)%n]; wrap=1 if t==n-1 else 0
            s=E[o2]-E[o1]; c=F(floors[t]-floors[(t+1)%n]+wrap)
            if s==0:
                if not (c<=ONE7): feasible=False; break
            elif s>0: hi=min(hi,(ONE7-c)/s)
            else: lo=max(lo,(ONE7-c)/s)
            if lo>=hi: feasible=False; break
        if feasible and lo<hi: good.append((lo,hi))
    good.sort(); out=[]
    for a,b in good:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out
def meas(iv): return sum((b-a for a,b in iv),F(0))
def netmeas(E): return meas(net_intervals(E))

NX=30000
XS=[(i+0.5)/NX for i in range(NX)]
def netmeas_grid(E):
    Ef=[float(e) for e in E]; th=1.0/7.0; cnt=0
    for x in XS:
        pts=sorted((e*x)%1.0 for e in Ef); mg=0.0; prev=pts[0]
        for t in range(1,len(pts)):
            d=pts[t]-prev
            if d>mg: mg=d
            prev=pts[t]
        d=pts[0]+1.0-pts[-1]
        if d>mg: mg=d
        if mg<=th: cnt+=1   # net = all gaps <= 1/7
    return cnt/NX

def ascend(k,spread,restarts,iters):
    """maximize netmeas (= minimize mu)."""
    best=-1.0; bestE=None
    for r in range(restarts):
        E=sorted(set([0]+random.sample(range(1,spread+1),k-1)))
        if len(E)<k: continue
        cur=netmeas_grid(E)
        for it in range(iters):
            idx=random.randrange(1,len(E)); old=E[idx]
            step=random.choice([-2,-1,1,2,-4,4,-7,7,-1,1])
            cand=old+step
            if cand<=0 or cand>spread or cand in E: continue
            E2=sorted(E[:idx]+[cand]+E[idx+1:]); v=netmeas_grid(E2)
            if v>=cur: E=E2; cur=v
        if cur>best: best=cur; bestE=E[:]
    return best,bestE

cap={8:F(2243,5880),9:F(1979,4004),10:F(55,91),11:F(66,91),12:F(6,7)}
if __name__=="__main__":
    print("MAXIMIZE meas(N(E)) = worst case for spread bound. Compare to cap_k.")
    for k in range(8,13):
        capk=cap[k]
        gbest=-1.0; gbestE=None
        spreads=[k, k+1, k+2, k+4, k+7, 2*k, 3*k]
        for sp in spreads:
            bv,bE=ascend(k,sp,restarts=25,iters=200)
            if bE and bv>gbest: gbest=bv; gbestE=bE
        # exact-verify the best found and the consecutive
        exC=netmeas(list(range(k)))
        exB=netmeas(gbestE) if gbestE else None
        print(f"\nk={k}: cap_k={float(capk):.5f}")
        print(f"   consec net (EXACT) = {exC} = {float(exC):.5f}")
        print(f"   descent best grid~{gbest:.5f}  E={gbestE}  EXACT={float(exB):.5f}({exB})")
        worst=max(exC,exB)
        print(f"   => worst net found = {float(worst):.5f}  vs cap_k={float(capk):.5f}  "
              f"{'OK (under cap, margin %.4f)'%float(capk-worst) if worst<=capk else '*** OVER CAP ***'}")
