#!/usr/bin/env python3
"""
lrc14_gaphole_kps-S6-wf.py
Test the SPREAD-CUTOFF lemma that powers the finite reduction:

  LEMMA(?):  if net(E) > 0  then every consecutive gap of sorted E is <= G0.
  (Equivalently: a large 'hole' in E forces net(E)=0.)
  If true with explicit G0, then for |E|=k fixed, net(E)>0 => max(E)<= (k-1)*G0,
  bounding the spread -> finite exhaustion.

We probe: among E with net(E)>0, what is the MAX consecutive gap observed?
Also test variants: max gap as a function of the OTHER gaps; and whether the
relevant cutoff is on max gap or on something subtler.

We scan exhaustively k=8 sets with bounded max(E) and record, for each E with
net>0, the maximum consecutive gap.  We seek the empirical G0 and a counterexample
search to see if net can stay positive with a large hole.
"""
import sys, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
from math import gcd
from functools import reduce
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
def maxgap_E(E):
    E=sorted(E); return max(E[i+1]-E[i] for i in range(len(E)-1))

if __name__=="__main__":
    print("[SPREAD-CUTOFF] among k=8 E with net>0, max consecutive gap of E?")
    # exhaustive over all 8-subsets of {0,...,W} containing 0, primitive
    for W in [10,12,14,16]:
        maxhole=0; maxholeE=None; npos=0
        for body in itertools.combinations(range(1,W+1),7):
            E=(0,)+body
            if E[-1]!=W: continue  # ensure spread exactly W (max element = W)
            if reduce(gcd,E)!=1: continue
            nm=netmeas(list(E))
            if nm>0:
                npos+=1
                h=maxgap_E(E)
                if h>maxhole: maxhole=h; maxholeE=E
        print(f"   W={W:2d}: #(net>0 primitive, max=W)={npos:5d}  max consec-gap among them={maxhole}  E={maxholeE}")

    print("\n[DIRECT hole test] does a single big hole kill net? scan E=run+hole+run:")
    # E = {0..a-1} U {a+H .. a+H+(8-a)-1}: two runs separated by hole H
    for a in range(1,8):
        b=8-a
        for H in range(1,12):
            E=list(range(a))+[a-1+H+i for i in range(b)]
            E=sorted(set(E))
            if len(E)!=8: continue
            if reduce(gcd,E)!=1:
                # skip non-primitive (scaling handles them)
                pass
            nm=netmeas(E)
            tag = 'NET>0' if nm>0 else '.'
            if nm>0:
                print(f"   runs {a}+{b}, hole H={H:2d}: E={E} net={float(nm):.5f} {tag}")
