#!/usr/bin/env python3
"""packet_density_kps_S128c35.py -- density of BONF5>0 among random 13-packets in [300,4000],
exact rational BONF5 per packet, stratified by blocker content (3APs, Sidon violations, small ratios)."""
import sys, random
from fractions import Fraction as F
from math import comb, gcd
sys.stdout.reconfigure(line_buffering=True)
random.seed(20260716)
def bonf5(speeds):
    ev=[]
    for v in speeds:
        for j in range(v):
            lo=F(14*j-1,14*v); hi=F(14*j+1,14*v)
            if lo<0:
                ev.append((F(0),1)); ev.append((hi,-1)); ev.append((lo+1,1)); ev.append((F(1),-1))
            else:
                ev.append((lo,1)); ev.append((hi,-1))
    ev.sort()
    mu=[F(0)]*14; depth=0; last=F(0)
    for x,d in ev:
        if x>last: mu[depth]+=x-last; last=x
        depth+=d
    if F(1)>last: mu[depth]+=F(1)-last
    b5=mu[0]-sum(comb(d-1,5)*mu[d] for d in range(6,14))
    return b5,mu
def blockers(xs):
    ap3=0; sid=0; rat=0
    ss={}
    for i in range(13):
        for j in range(i+1,13):
            s=xs[i]+xs[j]
            if s in ss: sid+=1
            ss[s]=1
    for i in range(13):
        for j in range(13):
            if i<j:
                for k in range(13):
                    if k!=i and k!=j and xs[i]+xs[j]==2*xs[k]: ap3+=1
    for i in range(13):
        for j in range(13):
            if i!=j:
                g=gcd(xs[i],xs[j])
                if xs[i]//g<=13 and xs[j]//g<=13: rat+=1
    return ap3,sid,rat//2
N=200
pos=0; results=[]
for trial in range(N):
    xs=sorted(random.sample(range(300,4001),13))
    b5,mu=bonf5(xs)
    a3,sd,rt=blockers(xs)
    ok=b5>0
    pos+=ok
    results.append((float(b5),a3,sd,rt))
    if trial%20==19: print("  %d/%d: BONF5>0 rate %.3f"%(trial+1,N,pos/(trial+1)),flush=True)
print("DENSITY: %d/%d = %.3f random 13-packets in [300,4000] have BONF5 > 0"%(pos,N,pos/N))
# stratify
import statistics
clean=[r for r in results if r[1]==0 and r[2]==0 and r[3]==0]
dirty=[r for r in results if not(r[1]==0 and r[2]==0 and r[3]==0)]
print("fully-clean (no 3AP, Sidon, no ratio<=13): %d packets, BONF5>0 rate %.3f, mean %.4f"%(
    len(clean),sum(1 for r in clean if r[0]>0)/max(1,len(clean)),statistics.mean(r[0] for r in clean) if clean else 0))
print("with blockers: %d packets, BONF5>0 rate %.3f, mean %.4f"%(
    len(dirty),sum(1 for r in dirty if r[0]>0)/max(1,len(dirty)),statistics.mean(r[0] for r in dirty) if dirty else 0))
for name,idx in [("3AP",1),("Sidon-viol",2),("ratio<=13",3)]:
    have=[r for r in results if r[idx]>0]
    if have:
        print("  packets with %s>0: %d, BONF5>0 rate %.3f, mean BONF5 %.4f"%(name,len(have),sum(1 for r in have if r[0]>0)/len(have),statistics.mean(r[0] for r in have)))
print("DONE")
