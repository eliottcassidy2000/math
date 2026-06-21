"""Fast flushed independent checks for OPEN-Q-108 synthesis critic. EXACT."""
import sys
from fractions import Fraction as F
import random
from math import gcd
from functools import reduce

def bp(E):
    s=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for j in range(7):
            m=0
            while True:
                xv=(F(j,7)+m)/e
                if xv>=1: break
                if xv>=0: s.add(xv)
                m+=1
    return sorted(b for b in s if 0<=b<1)

def p0p1(E):
    E=sorted(set(E)); B=bp(E); a=F(0);b=F(0)
    for lo,hi in zip(B,B[1:]+[F(1)]):
        if hi<=lo: continue
        mid=(lo+hi)/2; miss=set(range(1,7))-set(int((e*mid)%1*7) for e in E)
        if len(miss)==0:a+=hi-lo
        elif len(miss)==1:b+=hi-lo
    return a,b

def Vcount(E):
    E=sorted(set(E)); B=bp(E)
    comp=0; prev=False
    for lo,hi in zip(B,B[1:]+[F(1)]):
        if hi<=lo: continue
        mid=(lo+hi)/2; miss=len(set(range(1,7))-set(int((e*mid)%1*7) for e in E))
        cur=(miss==1)
        if cur and not prev: comp+=1
        prev=cur
    return comp

caps = {8:F(2243,5880),9:F(1979,4004),10:F(55,91),11:F(66,91),12:F(6,7)}
def P(*a): print(*a); sys.stdout.flush()

random.seed(7)
P("=== (1) chain identity exact (telescoping) ===")
fails=0
for _ in range(60):
    core=sorted(set([0])|set(random.sample(range(1,12),3)))
    fars=sorted(random.sample(range(15,120), 9-len(core)))
    Es=[core]; cur=list(core)
    for w in fars: cur=sorted(cur+[w]); Es.append(cur)
    p0E,_=p0p1(Es[-1]); acc,_=p0p1(core)
    for i in range(1,len(Es)):
        p0i,_=p0p1(Es[i]); _,p1im1=p0p1(Es[i-1])
        acc += F(1,7)*p1im1 + (p0i - p0p1(Es[i-1])[0] - F(1,7)*p1im1)
    fails += (acc!=p0E)
P(f"  chain fails: {fails}/60")

P("=== (2) V-growth: V<=C*k FALSE, V<=7*sigma rigorous ===")
maxVnz=F(0); maxVsig=F(0); v7=0
for _ in range(200):
    k=random.randint(8,12); E=set([0])
    while len(E)<k: E.add(random.randint(1,80))
    E=sorted(E); V=Vcount(E); sig=sum(E); nz=len([e for e in E if e>0])
    maxVnz=max(maxVnz,F(V,nz)); maxVsig=max(maxVsig,F(V,sig) if sig else F(0))
    if V>7*sig: v7+=1
P(f"  max V/#nonzero={float(maxVnz):.3f} (large => V<=C*k FALSE)")
P(f"  max V/sigma={float(maxVsig):.3f}")
P(f"  V>7*sigma violations: {v7}/200")

P("=== (3) k=9 separated two-far cutoff (core=consec7, w1, w2=2w1) ===")
for w1 in [40,45,48,49,50,55,60]:
    E=sorted(list(range(0,7))+[w1,2*w1])
    p0,_=p0p1(E)
    P(f"  w1={w1} p0={float(p0):.5f} cap9={float(caps[9]):.5f} {'OK' if p0<=caps[9] else 'VIOL'}")

P("=== (4) targeted wide hunt (resonant/AP/2-block), p0>cap? scale<=140 ===")
viol=0; tot=0; thin=F(1)
for _ in range(1200):
    k=random.randint(8,12); st=random.randint(0,3); E=set([0])
    if st==0:
        for i in range(1,k-1): E.add(i)
        E.add(random.randint(15,140))
    elif st==1:
        step=random.randint(2,7)
        for i in range(1,k): E.add(i*step)
    elif st==2:
        b=random.randint(20,90)
        for i in range(k//2): E.add(i+1)
        for i in range(k-k//2): E.add(b+i)
    else:
        while len(E)<k: E.add(random.randint(1,140))
    E=sorted(E)
    if len(E)<k or max(E)-min(E)<=14: continue
    if reduce(gcd,[e for e in E if e>0])!=1: continue
    p0,_=p0p1(E); tot+=1; m=caps[k]-p0
    if m<thin: thin=m
    if p0>caps[k]:
        viol+=1
        if viol<=5: P(f"  VIOL k={k} {E} p0={float(p0):.4f}")
P(f"  hunt: {viol}/{tot} violations, thinnest margin={float(thin):.5f}")
P("DONE")
