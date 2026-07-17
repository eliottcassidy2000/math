#!/usr/bin/env python3
"""depth_spectrum_leeyang_kps_S128c35.py -- kind-pasteur S128 cont.35.
THE DEPTH-SPECTRUM REFRAME: p(t) = sum_d mu_d (1-t)^d, mu_d = mu{exactly d bad sets active}.
S_k = sum_d C(d,k) mu_d. Lee-Yang admissibility = real-rootedness of ptilde(s)=sum mu_d s^d.
Exact Fraction sweep per packet; BONF5 exact; root geometry; Newton on mu_d."""
import sys
from fractions import Fraction as F
from math import comb, gcd
sys.stdout.reconfigure(line_buffering=True)
LAM=F(1,14)

def depth_spectrum(speeds):
    ev=[]
    for v in speeds:
        for j in range(v):
            lo=F(14*j-1,14*v); hi=F(14*j+1,14*v)
            if lo<0:
                ev.append((F(0),1)); ev.append((hi,-1))
                ev.append((lo+1,1)); ev.append((F(1),-1))
            else:
                ev.append((lo,1)); ev.append((hi,-1))
    ev.sort()
    n=len(speeds)
    mu=[F(0)]*(n+1)
    depth=0; last=F(0)
    for x,d in ev:
        if x>last: mu[depth]+=x-last; last=x
        depth+=d
    if F(1)>last: mu[depth]+=F(1)-last
    assert depth==0
    assert sum(mu)==1
    return mu

def analyze(name,speeds):
    n=len(speeds)
    mu=depth_spectrum(speeds)
    S=[sum(comb(d,k)*mu[d] for d in range(n+1)) for k in range(n+1)]
    b5=sum((-1)**k*S[k] for k in range(6))
    p1=sum((-1)**k*S[k] for k in range(n+1))
    # roots of ptilde(s) = sum mu_d s^d (nonneg coeffs; real roots must be <= 0; p roots t = 1-s)
    co=[float(m) for m in mu]
    while co and co[-1]==0: co.pop()
    deg=len(co)-1
    ws=[complex(0.4,0.9)**i*1.5 for i in range(deg)]
    for _ in range(3000):
        new=[]
        for i,wi in enumerate(ws):
            num=sum(co[j]*wi**j for j in range(deg+1)); den=co[-1]
            for j,wj in enumerate(ws):
                if j!=i: den*=(wi-wj)
            new.append(wi-num/den if den!=0 else wi)
        ws=new
    nreal=sum(1 for w in ws if abs(w.imag)<1e-7*max(1,abs(w)))
    # Newton on mu (normalized): b_d = mu_d / C(deg,d) on the trimmed support? use raw log-concavity + normalized Newton
    lcfail=sum(1 for d in range(1,n) if mu[d]*mu[d]<mu[d-1]*mu[d+1])
    nwfail=0
    for d in range(1,deg):
        a0=co[d]/comb(deg,d); am=co[d-1]/comb(deg,d-1); ap=co[d+1]/comb(deg,d+1)
        if a0*a0<am*ap: nwfail+=1
    print("%-14s n=%d  mu_0(good)=%.6f  BONF5=%+.6f  p(1)=%.2e | real roots %d/%d  Newton-fails %d  logconc-fails %d"%(
        name,n,float(mu[0]),float(b5),float(p1),nreal,deg,nwfail,lcfail))
    print("   mu_d: %s"%" ".join("%.4f"%float(m) for m in mu))
    print("   S1..S6: %s"%" ".join("%.4f"%float(S[k]) for k in range(1,7)))
    return mu,S,b5

tight=list(range(1,14))
opus30Z=[420,450,510,570,690,870,1230,1770,2370,3210,4170,5190,7230]
deepwell=list(range(1,13))+[182]
geo=[5*2**i for i in range(13)]
# greedy corrected-clean: Sidon + no 3-AP + no exact ratio with num,den <= 30 + no multiple of 13, from 300 upward, odd-ish steps
def clean_packet():
    xs=[]
    x=307
    import math
    while len(xs)<13:
        ok=True
        if x%13==0 or x%7==0: ok=False
        if ok:
            for a in xs:
                for b in xs:
                    if a<b and (a+x==2*b or b+x==2*a or x+x==a+b): ok=False; break
                if not ok: break
        if ok:
            sums={}
            for a in xs+[x]:
                for b in xs+[x]:
                    if a<b:
                        if a+b in sums: ok=False; break
                        sums[a+b]=1
                if not ok: break
        if ok:
            for a in xs:
                g=gcd(a,x)
                if a//g<=30 and x//g<=30: ok=False; break
        if ok: xs.append(x)
        x+=97 if len(xs)%2 else 89
        x+=x//19
    return xs
cc=clean_packet()
print("corrected-clean packet:",cc)
for name,pk in [("tight",tight),("deepwell",deepwell),("geometric",geo),("opus30Z",opus30Z),("corr-clean",cc)]:
    analyze(name,pk)
# cross-check tight S_k against cont.34 exact integers
mu,S,_=depth_spectrum(tight),None,None
print("cross-check tight: S_13 =",sum(comb(d,13)*depth_spectrum(tight)[d] for d in range(14)),"(expect 1/91)")
print("DONE")
