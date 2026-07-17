#!/usr/bin/env python3
"""bonf_leverage_identity_kps_S128c35.py -- (1) referee BONF_m = mu_0 + (-1)^m sum_{d>m} C(d-1,m) mu_d
exactly on all 5 packets, all m; (2) exact rationals for the certified packet incl. property audit;
(3) the tight system's atom decomposition."""
import sys
from fractions import Fraction as F
from math import comb, gcd
sys.stdout.reconfigure(line_buffering=True)
def depth_spectrum(speeds):
    ev=[]
    for v in speeds:
        for j in range(v):
            lo=F(14*j-1,14*v); hi=F(14*j+1,14*v)
            if lo<0:
                ev.append((F(0),1)); ev.append((hi,-1)); ev.append((lo+1,1)); ev.append((F(1),-1))
            else:
                ev.append((lo,1)); ev.append((hi,-1))
    ev.sort()
    mu=[F(0)]*(len(speeds)+1); depth=0; last=F(0)
    for x,d in ev:
        if x>last: mu[depth]+=x-last; last=x
        depth+=d
    if F(1)>last: mu[depth]+=F(1)-last
    return mu
tight=list(range(1,14))
opus30Z=[420,450,510,570,690,870,1230,1770,2370,3210,4170,5190,7230]
deepwell=list(range(1,13))+[182]
geo=[5*2**i for i in range(13)]
cc=[307,425,541,671,800,944,1087,1413,1943,2147,2570,3056,3310]
packs=[("tight",tight),("deepwell",deepwell),("geometric",geo),("opus30Z",opus30Z),("corr-clean",cc)]
print("== (1) identity referee: BONF_m == mu_0 + (-1)^m sum_{d>m} C(d-1,m) mu_d  (exact, all m=1..12) ==")
ok=True
for name,pk in packs:
    n=len(pk); mu=depth_spectrum(pk)
    S=[sum(comb(d,k)*mu[d] for d in range(n+1)) for k in range(n+1)]
    for m in range(1,n):
        lhs=sum((-1)**k*S[k] for k in range(m+1))
        rhs=mu[0]+(-1)**m*sum(comb(d-1,m)*mu[d] for d in range(m+1,n+1))
        if lhs!=rhs: ok=False; print("  FAIL",name,m)
print("  identity holds on all 5 packets, all m:",ok)
print("== (2) THE CERTIFIED PACKET (exact rationals) ==")
mu=depth_spectrum(cc); n=13
S=[sum(comb(d,k)*mu[d] for d in range(14)) for k in range(14)]
b5=sum((-1)**k*S[k] for k in range(6))
print("  speeds:",cc)
print("  BONF5 =",b5,"=",float(b5)," > 0:",b5>0)
print("  mu_0 (good measure) =",mu[0],"=",float(mu[0]))
print("  weighted tail sum_{d>=6} C(d-1,5) mu_d =",mu[0]-b5,"=",float(mu[0]-b5))
# property audit
sums={}; sidon=True; ap3=False; ratio=None
for i in range(13):
    for j in range(i+1,13):
        s=cc[i]+cc[j]
        if s in sums: sidon=False
        sums[s]=(i,j)
for i in range(13):
    for j in range(13):
        for k in range(13):
            if i!=j and j!=k and i!=k and cc[i]+cc[k]==2*cc[j]: ap3=True
minratio=(999,999)
for i in range(13):
    for j in range(13):
        if i==j: continue
        g=gcd(cc[i],cc[j]); p,q=cc[j]//g,cc[i]//g
        if max(p,q)<minratio[0]+minratio[1]: pass
        if q<=30 and p<=30: ratio=(cc[i],cc[j],p,q)
print("  audit: Sidon(sums distinct)=%s ; has-3AP=%s ; small ratio (p,q<=30)=%s ; 7|x or 13|x: %s"%(
    sidon,ap3,ratio,[x for x in cc if x%7==0 or x%13==0]))
print("== (3) tight-system atom decomposition ==")
mu=depth_spectrum(tight)
S5=sum(comb(d,5)*mu[d] for d in range(14))
atom=comb(13,5)*mu[13]
print("  S_5 = %s = %.4f ; origin-atom contribution C(13,5)*mu_13 = %s = %.4f (%.1f%%)"%(S5,float(S5),atom,float(atom),100*float(atom/S5)))
print("  BONF5 leverage of the origin atom: -C(12,5)*mu_13 =",-comb(12,5)*mu[13],"=",float(-comb(12,5)*mu[13]))
print("  kill threshold: an atom at depth 13 of mass > equid-BONF5/792 = %.2e breaks the certificate; deep well has 1/91 = %.2e (%.0fx over)"%(0.0821/792,1/91,(1/91)/(0.0821/792)))
print("DONE")
