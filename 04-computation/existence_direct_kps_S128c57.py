#!/usr/bin/env python3
"""existence_direct_kps_S128c57.py -- kind-pasteur S128 cont.57.
CLOSING THE EXISTENCE STEP -- no divisor count needed.

DIRECT CONSTRUCTION.  V = P u {v1,v2}, |P| = 11 distinct positive ints, mu = min P,
M = max P, v2 = v1 + d (1 <= d <= 5), v1 > 13M.  Put  q = v1 - mu.  Then
  v1 = q + mu          => v1 = mu (mod q),        e1 = mu
  v2 = q + mu + d      => v2 = mu + d (mod q),    e2 = mu + d
and since |P| = 11 distinct positive integers, M >= mu + 10 > mu + 5 >= mu + d, so
  E := P u {e1,e2}  is contained in [mu, M],  hence e_min = mu, e_max = M.
The band interval [q/(14 mu), 13q/(14 M)] has length  q(13 mu - M)/(14 M mu),
so it contains an integer whenever  q >= 14 M mu/(13 mu - M)  (needs M < 13 mu).
Since q = v1 - mu > 13M - mu, the construction succeeds whenever
        13M - mu  >=  14 M mu/(13 mu - M).                      (SIZE)
VERIFY: the identity, the containment, and (SIZE); enumerate the exceptional finite set."""
import sys, random
from fractions import Fraction as F
from math import ceil
sys.stdout.reconfigure(line_buffering=True)
def nd(x):
    fx=x-int(x); fx=fx+1 if fx<0 else fx; return min(fx,1-fx)
def certify(P,v1,v2):
    mu=min(P); M=max(P); d=v2-v1
    q=v1-mu
    if 13*mu-M<=0: return None,'M>=13mu'
    if mu+d>M: return None,'mu+d>M'
    # residues
    e1=v1%q; e2=v2%q
    if e1!=mu or e2!=(mu+d)%q: return None,'residue-mismatch'
    lo=F(q,14*mu); hi=F(13*q,14*M)
    if hi<lo: return None,'empty-band'
    a=ceil(lo)
    if a>hi: return None,'no-integer(SIZE fails)'
    V=sorted(list(P)+[v1,v2]); t=F(a,q)
    mn=min(nd(v*t) for v in V)
    return (q,a,mn), None
print("=== the 8 critical covering families (core {2..12}) ===")
P=list(range(2,13))
crit=[(168,169),(195,196),(208,210),(221,224),(234,238),(247,252),(294,299),(308,312)]
allok=True
for v1,v2 in crit:
    r,err=certify(P,v1,v2)
    if r is None:
        print("  %d,%d : FAIL (%s)"%(v1,v2,err)); allok=False
    else:
        q,a,mn=r; ok=mn>=F(1,14)
        if not ok: allok=False
        print("  %-9s q=%-5d a=%-4d min||v t||=%-8s >=1/14: %s"%("%d,%d"%(v1,v2),q,a,mn,ok))
print("  all 8:",allok)
print()
print("=== wide sweep: 11-subsets of {1..12}, v1 in (13M, 3000], d in 1..5 ===")
random.seed(57)
tot=0; okc=0; failmodes={}
for _ in range(600):
    drop=random.randint(1,12)
    P2=[x for x in range(1,13) if x!=drop]
    M=max(P2)
    v1=random.randint(13*M+1,3000); d=random.randint(1,5); v2=v1+d
    if len(set(P2+[v1,v2]))!=13: continue
    tot+=1
    r,err=certify(P2,v1,v2)
    if r is None:
        failmodes[err]=failmodes.get(err,0)+1; continue
    q,a,mn=r
    if mn>=F(1,14): okc+=1
    else: failmodes['bad-min']=failmodes.get('bad-min',0)+1
print("  certified %d / %d ; failure modes: %s"%(okc,tot,failmodes if failmodes else "none"))
print()
print("=== the exceptional set: which (mu,M) and v1 fail (SIZE)? ===")
print("  (SIZE): q = v1-mu >= 14*M*mu/(13mu-M)")
for P3 in [list(range(2,13)), [x for x in range(1,13) if x!=12], [x for x in range(1,13) if x!=2]]:
    mu,M=min(P3),max(P3)
    if 13*mu-M<=0:
        print("  P min %d max %d : M >= 13mu -- band unavailable"%(mu,M)); continue
    need=F(14*M*mu,13*mu-M)
    v1min=13*M+1
    print("  P min %-2d max %-2d : need q >= %-8s i.e. v1 >= %-8s ; smallest allowed v1 = %d -> %s"%(
        mu,M,need,need+mu,v1min,"OK" if v1min-mu>=need else "EXCEPTIONAL for v1 in [%d, %d]"%(v1min,int(need+mu))))
print("DONE")
