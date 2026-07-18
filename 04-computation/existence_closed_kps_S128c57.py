#!/usr/bin/env python3
"""existence_closed_kps_S128c57.py -- kind-pasteur S128 cont.57.
CLOSING THE EXISTENCE STEP OF THM-1018 -- unconditionally, with an EXPLICIT modulus.

THE CONSTRUCTION.  V = P u {v1,v2}, P = core of 11 distinct positive integers,
mu = min P, M = max P, v2 = v1 + d (1 <= d <= 5), v1 > 13M.  Put

        q = v1 + M ,      a = ceil( q / (14 mu) ) ,      t = a/q .

RESIDUES.  v1 < q so v1 = v1 (mod q) and q - v1 = M, hence e1 = M.
v2 = v1 + d < q so q - v2 = M - d, hence e2 = M - d.  Each core p <= M < q/2, so e_p = p.
Since 11 distinct positive integers force M - mu >= 10 > 5 >= d, we get M - d in [mu, M],
so E = P u {M, M-d} lies in [mu, M]:   e_min = mu,  e_max = M.   (RATIO INHERITED)

BAND (THM-1018 II).  Any integer a in [q/(14 e_min), 13q/(14 e_max)] certifies >= 1/14.
Length = q(13 mu - M)/(14 M mu) >= 1  iff  q >= 14 M mu/(13 mu - M).
Since q = v1 + M >= 14M + 1, this holds whenever  M <= 12 mu  (then 13mu - M >= mu).

  ==> CRITERION:  M <= 12*mu.   For any core P inside {1,...,12}: M <= 12 <= 12*mu.  ALWAYS.

General threshold 1/N (here N=14): criterion is M <= (N-2)*mu.
Compare the earlier q = v1 - mu construction: criterion M/mu <= (N-2+sqrt(N^2-4N))/2
= 6+sqrt(35) = 11.9161 at N=14 -- strictly weaker.  Both verified below."""
import sys, random, itertools
from fractions import Fraction as F
from math import ceil, isqrt
sys.stdout.reconfigure(line_buffering=True)
def nd(x):
    fx=x-int(x); fx=fx+1 if fx<0 else fx; return min(fx,1-fx)
def la(v,q):
    r=v%q; return min(r,q-r)
def cert(P,v1,v2,N=14):
    """the construction; returns (q,a,minnorm,e1,e2,emin,emax) or (None,reason)"""
    mu=min(P); M=max(P); d=v2-v1
    q=v1+M; a=ceil(F(q,N*mu))
    e1=la(v1,q); e2=la(v2,q)
    E=sorted(set(list(P)+[e1,e2])); emin,emax=min(E),max(E)
    hi=F((N-1)*q,N*emax)
    if a>hi: return None,('no-integer',q,a,float(hi))
    V=sorted(list(P)+[v1,v2]); t=F(a,q)
    mn=min(nd(v*t) for v in V)
    return (q,a,mn,e1,e2,emin,emax),None
print("### 1. THE IDENTITY: e1 = M, e2 = M-d, E in [mu,M] ###")
bad=0
for _ in range(4000):
    k=random.randint(1,12); P=[x for x in range(1,13) if x!=k]
    mu,M=min(P),max(P)
    v1=random.randint(13*M+1,10**6); d=random.randint(1,5); v2=v1+d
    q=v1+M
    if la(v1,q)!=M or la(v2,q)!=M-d: bad+=1; continue
    E=set(P)|{la(v1,q),la(v2,q)}
    if min(E)!=mu or max(E)!=M: bad+=1
print("  4000 random families: identity+containment failures =",bad)
print()
print("### 2. THE 8 CRITICAL COVERING FAMILIES (THM-1011 VII) ###")
P=list(range(2,13)); allok=True
for v1,v2 in [(168,169),(195,196),(208,210),(221,224),(234,238),(247,252),(294,299),(308,312)]:
    r,err=cert(P,v1,v2)
    if r is None: print("  %d,%d FAIL %s"%(v1,v2,err)); allok=False; continue
    q,a,mn,e1,e2,emin,emax=r; ok=mn>=F(1,14); allok&=ok
    print("  %-9s q=%-5d a=%-3d t=%-9s e=(%d,%d) min||vt||=%-8s >=1/14:%s"%(
        "%d,%d"%(v1,v2),q,a,F(a,q),e1,e2,mn,ok))
print("  all 8:",allok)
print()
print("### 3. THE FORMER EXCEPTIONAL SET (1,12 in core, v1 in [157,168]) ###")
tot=0; ok=0; worst=None
for k in range(2,12):
    P2=[x for x in range(1,13) if x!=k]
    for v1 in range(157,169):
        for d in range(1,6):
            v2=v1+d
            if len(set(P2+[v1,v2]))!=13: continue
            tot+=1
            r,err=cert(P2,v1,v2)
            if r and r[2]>=F(1,14):
                ok+=1
                if worst is None or r[2]<worst[0]: worst=(r[2],P2,v1,v2,r[0],r[1])
print("  %d/%d certified by q = v1+M  (previously ALL exceptional)"%(ok,tot))
if worst: print("  tightest: min||vt||=%s at core-drop core min%d max%d, v=(%d,%d), q=%d a=%d"%(
    worst[0],min(worst[1]),max(worst[1]),worst[2],worst[3],worst[4],worst[5]))
print()
print("### 4. EXHAUSTIVE over the whole stratum shape, v1 in (13M, 13M+4000] ###")
tot=0; ok=0; fails=[]
for k in range(1,13):
    P3=[x for x in range(1,13) if x!=k]; M=max(P3)
    for v1 in range(13*M+1,13*M+4001):
        for d in range(1,6):
            v2=v1+d
            if len(set(P3+[v1,v2]))!=13: continue
            tot+=1
            r,err=cert(P3,v1,v2)
            if r and r[2]>=F(1,14): ok+=1
            else: fails.append((k,v1,d,err if r is None else float(r[2])))
print("  certified %d / %d ; failures: %d"%(ok,tot,len(fails)))
if fails: print("  first fails:",fails[:5])
print()
print("### 5. GENERAL CORES (not inside {1..12}): criterion M <= 12*mu ###")
print("  ratio   n_tested  certified   (predicted: OK iff M <= 12 mu)")
for ratio_t in [(2,20),(2,24),(2,25),(3,36),(3,38),(5,60),(5,63),(1,12),(4,49)]:
    mu,M=ratio_t
    if M-mu<10: continue
    n=0; c=0
    for _ in range(400):
        mid=random.sample(range(mu+1,M),9)
        P4=sorted([mu]+mid+[M])
        if len(P4)!=11: continue
        v1=random.randint(13*M+1,13*M+5000); d=random.randint(1,5); v2=v1+d
        if len(set(P4+[v1,v2]))!=13: continue
        n+=1
        r,err=cert(P4,v1,v2)
        if r and r[2]>=F(1,14): c+=1
    print("  %2d/%-2d=%-6.3f  %4d      %4d        predicted %s"%(
        M,mu,M/mu,n,c,"ALL OK" if M<=12*mu else "FAIL (ratio>12)"))
print()
print("### 6. THE GENERAL-N LAW: criterion M <= (N-2) mu ###")
print("   N   (N-2)   old-construction bound (N-2+sqrt(N^2-4N))/2   new is better")
for N in [6,8,10,12,14,16,20]:
    old=(N-2+(N*N-4*N)**0.5)/2
    print("  %2d    %2d           %8.4f                        %s"%(N,N-2,old,N-2>=old))
print("DONE")
