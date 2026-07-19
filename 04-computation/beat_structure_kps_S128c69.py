#!/usr/bin/env python3
"""beat_structure_kps_S128c69.py -- kind-pasteur S128 cont.69.
THE BEAT STRUCTURE for the clustered majority (THM-1140's named next).

SETUP.  Two combs k_i, k_j differ in phase by frac((k_j-k_i) t): the relative phase drifts
with period 1/(k_j-k_i).  For r combs at a point t the configuration is the phase vector
(frac(k_1 t), ..., frac(k_r t)) mod 1.
  ALIGNED phases  -> teeth coincide -> gap ~ 6/(7k)  = 12/(14k)
  SPREAD phases   -> teeth interleave at spacing 1/(rk), width 1/(7k) each
                  -> gap ~ 1/(r k) - 1/(7k),  which at r=4 is 3/(28k)
and the four-comb theorem needs > 1/(7k) = 4/(28k).  So EQUALLY-SPREAD PHASES ARE FATAL
at r=4 (3/28 < 4/28) and the whole question is whether they can persist.

SWEEP TEST.  Across a core component of length ell the phase pair (i,j) sweeps ell*(k_j-k_i)
full cycles.  If that is >= 1 the phases cannot stay spread and alignment must occur
somewhere inside.  If it is << 1 the configuration is FROZEN and whichever configuration
holds, holds throughout.  Measure which regime the real clustered quadruples are in.
PRINT DATA ONLY."""
import sys, itertools, random
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
def nd(x):
    fx=x-int(x); fx=fx+1 if fx<0 else fx; return min(fx,1-fx)
def safe_set(P):
    bps={F(0),F(1)}
    for p in P:
        for j in range(p+1):
            for s in (F(1,14*p),-F(1,14*p)):
                v=F(j,p)+s
                if 0<=v<=1: bps.add(v)
    B=sorted(bps); out=[]
    for i in range(len(B)-1):
        a,b=B[i],B[i+1]
        if b<=a: continue
        if all(nd(p*((a+b)/2))>=F(1,14) for p in P): out.append((a,b))
    mg=[]
    for a,b in out:
        if mg and mg[-1][1]==a: mg[-1]=(mg[-1][0],b)
        else: mg.append((a,b))
    return mg
def rm(iv,k):
    out=[]
    for a,b in iv:
        jlo=int(a*k)-1; jhi=int(b*k)+1; cur=a
        for j in range(jlo,jhi+1):
            x=F(j,k)-F(1,14*k); y=F(j,k)+F(1,14*k)
            if y<=a or x>=b: continue
            x=max(x,a); y=min(y,b)
            if x>cur: out.append((cur,x))
            cur=max(cur,y)
        if cur<b: out.append((cur,b))
    return out
C8=[sorted(c) for c in itertools.combinations(range(1,13),8)]
print("### (1) SWEEP TEST: does the phase sweep inside a core component? ###")
print("  ell*(k_j-k_i) >= 1 means alignment is FORCED somewhere inside the component")
print("  core                       ell       killers            min sweep  max sweep  regime")
random.seed(69)
frozen=0; sweep=0
for tag,ks in [("worst of cont.68",(371,374,377,379)),("consecutive",(400,401,402,403)),
               ("step 5",(400,405,410,415)),("step 40",(400,440,480,520))]:
    for P in [[1,3,5,6,7,8,11,12],[2,3,4,5,6,7,8,9]]:
        iv=safe_set(P); ell=max(b-a for a,b in iv)
        ds=[ks[j]-ks[i] for i in range(4) for j in range(i+1,4)]
        sw=[float(ell)*d for d in ds]
        reg="FROZEN" if max(sw)<1 else ("partial" if min(sw)<1 else "SWEEPS")
        if max(sw)<1: frozen+=1
        else: sweep+=1
        print("  %-26s %-9.5f %-18s %-10.3f %-10.3f %s"%(str(P),float(ell),str(ks),min(sw),max(sw),reg))
print()
print("### (2) the worst case dissected: where does its surviving component sit? ###")
P=[1,3,5,6,7,8,11,12]; ks=(371,374,377,379)
iv0=safe_set(P); iv=iv0
for k in ks: iv=rm(iv,k)
L=max(b-a for a,b in iv)
seg=[s for s in iv if s[1]-s[0]==L][0]
mid=(seg[0]+seg[1])/2
print("  longest surviving component: [%s, %s], length %s = %.7f"%(seg[0],seg[1],L,float(L)))
print("  7*k4*L = %.4f"%float(7*ks[-1]*L))
print("  phases frac(k*t) at its midpoint (aligned means all near 0 or all near 1):")
for k in ks:
    ph=(k*mid)%1
    print("     k=%-5d frac = %-10.5f  dist to nearest integer = %.5f"%(k,float(ph),float(nd(k*mid))))
print("  gap predictions:  aligned 12/(14k4) = %.7f ;  4-spread 3/(28k4) = %.7f"%(
    12/(14*ks[-1]),3/(28*ks[-1])))
print("  needed for the theorem: 1/(7k4) = %.7f"%(1/(7*ks[-1])))
print()
print("### (3) across clustered quadruples: is the surviving gap ALIGNED-like or SPREAD-like? ###")
print("  classify by L*k4:  >=0.6 aligned-like, <=0.25 spread-like")
al=0; sp=0; mid_=0; worst=None
for _ in range(500):
    P=random.choice(C8); M=max(P); lo=13*M+1
    k1=random.randint(lo,lo+800)
    kk=[k1]
    for _ in range(3): kk.append(kk[-1]+random.randint(1,5))
    iv=safe_set(P)
    for k in kk: iv=rm(iv,k)
    if not iv: continue
    L=max(b-a for a,b in iv); v=float(L*kk[-1])
    if v>=0.6: al+=1
    elif v<=0.25: sp+=1
    else: mid_+=1
    if worst is None or v<worst[0]: worst=(v,tuple(P),tuple(kk))
print("  aligned-like (L*k4>=0.6): %d ; middle: %d ; spread-like (<=0.25): %d"%(al,mid_,sp))
print("  min L*k4 = %.4f at core %s killers %s"%(worst[0],list(worst[1]),worst[2]))
print("  theorem needs L*k4 > 1/7 = %.4f ; 4-spread prediction 3/28 = %.4f"%(1/7,3/28))
print("DONE")
