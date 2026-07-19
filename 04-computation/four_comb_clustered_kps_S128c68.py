#!/usr/bin/env python3
"""four_comb_clustered_kps_S128c68.py -- kind-pasteur S128 cont.68.
THE HARD HALF.  The spread half of the four-comb dichotomy holds with margin ~4.95.
73% of real quadruples are CLUSTERED at some step (ratio < 7/6), and the comb-union saving
there is only ~7% -- far less than hoped.  So measure the clustered case directly: what is
7*k4*L for fully clustered quadruples, where L is the longest component of
S(P) minus (D_k1 u D_k2 u D_k3 u D_k4)?  The four-comb theorem needs 7*k4*L > 1.
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
random.seed(681)
print("### fully CLUSTERED quadruples: 7*k4*L, theorem needs > 1 ###")
print("  regime                    tested  min 7k4L   median   #below 1")
for tag,spread in [("consecutive k,k+1,k+2,k+3",1),("step<=3",3),("step<=8",8),("step<=25",25)]:
    vals=[]; bad=0
    for _ in range(220):
        P=random.choice(C8); M=max(P); lo=13*M+1
        k1=random.randint(lo,lo+700)
        ks=[k1]
        for _ in range(3): ks.append(ks[-1]+random.randint(1,spread))
        iv=safe_set(P)
        for k in ks: iv=rm(iv,k)
        if not iv: bad+=1; vals.append(F(0)); continue
        L=max(b-a for a,b in iv)
        vals.append(7*ks[-1]*L)
    vals.sort()
    below=sum(1 for v in vals if v<=1)
    print("  %-25s %-7d %-10.4f %-8.4f %d"%(tag,len(vals),float(vals[0]),float(vals[len(vals)//2]),below))
print()
print("### the worst clustered cases: what do they look like? ###")
worst=None
for _ in range(900):
    P=random.choice(C8); M=max(P); lo=13*M+1
    k1=random.randint(lo,lo+900)
    ks=[k1]
    for _ in range(3): ks.append(ks[-1]+random.randint(1,4))
    iv=safe_set(P)
    for k in ks: iv=rm(iv,k)
    if not iv:
        worst=(F(0),tuple(P),tuple(ks)); break
    L=max(b-a for a,b in iv)
    v=7*ks[-1]*L
    if worst is None or v<worst[0]: worst=(v,tuple(P),tuple(ks))
v,P,ks=worst
print("  min 7*k4*L over 900 tight-clustered quadruples: %.4f"%float(v))
print("    at core %s, killers %s"%(list(P),ks))
print("    (theorem needs > 1 ; spread half measured 4.95)")
print()
print("### does S(P) minus four combs ever become EMPTY? ###")
emp=0; tot=0
for _ in range(500):
    P=random.choice(C8); M=max(P); lo=13*M+1
    ks=sorted(random.sample(range(lo,lo+1200),4))
    iv=safe_set(P)
    for k in ks: iv=rm(iv,k)
    tot+=1
    if not iv: emp+=1
print("  empty remainders: %d / %d"%(emp,tot))
print("DONE")
