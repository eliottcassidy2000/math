#!/usr/bin/env python3
"""four_comb_dichotomy_kps_S128c68.py -- kind-pasteur S128 cont.68.
TOWARD THE FOUR-COMB THEOREM (uniform r=5), extending codex THM-1097 (three combs, r=4).

WHY THE CRUDE ROUTE DIES AT FOUR COMBS.  With the sharp one-comb bound
|I cap D_k| <= |I|/7 + 6/(49k) and the component count <= 1 + sum(ell*k_i + 1), the
longest surviving component obeys
    L >= [ (7-r) ell/7 - (6/49) sum 1/k_i ] / [ (r+1) + ell sum k_i ] ,
and asking L > 1/(7 k_max) needs (7-r) > r, i.e. r < 3.5.  So three combs work and four
do not -- exactly where THM-1097 stops.

THE PROPOSED REPLACEMENT: a GAP RECURSION with a spread/clustered dichotomy.
Inside a component of length lam, removing D_k leaves a full gap of length 6/(7k) whenever
lam >= 1/k (the component contains a whole period).  Iterating,
    L_j >= 6/(7 k_j)   provided   L_{j-1} >= 1/k_j ,  i.e.  k_j >= (7/6) k_{j-1}.
So if every step is SPREAD (ratio >= 7/6) the four-comb bound holds with a factor 6 to
spare.  The residual is the CLUSTERED case, where some k_j < (7/6)k_{j-1} -- and there the
two combs nearly coincide, so their union should cost far less than two independent combs.
Measure both halves.  PRINT DATA ONLY."""
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
print("### (0) the atlas: longest component of S(P) over the 495 eight-speed cores ###")
ells=[]
for P in C8:
    iv=safe_set(P); ells.append((max(b-a for a,b in iv),tuple(P)))
ells.sort()
print("  min ell = %s = %.6f at core %s ; max ell = %.6f"%(
    ells[0][0],float(ells[0][0]),list(ells[0][1]),float(ells[-1][0])))
print("  min ell * (13*maxP+1) over cores: %.4f"%min(float(e)*(13*max(p)+1) for e,p in ells))
print()
print("### (1) SPREAD half: every step ratio >= 7/6 ==> L >= 6/(7 k4)? ###")
random.seed(68)
ok=0; tot=0; worst=None
for _ in range(300):
    P=random.choice(C8); M=max(P); lo=13*M+1
    k1=random.randint(lo,lo+400)
    ks=[k1]
    for _ in range(3):
        ks.append(int(ks[-1]*F(7,6))+random.randint(0,60))
    iv=safe_set(P)
    for k in ks: iv=rm(iv,k)
    if not iv: continue
    tot+=1
    L=max(b-a for a,b in iv)
    r=L*7*ks[-1]      # want > 1 for the theorem, >= 6 for the recursion claim
    if L>F(1,7*ks[-1]): ok+=1
    if worst is None or r<worst[0]: worst=(r,tuple(P),tuple(ks))
print("  spread quadruples: %d ; satisfying L > 1/(7k4): %d"%(tot,ok))
print("  worst 7*k4*L = %.4f (theorem needs > 1; recursion predicts >= 6) at core %s ks=%s"%(
    float(worst[0]),list(worst[1]),worst[2]))
print()
print("### (2) CLUSTERED half: some ratio < 7/6 -- is the comb union cheap? ###")
print("  measure |D_a cup D_b| within a unit window vs 2/7 (independent) for close a,b")
for a in [157,300,701]:
    for rat in ['1.00(a,a+1)','1.02','1.05','1.10','1.16(=7/6)','1.50','3.00']:
        if rat.startswith('1.00'): b=a+1
        else: b=int(a*float(rat.split('(')[0]))
        if b<=a: b=a+1
        P=[1,2,3,4,5,6,7,8]
        iv=[(F(0),F(1))]
        m1=sum(y-x for x,y in rm(iv,a))
        m2=sum(y-x for x,y in rm(rm(iv,a),b))
        cost1=1-m1; cost2=1-m2
        print("  a=%-5d b=%-6d  |D_a|=%.5f  |D_a u D_b|=%.5f  vs 2/7=%.5f  ratio=%.3f"%(
            a,b,float(cost1),float(cost2),2/7,float(cost2)/(2/7)))
    print()
print("### (3) how often is a real quadruple clustered at some step? ###")
cl=0; tot2=0
for P in C8[:60]:
    M=max(P); lo=13*M+1
    for _ in range(60):
        ks=sorted(random.sample(range(lo,lo+900),4))
        tot2+=1
        if any(ks[i+1]<F(7,6)*ks[i] for i in range(3)): cl+=1
print("  quadruples with at least one clustered step (ratio < 7/6): %d / %d = %.3f"%(cl,tot2,cl/tot2))
print("DONE")
