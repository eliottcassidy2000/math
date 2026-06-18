# LRC(14) — MECHANISM probe for the private-q7-clearing-crossing angle (kps-S4-wf)
#
# Goal: understand EXACTLY what forces j >= D/14 at the optimal crossing tau*=j/D,
# and whether the parked multiple-of-14 runner's private q-obligation is the driver.
#
# Focus on the principal single-drop tower {1..12, w} with w=14m (parked), the tightest
# family (codex's {1..12,182}, M=14/183). Prove a closed form and locate where the bound
# j>=D/14 comes from. Then test the GENERAL claim on mixed cores.
#
# Exact rationals throughout.

from fractions import Fraction as F
from math import gcd
from itertools import combinations

def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return r if r<=F(1,2) else 1-r

def cand(S):
    S=sorted(set(S)); C=set()
    for v in S:
        k=0
        while F(2*k+1,2*v)<=F(1,2): C.add(F(2*k+1,2*v)); k+=1
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            for d in (S[i]+S[j],S[j]-S[i]):
                if d>0:
                    k=1
                    while F(k,d)<=F(1,2): C.add(F(k,d)); k+=1
    C.add(F(1,2)); return C

def Mval(S):
    b=F(0); bt=None
    for t in cand(S):
        v=min(nrm(x*t) for x in S)
        if v>b: b=v; bt=t
    return b,bt

def is_cov(S): return all(any(v%q==0 for v in S) for q in range(2,15))
def is_prim(S):
    g=0
    for v in S: g=gcd(g,v)
    return g==1

def binding_set(S,t,M): return [v for v in S if nrm(v*t)==M]

print("="*72)
print("PART A: principal tower {1..12, w=14m}, m=1..40. Closed form for M, tau*.")
print("="*72)
print(f"{'m':>3} {'w':>5} {'M':>10} {'M*14':>8} {'tau*':>12} {'bindingpair/D/j':>20} {'j>=D/14?'}")
rows=[]
for m in range(1,41):
    w=14*m
    S=sorted(set(range(1,13))|{w})
    if len(S)!=13: continue
    cov=is_cov(S); prim=is_prim(S)
    M,t=Mval(S)
    B=binding_set(S,t,M)
    # express tau* = j/D for binding pair, prefer (1,w): D=1+w
    rep=None
    for va,vb in combinations(sorted(B),2):
        for D in (va+vb,abs(va-vb)):
            if D==0: continue
            if (t*D).denominator==1:
                j=int(t*D)
                rep=(va,vb,D,j);
                if va==1 and vb==w: break
        if rep and rep[0]==1 and rep[1]==w: break
    flag = (rep is not None and 14*rep[3]>=rep[2])
    print(f"{m:>3} {w:>5} {str(M):>10} {str(M*14):>8} {str(t):>12} {str(rep):>20} {flag} cov={cov}")
    rows.append((m,w,M,t,B,rep))

# Conjecture a closed form: from codex, m=13 -> {1..12,182} gives 14/183, D=183=1+182, j=14.
# Check M = 14/(w+1) when w>= something? and tau*=14/(w+1)?
print("\nClosed-form check: is M == 14/(w+1) and tau*==14/(w+1) for large m?")
for m,w,M,t,B,rep in rows:
    pred=F(14,w+1)
    print(f"  m={m} w={w}: M={M} pred14/(w+1)={pred} match={M==pred}  tau*={t} tpredmatch={t==pred}")

print("="*72)
print("PART B: WHY j>=D/14? At tau*=j/D with binding pair (1,w), w=D-1.")
print(" ||1*tau*||=j/D (the unit runner). ||w*tau*||=||(D-1)j/D||=||-j/D||=j/D. matched.")
print(" Need j/D>=1/14 i.e. 14j>=D=w+1. The OTHER 11 runners (2..12) must also clear j/D.")
print(" The binding LEVEL is set by unit runner (1) and parked runner (w) TOGETHER.")
print(" The covering obligation q=14 (only w is mult of 14) is what put w=14m in the set;")
print(" w=14m makes D=14m+1 and the question is whether j can be as small as ceil((w+1)/14)=m+1.")
print("="*72)

# Probe: for {1..12,w}, what is the EXACT j as function of w? Print j and m+1 (=ceil(D/14)).
print(f"{'m':>3} {'w':>5} {'D':>6} {'j':>4} {'ceil(D/14)=m+1':>14} {'j-(m+1)':>8} {'M=j/D':>10}")
for m,w,M,t,B,rep in rows:
    if rep is None: continue
    va,vb,D,j=rep
    import math
    cd=(D+13)//14
    print(f"{m:>3} {w:>5} {D:>6} {j:>4} {cd:>14} {j-cd:>8} {str(F(j,D)):>10}")
