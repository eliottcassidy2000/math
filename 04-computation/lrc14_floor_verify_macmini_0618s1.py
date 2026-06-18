#!/usr/bin/env python3
"""
lrc14_floor_verify_macmini_0618s1.py  (mac-mini-2026-06-18-S1)

Verify the c0=1/84 floor story:
 (1) finite-w0 good-period density rho_K -> 1/84 for the WORST shape P={1,2,3,12},
     consecutive cluster k=9; confirm > 0 at large finite w0.
 (2) a CONCRETE worst-shape covering S3 set (cluster contains a multiple of 14): exact M(S)>=1/14?
 (3) consecutive-extremality EXHAUSTIVE for k=6,7 (was k=4,5 done): 0 violations expected.
"""
import sys, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
H=F(1,14); TWO7=F(2,7)

def danger_arcs(u,h=H):
    iv=[]
    for j in range(u):
        c=F(j,u); a=(c-h/u)%1; b=(c+h/u)%1
        if a<b: iv.append((a,b))
        else: iv.append((a,F(1))); iv.append((F(0),b))
    return iv
def merge(iv):
    iv=sorted(iv); out=[]
    for a,b in iv:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out
def safe_set(A,h=H):
    if not A: return [(F(0),F(1))]
    dz=merge([iv for u in A for iv in danger_arcs(u,h)]); safe=[]; prev=F(0)
    for a,b in dz:
        if a>prev: safe.append((prev,a))
        prev=max(prev,b)
    if prev<1: safe.append((prev,F(1)))
    return safe
def meas(arcs): return sum((b-a for a,b in arcs),F(0))
def intersect(A,B):
    out=[]; i=j=0
    while i<len(A) and j<len(B):
        lo=max(A[i][0],B[j][0]); hi=min(A[i][1],B[j][1])
        if lo<hi: out.append((lo,hi))
        if A[i][1]<B[j][1]: i+=1
        else: j+=1
    return out
def good_set_exact(E):
    E=sorted(set(E)); k=len(E); diffs=set()
    for a in range(k):
        for b in range(a+1,k): diffs.add(E[b]-E[a])
    bps=set([F(0),F(1)])
    for d in diffs:
        for m in range(0,d+1): bps.add(F(m,d))
    bps=sorted(x for x in bps if 0<=x<=1); good=[]
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        pts=sorted(((e*xm)%1,e) for e in E); order=[e for _,e in pts]
        floors=[int((e*xm)//1) for e in order]
        for idx in range(k):
            e_cur=order[idx]; f_cur=floors[idx]
            if idx<k-1: e_nx=order[idx+1]; f_nx=floors[idx+1]; wrap=F(0)
            else: e_nx=order[0]; f_nx=floors[0]; wrap=F(1)
            A=F(e_nx-e_cur); Cc=F(f_cur-f_nx)+wrap
            if A==0:
                if Cc>TWO7: good.append((x0,x1))
                continue
            xb=(TWO7-Cc)/A
            lo,hi=(max(x0,xb),x1) if A>0 else (x0,min(x1,xb))
            if lo<hi: good.append((lo,hi))
    return merge(good)

# ---------- finite-w0 good-period density ----------
def rho_K_finite(P, Lspeeds):
    Vmax=max(Lspeeds); A=sorted(set(list(P)+[u for u in Lspeeds if u!=Vmax]))
    thresh=F(1,7*Vmax); safeA=safe_set(A); good=set()
    for (lo,hi) in safeA:
        jlo=int(lo*Vmax); jhi=int(hi*Vmax)+2
        for j in range(max(0,jlo),min(Vmax,jhi)):
            ilo=F(14*j+1,14*Vmax); ihi=F(14*j+13,14*Vmax)
            x0=max(lo,ilo); x1=min(hi,ihi)
            if x1-x0>thresh: good.add(j)
    return F(len(good),Vmax)

def Mval(S):
    """exact M(S)=max_tau min_v ||v tau||, via the safe-set at the achieved level (binary-ish):
       compute as the max over arc-structure: M(S)=sup{h: safe_set(S,h) nonempty}. We bisect on h."""
    # M is achieved at a binding crossing; compute max over candidate h = (a)/(2 b) crossing points.
    # Simpler exact: M(S) = max over tau of min_v ||v tau||; the max is at a point equidistant between
    # two nearest danger centers. Use the safe_set at level h and bisect h in rationals with denom<=2*maxS.
    lo,hi=F(0),F(1,2)
    # candidate levels: ||v tau|| values at crossing tau=k/(v±v'); gather, sort, test.
    S=sorted(S); cands=set([F(1,14)])
    for a in range(len(S)):
        for b in range(len(S)):
            if a==b: continue
            va,vb=S[a],S[b]
            for k1 in range(va+1):
                for k2 in range(vb+1):
                    den=va+vb
                    if den==0: continue
                    # crossing where ||va t||=||vb t|| ... too many; fallback: test h on a grid of rationals
    # fallback: bisection on h with safe_set nonemptiness (exact)
    losafe=F(0); hisafe=F(1,2);
    for _ in range(24):
        mid=(losafe+hisafe)/2
        if safe_set(S,mid): losafe=mid
        else: hisafe=mid
    return losafe  # lower bound on M (safe nonempty); good enough to confirm >=1/14

print("="*80)
print("(1) finite-w0 rho_K for WORST shape P={1,2,3,12}, consecutive cluster k=9 -> 1/84?")
print("="*80)
P=[1,2,3,12]
for w0 in (500,1500,4000,9000):
    L=[w0+d for d in range(9)]
    rk=rho_K_finite(P,L)
    print(f"   w0={w0:5d}: rho_K={float(rk):.6f}  (Vmax={max(L)})   target 1/84={1/84:.6f}")

print("\n"+"="*80)
print("(2) CONCRETE worst-shape covering S3 set: cluster has a multiple of 14, exact M>=1/14?")
print("="*80)
# choose w0 so the consecutive 9-cluster contains a multiple of 14: w0=2002 (2002=14*143) -> cluster 2002..2010 has 2002
for w0 in (2002, 5004, 9996):
    L=[w0+d for d in range(9)]
    S=sorted(set(P+L))
    has14=[u for u in S if u%14==0]
    rk=rho_K_finite(P,L)
    print(f"   S=P∪{{{w0}..{w0+8}}}: mult-of-14 in cluster={has14}, rho_K={float(rk):.5f} (#good>0 => M>=1/14)")

print("\n"+"="*80)
print("(3) consecutive-extremality EXHAUSTIVE k=6,7 (max offset <= 11)")
print("="*80)
for k in (6,7):
    cons=meas(good_set_exact(list(range(k)))); viol=0; checked=0; worst=(cons,None)
    for rest in itertools.combinations(range(1,12),k-1):
        E=[0]+list(rest); checked+=1
        m=meas(good_set_exact(E))
        if m<worst[0]: worst=(m,E)
        if m<cons: viol+=1
    print(f"   k={k}: consecutive mu={float(cons):.5f}; checked {checked}; violations={viol}; "
          f"min-mu shape={worst[1] if worst[1] else 'consecutive'} ({float(worst[0]):.5f})")
print("\nDONE.")
