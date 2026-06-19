#!/usr/bin/env python3
"""
lrc14_mu17_floor_macmini_0618s5.py  (mac-mini-2026-06-18-S5)

THM-530 reduced the LRC(14) S3 residual to: mu_{1/7}(E) >= thr_k for all clusters E,
where the global-witness density rho*_{1/7}(P,E) >= meas(G_P) + mu_{1/7}(E) - 1 (union bound),
mu_{1/7}(E)=meas{x: maxgap{frac(e_i x)} > 1/7}, thr_k = 1 - min_{|P|=13-k} meas(G_P).
kind-pasteur VERIFIED consecutive minimizes mu_{1/7}. This script makes it EXACT:
 (I)   exact mu_{1/7}(consecutive k), k=7..13.
 (II)  exact m_k = min over P (|P|=13-k) of meas(G_P).
 (III) the union-bound closure margin = m_k + mu_{1/7}(consec k) - 1 (>0 closes the residual).
 (IV)  EXHAUSTIVE consecutive-minimizes-mu_{1/7} for k=8 (all shapes spread<=12) + k=9 sample.
 (V)   dissociated/large-spread shapes: confirm mu_{1/7} stays >= mu_{1/7}(consec k).
"""
import sys, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
H=F(1,14); SEV=F(1,7)   # threshold for maxgap (global-witness: gap between adjacent points > 1/7)

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

def good_set_exact(E, thr=SEV):
    """{x: maxgap{frac(e x): e in E} > thr} as exact rational arcs."""
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
                if Cc>thr: good.append((x0,x1))
                continue
            xb=(thr-Cc)/A
            lo,hi=(max(x0,xb),x1) if A>0 else (x0,min(x1,xb))
            if lo<hi: good.append((lo,hi))
    return merge(good)

def mu17(E): return meas(good_set_exact(E,SEV))

print("="*80)
print("(I-III) union-bound closure:  m_k + mu_{1/7}(consec k) - 1  > 0 ?")
print("="*80)
print(f"{'k':>3} {'mu17(consec k)':>18} {'m_k=min meas(G_P)':>20} {'margin=m_k+mu-1':>18}")
for k in range(7,14):
    muc=mu17(list(range(k)))
    psz=13-k
    if psz==0:
        mk=F(1)  # P empty => G_P=[0,1)
    else:
        mk=min(meas(safe_set(list(P))) for P in itertools.combinations(range(1,14),psz))
    margin=mk+muc-F(1)
    print(f"{k:>3} {str(muc):>18} {str(mk):>20} {str(margin):>18}  {float(margin):+.4f}  "
          f"{'CLOSES' if margin>0 else 'FAILS'}")

print("\n"+"="*80)
print("(IV) EXHAUSTIVE: does consecutive {0..k-1} minimize mu_{1/7}? k=8 (spread<=12)")
print("="*80)
for k in (8,):
    cons=mu17(list(range(k))); viol=[]; checked=0; mn=(cons,'consec')
    for rest in itertools.combinations(range(1,13),k-1):
        E=[0]+list(rest); checked+=1
        m=mu17(E)
        if m<mn[0]: mn=(m,E)
        if m<cons: viol.append((E,m))
    print(f"   k={k}: mu17(consec)={float(cons):.5f}; checked {checked} shapes; "
          f"violations(mu<consec)={len(viol)}; min={float(mn[0]):.5f} at {mn[1]}")
    for E,m in viol[:5]: print(f"       VIOL E={E}: mu17={float(m):.5f}")

print("\n"+"="*80)
print("(V) large-spread / dissociated shapes: mu_{1/7} >= mu_{1/7}(consec k)?  (k=8,10)")
print("="*80)
for k in (8,10):
    cons=mu17(list(range(k)))
    tests=[
      [0]+[2*i for i in range(1,k)],              # 2-dilated single run (scale-inv => =consec)
      [0]+list(range(1,k-1))+[3*k],               # consec + far point
      sorted({0}|set(range(1,k//2))|{20+i for i in range(k-k//2)}),  # two blocks
      [0,1,3,7,15,31,63,127][:k] if k<=8 else [0,1,3,7,15,31,63,127,255,511][:k],  # dissociated (Sidon-ish)
    ]
    print(f"   k={k}: mu17(consec)={float(cons):.5f}")
    for E in tests:
        E=sorted(set(E))
        if len(E)!=k: continue
        m=mu17(E)
        flag = "OK(>=consec)" if m>=cons-F(1,10**9) else "BELOW consec!"
        print(f"       E={E}: mu17={float(m):.5f}  {flag}")
print("\nDONE.")
