#!/usr/bin/env python3
"""
lrc14_mu17_extremal_macmini_0618s5.py  (mac-mini-2026-06-18-S5)

The remaining LRC(14)-S3 crux: mu_{1/7}(E) >= mu_{1/7}(consecutive_k) for all integer E.
 (I)   AP-invariance sanity: mu17(AP {a,a+d,..}) == mu17(consec) for various a,d (PROVED; check).
 (II)  EXHAUSTIVE "consecutive minimizes mu17": k=9,10,11 over bounded spread; 0 violations?
 (III) tightest non-AP competitor (which shapes get closest to consec from above).
 (IV)  multi-block separation RATE: mu17(block1 ⊔ (g+block2)) as g->inf; stays >= consec? rate?
       (the remaining tail: non-AP large-spread shapes converge to a block-average >= the min.)
"""
import sys, itertools, random
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
random.seed(6185)
SEV=F(1,7)
def good_set_exact(E, thr=SEV):
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
        fl=[int((e*xm)//1) for e in order]
        for idx in range(k):
            ec=order[idx]; fc=fl[idx]
            if idx<k-1: en=order[idx+1]; fn=fl[idx+1]; wrap=F(0)
            else: en=order[0]; fn=fl[0]; wrap=F(1)
            A=F(en-ec); Cc=F(fc-fn)+wrap
            if A==0:
                if Cc>thr: good.append((x0,x1))
                continue
            xb=(thr-Cc)/A
            lo,hi=(max(x0,xb),x1) if A>0 else (x0,min(x1,xb))
            if lo<hi: good.append((lo,hi))
    good=sorted(good); out=[]
    for a,b in good:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out
def mu17(E): return float(sum(b-a for a,b in good_set_exact(E)))
def mu17F(E): return sum((b-a for a,b in good_set_exact(E)),F(0))

print("="*78)
print("(I) AP-invariance: mu17(consec k) vs mu17(AP a+d*{0..k-1})  (must be EQUAL)")
print("="*78)
for k in (8,10):
    base=mu17F(list(range(k)))
    print(f"   k={k}: mu17(consec)={base}")
    for a,d in [(0,1),(3,1),(0,2),(5,3),(2,7)]:
        E=[a+d*j for j in range(k)]
        m=mu17F(E)
        print(f"      AP a={a} d={d}: mu17={m}  {'EQUAL' if m==base else 'DIFFERENT!'}")

print("\n"+"="*78)
print("(II) EXHAUSTIVE: consecutive minimizes mu17? k=9,10,11 (bounded spread)")
print("="*78)
for k,maxsp in [(9,13),(10,14),(11,15)]:
    cons=mu17F(list(range(k))); viol=[]; checked=0; mn=(cons,'consec')
    for rest in itertools.combinations(range(1,maxsp+1), k-1):
        E=[0]+list(rest); checked+=1
        m=mu17F(E)
        if m<mn[0]: mn=(m,E)
        if m<cons: viol.append((E,float(m)))
    print(f"   k={k} (spread<={maxsp}): mu17(consec)={float(cons):.5f}; checked {checked}; "
          f"violations={len(viol)}; min={float(mn[0]):.5f} {'at consec' if mn[1]=='consec' else 'at '+str(mn[1])}")
    for E,m in viol[:4]: print(f"       VIOL E={E}: {m:.5f}")

print("\n"+"="*78)
print("(III) tightest non-AP competitors to consec (closest from above), k=9")
print("="*78)
k=9; cons=mu17F(list(range(k))); rows=[]
for rest in itertools.combinations(range(1,14),k-1):
    E=[0]+list(rest); m=mu17F(E)
    if m>cons: rows.append((m-cons,E,float(m)))
rows.sort()
for gap,E,m in rows[:6]:
    print(f"   E={E}: mu17={m:.5f}  (excess over consec = {float(gap):.5f})")

print("\n"+"="*78)
print("(IV) multi-block separation rate: mu17({0..3} ⊔ {g+0..g+3}) for k=8, growing g")
print("="*78)
b1=[0,1,2,3]
consk=mu17F(list(range(8)))
print(f"   mu17(consec 8) = {float(consk):.5f}  (the per-k minimum)")
for g in (4,8,16,40,100,400,2000):
    E=sorted(set(b1+[g+j for j in range(4)]))
    if len(E)!=8: continue
    m=mu17F(E)
    print(f"   g={g:5d}: mu17={float(m):.5f}  {'>=consec OK' if m>=consk else 'BELOW consec!'}")
print("   (single AP-block would stay == consec exactly; two separated blocks -> block-average >= consec)")
print("\nDONE.")
